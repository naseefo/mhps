

import cProfile
import click
import re
import data.defaults.params as dflt
import glob, os
import pandas as pd
from data.defaults.param_manager import *
from data.defaults.param_manager import default
from mhps.earthquake import get_earthquake_list, get_total_excitations
from mhps.fixed import read_const_param, read_ss_var_param, get_total_variable_parameter_set, superstructure_propxy, fixed_simulator
from progressbar import progressbar
# from tqdm import tqdm
# from tqdm import tqdm_gui
from mhps.postprocessor import result_viewer, createfolder
from pathlib import Path
import numpy as np

def profile(fnc):
    
    """A decorator that uses cProfile to profile a function"""
    
    def inner(*args, **kwargs):
        
        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())
        return retval

    return inner

# Load default values
default_values = get_default_param_values()


@click.group()
def cli1():
    pass


@cli1.command()
@click.option('--const_param', '-cp',type=str, default= 'ConstantParameters.txt', help= 'File name containing constant parameters for the structure')
@click.option('--var_param', '-vp',type=str, default= 'SuperstructureVariableParameters.csv', help= 'File name containing variable parameters for the structure')
@click.option('--earthquakes','-eq', type=(str), default= "Excitations.csv", help= 'Earthquakes')
@click.option('--knor', '-knor', type=int, default=1, help="Normalizing the superstructure for given time period. 1 for normalized and 0 for un-normalized.")
@click.option('--results_type', '-r', type=str, default="aa1", help="Choice to select output results")
@click.option('--lxy', '-lxy', type=int, default=0)
@click.option('--folder', '-f', type=str, default="Result", help="Folder name to store result")
@click.option('--outputunits', type=list, default=['m/s2', 'cm/s', 'cm', 'kn', 'j'])
def fixed(const_param, var_param, earthquakes, knor, results_type, lxy, folder, outputunits):


    """
    Example command:
    
    >> mhps fixed -cp Constant_Parameters.txt -vp Variable_Parameters.txt \
        -eq 1 Excitations.csv -res [1 1] -ndiv 120

    where, 

    (1) Constant_Parameters.txt is the file name containing constant
     parameters for simulating bi-directional fixed response. Alternatively, 
     file names can also be chosen.

    (2) Variable_Parameters.txt is the file name containing variable 
    parameters arrange row-wise in for parametric studies of the 
    fixed-base structure. Alternatively, file names can also be chosen.

    (3) Excitations.csv is the file name containing the list of earthquakes
    to be simulated for the fixed-base response of the structure. There 
    are multiple options of specifying earthquakes.

    Ex.1: Specify multiple bi-directional earthquakes as 
    
    >> ... -eq 2 [eqx1, eqy1, dt1, unit1, dur1], [eqx2, eqy2, dt2, unit2, dur2],...

    eqx1: File name of earthquake data arranged row-wise to be
    applied in x-direction (if not extension is specified,by 
    default file has a *.txt extension).

    eqy1: File name of earthquake data arranged row-wise to be
    applied in y-direction (if not extension is specified,by 
    default file has a *.txt extension).

    dt1: Time step at which data is sampled in eqx1 and eqy1.

    unit1: Unit of measurement of acceleration data in eqx1 and
    eqy1.

    dur1: Duration of earthquake to be considered. Specifying 0
    will select the entire duration of the earhquake

    If no option for earthquake or harmonic function is not
    specified, the program automatically applies 1940 Elcentro
    in both x- and y- direction from earthquake file available 
    in the database, in this case - Cent_acc_00.txt in x-
    direction and Cent_acc_90.txt in y-direction.

    """

    folder = createfolder(folder)
    os.makedirs('results\\' + folder + '\\Time History Response')
    Path('results\\' + folder + "\\SimulationInfo.csv").touch()

    if earthquakes:
        total_eq = get_total_excitations(earthquakes)
        earthquake_generator = get_earthquake_list(earthquakes)
    else:
        earthquake_generator = get_earthquake_list('db elc')

    total_param = get_total_variable_parameter_set(var_param)
    print("Total Parameters: " + str(total_param))
   
    maxnst, am, ak = read_const_param(const_param)

    # for i in progressbar(range(total_eq), redirect_stdout=True):
    peakvalues = None
    for i in range(total_eq):
        ref, xg, yg, dt, ndiv, ndt = next(earthquake_generator)
        # print(xg)
        # print(max(xg))
        # print(yg)
        # print("Excitation Ref. No. : " + str(ref))
        superstructure_param_generator = read_ss_var_param(var_param)
        # for j in progressbar(range(total_param), redirect_stdout=True):
        # for j in tqdm_gui(range(total_param)):
        for j in range(total_param):
            ijk, nst, tx1, zeta, rtytx = next(superstructure_param_generator)
            if ((i == 0) and (j==0)) or ((p_nst, p_tx1, p_rtytx, p_zeta != nst, tx1, rtytx, zeta)):
                smx, skx, cdx, smy, sky, cdy = superstructure_propxy(nst, tx1, rtytx, am, ak, zeta, knor)
                p_nst, p_tx1, p_rtytx, p_zeta = nst, tx1, rtytx, zeta
            result, model = fixed_simulator(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, nst, smx[0:nst,0:nst], skx[0:nst,0:nst], cdx[0:nst,0:nst], smy[0:nst,0:nst], sky[0:nst,0:nst], cdy[0:nst,0:nst])
            analysis_folder = 'results\\' + folder + '\\Time History Response\\' + 'ANA-EQ-' + str(result.eq_ref) + '-PARAM-' + str(result.ijk)
            os.makedirs(analysis_folder)
            peakvaluesparam, peakhead = result_viewer(result, model, results_type, analysis_folder)
         
            if peakvaluesparam is not None:
                if j == 0:
                    peakvalues = peakvaluesparam
                else:
                    peakvalues = np.vstack((peakvalues, peakvaluesparam))


        
        if peakvalues is not None:
            if i == 0:
                # Path('results\\' + folder + "\\Peak.csv").touch()
                peakmat = peakvalues
                peakmathead = [x + ", EQ-" + str(result.eq_ref) for x in peakhead]
            else:
                peakmat = np.hstack((peakmat, peakvalues))
                peakmathead += [x + ", EQ-" + str(result.eq_ref) for x in peakhead]
            if i == (total_eq - 1):
                peakmat = pd.DataFrame(peakmat)
                peakmat.columns = peakmathead
                peakmat.to_csv('results\\' + folder + "\\Peak.csv", mode='w', sep=',')

    
    return None


@cli1.command()
def spectral():
    pass

@cli1.command()
def fft():
    pass

@cli1.command()
def biso():
    print("I am in base isolated")
    pass




cli = click.CommandCollection(sources=[cli1, default])

if __name__ == '__main__':
    cli()


# TO DO: Find SA, SD, SV graphs for EQs

# TO DO: Identify modal frequencies for fixed base structure

# TO DO: Identify modal frequencies for base-isolated structure

# TO DO: Run time history analysis for fixed base symmetric structure

# TO DO: Run time history analysis for base-isolated symmetric structure

# TO DO: Run time history analysis for fixed base asymmetric structure

# TO DO: Run time history analysis for base-isolated asymmetric structure

# TO DO: Impact facility at isolator level


# TO DO: Isolator Types:
#   1.  OSBI
#   2.  PF
#   3.  CRR
#   4.  ERR
#   5.  FPS
#   6.  VFPI
#   7.  LRB
#   8.  HDR
#   9.  NZS
#   10. EDF
