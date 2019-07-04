#!/usr/bin/env python3

import cProfile
import pstats
import click
import re
import data.defaults.params as dflt
import glob, os
import pandas as pd
from data.defaults.param_manager import *
from data.defaults.param_manager import default
from mhps.earthquake import get_earthquake_list, get_total_excitations
from mhps.fixed import read_const_param, read_ss_var_param, get_total_variable_parameter_set, superstructure_propxy, fixed_simulator
from mhps.fixedtorsion import read_const_param_torsion, read_ss_torsion_var_param, get_total_variable_parameter_set_torsion, superstructure_torsion_propxy
from mhps.isolator_l import addlinear_iso, simulator_linear, read_iso_l_var_param
from mhps.isolator_pf import simulator_pf, read_iso_pf_var_param, IsoPFModel
from mhps.isolator_boucwen import simulator_boucwen, read_iso_boucwen_var_param, IsoBoucWenModel, addboucwen_iso
from mhps.isolator_osbi import simulator_osbi, read_iso_osbi_var_param, IsoOSBIModel
from mhps.isolator_osbix import simulator_osbix
from progressbar import progressbar
from mhps.postprocessor import result_viewer, createfolder
from pathlib import Path
import numpy as np
import io
import csv
from mhps.awsmanager import upload
import shutil
import matplotlib.pyplot as plt
import time
import datetime
from mhps.colorfn import prGreen, prCyan


np.seterr(all='warn')

def profile(fnc):
    
    """A decorator that uses cProfile to profile a function"""
    
    def inner(*args, **kwargs):
        
        pr = cProfile.Profile()
        pr.enable()
        retval = fnc(*args, **kwargs)
        pr.disable()
        s = io.StringIO()
        sortby = 'time'
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
@click.option('--var_param', '-vp',type=str, default= 'SS_VP.csv', help= 'File name containing variable parameters for the structure')
@click.option('--earthquakes','-eq', type=(str), default= "Excitations.csv", help= 'Earthquakes')
@click.option('--knor', '-knor', type=int, default=1, help="Normalizing the superstructure for given time period. 1 for normalized and 0 for un-normalized.")
@click.option('--results_type', '-r', type=str, default="g*, aa1, paa1", help="Choice to select output results")
@click.option('--lxy', '-lxy', type=int, default=0)
@click.option('--folder', '-f', type=str, default="result", help="Folder name to store result")
@click.option('--folder', '-f', type=str, default="result", help="Folder name to store result")
@click.option('--outputunits', type=list, default=['m/s2', 'cm/s', 'cm', 'kn', 'j'])
@click.option('--screen/--no-screen', default=False)
def fixed(const_param, var_param, earthquakes, knor, results_type, lxy, folder, outputunits, screen):
    fixed_fn(const_param, var_param, earthquakes, knor, results_type, lxy, folder, outputunits, screen)
    return None

# @profile
def fixed_fn(const_param, var_param, earthquakes, knor, results_type, lxy, folder, outputunits, screen):


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
    upload_preference = input('Do you wish to upload (default: no) ? [y/n] : ')

    # RESULT FOLDER SETUP
    folder = createfolder(folder)
    os.makedirs(os.path.join('results', folder, 'Time History Response'))
    # os.makedirs('results\\' + folder + '\\Time History Response')
    simulationdesc = input('Enter simulation description [None] : ')
    if simulationdesc:
        Path(os.path.join('results', folder, "SimulationInfo.csv")).touch()
        # Path('results\\' + folder + "\\SimulationInfo.csv").touch()
        with open(os.path.join('results', folder, "SimulationInfo.csv"), 'w') as f:
            for line in simulationdesc:
                f.write(line)
        f.close()

    # EARTHQUAKE GENERATOR SETUP
    if earthquakes:
        total_eq = get_total_excitations(earthquakes)
        earthquake_generator = get_earthquake_list(earthquakes)
    else:
        earthquake_generator = get_earthquake_list('db elc')

    total_param = get_total_variable_parameter_set(var_param)
    print("Total Parameters: " + str(total_param))
    
    # CONSTANT PARAMETER SETUP
    maxnst, am, ak = read_const_param(const_param)

    s_rate2 = 0.33
    s_rate = 1.5
    tot_time_eq2 = 0.0
    nparam = 1.0

    peakvalues = None
    start_t_eq = time.time()
    for i in range(total_eq):
        start_t = time.time()
        ref, xg, yg, zg, dt, ndiv, ndt = next(earthquake_generator)

        teq2 = (ndt-1)*dt*ndiv
        tot_time_eq2 = (tot_time_eq2 + teq2)/ref
        coeff1 = (0.0379*pow(ndiv, 1.5073))*pow(teq2, -0.6251*pow(ndiv, 0.125))
        time_current_eq = nparam*teq2*coeff1    #max(s_rate*total_param, s_rate2*teq2*total_param)
        time_remaining = (total_eq-ref+1)*(tot_time_eq2*coeff1*nparam)
        print('FACT: Duration of EQ = %8.4fs'%(teq2) + " | Start Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print('EXPECTATION: Analysis Duration = %8.2fs (%8.2f m) | Total Duration Remaining = %8.2fs (%6.2f m, %5.2fh)'%(time_current_eq, time_current_eq/60, time_remaining, time_remaining/60, time_remaining/60/60) + ' | ' +prGreen('Estimated Completion Time: ') + (datetime.datetime.now()+datetime.timedelta(seconds = time_current_eq)).strftime("%Y-%m-%d %H:%M:%S"))

        # SUPERSTRUCTURE VARIABLE PARAMETER SETUP
        superstructure_param_generator = read_ss_var_param(var_param)
        for j in range(total_param):
            ijk, nst, tx1, zeta, rtytx = next(superstructure_param_generator)

            if ((i == 0) and (j==0)) or ((p_nst, p_tx1, p_rtytx, p_zeta != nst, tx1, rtytx, zeta)):
                smx, skx, cdx, smy, sky, cdy = superstructure_propxy(nst, tx1, rtytx, am, ak, zeta, knor)
                print(smx)
                print(skx)
                p_nst, p_tx1, p_rtytx, p_zeta = nst, tx1, rtytx, zeta
                
            result, model = fixed_simulator(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, nst, smx[0:nst,0:nst], skx[0:nst,0:nst], cdx[0:nst,0:nst], smy[0:nst,0:nst], sky[0:nst,0:nst], cdy[0:nst,0:nst], screen)
            
            # analysis_folder = 'results\\' + folder + '\\Time History Response\\' + 'ANA-EQ-' + str(result.eq_ref) + '-PARAM-' + str(result.ijk)
            
            peakvaluesparam, peakvaluesparamhead = result_viewer(result, model, results_type, folder)
         
            if peakvaluesparam is not None:
                if j == 0:
                    peakvalues = peakvaluesparam
                else:
                    peakvalues = np.vstack((peakvalues, peakvaluesparam))


        
        if peakvalues is not None:
            if i == 0:
                # Path('results\\' + folder + "\\Peak.csv").touch()
                peakmat = peakvalues
                peakmathead = [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            else:
                peakmat = np.hstack((peakmat, peakvalues))
                peakmathead += [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            if i == (total_eq - 1):
                if total_param == 1:
                    peakmat = pd.DataFrame(peakmat.reshape(-1, len(peakmat)))
                else:
                    peakmat = pd.DataFrame(peakmat)
                # print(peakmat)
                # print(peakvaluesparamhead)
                peakmat.columns = peakmathead
                peakmat.to_csv(os.path.join("results", folder, "Peak.csv"), mode='w', sep=',', index=False)
                # peakmat.to_csv('results\\' + folder + "\\Peak.csv", mode='w', sep=',', index=False)
        
        end_t = time.time()
        eq_dur = end_t - start_t
        s_rate = eq_dur/total_param
        s_rate2 = eq_dur/teq2/total_param
        nparam = eq_dur/(teq2*coeff1)
        print("REALIZATION: Analysis Duration = %8.4fs | Parameter Speed Rate = %8.4fs | EQ Speed Rate = %8.4fs | Total Time Elapsed = %8.4f"%(eq_dur, s_rate, s_rate2, (end_t - start_t_eq)))
        print('STATUS: Earthquake #%d completed successfully!'%(ref) + " | End Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
    end_t_eq = time.time()
    print('STATUS OVERALL: Complete Duration = %8.4f'%(end_t_eq - start_t_eq))

    try:
        if upload_preference == 'y':
            data = pd.read_csv("accesskeys.csv")
            access_id = data['Access key ID'][0]
            access_secret = data['Secret access key'][0]
            upload(folder, access_id, access_secret)
            shutil.rmtree(os.path.join('results',folder))
    except:
        print('Result stored locally')
        chk = input('Do you wish to upload? (y/n) :')
        if chk == 'y':
            access_id = input('Enter Access ID: ')
            access_secret = input('Enter Access Secret Key: ')
            upload(folder, access_id, access_secret)
    
    

    return None

@cli1.command()
@click.option('--earthquakes','-eq', type=(str), default= "Excitations.csv", help= 'Earthquakes')
@click.option('--unit', '-u', type=(str), default='m/s2')
@click.option('--screen/--no-screen', default=False)
def seeq(earthquakes, unit, screen):
    if earthquakes:
        total_eq = get_total_excitations(earthquakes)
        earthquake_generator = get_earthquake_list(earthquakes)
    else:
        earthquake_generator = get_earthquake_list('db elc')
    
    if unit == "m/s2":
        scale = 1.0
    elif unit == "g":
        scale = 1.0/9.81
    elif unit == "cm/s2":
        scale = 1.0/981

    for i in range(total_eq):
        ref, xg, yg, zg, dt, ndiv, ndt = next(earthquake_generator)
        time = np.arange(0, (ndt-1)*dt*ndiv+dt, dt)
        # print(ndt, dt, ndt*dt, (ndt-1)*dt*ndiv, ndiv)
        # print(xg.size, yg.size, zg.size, time.size)

        plt.figure(i)

        # Plot X
        plt.subplot(311)
        plt.plot(time, xg*scale)
        plt.xlabel('Time (sec)')
        plt.ylabel('Acceleration (' + unit +')')
        plt.title('Acceleration in X-Dir')
        plt.grid(True)

        # Plot Y
        plt.subplot(312)
        plt.plot(time, yg*scale)
        plt.xlabel('Time (sec)')
        plt.ylabel('Acceleration (' + unit +')')
        plt.title('Acceleration in Y-Dir')
        plt.grid(True)

        # Plot Z
        plt.subplot(313)
        plt.plot(time, zg*scale)
        plt.xlabel('Time (sec)')
        plt.ylabel('Acceleration (' + unit +')')
        plt.title('Acceleration in Z-Dir')
        plt.grid(True)
        
        plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25, wspace=0.35)
        plt.show()



@cli1.command()
@click.option('--earthquakes','-eq', type=(str), default= "Excitations.csv", help= 'Earthquakes')
@click.option('--results_type', '-r', type=str, default="paa1, pv1, pd1", help="Choice to select output results")
@click.option('--folder', '-f', type=str, default="Result", help="Folder name to store result")
@click.option('--damping', '-d', type=float, default=0.05, help="Critical damping ratio")
@click.option('--screen/--no-screen', default=False)
def spectral(earthquakes, results_type, folder, damping, screen):
    np.seterr(all='warn')
    """

    SAMPLE COMMAND: py -3.6-32 mhps.py spectral

    """
    
    var_param = os.path.join("data", "defaults", "SpectralVariableParameters.csv")
    const_param = os.path.join("data", "defaults", "SpectralConstantParameters.txt")
    knor = 1
    lxy = 1
    
    upload_preference = input('Do you wish to upload (default: no) ? [y/n] : ')

    # RESULT FOLDER SETUP
    folder = createfolder(folder)
    os.makedirs(os.path.join('results', folder, 'Time History Response'))
    # os.makedirs('results\\' + folder + '\\Time History Response')
    simulationdesc = input('Enter simulation description [None] : ')
    if simulationdesc:
        Path(os.path.join('results', folder, "SimulationInfo.csv")).touch()
        # Path('results\\' + folder + "\\SimulationInfo.csv").touch()
        with open(os.path.join('results', folder, "SimulationInfo.csv"), 'w') as f:
            for line in simulationdesc:
                f.write(line)
        f.close()

    # EARTHQUAKE GENERATOR SETUP
    if earthquakes:
        total_eq = get_total_excitations(earthquakes)
        earthquake_generator = get_earthquake_list(earthquakes)
    else:
        earthquake_generator = get_earthquake_list('db elc')

    total_param = get_total_variable_parameter_set(var_param)
    print("Total Parameters: " + str(total_param))
    
    # CONSTANT PARAMETER SETUP
    maxnst, am, ak = read_const_param(const_param)

    s_rate2 = 0.33
    s_rate = 1.5
    tot_time_eq2 = 0.0
    nparam = 1.0
    start_t_eq = time.time()

    peakvalues = None
    for i in range(total_eq):
        start_t = time.time()
        ref, xg, yg, zg, dt, ndiv, ndt = next(earthquake_generator)

        teq2 = (ndt-1)*dt*ndiv
        tot_time_eq2 = (tot_time_eq2 + teq2)/ref
        coeff1 = (0.0379*pow(ndiv, 1.5073))*pow(teq2, -0.6251*pow(ndiv, 0.125))
        time_current_eq = nparam*teq2*coeff1    #max(s_rate*total_param, s_rate2*teq2*total_param)
        time_remaining = (total_eq-ref+1)*(tot_time_eq2*coeff1*nparam)
        print('FACT: Duration of EQ = %8.4fs'%(teq2) + " | Start Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print('EXPECTATION: Analysis Duration = %8.2fs (%8.2f m) | Total Duration Remaining = %8.2fs (%6.2f m, %5.2fh)'%(time_current_eq, time_current_eq/60, time_remaining, time_remaining/60, time_remaining/60/60) + ' | Estimated Completion Time: ' + (datetime.datetime.now()+datetime.timedelta(seconds = time_current_eq)).strftime("%Y-%m-%d %H:%M:%S"))
        
        # SUPERSTRUCTURE VARIABLE PARAMETER SETUP
        superstructure_param_generator = read_ss_var_param(var_param)
        for j in range(total_param):
            ijk, nst, tx1, zeta, rtytx = next(superstructure_param_generator)
            zeta = damping
            if ((i == 0) and (j==0)) or ((p_nst, p_tx1, p_rtytx, p_zeta != nst, tx1, rtytx, zeta)):
                smx, skx, cdx, smy, sky, cdy = superstructure_propxy(nst, tx1, rtytx, am, ak, zeta, knor)
                p_nst, p_tx1, p_rtytx, p_zeta = nst, tx1, rtytx, zeta
            result, model = fixed_simulator(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, nst, smx[0:nst,0:nst], skx[0:nst,0:nst], cdx[0:nst,0:nst], smy[0:nst,0:nst], sky[0:nst,0:nst], cdy[0:nst,0:nst], screen)
            # analysis_folder = 'results\\' + folder + '\\Time History Response\\' + 'ANA-EQ-' + str(result.eq_ref) + '-PARAM-' + str(result.ijk)
            # analysis_folder = os.path.join('results', folder, 'Time History Response', 'ANA-EQ-' + str(result.eq_ref) + '-PARAM-' + str(result.ijk))
            # os.makedirs(analysis_folder)
            peakvaluesparam, peakvaluesparamhead = result_viewer(result, model, results_type, folder)
         
            if peakvaluesparam is not None:
                if j == 0:
                    peakvalues = peakvaluesparam
                else:
                    peakvalues = np.vstack((peakvalues, peakvaluesparam))
        print('EQ' + str(ref) + 'Successful!')


        
        if peakvalues is not None:
            if i == 0:
                # Path('results\\' + folder + "\\Peak.csv").touch()
                dtype1 = np.dtype([('IJK', 'i4'), ('NST', 'i4'), ('TX1', 'f4'), ('ZETA','f4'), ('RTYTX','f4')])
                ijk, nst, tx1, zeta, rtytx = np.loadtxt(var_param, delimiter=',', usecols=range(5), skiprows=1, unpack = True, dtype=dtype1)
                vector_length = len(ijk)
                ijk = ijk.reshape(vector_length,1)
                nst = nst.reshape(vector_length,1)
                tx1 = tx1.reshape(vector_length,1)
                zeta = zeta.reshape(vector_length,1)
                rtytx = rtytx.reshape(vector_length,1)

                peakmat = np.hstack((ijk, nst, tx1, zeta, rtytx, peakvalues))
                peakmathead = ['#', 'N', 'T', 'ZETA', 'RTYTX'] +  [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            else:
                peakmat = np.hstack((peakmat, peakvalues))
                peakmathead += [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            if i == (total_eq - 1):
                if total_param == 1:
                    peakmat = pd.DataFrame(peakmat.reshape(-1, len(peakmat)))
                else:
                    peakmat = pd.DataFrame(peakmat)
                # print(peakmat)
                # print(peakvaluesparamhead)
                peakmat.columns = peakmathead
                # peakmat.to_csv('results\\' + folder + "\\Peak.csv", mode='w', sep=',', index=False)
                peakmat.to_csv(os.path.join("results", folder, "Peak.csv"), mode='w', sep=',', index=False)
        end_t = time.time()
        eq_dur = end_t - start_t
        s_rate = eq_dur/total_param
        s_rate2 = eq_dur/teq2/total_param
        nparam = eq_dur/(teq2*coeff1)
        print("REALIZATION: Analysis Duration = %8.4fs | Parameter Speed Rate = %8.4fs | EQ Speed Rate = %8.4fs | Total Time Elapsed = %8.4f"%(eq_dur, s_rate, s_rate2, (end_t - start_t_eq)))
        print('STATUS: Earthquake #%d completed successfully!'%(ref) + " | End Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
    end_t_eq = time.time()
    print('STATUS OVERALL: Complete Duration = %8.4f'%(end_t_eq - start_t_eq))
        
    try:
        if upload_preference == 'y':
            data = pd.read_csv("accesskeys.csv")
            access_id = data['Access key ID'][0]
            access_secret = data['Secret access key'][0]
            upload(folder, access_id, access_secret)
            shutil.rmtree(os.path.join('results',folder))
    except:
        print('Result stored locally')
        chk = input('Do you wish to upload? (y/n) :')
        if chk == 'y':
            access_id = input('Enter Access ID: ')
            access_secret = input('Enter Access Secret Key: ')
            upload(folder, access_id, access_secret)

    
    return None

@cli1.command()
def fft():
    pass

@cli1.command()
@click.option('--const_param', '-cp',type=str, default= 'ConstantParameters.txt', help= 'File name containing constant parameters for the structure')
@click.option('--var_param', '-vp',type=str, default= 'SS_VP.csv', help= 'File name containing variable parameters for the structure')
@click.option('--iso_param', '-vp',type=str, default= 'L_VP.csv', help= 'File name containing variable parameters for linear base isolator')
@click.option('--earthquakes','-eq', type=(str), default= "Excitations.csv", help= 'Earthquakes')
@click.option('--knor', '-knor', type=int, default=1, help="Normalizing the superstructure for given time period. 1 for normalized and 0 for un-normalized.")
@click.option('--results_type', '-r', type=str, default="g*, aa1, db, fb, paa1, pdb", help="Choice to select output results")
@click.option('--folder', '-f', type=str, default="Result", help="Folder name to store result")
@click.option('--outputunits', type=list, default=['m/s2', 'cm/s', 'cm', 'kn', 'j'])
@click.option('--lxy', '-lxy', type=int, default=0)
@click.option('--screen/--no-screen', default=False)
def biso_linear(const_param, var_param, iso_param, earthquakes, knor, results_type, folder, outputunits, lxy, screen):
    print("I am in linear isolator module")

    upload_preference = input('Do you wish to upload (default: no) ? [y/n] : ')

    # RESULT FOLDER SETUP
    folder = createfolder(folder)
    os.makedirs(os.path.join('results', folder, 'Time History Response'))
    # os.makedirs('results\\' + folder + '\\Time History Response')
    simulationdesc = input('Enter simulation description [None] : ')
    if simulationdesc:
        Path(os.path.join('results', folder, "SimulationInfo.csv")).touch()
        # Path('results\\' + folder + "\\SimulationInfo.csv").touch()
        with open(os.path.join('results', folder, "SimulationInfo.csv"), 'w') as f:
            for line in simulationdesc:
                f.write(line)
        f.close()

    # EARTHQUAKE GENERATOR SETUP
    if earthquakes:
        total_eq = get_total_excitations(earthquakes)
        earthquake_generator = get_earthquake_list(earthquakes)
    else:
        earthquake_generator = get_earthquake_list('db elc')

    total_param = get_total_variable_parameter_set(var_param)
    print("Total Parameters: " + str(total_param))
    
    # CONSTANT PARAMETER SETUP
    maxnst, am, ak = read_const_param(const_param)

    s_rate2 = 0.33
    s_rate = 1.5
    tot_time_eq2 = 0.0
    nparam = 1.0
    start_t_eq = time.time()

    peakvalues = None
    for i in range(total_eq):
        start_t = time.time()
        ref, xg, yg, zg, dt, ndiv, ndt = next(earthquake_generator)

        teq2 = (ndt-1)*dt*ndiv
        tot_time_eq2 = (tot_time_eq2 + teq2)/ref
        coeff1 = (0.0379*pow(ndiv, 1.5073))*pow(teq2, -0.6251*pow(ndiv, 0.125))
        time_current_eq = nparam*teq2*coeff1    #max(s_rate*total_param, s_rate2*teq2*total_param)
        time_remaining = (total_eq-ref+1)*(tot_time_eq2*coeff1*nparam)
        print('FACT: Duration of EQ = %8.4fs'%(teq2) + " | Start Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print('EXPECTATION: Analysis Duration = %8.2fs (%8.2f m) | Total Duration Remaining = %8.2fs (%6.2f m, %5.2fh)'%(time_current_eq, time_current_eq/60, time_remaining, time_remaining/60, time_remaining/60/60) + ' | Estimated Completion Time: ' + (datetime.datetime.now()+datetime.timedelta(seconds = time_current_eq)).strftime("%Y-%m-%d %H:%M:%S"))
        
        # SUPERSTRUCTURE VARIABLE PARAMETER SETUP
        superstructure_param_generator = read_ss_var_param(var_param)
        
        ##? # ISOLATOR VARIABLE PARAMETER SETUP
        ##? 
        isolator_param_generator = read_iso_l_var_param(iso_param)

        for j in range(total_param):
            ijk, nst, tx1, zeta, rtytx = next(superstructure_param_generator)
            ##? 
            ijk, rmbm, tbx, zetabx, rtytxb, rzxzyb = next(isolator_param_generator)

            ##? 
            # if ((i == 0) and (j==0)) or ((p_nst, p_tx1, p_rtytx, p_zeta != nst, tx1, rtytx, zeta)):
            if ((i == 0) and (j==0)) or ((p_nst, p_tx1, p_rtytx, p_zeta, p_rmbm, p_tbx, p_zetabx, p_rtytxb, p_rzxzyb != nst, tx1, rtytx, zeta, rmbm, tbx, zetabx, rtytxb, rzxzyb)):
                smx, skx, cdx, smy, sky, cdy = superstructure_propxy(nst, tx1, rtytx, am, ak, zeta, knor)
                ##? 
                smx, skx, cdx, smy, sky, cdy = addlinear_iso(smx, skx, cdx, smy, sky, cdy, nst, rmbm, tbx, zetabx, rtytxb, rzxzyb)
                # p_nst, p_tx1, p_rtytx, p_zeta = nst, tx1, rtytx, zeta
                ##? 
                p_nst, p_tx1, p_rtytx, p_zeta, p_rmbm, p_tbx, p_zetabx, p_rtytxb, p_rzxzyb  = nst, tx1, rtytx, zeta, rmbm, tbx, zetabx, rtytxb, rzxzyb
            
            result, model = simulator_linear(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, nst+1, smx, skx, cdx, smy, sky, cdy, screen)    
            #result, model = fixed_simulator(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, nst, smx[0:nst,0:nst], skx[0:nst,0:nst], cdx[0:nst,0:nst], smy[0:nst,0:nst], sky[0:nst,0:nst], cdy[0:nst,0:nst])
            # analysis_folder = 'results\\' + folder + '\\Time History Response\\' + 'ANA-EQ-' + str(result.eq_ref) + '-PARAM-' + str(result.ijk)
            # os.makedirs(analysis_folder)
            peakvaluesparam, peakvaluesparamhead = result_viewer(result, model, results_type, folder)
         
            if peakvaluesparam is not None:
                if j == 0:
                    peakvalues = peakvaluesparam
                else:
                    peakvalues = np.vstack((peakvalues, peakvaluesparam))


        
        if peakvalues is not None:
            if i == 0:
                # Path('results\\' + folder + "\\Peak.csv").touch()
                peakmat = peakvalues
                peakmathead = [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            else:
                peakmat = np.hstack((peakmat, peakvalues))
                peakmathead += [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            if i == (total_eq - 1):
                if total_param == 1:
                    peakmat = pd.DataFrame(peakmat.reshape(-1, len(peakmat)))
                else:
                    peakmat = pd.DataFrame(peakmat)
                # print(peakmat)
                # print(peakvaluesparamhead)
                peakmat.columns = peakmathead
                peakmat.to_csv(os.path.join("results", folder, "Peak.csv"), mode='w', sep=',', index=False)
        end_t = time.time()
        eq_dur = end_t - start_t
        s_rate = eq_dur/total_param
        s_rate2 = eq_dur/teq2/total_param
        nparam = eq_dur/(teq2*coeff1)
        print("REALIZATION: Analysis Duration = %8.4fs | Parameter Speed Rate = %8.4fs | EQ Speed Rate = %8.4fs | Total Time Elapsed = %8.4f"%(eq_dur, s_rate, s_rate2, (end_t - start_t_eq)))
        print('STATUS: Earthquake #%d completed successfully!'%(ref) + " | End Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
    end_t_eq = time.time()
    print('STATUS OVERALL: Complete Duration = %8.4f'%(end_t_eq - start_t_eq))

    try:
        if upload_preference == 'y':
            data = pd.read_csv("accesskeys.csv")
            access_id = data['Access key ID'][0]
            access_secret = data['Secret access key'][0]
            upload(folder, access_id, access_secret)
            shutil.rmtree(os.path.join('results',folder))
    except:
        print('Result stored locally')
        chk = input('Do you wish to upload? (y/n) :')
        if chk == 'y':
            access_id = input('Enter Access ID: ')
            access_secret = input('Enter Access Secret Key: ')
            upload(folder, access_id, access_secret)
    
    return None


@cli1.command()
@click.option('--const_param', '-cp',type=str, default= 'ConstantParameters.txt', help= 'File name containing constant parameters for the structure')
@click.option('--var_param', '-vp',type=str, default= 'SS_VP.csv', help= 'File name containing variable parameters for the structure')
@click.option('--iso_param', '-vp',type=str, default= 'PF_VP.csv', help= 'File name containing variable parameters for linear base isolator')
@click.option('--earthquakes','-eq', type=(str), default= "Excitations.csv", help= 'Earthquakes')
@click.option('--knor', '-knor', type=int, default=1, help="Normalizing the superstructure for given time period. 1 for normalized and 0 for un-normalized.")
@click.option('--results_type', '-r', type=str, default="g*, aa1, db, f*, paa1, pdb", help="Choice to select output results")
@click.option('--folder', '-f', type=str, default="Result", help="Folder name to store result")
@click.option('--outputunits', type=list, default=['m/s2', 'cm/s', 'cm', 'kn', 'j'])
@click.option('--lxy', '-lxy', type=int, default=0)
@click.option('--screen/--no-screen', default=False)
def biso_pf(const_param, var_param, iso_param, earthquakes, knor, results_type, folder, outputunits, lxy, screen):

    upload_preference = input('Do you wish to upload (default: no) ? [y/n] : ')

    # RESULT FOLDER SETUP
    folder = createfolder(folder)
    os.makedirs(os.path.join('results', folder, 'Time History Response'))
    # os.makedirs('results\\' + folder + '\\Time History Response')
    simulationdesc = input('Enter simulation description [None] : ')
    if simulationdesc:
        Path(os.path.join('results', folder, "SimulationInfo.csv")).touch()
        # Path('results\\' + folder + "\\SimulationInfo.csv").touch()
        with open(os.path.join('results', folder, "SimulationInfo.csv"), 'w') as f:
            for line in simulationdesc:
                f.write(line)
        f.close()

    # EARTHQUAKE GENERATOR SETUP
    if earthquakes:
        total_eq = get_total_excitations(earthquakes)
        earthquake_generator = get_earthquake_list(earthquakes)
    else:
        earthquake_generator = get_earthquake_list('db elc')

    total_param = get_total_variable_parameter_set(var_param)
    print("Total Parameters: " + str(total_param))
    
    # CONSTANT PARAMETER SETUP
    maxnst, am, ak = read_const_param(const_param)

    s_rate2 = 0.33
    s_rate = 1.5
    tot_time_eq2 = 0.0
    nparam = 1.0
    start_t_eq = time.time()

    peakvalues = None
    for i in range(total_eq):
        start_t = time.time()
        ref, xg, yg, zg, dt, ndiv, ndt = next(earthquake_generator)

        teq2 = (ndt-1)*dt*ndiv
        tot_time_eq2 = (tot_time_eq2 + teq2)/ref
        coeff1 = (0.0379*pow(ndiv, 1.5073))*pow(teq2, -0.6251*pow(ndiv, 0.125))
        time_current_eq = nparam*teq2*coeff1    #max(s_rate*total_param, s_rate2*teq2*total_param)
        time_remaining = (total_eq-ref+1)*(tot_time_eq2*coeff1*nparam)
        print('FACT: Duration of EQ = %8.4fs'%(teq2) + " | Start Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print('EXPECTATION: Analysis Duration = %8.2fs (%8.2f m) | Total Duration Remaining = %8.2fs (%6.2f m, %5.2fh)'%(time_current_eq, time_current_eq/60, time_remaining, time_remaining/60, time_remaining/60/60) + ' | Estimated Completion Time: ' + (datetime.datetime.now()+datetime.timedelta(seconds = time_current_eq)).strftime("%Y-%m-%d %H:%M:%S"))
        
        # SUPERSTRUCTURE VARIABLE PARAMETER SETUP
        superstructure_param_generator = read_ss_var_param(var_param)
        
        # ISOLATOR VARIABLE PARAMETER SETUP
        isolator_param_generator = read_iso_pf_var_param(iso_param)

        for j in range(total_param):
            ijk, nst, tx1, zeta, rtytx = next(superstructure_param_generator)
            iso = next(isolator_param_generator)

            if ((i == 0) and (j==0)) or ((p_nst, p_tx1, p_rtytx, p_zeta, p_iso != nst, tx1, rtytx, zeta, iso)):
                smx, skx, cdx, smy, sky, cdy = superstructure_propxy(nst, tx1, rtytx, am, ak, zeta, knor)
                smx, skx, cdx, smy, sky, cdy = addlinear_iso(smx, skx, cdx, smy, sky, cdy, nst, iso.rmbm, iso.tbx, iso.zetabx, iso.rtytxb, iso.rzxzyb)
                p_nst, p_tx1, p_rtytx, p_zeta  = nst, tx1, rtytx, zeta
                p_iso = iso
            
            result, model = simulator_pf(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, nst+1, smx, skx, cdx, smy, sky, cdy, iso, screen)    
            # analysis_folder = 'results\\' + folder + '\\Time History Response\\' + 'ANA-EQ-' + str(result.eq_ref) + '-PARAM-' + str(result.ijk)
            # os.makedirs(analysis_folder)
            peakvaluesparam, peakvaluesparamhead = result_viewer(result, model, results_type, folder)
         
            if peakvaluesparam is not None:
                if j == 0:
                    peakvalues = peakvaluesparam
                else:
                    peakvalues = np.vstack((peakvalues, peakvaluesparam))

        if peakvalues is not None:
            if i == 0:
                if total_param == 1:
                    peakmat = pd.DataFrame(peakmat.reshape(-1, len(peakmat)))
                else:
                    peakmat = pd.DataFrame(peakmat)
                peakmathead = [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            else:
                peakmat = np.hstack((peakmat, peakvalues))
                peakmathead += [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            if i == (total_eq - 1):
                peakmat = pd.DataFrame(peakmat)
                peakmat.columns = peakmathead
                peakmat.to_csv(os.path.join("results", folder, "Peak.csv"), mode='w', sep=',', index=False)
        end_t = time.time()
        eq_dur = end_t - start_t
        s_rate = eq_dur/total_param
        s_rate2 = eq_dur/teq2/total_param
        nparam = eq_dur/(teq2*coeff1)
        print("REALIZATION: Analysis Duration = %8.4fs | Parameter Speed Rate = %8.4fs | EQ Speed Rate = %8.4fs | Total Time Elapsed = %8.4f"%(eq_dur, s_rate, s_rate2, (end_t - start_t_eq)))
        print('STATUS: Earthquake #%d completed successfully!'%(ref) + " | End Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
    end_t_eq = time.time()
    print('STATUS OVERALL: Complete Duration = %8.4f'%(end_t_eq - start_t_eq))

    try:
        if upload_preference == 'y':
            data = pd.read_csv("accesskeys.csv")
            access_id = data['Access key ID'][0]
            access_secret = data['Secret access key'][0]
            upload(folder, access_id, access_secret)
            shutil.rmtree(os.path.join('results',folder))
    except:
        print('Result stored locally')
        chk = input('Do you wish to upload? (y/n) :')
        if chk == 'y':
            access_id = input('Enter Access ID: ')
            access_secret = input('Enter Access Secret Key: ')
            upload(folder, access_id, access_secret)
    
    return None


@cli1.command()
@click.option('--const_param', '-cp',type=str, default= 'ConstantParameters.txt', help= 'File name containing constant parameters for the structure')
@click.option('--var_param', '-vp',type=str, default= 'SS_VP.csv', help= 'File name containing variable parameters for the structure')
@click.option('--iso_param', '-vp',type=str, default= 'OSBI_VP.csv', help= 'File name containing variable parameters for linear base isolator')
@click.option('--earthquakes','-eq', type=(str), default= "Excitations.csv", help= 'Earthquakes')
@click.option('--knor', '-knor', type=int, default=1, help="Normalizing the superstructure for given time period. 1 for normalized and 0 for un-normalized.")
@click.option('--results_type', '-r', type=str, default="g*, aa1, db, f*, ro*, sm*, paa1, pdb", help="Choice to select output results")
@click.option('--folder', '-f', type=str, default="result", help="Folder name to store result")
@click.option('--outputunits', type=list, default=['m/s2', 'cm/s', 'cm', 'kn', 'j'])
@click.option('--lxy', '-lxy', type=int, default=0)
@click.option('--screen/--no-screen', default=False)
def biso_osbi(const_param, var_param, iso_param, earthquakes, knor, results_type, folder, outputunits, lxy, screen):

    upload_preference = input('Do you wish to upload (default: no) ? [y/n] : ')

    # RESULT FOLDER SETUP
    folder = createfolder(folder)
    os.makedirs(os.path.join('results', folder, 'Time History Response'))
    simulationdesc = input('Enter simulation description [None] : ')
    if simulationdesc:
        Path(os.path.join('results', folder, "SimulationInfo.csv")).touch()
        with open(os.path.join('results', folder, "SimulationInfo.csv"), 'w') as f:
            for line in simulationdesc:
                f.write(line)
        f.close()

    # EARTHQUAKE GENERATOR SETUP
    if earthquakes:
        total_eq = get_total_excitations(earthquakes)
        earthquake_generator = get_earthquake_list(earthquakes)
    else:
        earthquake_generator = get_earthquake_list('db elc')

    total_param = get_total_variable_parameter_set(var_param)
    print("\nTotal Parameters: " + str(total_param))
    print("Total EQs: " + str(total_eq) +'\n\n')
    
    # CONSTANT PARAMETER SETUP
    maxnst, am, ak = read_const_param(const_param)

    s_rate2 = 0.33
    s_rate = 1.5
    tot_time_eq2 = 0.0
    nparam = 1.0
    start_t_eq = time.time()

    peakvalues = None
    for i in range(total_eq):
        start_t = time.time()
        ref, xg, yg, zg, dt, ndiv, ndt = next(earthquake_generator)

        teq2 = (ndt-1)*dt*ndiv
        tot_time_eq2 = (tot_time_eq2 + teq2)/ref
        coeff1 = (0.0379*pow(ndiv, 1.5073))*pow(teq2, -0.6251*pow(ndiv, 0.125))
        time_current_eq = nparam*teq2*coeff1    #max(s_rate*total_param, s_rate2*teq2*total_param)
        time_remaining = (total_eq-ref+1)*(tot_time_eq2*coeff1*nparam)
        print('FACT: Duration of EQ = %8.4fs'%(teq2) + " | Start Time: " + prCyan(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        print('EXPECTATION: Analysis Duration = %8.2fs (%8.2f m) | Total Duration Remaining = %8.2fs (%6.2f m, %5.2fh)'%(time_current_eq, time_current_eq/60, time_remaining, time_remaining/60, time_remaining/60/60) + ' | Estimated Completion Time: ' + (datetime.datetime.now()+datetime.timedelta(seconds = time_current_eq)).strftime("%Y-%m-%d %H:%M:%S"))
        
        # SUPERSTRUCTURE VARIABLE PARAMETER SETUP
        superstructure_param_generator = read_ss_var_param(var_param)
        
        # ISOLATOR VARIABLE PARAMETER SETUP
        isolator_param_generator = read_iso_osbi_var_param(iso_param, am[0], 4)

        for j in range(total_param):
            ijk, nst, tx1, zeta, rtytx = next(superstructure_param_generator)
            iso = next(isolator_param_generator)
            if screen == True:
                print(iso)

            if ((i == 0) and (j==0)) or ((p_nst, p_tx1, p_rtytx, p_zeta, p_iso != nst, tx1, rtytx, zeta, iso)):
                smx, skx, cdx, smy, sky, cdy = superstructure_propxy(nst, tx1, rtytx, am, ak, zeta, knor)
                smx, skx, cdx, smy, sky, cdy = addlinear_iso(smx, skx, cdx, smy, sky, cdy, nst, iso.rmbm, iso.tbx, iso.zetabx, iso.rtytxb, iso.rzyzxb)
                p_nst, p_tx1, p_rtytx, p_zeta  = nst, tx1, rtytx, zeta
                p_iso = iso
            
            result, model = simulator_osbi(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, nst+1, smx, skx, cdx, smy, sky, cdy, iso, screen)    
            # analysis_folder = os.path.join('results', folder, 'Time History Response', 'ANA-EQ-' + str(result.eq_ref) + '-PARAM-' + str(result.ijk))
            # os.makedirs(analysis_folder)
            peakvaluesparam, peakvaluesparamhead = result_viewer(result, model, results_type, folder)
         
            if peakvaluesparam is not None:
                if j == 0:
                    peakvalues = peakvaluesparam
                else:
                    peakvalues = np.vstack((peakvalues, peakvaluesparam))


        if peakvalues is not None:
            if i == 0:

                #

                # This place is dedicated to add model info in peak results as in spectral

                # 
                peakmat = peakvalues
                peakmathead = [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            else:
                peakmat = np.hstack((peakmat, peakvalues))
                peakmathead += [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            if i == (total_eq - 1):
                if total_param == 1:
                    peakmat = pd.DataFrame(peakmat.reshape(-1, len(peakmat)))
                else:
                    peakmat = pd.DataFrame(peakmat)
                peakmat.columns = peakmathead
                peakmat.to_csv(os.path.join("results", folder, "Peak.csv"), mode='w', sep=',', index=False)
        end_t = time.time()
        eq_dur = end_t - start_t
        s_rate = eq_dur/total_param
        s_rate2 = eq_dur/teq2/total_param
        nparam = eq_dur/(teq2*coeff1)
        print("REALIZATION: Analysis Duration = %8.4fs | Parameter Speed Rate = %8.4fs | EQ Speed Rate = %8.4fs | Total Time Elapsed = %8.4f"%(eq_dur, s_rate, s_rate2, (end_t - start_t_eq)))
        print('STATUS: Earthquake #%d completed successfully!'%(ref) + " | End Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
    end_t_eq = time.time()
    print('STATUS OVERALL: Complete Duration = %8.4f'%(end_t_eq - start_t_eq))

    try:
        if upload_preference == 'y':
            data = pd.read_csv("accesskeys.csv")
            access_id = data['Access key ID'][0]
            access_secret = data['Secret access key'][0]
            upload(folder, access_id, access_secret)
            shutil.rmtree(os.path.join('results',folder))
    except:
        print('Result stored locally')
        chk = input('Do you wish to upload? (y/n) :')
        if chk == 'y':
            access_id = input('Enter Access ID: ')
            access_secret = input('Enter Access Secret Key: ')
            upload(folder, access_id, access_secret)
    
    return None

@cli1.command()
@click.option('--const_param', '-cp',type=str, default= 'ConstantParameters.txt', help= 'File name containing constant parameters for the structure')
@click.option('--var_param', '-vp',type=str, default= 'SS_VP.csv', help= 'File name containing variable parameters for the structure')
@click.option('--iso_param', '-vp',type=str, default= 'OSBI_VP.csv', help= 'File name containing variable parameters for linear base isolator')
@click.option('--earthquakes','-eq', type=(str), default= "Excitations.csv", help= 'Earthquakes')
@click.option('--knor', '-knor', type=int, default=1, help="Normalizing the superstructure for given time period. 1 for normalized and 0 for un-normalized.")
@click.option('--results_type', '-r', type=str, default="g*, aa1, db, f*, paa1, pdb", help="Choice to select output results")
@click.option('--folder', '-f', type=str, default="Result", help="Folder name to store result")
@click.option('--outputunits', type=list, default=['m/s2', 'cm/s', 'cm', 'kn', 'j'])
@click.option('--lxy', '-lxy', type=int, default=0)
@click.option('--screen/--no-screen', default=False)
def biso_boucwen(const_param, var_param, iso_param, earthquakes, knor, results_type, folder, outputunits, lxy, screen):

    upload_preference = input('Do you wish to upload (default: no) ? [y/n] : ')

    # RESULT FOLDER SETUP
    folder = createfolder(folder)
    os.makedirs(os.path.join('results', folder, 'Time History Response'))
    simulationdesc = input('Enter simulation description [None] : ')
    if simulationdesc:
        Path(os.path.join('results', folder, "SimulationInfo.csv")).touch()
        with open(os.path.join('results', folder, "SimulationInfo.csv"), 'w') as f:
            for line in simulationdesc:
                f.write(line)
        f.close()

    # EARTHQUAKE GENERATOR SETUP
    if earthquakes:
        total_eq = get_total_excitations(earthquakes)
        earthquake_generator = get_earthquake_list(earthquakes)
    else:
        earthquake_generator = get_earthquake_list('db elc')

    total_param = get_total_variable_parameter_set(var_param)
    print("Total Parameters: " + str(total_param))
    
    # CONSTANT PARAMETER SETUP
    maxnst, am, ak = read_const_param(const_param)

    s_rate2 = 0.33
    s_rate = 1.5
    tot_time_eq2 = 0.0
    nparam = 1.0
    start_t_eq = time.time()

    peakvalues = None
    for i in range(total_eq):
        start_t = time.time()
        ref, xg, yg, zg, dt, ndiv, ndt = next(earthquake_generator)

        teq2 = (ndt-1)*dt*ndiv
        tot_time_eq2 = (tot_time_eq2 + teq2)/ref
        coeff1 = (0.0379*pow(ndiv, 1.5073))*pow(teq2, -0.6251*pow(ndiv, 0.125))
        time_current_eq = nparam*teq2*coeff1    #max(s_rate*total_param, s_rate2*teq2*total_param)
        time_remaining = (total_eq-ref+1)*(tot_time_eq2*coeff1*nparam)
        print('FACT: Duration of EQ = %8.4fs'%(teq2) + " | Start Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        print('EXPECTATION: Analysis Duration = %8.2fs (%8.2f m) | Total Duration Remaining = %8.2fs (%6.2f m, %5.2fh)'%(time_current_eq, time_current_eq/60, time_remaining, time_remaining/60, time_remaining/60/60) + ' | Estimated Completion Time: ' + (datetime.datetime.now()+datetime.timedelta(seconds = time_current_eq)).strftime("%Y-%m-%d %H:%M:%S"))
        
        # SUPERSTRUCTURE VARIABLE PARAMETER SETUP
        superstructure_param_generator = read_ss_var_param(var_param)
        
        # ISOLATOR VARIABLE PARAMETER SETUP
        isolator_param_generator = read_iso_boucwen_var_param(iso_param)

        for j in range(total_param):
            ijk, nst, tx1, zeta, rtytx = next(superstructure_param_generator)
            iso = next(isolator_param_generator)

            if ((i == 0) and (j==0)) or ((p_nst, p_tx1, p_rtytx, p_zeta, p_iso != nst, tx1, rtytx, zeta, iso)):
                smx, skx, cdx, smy, sky, cdy = superstructure_propxy(nst, tx1, rtytx, am, ak, zeta, knor)
                smx, skx, cdx, smy, sky, cdy = addboucwen_iso(smx, skx, cdx, smy, sky, cdy, nst, iso.rmbm, iso.tbx, iso.zetabx, iso.rtytxb, iso.rzyzxb)
                p_nst, p_tx1, p_rtytx, p_zeta  = nst, tx1, rtytx, zeta
                p_iso = iso
            
            result, model = simulator_boucwen(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, nst+1, smx, skx, cdx, smy, sky, cdy, iso, screen)    
            # analysis_folder = 'results\\' + folder + '\\Time History Response\\' + 'ANA-EQ-' + str(result.eq_ref) + '-PARAM-' + str(result.ijk)
            # os.makedirs(analysis_folder)
            peakvaluesparam, peakvaluesparamhead = result_viewer(result, model, results_type, folder)
         
            if peakvaluesparam is not None:
                if j == 0:
                    peakvalues = peakvaluesparam
                else:
                    peakvalues = np.vstack((peakvalues, peakvaluesparam))

        if peakvalues is not None:
            if i == 0:
                peakmat = peakvalues
                peakmathead = [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            else:
                peakmat = np.hstack((peakmat, peakvalues))
                peakmathead += [s+'-EQ-'+str(result.eq_ref) for s in peakvaluesparamhead]
            if i == (total_eq - 1):
                if total_param == 1:
                    peakmat = pd.DataFrame(peakmat.reshape(-1, len(peakmat)))
                else:
                    peakmat = pd.DataFrame(peakmat)
                peakmat.columns = peakmathead
                peakmat.to_csv(os.path.join("results", folder, "Peak.csv"), mode='w', sep=',', index=False)
        end_t = time.time()
        eq_dur = end_t - start_t
        s_rate = eq_dur/total_param
        s_rate2 = eq_dur/teq2/total_param
        nparam = eq_dur/(teq2*coeff1)
        print("REALIZATION: Analysis Duration = %8.4fs | Parameter Speed Rate = %8.4fs | EQ Speed Rate = %8.4fs | Total Time Elapsed = %8.4f"%(eq_dur, s_rate, s_rate2, (end_t - start_t_eq)))
        print('STATUS: Earthquake #%d completed successfully!'%(ref) + " | End Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n\n")
    end_t_eq = time.time()
    print('STATUS OVERALL: Complete Duration = %8.4f'%(end_t_eq - start_t_eq))

    try:
        if upload_preference == 'y':
            data = pd.read_csv("accesskeys.csv")
            access_id = data['Access key ID'][0]
            access_secret = data['Secret access key'][0]
            upload(folder, access_id, access_secret)
            shutil.rmtree(os.path.join('results',folder))
    except:
        print('Result stored locally')
        chk = input('Do you wish to upload? (y/n) :')
        if chk == 'y':
            access_id = input('Enter Access ID: ')
            access_secret = input('Enter Access Secret Key: ')
            upload(folder, access_id, access_secret)
    
    return None


##########################################################################################
#                                                                                        #
#                                   TORSIONAL CODE                                       #
#                                                                                        #
##########################################################################################




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
#   2.  PF (Done)
#   3.  CRR
#   4.  ERR
#   5.  FPS (Done)
#   6.  VFPI
#   7.  LRB (Done)
#   8.  HDR
#   9.  NZS
#   10. EDF
