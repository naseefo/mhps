import os
import sys
import re
import numpy as np
from pathlib import Path
import pandas as pd
import csv
import cProfile
import io
import pstats

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


class ResultFixedXY:
    def __init__(self, eq_refi= 0.0, ijki= 0.0, timei= 0.0, gxi= 0.0, dxi= 0.0, vxi= 0.0, axi= 0.0, aaxi= 0.0, gyi= 0.0, dyi= 0.0, vyi= 0.0, ayi= 0.0, aayi= 0.0, fxi= 0.0, fyi= 0.0, eki= 0.0, edi= 0.0, esi= 0.0, eii= 0.0, errori= 0.0, smxi= 0.0, skxi= 0.0, cdxi= 0.0, smyi= 0.0, skyi= 0.0, cdyi= 0.0, roll = 0, theta_0 = 0.0, theta_r= 0.0, theta_r_dot2= 0.0, zbd= 0.0, zbd_dot2= 0.0, Fs1x= 0.0, Fs1y= 0.0, Fs2x= 0.0, Fs2y= 0.0, Fbx= 0.0, Fby= 0.0, F_axial= 0.0, Dc_axial= 0.0, Strain_axial=0.0, t_si =0.0, t_bi = 0.0, f_bi =0.0, t_sci =0.0, t_bci = 0.0, f_bci =0.0): # del
        self.eq_ref = eq_refi
        self.ijk = ijki
        self.time = timei
        self.gx = gxi
        self.dx = dxi
        self.vx = vxi
        self.ax = axi
        self.aax = aaxi
        self.gy = gyi
        self.dy = dyi
        self.vy = vyi
        self.ay = ayi
        self.aay = aayi
        self.fx = fxi
        self.fy = fyi
        self.ek = eki
        self.ed = edi
        self.es = esi
        self.ei = eii
        self.error = errori
        self.smx = smxi
        self.skx = skxi
        self.cdx = cdxi
        self.smy = smyi
        self.sky = skyi
        self.cdy = cdyi
        self.roll = roll
        self.theta_0 = theta_0
        self.theta_r = theta_r
        self.theta_r_dot2 = theta_r_dot2
        self.zbd = zbd
        self.zbd_dot2 = zbd_dot2
        self.Fs1x = Fs1x
        self.Fs1y = Fs1y
        self.Fs2x = Fs2x
        self.Fs2y = Fs2y
        self.Fbx = Fbx
        self.Fby = Fby
        self.F_axial = F_axial
        self.Dc_axial = Dc_axial
        self.Strain_axial = Strain_axial

        self.t_s = t_si
        self.t_sc = t_sci
        self.t_b = t_bi
        self.t_bc = t_bci
        self.f_b = f_bi
        self.f_bc = f_bci
      

class ModelInfo:
    def __init__(self, nsti):
        self.nst = nsti

def get_result(result, responsevariable, floorstart, floorend, peaktype, dirn):
    if responsevariable == 'g':
        vector1 = result.gx
        vector2 = result.gy
        vector = np.hstack((vector1, vector2))
        vectorhead = ['GX', 'GY']
    elif responsevariable == 'd':
        vector1 = result.dx[:, floorstart:floorend+1]
        vector2 = result.dy[:, floorstart:floorend+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("RDX-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("RDY-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    elif responsevariable == 'v':
        vector1 = result.vx[:, floorstart:floorend+1]
        vector2 = result.vy[:, floorstart:floorend+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("RVX-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("RVY-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    elif responsevariable == 'a':
        vector1 = result.ax[:, floorstart:floorend+1]
        vector2 = result.ay[:, floorstart:floorend+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("RAX-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("RAY-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    elif responsevariable == 'aa':
        vector1 = result.aax[:, floorstart:floorend+1]
        vector2 = result.aay[:, floorstart:floorend+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("AAX-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("AAY"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    elif responsevariable == 'f':
        vector1 = result.fx[:, floorstart:floorend+1]
        vector2 = result.fy[:, floorstart:floorend+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("FX-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("FY-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    elif responsevariable == 'en':
        vector1 = result.ek
        vector2 = result.ed
        vector3 = result.es
        vector4 = result.ei
        vector = np.hstack((vector1, vector2, vector3, vector4))
        vectorhead = ['KE', 'DE', 'SE', 'IE']
    elif responsevariable == 'er':
        vector = result.error
        vectorhead = ['Relative Error']
    elif responsevariable == 'ro':
        vector = result.roll
        vectorhead = ['ROLL']
    elif responsevariable == 'sm':
        vector1 = result.theta_0
        vector2 = result.theta_r
        vector3 = result.theta_r_dot2
        vector4 = result.zbd
        vector5 = result.zbd_dot2
        vector6 = result.Fs1x
        vector7 = result.Fs1y
        vector8 = result.Fs2x
        vector9 = result.Fs2y
        vector10 = result.Fbx
        vector11 = result.Fby
        vector12 = result.F_axial
        vector13 = result.Dc_axial
        vector14 = result.Strain_axial
        vector = np.hstack((vector1, vector2, vector3, vector4, vector5, vector6, vector7, vector8, vector9, vector10, vector11, vector12, vector13, vector14))
        vectorhead = ['theta_0', 'theta_r', 'theta_r_dot2', 'zbd', 'zbd_dot2', 'Fs1x', 'Fs1y', 'Fs2x', 'Fs2y', 'Fbx', 'Fby', 'F_Axial', 'Dc_Axial', 'Strain_Axial']
    
    elif responsevariable == 'tsd':
        si = 0
        floorstart = floorstart + si
        floorend = floorend + si
        vector = result.t_s[:, floorstart:floorend+1]
        vectorhead = [("SLAB-D-"+str(x+1)) for x in range(floorstart,floorend+1)]
    elif responsevariable == 'tbd':
        si = 0
        floorstart = floorstart + si
        floorend = floorend + si
        vector = result.t_b[:, floorstart:floorend+1]
        vectorhead = [("ISO-D-"+str(x+1)) for x in range(floorstart,floorend+1)]

    elif responsevariable == 'tsv':
        si = 3
        floorstart = floorstart + si
        floorend = floorend + si
        vector = result.t_s[:, floorstart:floorend+1]
        vectorhead = [("SLAB-V-"+str(x+1)) for x in range(floorstart,floorend+1)]
    elif responsevariable == 'tbv':
        si = 3
        floorstart = floorstart + si
        floorend = floorend + si
        vector = result.t_b[:, floorstart:floorend+1]
        vectorhead = [("ISO-V-"+str(x+1)) for x in range(floorstart,floorend+1)]

    elif responsevariable == 'tsa':
        si = 6
        floorstart = floorstart + si
        floorend = floorend + si
        vector = result.t_s[:, floorstart:floorend+1]
        vectorhead = [("SLAB-A-"+str(x+1)) for x in range(floorstart,floorend+1)]
    elif responsevariable == 'tba':
        si = 6
        floorstart = floorstart + si
        floorend = floorend + si
        vector = result.t_b[:, floorstart:floorend+1]
        vectorhead = [("ISO-A-"+str(x+1)) for x in range(floorstart,floorend+1)]

    elif responsevariable == 'tsaa':
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2
        si = 9
        floorstart = floorstart + si
        floorend = floorend + si
        vector = result.t_s[:, floorstart:floorend+1]
        vectorhead = [("SLAB-AA-"+str(x+1)) for x in range(floorstart,floorend+1)]
    elif responsevariable == 'tbaa':
        si = 9
        floorstart = floorstart + si
        floorend = floorend + si
        vector = result.t_b[:, floorstart:floorend+1]
        vectorhead = [("ISO-AA-"+str(x+1)) for x in range(floorstart,floorend+1)]
    
    elif responsevariable == 'tf':
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2
        si = 0
        floorstart = floorstart + si
        floorend = floorend + si
        vector = result.f_b[:, floorstart:floorend+1]
        vectorhead = [("BASESHEAR-"+str(x+1)) for x in range(floorstart,floorend+1)]
    
    elif responsevariable == 'tscd':
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2
        if (floorend - floorstart != 0):
            floorend = floorend+1
        si = 0
        br = si + 4
        floorstartX = floorstart + si
        floorendX = floorend + si
        vector1 = result.t_sc[:, floorstartX:floorendX+1]
        floorstartY = floorstart + br
        floorendY = floorend + br
        vector2 = result.t_sc[:, floorstartY:floorendY+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("SLAB-DX-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("SLAB-DY-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    
    elif responsevariable == 'tbcd':
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2

        if (floorend - floorstart != 0):
            floorend = floorend+1
        si = 0
        br = si + 4
        floorstartX = floorstart + si
        floorendX = floorend + si
        vector1 = result.t_bc[:, floorstartX:floorendX+1]
        floorstartY = floorstart + br
        floorendY = floorend + br
        vector2 = result.t_bc[:, floorstartY:floorendY+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("ISO-DX-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("ISO-DY-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
        print(vectorhead)

    elif responsevariable == 'tscv':
        print(floorstart, floorend)
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2

        if (floorend - floorstart != 0):
            floorend = floorend+1
        si = 8
        br = si + 4
        floorstartX = floorstart + si
        floorendX = floorend + si
        vector1 = result.t_sc[:, floorstartX:floorendX+1]
        floorstartY = floorstart + br
        floorendY = floorend + br
        vector2 = result.t_sc[:, floorstartY:floorendY+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("SLAB-VX-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("SLAB-VY-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    
    elif responsevariable == 'tbcv':
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2

        if (floorend - floorstart != 0):
            floorend = floorend+1
        si = 8
        br = si + 3
        floorstartX = floorstart + si
        floorendX = floorend + si
        vector1 = result.t_bc[:, floorstartX:floorendX+1]
        floorstartY = floorstart + br
        floorendY = floorend + br
        vector2 = result.t_bc[:, floorstartY:floorendY+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("ISO-VX-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("ISO-VY-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    
    elif responsevariable == 'tsca':
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2

        if (floorend - floorstart != 0):
            floorend = floorend+1
        si = 16
        br = si + 4
        floorstartX = floorstart + si
        floorendX = floorend + si
        vector1 = result.t_sc[:, floorstartX:floorendX+1]
        floorstartY = floorstart + br
        floorendY = floorend + br
        vector2 = result.t_sc[:, floorstartY:floorendY+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("SLAB-AX-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("SLAB-AY-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady

    elif responsevariable == 'tbca':
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2

        if (floorend - floorstart != 0):
            floorend = floorend+1
        si = 16
        br = si + 4
        floorstartX = floorstart + si
        floorendX = floorend + si
        vector1 = result.t_bc[:, floorstartX:floorendX+1]
        floorstartY = floorstart + br
        floorendY = floorend + br
        vector2 = result.t_bc[:, floorstartY:floorendY+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("ISO-AX-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("ISO-AY-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    
    elif responsevariable == 'tscaa':
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2
        if (floorend - floorstart != 0):
            floorend = floorend+1
        si = 24
        br = si + 4
        floorstartX = floorstart + si
        floorendX = floorend + si
        vector1 = result.t_sc[:, floorstartX:floorendX+1]
        floorstartY = floorstart + br
        floorendY = floorend + br
        vector2 = result.t_sc[:, floorstartY:floorendY+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("SLAB-AAX-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("SLAB-AAY-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady

    elif responsevariable == 'tbcaa':
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2

        if (floorend - floorstart != 0):
            floorend = floorend+1
        si = 24
        br = si + 4
        floorstartX = floorstart + si
        floorendX = floorend + si
        vector1 = result.t_bc[:, floorstartX:floorendX+1]
        floorstartY = floorstart + br
        floorendY = floorend + br
        vector2 = result.t_bc[:, floorstartY:floorendY+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("ISO-AAX-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("ISO-AAY-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady

    elif responsevariable == 'tfc':
        if (floorend - floorstart == 5):
            floorstart = 0
            floorend = 2
        if (floorend - floorstart != 0):
            floorend = floorend+1
        si = 0
        br = si + 4
        floorstartX = floorstart + si
        floorendX = floorend + si
        vector1 = result.f_bc[:, floorstartX:floorendX+1]
        floorstartY = floorstart + br
        floorendY = floorend + br
        vector2 = result.f_bc[:, floorstartY:floorendY+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("BASESHEAR-X-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("BASESHEAR-Y-C-"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    

    if peaktype == 1:
        responsevalues = np.absolute(vector).max(axis=0)
    elif peaktype == 0:
        responsevalues = vector
    return responsevalues, vectorhead


def result_viewer(result, model, results_type, folder):
    results_type = results_type.split(',')
    i, j = 0, 0
    peakvalues = None
    timeresponses = None
    peakhead = None
    for rpattern in results_type:
        responsevariable, floorstart, floorend, peaktype, dirn = pattern_reader(rpattern.strip(), model.nst)
        # print(responsevariable, floorstart, floorend, peaktype, dirn)
        responsevalues, vectorhead = get_result(result, responsevariable, floorstart, floorend, peaktype, dirn)
        if peaktype == 1:
            if i == 0:
                peakvalues = responsevalues
                peakhead = vectorhead
                i = 1
            else:
                peakvalues = np.hstack((peakvalues, responsevalues))
                peakhead += vectorhead
        else:
            if j == 0:
                timeresponses = np.hstack((result.time, responsevalues))
                timeresponsehead = ['Time'] + vectorhead
                j = 1
            else:
                timeresponses = np.hstack((timeresponses, responsevalues))
                timeresponsehead += vectorhead
    
    if timeresponses is not None:
        analysis_folder = os.path.join('results', folder,'Time History Response','ANA-EQ-' + str(result.eq_ref) + '-PARAM-' + str(result.ijk))
        if not os.path.isdir(analysis_folder):
            os.makedirs(analysis_folder)
        timeresponses = pd.DataFrame(timeresponses)
        timeresponses.columns = timeresponsehead
        superstructure = pd.DataFrame(result.skx)
        timeresponses.to_csv(os.path.join(analysis_folder, "TimeDomainResponses.csv"), mode='w', sep=',', index=False)
        superstructure.to_csv(os.path.join(analysis_folder, "SimulationInfo.csv"), mode='w', index=False)
    return peakvalues, peakhead

def pattern_reader(rpattern, nst):
    pattern_res = '^([p]|)([a-z]+|[f,e,r]+)([0-9\-\*\s]+|b$)'
    flag = 0
    for i in range(0,23):
        pattern_check = re.findall(pattern_res, rpattern)
        # print(pattern_check)
        if pattern_check[0]:
            if pattern_check[0][0] == 'p':
                peaktype = 1
            else:
                peaktype = 0
            
            responsevariable = pattern_check[0][1]

            chk1 = re.findall('([\*])', pattern_check[0][2])
            chk2 = re.findall('([0-9]+)\-([0-9]+)', pattern_check[0][2])
            chk3 = re.findall('([0-9]+)', pattern_check[0][2])
            chk4 = re.findall('([b])', pattern_check[0][2])

            if chk1:
                floorstart = 0
                floorend = nst-1
                flag = 1
                break
            elif chk2:
                floorstart = int(chk2[0][0]) - 1
                floorend = int(chk2[0][1]) - 1
                flag = 1
                break
            elif chk3:
                floorstart = int(chk3[0][0]) - 1
                floorend = int(chk3[0][0]) - 1
                flag = 1
                break
            elif chk4:
                floorstart = nst-1
                floorend = nst-1
                flag = 1
                break
            else:
                responsevariable += 'b'
                flag = 0
                break
            
    if flag == 0:
        print("ERROR: Discrepancy in response variable keys.")
        sys.exit()
    dirn = 0
    return responsevariable, floorstart, floorend, peaktype, dirn





def createfolder(folder):
    try:      
        if os.path.exists(os.path.join('results', folder)):
            for i in range(1,1000):
                if not os.path.exists(os.path.join('results', folder + '-' + str(i))):
                    folder = folder + '-' + str(i)
                    folder_input = input('Folder already exists! Enter new folder name ['+ folder + '] :')
                    if folder_input == "":
                        os.makedirs(os.path.join('results', folder))
                    else:
                        folder = folder_input
                        os.makedirs(os.path.join('results', folder))
                    break
        else:
            os.makedirs(os.path.join('results', folder))
    except OSError:
        print('Error creating directory')
    
    return folder





            

