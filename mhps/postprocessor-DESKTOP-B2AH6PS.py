import os
import sys
import re
import numpy as np
from pathlib import Path
import pandas as pd

class ResultFixedXY:
    def __init__(self, eq_refi, ijki, timei, gxi, dxi, vxi, axi, aaxi, gyi, dyi, vyi, ayi, aayi, fxi, fyi, eki, edi, esi, eii, errori, smx, skx, cdx, smy, sky, cdy):
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
        self.smx = smx 
        self.skx = skx
        self.cdx = cdx
        self.smy = smy
        self.sky = sky
        self.cdy = cdy

class ModelInfo:
    def __init__(self, nsti):
        self.nst = nsti

def get_result(result, responsevariable, floorstart, floorend, peaktype):
    if responsevariable == 'g':
        vector1 = result.gx
        vector2 = result.gy
        vector = np.hstack((vector1, vector2))
        vectorhead = ['Ground Acceleration in X', 'Ground Acceleration in Y']
    elif responsevariable == 'd':
        vector1 = result.dx[:, floorstart:floorend+1]
        vector2 = result.dy[:, floorstart:floorend+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("Relative Displacement in X DOF No #"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("Relative Displacement in Y DOF No #"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    elif responsevariable == 'v':
        vector1 = result.vx[:, floorstart:floorend+1]
        vector2 = result.vy[:, floorstart:floorend+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("Relative Velocity in X DOF No #"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("Relative Velocity in Y DOF No #"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    elif responsevariable == 'a':
        vector1 = result.ax[:, floorstart:floorend+1]
        vector2 = result.ay[:, floorstart:floorend+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("Relative Acceleration in X DOF No #"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("Relative Acceleration in Y DOF No #"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    elif responsevariable == 'aa':
        vector1 = result.aax[:, floorstart:floorend+1]
        vector2 = result.aay[:, floorstart:floorend+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("Absolute Acceleration in X DOF No #"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("Absolute Acceleration in Y DOF No #"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    elif responsevariable == 'f':
        vector1 = result.fx[:, floorstart:floorend+1]
        vector2 = result.fy[:, floorstart:floorend+1]
        vector = np.hstack((vector1, vector2))
        vectorheadx = [("Base Shear in X DOF No #"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorheady = [("Base Shear in Y DOF No #"+str(x+1)) for x in range(floorstart,floorend+1)]
        vectorhead = vectorheadx + vectorheady
    elif responsevariable == 'en':
        vector1 = result.eki
        vector2 = result.edi
        vector3 = result.esi
        vector4 = result.eii
        vector = np.hstack((vector1, vector2, vector3, vector4))
        vectorhead = ['Kinetic Energy', 'Dissipative Energy', 'Strain Energy', 'Input Energy']
    elif responsevariable == 'er':
        vector = result.errori
        vectorhead = ['Relative Error']

    if peaktype == 1:
        responsevalues = vector.max(axis=0)
    elif peaktype == 0:
        responsevalues = vector
    return responsevalues, vectorhead

def result_viewer(result, model, results_type, analysis_folder):
    results_type = results_type.split(',')
    i, j = 0, 0
    peakvalues = None
    timeresponses = None
    peakhead = None
    for rpattern in results_type:
        responsevariable, floorstart, floorend, peaktype = pattern_reader(rpattern.strip(), model.nst)
        responsevalues, vectorhead = get_result(result, responsevariable, floorstart, floorend, peaktype)
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
        timeresponses = pd.DataFrame(timeresponses)
        print(timeresponsehead)
        timeresponses.columns = timeresponsehead
        timeresponses.to_csv(analysis_folder + "\\TimeDomainResponses.csv", mode='w', sep=',')
    return peakvalues, peakhead

def pattern_reader(rpattern, nst):
    pattern_res = '^([p]|)([a-z]+|[f,e,r]+)([0-9\-\*\s]+|b$)'
    flag = 0
    for i in range(0,23):
        pattern_check = re.findall(pattern_res, rpattern)
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
                floorstart = nst
                floorend = nst
                flag = 1
                break
            else:
                responsevariable += 'b'
                flag = 0
                break
            
    if flag == 0:
        print("ERROR: Discrepancy in response variable keys.")
        sys.exit()
    
    return responsevariable, floorstart, floorend, peaktype

def createfolder(folder):
    try:      
        if os.path.exists('results\\' + folder):
            for i in range(1,1000):
                if not os.path.exists('results\\' + folder + '-' + str(i)):
                    folder = folder + '-' + str(i)
                    folder_input = input('Folder already exists! Enter new folder name ['+ folder + '] :')
                    if folder_input == "":
                        os.makedirs('results\\' + folder)
                    else:
                        folder = folder_input
                        os.makedirs('results\\' + folder)
                    break
        else:
            os.makedirs('results\\' + folder)
    except OSError:
        print('Error creating directory')
    
    return folder





            

