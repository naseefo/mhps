

from data.earthquakes.eqs import find_eq_from_db, add_eq_to_db, del_eq_from_db, lis_eq_from_db
import numpy as np
import re
import pandas as pd
import numexpr as ne
import os
import math


def sample():
    print("I am sample")

def get_total_excitations(eq):
    count = 1
    if eq == "Excitations.csv":
        f1 = open('Excitations.csv', 'r')
        for eq_line in f1:
            if eq_line[0] != "#":
                if eq_line.strip():
                    count, commands = pattern_reader(eq_line.strip(), count, True)
        f1.close()
    else:
        count, command = pattern_reader(eq.strip(), count, True)
    return count-1

def get_earthquake_list(eq):

    # pattern_eq_1 = "EQ"
    # pattern_check_1 = re.findall(pattern_eq_1, eq_string)

    ref_count = 1
    if eq == "Excitations.csv":
        f1 = open('Excitations.csv', 'r')
        for eq_line in f1:
            if eq_line[0] != "#":
                if eq_line.strip():
                    ref_count, commands = pattern_reader(eq_line.strip(), ref_count, False)
                    # print(commands)
                    for command in commands:
                        if command[0] == 1:
                            ref, xg, yg, zg, dt, ndiv, ndt = eq_finder(command)
                            yield ref, xg, yg, zg, dt, ndiv, ndt
                        else:
                            ref, xg, yg, dt, ndiv, ndt = harmonic_finder(command)
                            zg = np.zeros(xg.size, dtype='float')
                            yield ref, xg, yg, zg, dt, ndiv, ndt
        f1.close()
    else:
        ref_count, commands = pattern_reader(eq.strip(), ref_count, False)
        for command in commands:
            if command[0] == 1:
                ref, xg, yg, zg, dt, ndiv, ndt = eq_finder(command)
                yield ref, xg, yg, zg, dt, ndiv, ndt
            else:
                ref, xg, yg, dt, ndiv, ndt = harmonic_finder(command)
                zg = np.zeros(xg.size, dtype='float')
                yield ref, xg, yg, zg, dt, ndiv, ndt
    
def eq_path_finder(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)    

def eq_finder(command):
    # eq_type, ref_count, xg_filename, yg_filename, dt, unit, scale, dur, ndiv
    # [1, ref_count, xg_file, yg_file, zg_file, dt, unit, scale, dur, ndiv, rot]
    # eq_type, ref_count, xg_filename, yg_filename, zg_filename, dt, unit, scale, dur, ndiv, alpha1, alpha2, strike, rot
    # eq_type, ref_count, xg_filename, yg_filename, zg_filename, dt, unit, scale, dur, ndiv, alpha1, alpha2, strike, rotx, roty

    alpha1 = command[10]
    alpha2 = command[11]
    strike = command[12]
    rotx = command[13]
    roty = command[14]
    dirn = command[15]
    rot_flag = command[16]

    

    path = os.path.join('data', 'earthquakes', 'not_classified')
    ndt = int(command[8]/command[5]) + 1
    scale = command[7]
   

    if command[2] != '' and command[3] =='':
        # xg = pd.read_csv(path + command[2], header=None).values
        xg = pd.read_csv(eq_path_finder(command[2], os.path.join("data", "earthquakes")), header=None).values
        xg = xg*scale*scale_finder(command[6])
        yg = np.zeros(xg.size, dtype='float')
        # xg = np.insert(xg, 0, 0.0, axis=0)
        # yg = np.insert(yg, 0, 0.0, axis=0)
        len_xg = xg.size
        len_yg = yg.size
        ndt = min(ndt, len_xg, len_yg)
        if len_xg >= ndt:
            xg = xg[:ndt]
        if len_yg >= ndt:
            yg = yg[:ndt]
    elif command[3] != '' and command[2] =='':
        # yg = pd.read_csv(path + command[3], header=None).values
        yg = pd.read_csv(eq_path_finder(command[3], os.path.join("data", "earthquakes")), header=None).values
        yg = yg*scale*scale_finder(command[6])
        xg = np.zeros(yg.size, dtype='float')
        # xg = np.insert(xg, 0, 0.0, axis=0)
        # yg = np.insert(yg, 0, 0.0, axis=0)
        len_xg = xg.size
        len_yg = yg.size
        if len_xg >= ndt:
            xg = xg[:ndt]
        if len_yg >= ndt:
            yg = yg[:ndt]
    else:
        # xg = pd.read_csv(path + command[2], header=None).values
        h1 = pd.read_csv(eq_path_finder(command[2], os.path.join("data", "earthquakes")), header=None).values
        # yg = pd.read_csv(path + command[3], header=None).values
        h2 = pd.read_csv(eq_path_finder(command[3], os.path.join("data", "earthquakes")), header=None).values
        h1 = h1*scale*scale_finder(command[6])
        h2 = h2*scale*scale_finder(command[6])
        
        
        # xg = np.insert(xg, 0, 0.0, axis=0)
        # yg = np.insert(yg, 0, 0.0, axis=0)
        len_xg = h1.size
        len_yg = h2.size
        ndt = min(ndt, len_xg, len_yg)
        # print(ndt)
        if len_xg >= ndt:
            h1 = h1[:ndt]
        if len_yg >= ndt:
            h2 = h2[:ndt]

        if rot_flag == True:
            # print(alpha1, alpha2, strike, rotx, roty)
            # print(strike - alpha1 - rotx)
            # print(strike - alpha2 - roty)
            beta1 = (strike - alpha1 - rotx)*math.pi/180.0
            beta2 = (strike - alpha2 - roty)*math.pi/180.0
            # print(beta1)
            # print(beta2)
            # Page 25/D-17 AK Chopra: Selecting and Scaling Earthquake Ground Motions for Performing Response-History Analyses 
            yg = h1*math.cos(beta1) - h2*math.sin(beta2)  # Fault-Normal component
            xg = h1*math.sin(beta1) + h2*math.cos(beta2)  # Fault-Parallel component
        elif rot_flag == False:
            xg = h1
            yg = h2

        if dirn == "y":
            xg = np.zeros(xg.size, dtype='float')
        elif dirn == "x":
            yg = np.zeros(yg.size, dtype='float')

    if command[4] != '':
        zg = pd.read_csv(eq_path_finder(command[4], os.path.join("data", "earthquakes")), header=None).values
        zg = zg*scale*scale_finder(command[6])
        len_xg = xg.size
        len_yg = yg.size
        len_zg = zg.size
        ndt = min(ndt, len_xg, len_yg, len_zg)
        # print(len_zg)
        # print(ndt)
        if len_zg >= ndt:
            zg = zg[:ndt]
        if len_zg <= ndt:
            xg = xg[:ndt]
            yg = yg[:ndt]
        # print(xg.size)
        # print(xg/9.81)
        # print(yg.size)
        # print(yg/9.81)
        # print(zg.size)
        # print(zg/9.81)
    else:
        zg = np.zeros(xg.size, dtype='float')
    
    



    xg = excitation_divide(xg, command[9])
    yg = excitation_divide(yg, command[9])
    zg = excitation_divide(zg, command[9])
    dt = command[5]/command[9]
    return command[1], xg, yg, zg, dt, command[9], ndt

def excitation_divide(eq_data,ndiv):

    res = np.empty(ndiv*(eq_data.size-1)+1,dtype='float')
    start = eq_data[:-1]
    stop = np.concatenate((eq_data[1:-1],eq_data[-1,None]))
    s0 = start[:,None]
    s1 = stop[:,None]
    r = np.arange(ndiv)
    mform = ne.evaluate('((1.0/ndiv)*(s1-s0))*r+s0') 
    res[:-1] = mform.reshape(ndiv*(eq_data.size-1))
    res[-1] = eq_data[-1]

    return res

def harmonic_finder(command):
    # eq_type, ref_count, ax, tx, phx, ay, ty, phy, dt, unit, dur, ndiv
    
    
    ndt = command[10]/command[8]
    
    
    if command[2] != '' and command[5] =='':
        
        t = np.arange(0, command[10], command[8])
        xg = command[2]*np.sin(2.0*np.pi/command[3]*t + command[4])
        
        xg = xg*scale_finder(command[9])
        yg = np.zeros(xg.size, dtype='float')
        # len_xg = xg.size
        # len_yg = yg.size
        
        # ndt = min(ndt, len_xg, len_yg)
        # print(xg)
        # print(yg)
        # if len_xg >= ndt:
        #     xg = xg[:ndt]
            
        # if len_yg >= ndt:
        #     yg = yg[:ndt]
        #     print(command)    
        
    elif command[5] != '' and command[2] =='':
        t = np.arange(0, command[10], command[8])
        yg = command[5]*np.sin(2.0*np.pi/command[6]*t + command[7])
        yg = yg*scale_finder(command[9])
        xg = np.zeros(yg.size, dtype='float')
        # len_xg = xg.size
        # len_yg = yg.size
        # if len_xg >= ndt:
        #     xg = xg[:ndt]
        # if len_yg >= ndt:
        #     yg = yg[:ndt]
    else:
        
        t = np.arange(0, command[10], command[8])
        xg = command[2]*np.sin(2.0*np.pi/command[3]*t + command[4])
        yg = command[5]*np.sin(2.0*np.pi/command[6]*t + command[7])
        xg = xg*scale_finder(command[9])
        yg = yg*scale_finder(command[9])
        # len_xg = xg.size
        # len_yg = yg.size
        # ndt = min(ndt, len_xg, len_yg)
        # if len_xg >= ndt:
        #     xg = xg[:ndt]
        # if len_yg >= ndt:
        #     yg = yg[:ndt]
    
    xg = excitation_divide(xg, command[11])
    yg = excitation_divide(yg, command[11])
    dt = command[8]/command[11]

    return command[1], xg, yg, dt, command[11], len(t)

def pattern_reader(eq_string, ref_count, count_flag):
    # ^([\w\.\s]+),([\w\.\s]+),([\s0-9\.]+),([\sa-z0-9/]+),([\s0-9\.]+),([\s0-9\.]+)(,[\s0-9]+|)\;
    # ^([\w\.\s]+),([\w\.\s]+),([\s0-9\.]+),([\sa-z0-9/]+),([\s0-9\.]+),([\s0-9\.]+),([\s0-9]+)(,[\s0-9\.]+|)\;
    # ^([\w\.\s]+),([\w\.\s]+),([\s0-9\.]+),([\sa-z0-9/]+),([\s0-9\.]+),([\s0-9\.]+),([\s0-9]+)(,[\s0-9\.]+|)\;
    # ^([\w\.\s\-]+),([\w\.\s\-]+),([\w\.\s\-]+),([\s0-9\.]+),([\sa-z0-9/]+),([\s0-9\.]+),([\s0-9\.]+),([\s0-9]+)(,[\s0-9\.]+|)\;
    # ^([\w\.\s\-]+),([\w\.\s\-]+),([\w\.\s\-]+),([\s0-9\.]+),([\sa-z0-9/]+),([\s0-9\.]+),([\s0-9\.]+),([\s0-9]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+)(,[\s0-9\.]+|)\;
    # ^([\w\.\s\-]+),([\w\.\s\-]+),([\w\.\s\-]+),([\s0-9\.]+),([\sa-z0-9/]+),([\s0-9\.]+),([\s0-9\.]+),([\s0-9]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\w])(,[\w]+|)\;
    pattern_eq_1 = "^([\w\.\s]+),([\w\.\s]+),([\s0-9\.]+),([\sa-z0-9/]+),([\s0-9\.]+),([\s0-9\.]+),([\s0-9]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s\w]+)(,[\s\w]+|)\;"
    pattern_eq_2 = "^([\w\.\s]+)[,\s]+([s0-9\.]+)[,\s]+([a-z0-9/]+)[,\s]([\s0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([\w])[,\s]+([0-9]+);"
    pattern_eq_3 = "^([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([\w/]+)[,\s]+([0-9\.]+)[,\s]+([0-9\.]+)\;"
    pattern_eq_4 = "([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([\w/]+)[,\s]+([0-9\.]+)[,\s]+([\w])[,\s]+([0-9\.]+)\;"
    pattern_eq_5 = "^([0-9\.]+)[,\s]+([0-9\.]+):([0-9\.]+):([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([0-9\.]+)[,\s]+([\w/]+)[,\s]+([0-9\.]+)[,\s]+([\w])[,\s]+([0-9\.]+)\;"
    pattern_eq_6 = "^([\w\.\s\-]+),([\w\.\s\-]+),([\w\.\s\-]+),([\s0-9\.]+),([\sa-z0-9/]+),([\s0-9\.]+),([\s0-9\.]+),([\s0-9]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s0-9\.\-]+),([\s\w]+)(,[\s\w]+|)\;"

    pattern_check_1 = re.findall(pattern_eq_1, eq_string)
    pattern_check_2 = re.findall(pattern_eq_2, eq_string)
    pattern_check_3 = re.findall(pattern_eq_3, eq_string)
    pattern_check_4 = re.findall(pattern_eq_4, eq_string)
    pattern_check_5 = re.findall(pattern_eq_5, eq_string)
    pattern_check_6 = re.findall(pattern_eq_6, eq_string)

    if pattern_check_1:
        ref_count, command = pattern_eq_1_handler(pattern_check_1[0], ref_count, count_flag) # [1, ref_count, xg_file, yg_file, zg_file, dt, unit, scale, dur, ndiv, alpha1, alpha2, strike, rot]
    elif pattern_check_2:
        ref_count, command = pattern_eq_2_handler(pattern_check_2[0], ref_count, count_flag) # [1, ref_count, xg_file, yg_file, zg_file, dt, scale, dur, ndiv, alpha1, alpha2, strike, rot]
    elif pattern_check_3:
        ref_count, command = pattern_eq_3_handler(pattern_check_3[0], ref_count, count_flag) # [2, ref_count, AXI, TX1, PHX1, AY1, TY1, PHY1, DT, UNIT, DUR, NDIV]
    elif pattern_check_4:
        ref_count, command = pattern_eq_4_handler(pattern_check_4[0], ref_count, count_flag) # [2, ref_count, AXI, TX1, PHX1, AY1, TY1, PHY1, DT, UNIT, DUR, NDIV]
    elif pattern_check_5:
        ref_count, command = pattern_eq_5_handler(pattern_check_5[0], ref_count, count_flag) # [2, ref_count, AXI, TX1, PHX1, AY1, TY1, PHY1, DT, UNIT, DUR, NDIV]
    elif pattern_check_6:
        ref_count, command = pattern_eq_6_handler(pattern_check_6[0], ref_count, count_flag) # [1, ref_count, xg_file, yg_file, zg_file, dt, unit, scale, dur, ndiv, alpha1, alpha2, strike, rot]

    if command is None:
        print("ERROR: Excitation file in wrong format!")
    
    return ref_count, command

def pattern_eq_1_handler(eq_string, ref_count, count_flag):
    # Format: ax, ay, dt, unit, scale, dur, ndiv
    # e.g. Cent_acc_00.txt, Cent_acc_00.txt, 0.02, cm/s2, 1.0, 20, 120;
    # print("I am in pattern 1 handler")

    if count_flag:
        ref_count += 1
        command = []
        return ref_count, command

    eq_type = 1
    xg_filename = eq_string[0].strip(', ')
    yg_filename = eq_string[1].strip(', ')
    zg_filename = ''
    dt = float(eq_string[2].strip(', '))
    unit = eq_string[3].strip(', ')
    scale = float(eq_string[4].strip(', '))
    dur = float(eq_string[5].strip(', '))
    ndiv = int(eq_string[6].strip(', '))
    alpha1 = float(eq_string[7].strip(', '))
    alpha2 = float(eq_string[8].strip(', '))
    strike = float(eq_string[9].strip(', '))
    rotx = float(eq_string[10].strip(', '))
    roty = float(eq_string[11].strip(', '))
    dirn = eq_string[12].strip(', ')
    rot_flag = eq_string[13].strip(', ')
    if rot_flag == 'ry':
        rot_flag = True
    elif rot_flag == 'rn':
        rot_flag = False
    

    command = [(eq_type, ref_count, xg_filename, yg_filename, zg_filename, dt, unit, scale, dur, ndiv, alpha1, alpha2, strike, rotx, roty, dirn, rot_flag)]
    ref_count += 1
    return ref_count, command

def pattern_eq_2_handler(eq_string, ref_count, count_flag):
    # Format: a, dt, unit, dur, dirn, ndiv;
    # e.g. Cent_acc_00.csv, 0.02, g, 1.0, 10, y, 120;
    # print("I am in pattern 2 handler")

    if count_flag:
        ref_count += 1
        command = []
        return ref_count, command

    eq_type = 1

    filename = eq_string[0].strip(', ')
    dt = float(eq_string[1].strip(', '))
    unit = eq_string[2].strip(', ')
    scale = float(eq_string[3].strip(', '))
    dur = float(eq_string[4].strip(', '))
    dirn = eq_string[5].strip(', ')
    ndiv = int(eq_string[6].strip(', '))
    
    if dirn == 'x':
        xg_filename = filename
        yg_filename = ''
    else:
        xg_filename = ''
        yg_filename = filename
    zg_filename = ''

    alpha1 = 0.0
    alpha2 = 0.0
    strike = 0.0
    rotx = 0.0
    roty = 0.0
    rot_flag = False

    command = [(eq_type, ref_count, xg_filename, yg_filename, zg_filename, dt, unit, scale, dur, ndiv, alpha1, alpha2, strike, rotx, roty, dirn, rot_flag)]
    ref_count += 1
    return ref_count, command

def pattern_eq_3_handler(eq_string, ref_count, count_flag):
    # Format: AXI, TX1, PHX1, AY1, TY1, PHY1, DT, UNIT, DUR, NDIV
    # e.g. 3.0, 0.1, 1.5, 2.0, 0.2, 1.0, 0.02, cm/s2, 20, 120;
    # print("I am in pattern 3 handler")
    if count_flag:
        ref_count += 1
        command = []
        return ref_count, command

    eq_type = 2
    
    ax = float(eq_string[0].strip(', '))
    tx = float(eq_string[1].strip(', '))
    phx = float(eq_string[2].strip(', '))
    ay = float(eq_string[3].strip(', '))
    ty = float(eq_string[4].strip(', '))
    phy = float(eq_string[5].strip(', '))
    dt = float(eq_string[6].strip(', '))
    unit = eq_string[7].strip(', ')
    dur = float(eq_string[8].strip(', '))
    ndiv = int(eq_string[9].strip(', '))

    command = [(eq_type, ref_count, ax, tx, phx, ay, ty, phy, dt, unit, dur, ndiv)]
    ref_count += 1
    return ref_count, command

def pattern_eq_4_handler(eq_string, ref_count, count_flag):
    # Format: A, T, PH, DT, UNIT, DUR, DIRN, NDIV
    # e.g. 3.0, 0.1, 1.5, 0.02, cm/s2, 20, x, 120;
    # print("I am in pattern 4 handler")
    if count_flag:
        ref_count += 1
        command = []
        return ref_count, command

    eq_type = 2
    
    if eq_string[6].strip(', ') == 'x':
        ax = float(eq_string[0].strip(', '))
        tx = float(eq_string[1].strip(', '))
        phx = float(eq_string[2].strip(', '))
        ay = ''
        ty = ''
        phy = ''
    else:
        ax = ''
        tx = ''
        phx = ''
        ay = float(eq_string[0].strip(', '))
        ty = float(eq_string[1].strip(', '))
        phy = float(eq_string[2].strip(', '))
    
    dt = float(eq_string[3].strip(', '))
    unit = eq_string[4].strip(', ')
    dur = float(eq_string[5].strip(', '))
    ndiv = int(eq_string[7].strip(', '))

    command = [(eq_type, ref_count, ax, tx, phx, ay, ty, phy, dt, unit, dur, ndiv)]
    ref_count += 1
    return ref_count, command

def pattern_eq_5_handler(eq_string, ref_count, count_flag):
    # Format: A, T0:DTT:TN, PH, DT, UNIT, DUR, DIR, NDIV
    # e.g. 3.0, 0.1:0.01:5.0, 1.5, 0.02, cm/s2, 20, x, 120;
    # print("I am in pattern 5 handler")

    eq_type = 2

    t0 = float(eq_string[1].strip(', '))
    dtt = float(eq_string[2].strip(', '))
    tn = float(eq_string[3].strip(', '))

    if count_flag:
        ref_count += len(np.arange(t0, tn+dtt, dtt))
        command = []
        return ref_count, command

    command = []
    for t in np.arange(t0, tn+dtt, dtt):
        if eq_string[8].strip(', ') == 'x':
            ax = float(eq_string[0].strip(', '))
            tx = t
            phx = float(eq_string[4].strip(', '))
            ay = ''
            ty = ''
            phy = ''
        else:
            ax = ''
            tx = ''
            phx = ''
            ay = float(eq_string[0].strip(', '))
            ty = t
            phy = float(eq_string[4].strip(', '))
            
        dt = float(eq_string[5].strip(', '))
        unit = eq_string[6].strip(', ')
        dur = float(eq_string[7].strip(', '))
        ndiv = int(eq_string[9].strip(', '))

        command.append((eq_type, ref_count, ax, tx, phx, ay, ty, phy, dt, unit, dur, ndiv))
        ref_count += 1
    return ref_count, command

def pattern_eq_6_handler(eq_string, ref_count, count_flag):
    # Format: ax, ay, dt, unit, scale, dur, ndiv
    # e.g. Cent_acc_00.txt, Cent_acc_00.txt, Cent_acc_UP.txt, 0.02, cm/s2, 1.0, 20, 120;
    # print("I am in pattern 1 handler")
    # ref_count, command = pattern_eq_6_handler(pattern_check_6[0], ref_count, count_flag)

    if count_flag:
        ref_count += 1
        command = []
        return ref_count, command

    eq_type = 1
    xg_filename = eq_string[0].strip(', ')
    yg_filename = eq_string[1].strip(', ')
    zg_filename = eq_string[2].strip(', ')

    dt = float(eq_string[3].strip(', '))
    unit = eq_string[4].strip(', ')
    scale = float(eq_string[5].strip(', '))
    dur = float(eq_string[6].strip(', '))
    ndiv = int(eq_string[7].strip(', '))
    alpha1 = float(eq_string[8].strip(', '))
    alpha2 = float(eq_string[9].strip(', '))
    strike = float(eq_string[10].strip(', '))
    rotx = float(eq_string[11].strip(', '))
    roty = float(eq_string[12].strip(', '))
    dirn = eq_string[13].strip(', ')
    rot_flag = eq_string[14].strip(', ')
    if rot_flag == 'ry':
        rot_flag = True
    elif rot_flag == 'rn':
        rot_flag = False
    
    command = [(eq_type, ref_count, xg_filename, yg_filename, zg_filename, dt, unit, scale, dur, ndiv, alpha1, alpha2, strike, rotx, roty, dirn, rot_flag)]
    ref_count += 1
    return ref_count, command

def get_param_set_list(var_param):
    pass


def scale_finder(units):
    units = units.lower()
    if units == "cm/s2":
        scale = float(0.01)
    elif units == "m/s2":
        scale = float(1.0)
    elif units == "g":
        scale = float(9.81)
    else:
        scale = "undefined"
    return scale

def search_eq_working_folder(filename):
    file = glob.glob(filename)
    if file:
        data = pd.read_csv(filename, header=None).values
    else:
        return None
    return data