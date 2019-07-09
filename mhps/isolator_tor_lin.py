
import numpy as np
import math
from data.defaults.param_manager import default
from data.defaults.param_manager import get_default_param_values
import cProfile
import io
import pstats
import pandas as pd
from mhps.postprocessor import ResultFixedXY, ModelInfo
from math import sin, cos, tan, atan, pow, exp, sqrt, pi
import matplotlib.pyplot as plt


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

default_values = get_default_param_values()

def read_const_param_torsion(const_param_file):

    """
    This function reads all the constant set of parameters for the
    program which does not require to be read in iteration of the
    parametric study.
    """
    f1  = open(const_param_file, 'r')
    
    fm = float(f1.readline().split(':')[1])
    nb = int(f1.readline().split(':')[1])
    x = np.array([float(x.strip()) for x in f1.readline().split(':')[1:][0].strip().split(',')])
    y = np.array([float(x.strip()) for x in f1.readline().split(':')[1:][0].strip().split(',')])
    xb = np.array([float(x.strip()) for x in f1.readline().split(':')[1:][0].strip().split(',')])
    yb = np.array([float(x.strip()) for x in f1.readline().split(':')[1:][0].strip().split(',')])
    
    f1.close()
    
    return fm, nb, x, y, xb, yb

def read_ss_torsion_var_param(var_param_file):

    """
    This function reads all the set of parameters for the parametric
    studies and stores it in an array. It will be a one-time allocation
    to increase speed.
    """

    dtype1 = np.dtype([('IJK', 'i4'), ('TX1', 'f4'), ('ZETA','f4'), ('EXD','f4'), ('WRWX','f4')])

    ijk, tx1, zeta, exd, wrwx = np.loadtxt \
    (var_param_file, delimiter=',', usecols=range(5),\
     skiprows=1, unpack = True, dtype=dtype1)
    
    try:
        total_param = len(ijk)
    except:
        total_param = 1


    for i in range(0, total_param):
        try:
            yield ijk[i], tx1[i], zeta[i], exd[i], wrwx[i]
        except:
            yield ijk, tx1, zeta, exd, wrwx


def get_total_variable_parameter_set_torsion(var_param_file):
 
    dtype1 = np.dtype([('IJK', 'i4'), ('TX1', 'f4'), ('ZETA','f4'), ('EXD','f4'), ('WRWX','f4')])

    ijk, tx1, zeta, exd, wrwx = np.loadtxt(var_param_file, delimiter=',', usecols=range(5), skiprows=1, unpack = True, dtype=dtype1)
    try:
        return len(ijk)
    except:
        return 1


def superstructure_propxy_t(tx1, zeta, exd, wrwx, fm, nb, x, y, xb, yb):

    x = np.asarray(x, dtype=np.dtype('d'), order='F')
    y = np.asarray(y, dtype=np.dtype('d'), order='F')
    
    """
    Calculation of superstructure mass, stiffness and damping matrices
    in X- and Y- direction.
    """
    
    zeta = np.ones(6)*zeta
    wx = 2.0*pi/tx1
    wr = wx*wrwx
    am = [fm, fm, fm]

    # Calculation of msss matrix in x-direction
    sm = np.zeros(shape=(6,6), dtype=np.dtype('d'), order='F')
    sm = np.diag(am)    # ------> Mass matrix in x-direction

    # Calculation of stiffness matrix in x- and y-direction
    
    sk = np.zeros(shape=(6, 6), dtype=np.dtype('d'), order='F')
    ckx = np.zeros(shape=(4, ), dtype=np.dtype('d'), order='F')
    cky = np.zeros(shape=(4, ), dtype=np.dtype('d'), order='F')

    B = x[0] - x[1]

    akx = fm*pow(wx, 2.0)

    print("Width = %8.4f m, Time period = %8.4f s, Angular Frequency = %8.4f, Stiffness = %8.3f N/m"%(B, tx1, wx, akx))
      
    ckx[0] = 0.25*akx*(1.0 + 2.0*exd)
    ckx[1] = 0.25*akx*(1.0 - 2.0*exd)
    ckx[2] = 0.25*akx*(1.0 - 2.0*exd)
    ckx[3] = 0.25*akx*(1.0 + 2.0*exd)

    cky = ckx

    kxxs = np.sum(ckx)
    kxys = 0.0
    kxts = -1.0*np.sum(ckx*y)
    kyxs = 0.0
    kyys = np.sum(cky)
    kyts = np.sum(cky*x)
    ktxs = kxts
    ktys = kyts
    ktts = np.sum(ckx*np.square(y)) + np.sum(cky*np.square(x))
    sk[0,0] = kxxs
    sk[0,1] = kxys
    sk[0,2] = kxts
    sk[1,0] = kyxs
    sk[1,1] = kyys
    sk[1,2] = kyts
    sk[2,0] = ktxs
    sk[2,1] = ktys
    sk[2,2] = ktts

    fmr = sk[2,2]/pow(wr, 2.0)
    sm[2,2] = fmr

    # Calculation of damping matrix in x-direction
    D = np.dot(sk[0:3,0:3],np.linalg.inv(sm[0:3,0:3]))
    e,v = np.linalg.eig(D)
    idx = e.argsort()[::1]
    e = e[idx]
    #v = v[:,idx]

    T = 2.0*math.pi/np.sqrt(e)
    del e, idx, v

    nst = 3
    W = 2.0*math.pi/T
    c_w = (np.ones((nst,nst))*W[None,:]).T
    c_w = np.power(c_w,np.arange(-1,2*nst-2,2))
    a_i = np.dot(np.linalg.inv(c_w),2.0*zeta[0:nst].reshape(nst,1))
    cd = np.zeros(shape=(nst,nst),dtype=np.dtype('d'), order='F')
    for i in range(nst):
        cd[0:nst,0:nst] += np.dot(np.linalg.matrix_power(D,i), sm[0:nst,0:nst])*a_i[i]  # Damping matrix in x-direction
    del W, c_w, a_i, D, T

    return sm, sk, cd



def addlinear_iso_tor(sm, sk, cd, rmbm, tbx, zetabx, ebxd, wrwxb, xb, yb, nb):
                 
    """
    Calculation of superstructure mass, stiffness and damping matrices
    in X- and Y- direction with linear isolator.
    """
    print("RMBM = %3.2f, TBX = %3.2f, ZETABX = %3.2f, EBXD = %3.2f, WRWXB = %3.2f" %(rmbm, tbx, zetabx, ebxd, wrwxb))

    xb = np.asarray(xb, dtype=np.dtype('d'), order='F')
    yb = np.asarray(yb, dtype=np.dtype('d'), order='F')
    
    wbx = 2.0*pi/tbx
    fm = sm[0,0]
    bm = sm[0,0]*rmbm
    tm = fm + bm

    ckabx = (fm + bm)*wbx*wbx

    bkx = np.zeros(shape=(4, ), dtype=np.dtype('d'), order='F')
    bky = np.zeros(shape=(4, ), dtype=np.dtype('d'), order='F')

    bkx[0] = 0.25*ckabx*(1.0 + 2.0*ebxd)
    bkx[1] = 0.25*ckabx*(1.0 - 2.0*ebxd)
    bkx[2] = 0.25*ckabx*(1.0 - 2.0*ebxd)
    bkx[3] = 0.25*ckabx*(1.0 + 2.0*ebxd)

    bky = bkx

    kxxb = np.sum(bkx)
    kxyb = 0.0
    kxtb = -1.0*np.sum(bkx*yb)
    kyxb = 0.0
    kyyb = np.sum(bky)
    kytb = np.sum(bky*xb)
    ktxb = kxtb
    ktyb = kytb
    kttb = np.sum(bkx*np.square(yb)) + np.sum(bky*np.square(xb))
    sk[3,3] = kxxb
    sk[3,4] = kxyb
    sk[3,5] = kxtb
    sk[4,3] = kyxb
    sk[4,4] = kyyb
    sk[4,5] = kytb
    sk[5,3] = ktxb
    sk[5,4] = ktyb
    sk[5,5] = kttb

    smb = np.zeros(shape=(3, 3), dtype=np.dtype('d'), order='F')
    wrb = wbx*wrwxb
    bmr = sk[5,5]/pow(wrb, 2.0) - sm[2,2]
    smb[0,0] = bm + fm
    smb[1,1] = bm + fm
    smb[2,2] = bmr + sm[2,2]
    sm[0:3,3:6] = sm[0:3,0:3]
    sm[3:6,0:3] = sm[0:3,0:3]
    sm[3:6,3:6] = smb
    del smb

    bcx = np.zeros(shape=(4, ), dtype=np.dtype('d'), order='F')
    bcy = np.zeros(shape=(4, ), dtype=np.dtype('d'), order='F')

    cdabx = 2.0*zetabx*tm*wbx
    bcx[0] = 0.25*cdabx
    bcx[1] = 0.25*cdabx
    bcx[2] = 0.25*cdabx
    bcx[3] = 0.25*cdabx

    bcy = bcx

    cxxb = np.sum(bcx)
    cxyb = 0.0
    cxtb = -1.0*np.sum(bcx*yb)
    cyxb = 0.0
    cyyb = np.sum(bcy)
    cytb = np.sum(bcy*xb)
    ctxb = cxtb
    ctyb = cytb
    cttb = np.sum(bcx*np.square(yb)) + np.sum(bcy*np.square(xb))
    cd[3,3] = cxxb
    cd[3,4] = cxyb
    cd[3,5] = cxtb
    cd[4,3] = cyxb
    cd[4,4] = cyyb
    cd[4,5] = cytb
    cd[5,3] = ctxb
    cd[5,4] = ctyb
    cd[5,5] = cttb
    
    return sm, sk, cd



#@profile
def simulator_linear_tor(ref, xg, yg, dt, ndiv, ndt, ijk, sm, sk, cd, x, y, xb, yb, nb):

    gamma = 0.5
    beta = 1/6

    sm_inv = np.linalg.inv(sm)
    print(sm_inv)
    # sm_diag = np.diag(-1.0*sm).reshape(nst,1)
    r = np.zeros((6,3), dtype=np.dtype('d'), order='F')
    r[3:6,0:3] = np.diag([1.0, 1.0, 1.0])

    time = np.zeros((1, ndt), dtype=np.dtype('d'), order='F')
    
    d = np.zeros((3, ndt), dtype=np.dtype('d'), order='F')
    v = np.zeros((3, ndt), dtype=np.dtype('d'), order='F')
    a = np.zeros((3, ndt), dtype=np.dtype('d'), order='F')
    aa = np.zeros((3, ndt), dtype=np.dtype('d'), order='F')

    db = np.zeros((3, ndt), dtype=np.dtype('d'), order='F')
    vb = np.zeros((3, ndt), dtype=np.dtype('d'), order='F')
    ab = np.zeros((3, ndt), dtype=np.dtype('d'), order='F')
    aab = np.zeros((3, ndt), dtype=np.dtype('d'), order='F')

    gx = np.zeros((1, ndt), dtype=np.dtype('d'), order='F')
    gy = np.zeros((1, ndt), dtype=np.dtype('d'), order='F')
    f = np.zeros((ndt, 3), dtype=np.dtype('d'), order='F')

    dd = np.zeros((6, 1), dtype=np.dtype('d'), order='F')
    dv = np.zeros((6, 1), dtype=np.dtype('d'), order='F')
    da = np.zeros((6, 1), dtype=np.dtype('d'), order='F')
    dp = np.zeros((3, 1), dtype=np.dtype('d'), order='F')

    dcx = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    dcy = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    vcx = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    vcy = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    acx = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    acy = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    aacx = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    aacy = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')

    dbcx = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    dbcy = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    vbcx = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    vbcy = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    abcx = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    abcy = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    aabcx = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    aabcy = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')

    fcx = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')
    fcy = np.zeros((ndt, 4), dtype=np.dtype('d'), order='F')

    gx[0,0] = xg[0]
    gy[0,0] = yg[0]

    ek = np.zeros((ndt, 1), dtype=np.dtype('d'), order='F')
    ed = np.zeros((ndt, 1), dtype=np.dtype('d'), order='F')
    es = np.zeros((ndt, 1), dtype=np.dtype('d'), order='F')
    ei = np.zeros((ndt, 1), dtype=np.dtype('d'), order='F')
    error = np.zeros((ndt, 1), dtype=np.dtype('d'), order='F')

    d1 = np.ones((6, 1), dtype=np.dtype('d'), order='F')*0.0
    v1 = np.ones((6, 1), dtype=np.dtype('d'), order='F')*0.0
    a1 = np.ones((6, 1), dtype=np.dtype('d'), order='F')*0.0
    p1 = np.ones((6, 1), dtype=np.dtype('d'), order='F')*0.0
    ug1 = np.array([[xg[0]],[yg[0]], [0.0]])
    p1 = -1.0*np.dot(np.dot(sm, r), ug1)
    # print(p1)
    a1 = np.dot(sm_inv, p1 - np.dot(cd, v1) - np.dot(sk, d1))

    d2 = np.zeros((6, 1), dtype=np.dtype('d'), order='F')
    v2 = np.zeros((6, 1), dtype=np.dtype('d'), order='F')
    p2 = np.zeros((6, 1), dtype=np.dtype('d'), order='F')
    a2 = np.zeros((6, 1), dtype=np.dtype('d'), order='F')

    # na1x = (1.0/beta/np.power(dt,2))*sm + (gamma/beta/dt)*cd
    # na2x = (1.0/beta/dt)*sm + (gamma/beta - 1.0)*cd
    # na3x = (1.0/2.0/beta - 1.0)*sm + (gamma*dt/2.0/beta - dt)*cd

    na1x = (1.0/beta/np.power(dt,2))*sm + (gamma/beta/dt)*cd
    na2x = (1.0/beta/dt)*sm + gamma/beta*cd
    na3x = 1.0/2.0/beta*sm + (gamma*dt/2.0/beta - dt)*cd

    knx = sk + na1x
    knx_inv = np.linalg.inv(knx)

    index = 0

    t = 0.0
    time[0,index] = 0.0
    d[:,index] = d1[0:3,0]
    v[:,index] = v1[0:3,0]
    a[:,index] = a1[0:3,0]
    aa[0,index] = a1[0,0] + xg[0]
    aa[1,index] = a1[1,0] + yg[0]
    aa[2,index] = a1[2,0] + 0.0

    db[:,index] = d1[3:6,0]
    vb[:,index] = v1[3:6,0]
    ab[:,index] = a1[3:6,0]
    aab[0,index] = a1[3,0] + xg[0]
    aab[1,index] = a1[4,0] + yg[0]
    aab[2,index] = a1[5,0] + 0.0

    eki = 0.0
    edi = 0.0
    esi = 0.0
    eii = 0.0

    print("Linear-Mass Matrix")
    print(sm)
    print("Linear-Stiffness Matrix")
    print(sk)
    print("Linear-Damping Matrix")
    print(cd)

    for i in range(1,len(xg)):

        t += dt
        
        ug2 = np.array([[xg[i]],[yg[i]], [0.0]])
        p2 = -1.0*np.dot(np.dot(sm, r), ug2)
        dp = p2 - p1
        
        pcx1 = dp + np.dot(na2x, v1) + np.dot(na3x, a1)
        # pcx1[3:6,0] = pcx1[3:6,0] - 0.0
        dd = np.dot(knx_inv, pcx1)
        d2 = d1 + dd
        dv = (gamma/beta/dt)*dd - gamma/beta*v1 + dt*(1.0 - gamma/2.0/beta)*a1
        v2 = v1 + dv

        ep = p2 - np.dot(cd, v2) - np.dot(sk, d2)
        # ep[3:6,0] = ep[3:6,0] - 0.0
        a2 = np.dot(sm_inv, ep)

        eki = eki + 0.5*dt*(np.dot(np.dot(v2.T, sm),a2) + np.dot(np.dot(v1.T, sm),a1))
        edi = edi + 0.5*dt*(np.dot(np.dot(v2.T, cd),v2) + np.dot(np.dot(v1.T, cd),v1))
        esi = esi + 0.5*dt*(np.dot(np.dot(v2.T, sk),d2) + np.dot(np.dot(v1.T, sk),d1))
        eii = eii - 0.5*dt*(np.dot(np.dot(v2.T, sm), np.array([[xg[i]], [yg[i]], [0.0], [xg[i]], [yg[i]], [0.0]])) + np.dot(np.dot(v1.T, sm), np.array([[xg[i-1]], [yg[i-1]], [0.0], [xg[i-1]], [yg[i-1]], [0.0]])))

        d1, v1, p1, a1 = d2, v2, p2, a2 

        if not i%(ndiv):
            index += 1
            time[0,index] += t
            gx[0, index] = xg[i]
            gy[0, index] = yg[i]
            d[:, index] = d2[0:3,0]
            v[:, index] = v2[0:3,0]
            a[:, index] = a2[0:3,0]
            aa[0,index] = a2[0,0] + a2[3] + xg[i]
            aa[1,index] = a2[1,0] + a2[4] + yg[i]
            aa[2,index] = a2[2,0] + a2[5] +0.0

            db[:, index] = d2[3:6,0]
            vb[:, index] = v2[3:6,0]
            ab[:, index] = a2[3:6,0]
            aab[0,index] = a2[3,0] + xg[i]
            aab[1,index] = a2[4,0] + yg[i]
            aab[2,index] = a2[5,0] + 0.0

            for i in range(3):
                f[index, i] = sm[i,i]*aa[i,index] + sm[i+3, i+3]*aab[i, index]
      
            # Corner calculations
            for nc in range(0,4):
                dcx[index, nc] = d2[0,0] - y[nc]*d2[2,0]
                dcy[index, nc] = d2[1,0] + x[nc]*d2[2,0]

                vcx[index, nc] = v2[0,0] - y[nc]*v2[2,0]
                vcy[index, nc] = v2[1,0] + x[nc]*v2[2,0]

                acx[index, nc] = a2[0,0] - y[nc]*a2[2,0]
                acy[index, nc] = a2[1,0] + x[nc]*a2[2,0]

                aacx[index, nc] = aa[0,index] - y[nc]*aa[2,index]
                aacy[index, nc] = aa[1,index] + x[nc]*aa[2,index]

                dbcx[index, nc] = d2[3,0] - yb[nc]*d2[5,0]
                dbcy[index, nc] = d2[4,0] + xb[nc]*d2[5,0]

                vbcx[index, nc] = v2[3,0] - yb[nc]*v2[5,0]
                vbcy[index, nc] = v2[4,0] + xb[nc]*v2[5,0]

                abcx[index, nc] = a2[3,0] - yb[nc]*a2[5,0]
                abcy[index, nc] = a2[4,0] + xb[nc]*a2[5,0]

                aabcx[index, nc] = aab[0,index] - yb[nc]*aab[2,index]
                aabcy[index, nc] = aab[1,index] + xb[nc]*aab[2,index]

                fcx[index, nc] = f[index, 0] - yb[nc]*f[index, 2]
                fcy[index, nc] = f[index, 1] + xb[nc]*f[index, 2]


            
            ek[index, 0] = eki
            ed[index, 0] = edi
            es[index, 0] = esi
            ei[index, 0] = eii

            error[index, 0] = (edi + esi + eki - eii)/(abs(edi + esi) + eki + abs(eii))
    
    
    peakerror = max(abs(error))
    sumerror = sum(abs(error))
    peaktopaccx = max(abs(aa[0,:]))
    peaktopaccy = max(abs(aa[1,:]))
    peaktopacctheta = max(abs(aa[2,:]))
    peaktopdispx = max(abs(d[0,:]))
    peaktopdispy = max(abs(d[1,:]))
    peaktopdisptheta = max(abs(d[2,:]))
    peakbasedispx = max(abs(db[0,:]))
    peakbasedispy = max(abs(db[1,:]))
    peakbasedisptheta = max(abs(db[2,:]))

    print(" ")
    print("Simulation" + "\033[91m" + " SET%d-%d" %(ref, ijk) + "\033[0m" + ": Earthquake #: %d, Parameter #: %d" %(ref, ijk))
    print("Peak Error: % 8.6f" %(peakerror))
    print("Absolute Sum of Errors: % 8.6f" %(sumerror))
    print("Peak Top Floor Absolute Acceleration in X-Direction: % 8.6f m/s2" %(peaktopaccx))
    print("Peak Top Floor Absolute Acceleration in Y-Direction: % 8.6f m/s2" %(peaktopaccy))
    print("Peak Top Floor Absolute Acceleration in Theta-Direction: % 8.6f rad/s2" %(peaktopacctheta))
    print("Peak Top Floor Relative Displacement in X-Direction: % 8.6f cm" %(peaktopdispx*100.0))
    print("Peak Top Floor Relative Displacement in Y-Direction: % 8.6f cm" %(peaktopdispy*100.0))
    print("Peak Top Floor Relative Displacement in Theta-Direction: % 8.6f rad" %(peaktopdisptheta))
    print("Peak Isolator Displacement in X-Direction: % 8.6f cm" %(peakbasedispx*100.0))
    print("Peak Isolator Displacement in Y-Direction: % 8.6f cm" %(peakbasedispy*100.0))
    print("Peak Isolator Displacement in Theta-Direction: % 8.6f rad" %(peakbasedisptheta))
    
    t_s = np.hstack((d.T, v.T, a.T, aa.T))
    t_sc = np.hstack((dcx, dcy, vcx, vcy, acx, acy, aacx, aacy))
    t_b = np.hstack((db.T, vb.T, ab.T, aab.T))
    t_bc = np.hstack((dbcx, dbcy, vbcx, vbcy, abcx, abcy, aabcx, aabcy))
    f_b = f
    f_bc = np.hstack((fcx, fcy))
    

    result = ResultFixedXY(ref, ijk, time.T, gx.T, gyi = gy.T, eki = ek, edi = ed, esi = es, eii = ei, errori = error, t_si = t_s, t_bi = t_b, f_bi = f_b,t_sci = t_sc, t_bci = t_bc, f_bci = f_bc, smxi = sm, skxi = sk, cdxi = cd)
    model = ModelInfo(6)
 
    # plt.plot(np.transpose(time), d[1,:]*100)
    # plt.xlabel('Time (sec)')
    # plt.ylabel('Displacement (cm)')
    # plt.show()
    # print(d)
    # print("I am done")
    return result, model

