
import numpy as np
from scipy import linalg
import math
from data.defaults.param_manager import default
from data.defaults.param_manager import get_default_param_values
import cProfile
import io
import pstats
import pandas as pd
from mhps.postprocessor import ResultFixedXY, ModelInfo


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

def read_const_param(const_param_file):

    """
    This function reads all the constant set of parameters for the
    program which does not require to be read in iteration of the
    parametric study.
    """
    f1  = open(const_param_file, 'r')
    
    maxnst = int(f1.readline().split(':')[1])
    am = np.array([float(x.strip()) for x in f1.readline().split(':')[1:][0].strip().split(',')])
    ak = np.array([float(x.strip()) for x in f1.readline().split(':')[1:][0].strip().split(',')])
    
    f1.close()
    
    return maxnst, am, ak

def read_ss_var_param(var_param_file):

    """
    This function reads all the set of parameters for the parametric
    studies and stores it in an array. It will be a one-time allocation
    to increase speed.
    """
    
    dtype1 = np.dtype([('IJK', 'i4'), ('NST', 'i4'), ('TX1', 'f4'), ('ZETA','f4'), ('RTYTX','f4')])

    ijk, nst, tx1, zeta, rtytx = np.loadtxt \
    (var_param_file, delimiter=',', usecols=range(5),\
     skiprows=1, unpack = True, dtype=dtype1)
    
    try:
        total_param = len(ijk)
    except:
        total_param = 1


    for i in range(0, total_param):
        try:
            yield ijk[i], nst[i], tx1[i], zeta[i], rtytx[i]
        except:
            yield ijk, nst, tx1, zeta, rtytx


def get_total_variable_parameter_set(var_param_file):
 
    dtype1 = np.dtype([('IJK', 'i4'), ('NST', 'i4'), ('TX1', 'f4'), ('ZETA','f4')])

    ijk, nst, tx1, zeta = np.loadtxt(var_param_file, delimiter=',', usecols=range(4), skiprows=1, unpack = True, dtype=dtype1)
    try:
        return len(ijk)
    except:
        return 1

def superstructure_propxy(nst, tx1, rtytx, am, ak, zeta, knor):
    
    """
    Calculation of superstructure mass, stiffness and damping matrices
    in X- and Y- direction.
    """
    
    zeta = np.ones(nst)*zeta

    # Calculation of msss matrix in x-direction
    smx = np.zeros(shape=(nst+1,nst+1), dtype=np.dtype('d'), order='C')
    smx[0:nst,0:nst] = np.diag(am[0:nst])   # ------> Mass matrix in x-direction

    # Calculation of msss matrix in y-direction
    smy = np.zeros(shape=(nst+1,nst+1), dtype=np.dtype('d'), order='C')
    smy = smx                  # ------> Mass matrix in y-direction

    # Calculation of stiffness matrix in x- and y-direction
    temp = np.array([[1.,-1.],[-1.,1.]], dtype=np.dtype('d'))
    skx = np.zeros(shape=(nst+1,nst+1), dtype=np.dtype('d'), order='C')
    for i in range(nst-1):
        skx[i:i+2,i:i+2] += temp*ak[i]
    skx[-2,-2] += ak[-1]
    sky = np.zeros(shape=(nst+1,nst+1), dtype=np.dtype('d'), order='C')
    sky = skx
    del temp
    
    D = np.dot(skx[0:nst,0:nst],np.linalg.inv(smx[0:nst,0:nst]))
    e,v = np.linalg.eig(D)

    if knor == 1:
        wx1 = 2.0*np.pi/tx1
        RT = (wx1**2.0)/min(e)
        skx = RT*skx                # -------> Stiffness matrix in x-direction
        
        wy1 = wx1/rtytx
        RT = (wy1**2.0)/min(e)
        sky = RT*sky                # -------> Stiffness matrix in y-direction
        del RT
    del e, v

    # Calculation of damping matrix in x-direction
    D = np.dot(skx[0:nst,0:nst],np.linalg.inv(smx[0:nst,0:nst]))
    e,v = np.linalg.eig(D)
    idx = e.argsort()[::1]
    e = e[idx]
    #v = v[:,idx]

    T = 2.0*math.pi/np.sqrt(e)
    del e, idx, v

    W = 2.0*math.pi/T
    c_w = (np.ones((nst,nst))*W[None,:]).T
    c_w = np.power(c_w,np.arange(-1,2*nst-2,2))
    a_i = np.dot(np.linalg.inv(c_w),2.0*zeta[0:nst].reshape(nst,1))
    cdx = np.zeros(shape=(nst+1,nst+1),dtype=np.dtype('d'))
    for i in range(nst):
        cdx[0:nst,0:nst] += np.dot(np.linalg.matrix_power(D,i),smx[0:nst,0:nst])*a_i[i]  # Damping matrix in x-direction
    del W, c_w, a_i, D, T

    # Calculation of damping matrix in y-direction
    D = np.dot(sky[0:nst,0:nst],np.linalg.inv(smx[0:nst,0:nst]))
    e,v = np.linalg.eig(D)
    idx = e.argsort()[::1]
    e = e[idx]
    #v = v[:,idx]

    T = 2.0*math.pi/np.sqrt(e)
    del e, idx, v

    W = 2.0*math.pi/T
    c_w = (np.ones((nst,nst))*W[None,:]).T
    c_w = np.power(c_w,np.arange(-1,2*nst-2,2))
    a_i = np.dot(np.linalg.inv(c_w),2.0*zeta[0:nst].reshape(nst,1))
    cdy = np.zeros(shape=(nst+1,nst+1),dtype=np.dtype('d'))
    for i in range(nst):
        cdy[0:nst,0:nst] += np.dot(np.linalg.matrix_power(D,i),smx[0:nst,0:nst])*a_i[i]  # Damping matrix in y-direction
    del W, c_w, a_i, D, T

    return smx, skx, cdx, smy, sky, cdy

@profile
def fixed_simulator(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, nst, smx, skx, cdx, smy, sky, cdy, screen_on):
    
    # print(ndt)
    gamma = 0.5
    beta = 1/6

    # knx = skx + (gamma/beta/dt)*cdx + (1.0/beta/np.power(dt,2))*smx

    smx_inv = np.linalg.inv(smx)
    smx_diag = np.diag(-1.0*smx).reshape(nst,1)
    smy_inv = np.linalg.inv(smy)
    smy_diag = np.diag(-1.0*smy).reshape(nst,1)

    time = np.zeros((1, ndt), dtype=float)
    
    dx = np.zeros((nst, ndt), dtype=float)
    vx = np.zeros((nst, ndt), dtype=float)
    ax = np.zeros((nst, ndt), dtype=float)
    aax = np.zeros((nst, ndt), dtype=float)
    gx = np.zeros((1, ndt), dtype=float)
    dy = np.zeros((nst, ndt), dtype=float)
    vy = np.zeros((nst, ndt), dtype=float)
    ay = np.zeros((nst, ndt), dtype=float)
    aay = np.zeros((nst, ndt), dtype=float)
    gy = np.zeros((1, ndt), dtype=float)

    gx[0,0] = xg[0]
    gy[0,0] = yg[0]

    fx = np.zeros((ndt, nst), dtype=float)
    fy = np.zeros((ndt, nst), dtype=float)
    ek = np.zeros((ndt, 1), dtype=float)
    ed = np.zeros((ndt, 1), dtype=float)
    es = np.zeros((ndt, 1), dtype=float)
    ei = np.zeros((ndt, 1), dtype=float)
    error = np.zeros((ndt, 1), dtype=float)


    dx1 = np.ones((nst, 1), dtype=float)*0.0
    vx1 = np.ones((nst, 1), dtype=float)*0.0
    ax1 = np.ones((nst, 1), dtype=float)*0.0
    px1 = smx_diag*xg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ax1 = np.dot(smx_inv, px1 - np.dot(cdx, vx1) - np.dot(skx, dx1))
    dy1 = np.ones((nst, 1), dtype=float)*0.0
    vy1 = np.ones((nst, 1), dtype=float)*0.0
    ay1 = np.ones((nst, 1), dtype=float)*0.0
    py1 = smy_diag*yg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ay1 = np.dot(smy_inv, py1 - np.dot(cdy, vy1) - np.dot(sky, dy1))

    # I = np.ones((nst,1), dtype=float)

    dx2 = np.zeros((nst, 1), dtype=float)
    vx2 = np.zeros((nst, 1), dtype=float)
    px2 = np.zeros((nst, 1), dtype=float)
    ax2 = np.zeros((nst, 1), dtype=float)
    dy2 = np.zeros((nst, 1), dtype=float)
    vy2 = np.zeros((nst, 1), dtype=float)
    py2 = np.zeros((nst, 1), dtype=float)
    ay2 = np.zeros((nst, 1), dtype=float)
    

    na1x = (1.0/beta/np.power(dt,2))*smx + (gamma/beta/dt)*cdx
    na2x = (1.0/beta/dt)*smx + (gamma/beta - 1.0)*cdx
    na3x = (1.0/2.0/beta - 1.0)*smx + (gamma*dt/2.0/beta - dt)*cdx
    na1y = (1.0/beta/np.power(dt,2))*smy + (gamma/beta/dt)*cdy
    na2y = (1.0/beta/dt)*smy + (gamma/beta - 1.0)*cdy
    na3y = (1.0/2.0/beta - 1.0)*smy + (gamma*dt/2.0/beta - dt)*cdy

    knx = skx + na1x
    knx_inv = np.linalg.inv(knx)
    kny = sky + na1y
    kny_inv = np.linalg.inv(kny)

    index = 0

    t = 0.0
    time[0,index] = 0.0
    dx[:,index] = dx1[:,0]
    vx[:,index] = vx1[:,0]
    ax[:,index] = ax1[:,0]
    aax[:,index] = ax1[:,0] + xg[0]
    dy[:,index] = dy1[:,0]
    vy[:,index] = vy1[:,0]
    ay[:,index] = ay1[:,0]
    aay[:,index] = ay1[:,0] + yg[0]

    eki = 0.0
    edi = 0.0
    esi = 0.0
    eii = 0.0

    r = np.ones((nst, 1), dtype=float)

    for i in range(1,len(xg)):

        t += dt
        px2 = smx_diag*xg[i]
        pcx1 = px2 + np.dot(na1x, dx1) + np.dot(na2x, vx1) + np.dot(na3x, ax1)

        # if i == 1:
        #     print('I am in flags')
        #     print(knx_inv.flags)
        #     print(pcx1.flags)

        dx2 = np.dot(knx_inv, pcx1)
        vx2 = (gamma/beta/dt)*(dx2 - dx1) + (1.0 - gamma/beta)*vx1 + dt*(1.0 - gamma/2.0/beta)*ax1
        ax2 = np.dot(smx_inv, px2 - np.dot(cdx, vx2) - np.dot(skx, dx2))
        dx1, vx1, px1, ax1 = dx2, vx2, px2, ax2 
        
        py2 = smy_diag*yg[i]
        pcy1 = py2 + np.dot(na1y, dy1) + np.dot(na2y, vy1) + np.dot(na3y, ay1)
        dy2 = np.dot(kny_inv, pcy1)
        vy2 = (gamma/beta/dt)*(dy2 - dy1) + (1.0 - gamma/beta)*vy1 + dt*(1.0 - gamma/2.0/beta)*ay1
        ay2 = np.dot(smy_inv, py2 - np.dot(cdy, vy2) - np.dot(sky, dy2))

        eki = eki + 0.5*dt*(np.dot(np.dot(vx2.T, smx),ax2) + np.dot(np.dot(vx1.T, smx),ax1)) + 0.5*dt*(np.dot(np.dot(vy2.T, smy),ay2) + np.dot(np.dot(vy1.T, smy),ay1))
        edi = edi + 0.5*dt*(np.dot(np.dot(vx2.T, cdx),vx2) + np.dot(np.dot(vx1.T, cdx),vx1)) + 0.5*dt*(np.dot(np.dot(vy2.T, cdy),vy2) + np.dot(np.dot(vy1.T, cdy),vy1))
        esi = esi + 0.5*dt*(np.dot(np.dot(vx2.T, skx),dx2) + np.dot(np.dot(vx1.T, skx),dx1)) + 0.5*dt*(np.dot(np.dot(vy2.T, sky),dy2) + np.dot(np.dot(vy1.T, sky),dy1))
        # eii = eii - 0.5*dt*(np.dot(np.dot(vx2.T, smx),np.dot(r, xg[i])) + np.dot(np.dot(vx1.T, smx),np.dot(r, xg[i-1]))) - 0.5*dt*(np.dot(np.dot(vy2.T, smy),np.dot(r, yg[i])) + np.dot(np.dot(vy1.T, smy),np.dot(r, yg[i-1])))
        eii = eii - 0.5*dt*(np.dot(np.dot(vx2.T, smx),r*xg[i]) + np.dot(np.dot(vx1.T, smx),r*xg[i-1])) - 0.5*dt*(np.dot(np.dot(vy2.T, smy), r*yg[i]) + np.dot(np.dot(vy1.T, smy), r*yg[i-1]))
        

        dy1, vy1, py1, ay1 = dy2, vy2, py2, ay2


        if not i%(ndiv):
            index += 1
            time[0,index] += t
            gx[0,index] = xg[i]
            dx[:,index] = dx2[:,0]
            vx[:,index] = vx2[:,0]
            ax[:,index] = ax2[:,0]
            aax[:,index] = ax2[:,0] + xg[i]

            gy[0,index] = yg[i]
            dy[:,index] = dy2[:,0]
            vy[:,index] = vy2[:,0]
            ay[:,index] = ay2[:,0]
            aay[:,index] = ay2[:,0] + yg[i]

            for j in range(nst):
                # print(1.0*np.dot(smx_diag[0:j+1].T, aax[0:j+1,index]))
                fx[index, j] = 1.0*np.dot(smx_diag[0:j+1].T, aax[0:j+1,index])
                # print(fx[index, j])    
                fy[index, j] = 1.0*np.dot(smy_diag[0:j+1].T, aay[0:j+1,index])    
            # fx[index, 0] = 1.0*np.dot(smx_diag.T, aax[:,index])
            # fy[index, 0] = 1.0*np.dot(smy_diag.T, aay[:,index])
            ek[index, 0] = eki
            ed[index, 0] = edi
            es[index, 0] = esi
            ei[index, 0] = eii

            error[index, 0] = (edi + esi + eki - eii)/(abs(edi + esi) + eki + abs(eii))
    
    # for in 
    peakerror = max(abs(error))
    sumerror = sum(abs(error))
    
    peaktopaccx = max(abs(aax[0,:]))
    peaktopaccy = max(abs(aay[0,:]))
    peaktopdispx = max(abs(dx[0,:]))
    peaktopdispy = max(abs(dy[0,:]))

   

    if screen_on == True:
        print(" ")
        print("Simulation" + "\033[91m" + " SET%d-%d" %(ref, ijk) + "\033[0m" + ": Earthquake #: %d, Parameter #: %d" %(ref, ijk))
        print("Peak Error: % 8.6f" %(peakerror))
        print("Absolute Sum of Errors: % 8.6f" %(sumerror))
        print("Peak Top Floor Absolute Acceleration in X-Direction: % 8.6f m/s2" %(peaktopaccx))
        print("Peak Top Floor Absolute Acceleration in Y-Direction: % 8.6f m/s2" %(peaktopaccy))
        print("Peak Top Floor Relative Displacement in X-Direction: % 8.6f cm" %(peaktopdispx*100.0))
        print("Peak Top Floor Relative Displacement in Y-Direction: % 8.6f cm" %(peaktopdispy*100.0))
    
    result = ResultFixedXY(ref, ijk, time.T, gx.T, dx.T, vx.T, ax.T, aax.T, gy.T, dy.T, vy.T, ay.T, aay.T, fx, fy, ek, ed, es, ei, error, smx, skx, cdx, smy, sky, cdy)
    model = ModelInfo(nst)

    # acc = pd.DataFrame(np.hstack((time.T, gx.T)))

    # acc.to_csv("results\\Result.csv", mode='w', sep=',')

    return result, model

