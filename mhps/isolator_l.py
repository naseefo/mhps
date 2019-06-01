
import numpy as np
import math
from data.defaults.param_manager import default
from data.defaults.param_manager import get_default_param_values
import cProfile
import io
import pstats
import pandas as pd
from mhps.postprocessor import ResultFixedXY, ModelInfo
import math


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
    am = np.array([np.double(x.strip()) for x in f1.readline().split(':')[1:][0].strip().split(',')])
    ak = np.array([np.double(x.strip()) for x in f1.readline().split(':')[1:][0].strip().split(',')])
    
    f1.close()
    
    return maxnst, am, ak

def read_iso_l_var_param(var_param_file):

    """
    This function reads all the set of parameters for the parametric
    studies and stores it in an array. It will be a one-time allocation
    to increase speed.
    """
    
    dtype1 = np.dtype([('IJK', 'i4'), ('RMBM', 'd'), ('TBX', 'd'), ('ZETABX','d'), ('RTYTXB','d'), ('RZXZYB', 'd')])

    ijk, rmbm, tbx, zetabx, rtytxb, rzyzxb = np.loadtxt \
    (var_param_file, delimiter=',', usecols=range(6),\
     skiprows=1, unpack = True, dtype=dtype1)
    
    try:
        total_param = len(ijk)
    except:
        total_param = 1


    for i in range(0, total_param):
        try:
            yield ijk[i], rmbm[i], tbx[i], zetabx[i], rtytxb[i], rzyzxb[i]
        except:
            yield ijk, rmbm, tbx, zetabx, rtytxb, rzyzxb


def get_total_variable_parameter_set(var_param_file):
 
    dtype1 = np.dtype([('IJK', 'i4'), ('NST', 'i4'), ('TX1', 'd'), ('ZETA','d')])

    ijk, nst, tx1, zeta = np.loadtxt(var_param_file, delimiter=',', usecols=range(4), skiprows=1, unpack = True, dtype=dtype1)
    try:
        return len(ijk)
    except:
        return 1

def addlinear_iso(smx, skx, cdx, smy, sky, cdy, nst, rmbm, tbx, zetabx, rtytxb, rzyzxb):
                 
    """
    Calculation of superstructure mass, stiffness and damping matrices
    in X- and Y- direction with linear isolator.
    """
    print("NST = %d, RMBM = %3.2f, TBX = %3.2f, ZETABX = %3.2f, RTYTXB = %3.2f, RZXZYB = %3.2f" %(nst, rmbm, tbx, zetabx, rtytxb, rzyzxb))
    bm = smx[0,0]*rmbm
    tm = np.sum(np.diag(smx[0:nst,0:nst])) + bm
    smx[0:nst, nst] = np.diag(smx[0:nst,0:nst]) 
    smx[nst, nst] = bm
    # print(smx)

    smy[0:nst, nst] = np.diag(smy[0:nst,0:nst]) 
    smy[nst, nst] = bm

    if zetabx == 0:
        cdabx = 0.0
        cdaby = 0.0
    else:
        cdabx = 2.0*zetabx*(2*math.pi/tbx)*tm
        cdaby = 2.0*zetabx*rzyzxb*(2*math.pi/(tbx*rtytxb))*tm

    cdx[nst, nst-1] = -1.0*(cdx[nst-1,nst-2] + cdx[nst-1,nst-1])
    cdx[nst, nst] = cdabx
    cdy[nst, nst-1] = -1.0*(cdy[nst-1,nst-2] + cdy[nst-1,nst-1])
    cdy[nst, nst] = cdaby
    # print(cdx)

    

    if tbx > 49:
        ckabx = 0.0
        ckaby = 0.0
    else:
        ckabx = math.pow(2.0*math.pi/tbx, 2.0)*tm
        ckaby = math.pow(2.0*math.pi/(tbx*rtytxb), 2.0)*tm

    skx[nst, nst-1] = -1.0*(skx[nst-1,nst-2] + skx[nst-1,nst-1])
    skx[nst, nst] = ckabx
    sky[nst, nst-1] = -1.0*(sky[nst-1,nst-2] + sky[nst-1,nst-1])
    sky[nst, nst] = ckaby



    return smx, skx, cdx, smy, sky, cdy

#@profile
def simulator_linear1(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, ndof, smx, skx, cdx, smy, sky, cdy):
    

    gamma = 0.5
    beta = 1/6

    # knx = skx + (gamma/beta/dt)*cdx + (1.0/beta/np.power(dt,2))*smx

    smx_inv = np.linalg.inv(smx)
    smx_diag = np.diag(-1.0*smx).reshape(ndof,1)
    smy_inv = np.linalg.inv(smy)
    smy_diag = np.diag(-1.0*smy).reshape(ndof,1)

    time = np.zeros((1, ndt), dtype=np.dtype('d'), order ='F')
    
    dx = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    vx = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    ax = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    aax = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    gx = np.zeros((1, ndt), dtype=np.dtype('d'), order ='F')
    dy = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    vy = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    ay = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    aay = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    gy = np.zeros((1, ndt), dtype=np.dtype('d'), order ='F')

    gx[0,0] = xg[0]
    gy[0,0] = yg[0]

    fx = np.zeros((ndt, ndof), dtype=np.dtype('d'), order ='F')
    fy = np.zeros((ndt, ndof), dtype=np.dtype('d'), order ='F')
    ek = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F')
    ed = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F')
    es = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F')
    ei = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F')
    error = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F')


    dx1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    vx1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    ax1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    px1 = smx_diag*xg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ax1 = np.dot(smx_inv, px1 - np.dot(cdx, vx1) - np.dot(skx, dx1))
    dy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    vy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    ay1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    py1 = smy_diag*yg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ay1 = np.dot(smy_inv, py1 - np.dot(cdy, vy1) - np.dot(sky, dy1))

    # I = np.ones((ndof,1), dtype=np.dtype('d'), order ='F')

    dx2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    vx2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    px2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    ax2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    dy2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    vy2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    py2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    ay2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    

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

    smxi = np.diag(np.diag(smx))   # The right hand side matrix for mass has changes compared to the left hand side. 
    smyi = np.diag(np.diag(smy))
    eki = 0.0
    edi = 0.0
    esi = 0.0
    eii = 0.0

    r = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')
    for i in range(1,len(xg)):

        t += dt
        px2 = smx_diag*xg[i]
        pcx1 = px2 + np.dot(na1x, dx1) + np.dot(na2x, vx1) + np.dot(na3x, ax1)
        dx2 = np.dot(knx_inv, pcx1)
        vx2 = (gamma/beta/dt)*(dx2 - dx1) + (1.0 - gamma/beta)*vx1 + dt*(1.0 - gamma/2.0/beta)*ax1
        ax2 = np.dot(smx_inv, px2 - np.dot(cdx, vx2) - np.dot(skx, dx2))
        
        py2 = smy_diag*yg[i]
        pcy1 = py2 + np.dot(na1y, dy1) + np.dot(na2y, vy1) + np.dot(na3y, ay1)
        dy2 = np.dot(kny_inv, pcy1)
        vy2 = (gamma/beta/dt)*(dy2 - dy1) + (1.0 - gamma/beta)*vy1 + dt*(1.0 - gamma/2.0/beta)*ay1
        ay2 = np.dot(smy_inv, py2 - np.dot(cdy, vy2) - np.dot(sky, dy2))

        eki = eki + 0.5*dt*(np.dot(np.dot(vx2.T, smx),ax2) + np.dot(np.dot(vx1.T, smx),ax1)) + 0.5*dt*(np.dot(np.dot(vy2.T, smy),ay2) + np.dot(np.dot(vy1.T, smy),ay1))
        edi = edi + 0.5*dt*(np.dot(np.dot(vx2.T, cdx),vx2) + np.dot(np.dot(vx1.T, cdx),vx1)) + 0.5*dt*(np.dot(np.dot(vy2.T, cdy),vy2) + np.dot(np.dot(vy1.T, cdy),vy1))
        esi = esi + 0.5*dt*(np.dot(np.dot(vx2.T, skx),dx2) + np.dot(np.dot(vx1.T, skx),dx1)) + 0.5*dt*(np.dot(np.dot(vy2.T, sky),dy2) + np.dot(np.dot(vy1.T, sky),dy1))
        eii = eii - 0.5*dt*(np.dot(np.dot(vx2.T, smxi),np.dot(r, xg[i])) + np.dot(np.dot(vx1.T, smxi),np.dot(r, xg[i-1]))) - 0.5*dt*(np.dot(np.dot(vy2.T, smyi),np.dot(r, yg[i])) + np.dot(np.dot(vy1.T, smyi),np.dot(r, yg[i-1])))
        
        dx1, vx1, px1, ax1 = dx2, vx2, px2, ax2 
        dy1, vy1, py1, ay1 = dy2, vy2, py2, ay2 
        
        # RESULT STORAGE SECTION
        if not i%(ndiv):
            index += 1
            time[0,index] += t
            gx[0,index] = xg[i]
            dx[:,index] = dx2[:,0]
            vx[:,index] = vx2[:,0]
            ax[:,index] = ax2[:,0]
            aax[0:ndof,index] = ax2[:,0] + ax2[ndof-1,0] + xg[i]
            aax[ndof-1,index] = ax2[ndof-1,0] + xg[i]

            gy[0,index] = yg[i]
            dy[:,index] = dy2[:,0]
            vy[:,index] = vy2[:,0]
            ay[:,index] = ay2[:,0]
            aay[:,index] = ay2[:,0] + ay2[ndof-1,0] + yg[i]
            aay[ndof-1,index] = ay2[ndof-1,0] + yg[i]

            for j in range(ndof):
                fx[index, j] = 1.0*np.dot(smx_diag[0:j+1].T, aax[0:j+1,index])   
                fy[index, j] = 1.0*np.dot(smy_diag[0:j+1].T, aay[0:j+1,index])    
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

    print(" ")
    print("Simulation" + "\033[91m" + " SET%d-%d" %(ref, ijk) + "\033[0m" + ": Earthquake #: %d, Parameter #: %d" %(ref, ijk))
    print("Peak Error: % 8.6f" %(peakerror))
    print("Absolute Sum of Errors: % 8.6f" %(sumerror))
    print("Peak Top Floor Absolute Acceleration in X-Direction: % 8.6f m/s2" %(peaktopaccx))
    print("Peak Top Floor Absolute Acceleration in Y-Direction: % 8.6f m/s2" %(peaktopaccy))
    print("Peak Top Floor Relative Displacement in X-Direction: % 8.6f cm" %(peaktopdispx*100.0))
    print("Peak Top Floor Relative Displacement in Y-Direction: % 8.6f cm" %(peaktopdispy*100.0))
    
    result = ResultFixedXY(ref, ijk, time.T, gx.T, dx.T, vx.T, ax.T, aax.T, gy.T, dy.T, vy.T, ay.T, aay.T, fx, fy, ek, ed, es, ei, error, smx, skx, cdx, smy, sky, cdy)
    model = ModelInfo(ndof)

    # acc = pd.DataFrame(np.hstack((time.T, gx.T)))

    # acc.to_csv("results\\Result.csv", mode='w', sep=',')

    return result, model


def simulator_linear(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, ndof, smx, skx, cdx, smy, sky, cdy, screen_on):
    
    gamma = 0.5
    beta = 1/6

    # knx = skx + (gamma/beta/dt)*cdx + (1.0/beta/np.power(dt,2))*smx

    smx_inv = np.linalg.inv(smx)
    smx_diag = np.diag(-1.0*smx).reshape(ndof,1)
    smy_inv = np.linalg.inv(smy)
    smy_diag = np.diag(-1.0*smy).reshape(ndof,1)

    time = np.zeros((1, ndt), dtype=np.dtype('d'), order ='F')
    
    dx = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    vx = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    ax = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    aax = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    gx = np.zeros((1, ndt), dtype=np.dtype('d'), order ='F')
    ddx = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    dvx = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')

    dy = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    vy = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    ay = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    aay = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F')
    gy = np.zeros((1, ndt), dtype=np.dtype('d'), order ='F')
    ddy = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    dvy = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')

    gx[0,0] = xg[0]
    gy[0,0] = yg[0]

    fx = np.zeros((ndt, ndof), dtype=np.dtype('d'), order ='F')
    fy = np.zeros((ndt, ndof), dtype=np.dtype('d'), order ='F')
    ek = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F')
    ed = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F')
    es = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F')
    ei = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F')
    error = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F')

    dx1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    vx1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    ax1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    px1 = smx_diag*xg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ax1 = np.dot(smx_inv, px1 - np.dot(cdx, vx1) - np.dot(skx, dx1))
    dy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    vy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    ay1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    py1 = smy_diag*yg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ay1 = np.dot(smy_inv, py1 - np.dot(cdy, vy1) - np.dot(sky, dy1))

    # I = np.ones((ndof,1), dtype=np.dtype('d'), order ='F')
    dx2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    vx2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    px2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    ax2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    dy2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    vy2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    py2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    ay2 = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F')
    
    na1x = (1.0/beta/np.power(dt,2))*smx + (gamma/beta/dt)*cdx
    na2x = (1.0/beta/dt)*smx + gamma/beta*cdx
    na3x = 1.0/2.0/beta*smx + (gamma*dt/2.0/beta - dt)*cdx
    na1y = (1.0/beta/np.power(dt,2))*smy + (gamma/beta/dt)*cdy
    na2y = (1.0/beta/dt)*smy + gamma/beta*cdy
    na3y = 1.0/2.0/beta*smy + (gamma*dt/2.0/beta - dt)*cdy

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

    smxi = np.diag(np.diag(smx))   # The right hand side matrix for mass has changes compared to the left hand side. 
    smyi = np.diag(np.diag(smy))
    eki = 0.0
    edi = 0.0
    esi = 0.0
    eii = 0.0

    r = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')
    dpx = 0.0
    dpy = 0.0
    for i in range(1,len(xg)):

        t += dt
        px2 = smx_diag*xg[i]
        dpx = (px2-px1)
        pcx1 = dpx + np.dot(na2x, vx1) + np.dot(na3x, ax1)
        ddx = np.dot(knx_inv, pcx1)
        # print(ddx)
        dx2 = dx1 + ddx
        dvx = (gamma/beta/dt)*ddx - gamma/beta*vx1 + dt*(1.0 - gamma/2.0/beta)*ax1
        vx2 = vx1 + dvx
        ax2 = np.dot(smx_inv, px2 - np.dot(cdx, vx2) - np.dot(skx, dx2))
        
        py2 = smy_diag*yg[i]
        dpy = (py2-py1)
        pcy1 = dpy + np.dot(na2y, vy1) + np.dot(na3y, ay1)
        ddy = np.dot(kny_inv, pcy1)
        dy2 = dy1 + ddy
        dvy = (gamma/beta/dt)*ddy - gamma/beta*vy1 + dt*(1.0 - gamma/2.0/beta)*ay1
        vy2 = vy1 + dvy
        ay2 = np.dot(smy_inv, py2 - np.dot(cdy, vy2) - np.dot(sky, dy2))

        eki = eki + 0.5*dt*(np.dot(np.dot(vx2.T, smx),ax2) + np.dot(np.dot(vx1.T, smx),ax1)) + 0.5*dt*(np.dot(np.dot(vy2.T, smy),ay2) + np.dot(np.dot(vy1.T, smy),ay1))
        edi = edi + 0.5*dt*(np.dot(np.dot(vx2.T, cdx),vx2) + np.dot(np.dot(vx1.T, cdx),vx1)) + 0.5*dt*(np.dot(np.dot(vy2.T, cdy),vy2) + np.dot(np.dot(vy1.T, cdy),vy1))
        esi = esi + 0.5*dt*(np.dot(np.dot(vx2.T, skx),dx2) + np.dot(np.dot(vx1.T, skx),dx1)) + 0.5*dt*(np.dot(np.dot(vy2.T, sky),dy2) + np.dot(np.dot(vy1.T, sky),dy1))
        eii = eii - 0.5*dt*(np.dot(np.dot(vx2.T, smxi),np.dot(r, xg[i])) + np.dot(np.dot(vx1.T, smxi),np.dot(r, xg[i-1]))) - 0.5*dt*(np.dot(np.dot(vy2.T, smyi),np.dot(r, yg[i])) + np.dot(np.dot(vy1.T, smyi),np.dot(r, yg[i-1])))
        
        dx1, vx1, px1, ax1 = dx2, vx2, px2, ax2 
        dy1, vy1, py1, ay1 = dy2, vy2, py2, ay2 
         
        # RESULT STORAGE SECTION
        if not i%(ndiv):
            index += 1
            time[0,index] += t
            gx[0,index] = xg[i]
            dx[:,index] = dx2[:,0]
            vx[:,index] = vx2[:,0]
            ax[:,index] = ax2[:,0]
            aax[0:ndof,index] = ax2[:,0] + ax2[ndof-1,0] + xg[i]
            aax[ndof-1,index] = ax2[ndof-1,0] + xg[i]

            gy[0,index] = yg[i]
            dy[:,index] = dy2[:,0]
            vy[:,index] = vy2[:,0]
            ay[:,index] = ay2[:,0]
            aay[:,index] = ay2[:,0] + ay2[ndof-1,0] + yg[i]
            aay[ndof-1,index] = ay2[ndof-1,0] + yg[i]

            for j in range(ndof):
                fx[index, j] = 1.0*np.dot(smx_diag[0:j+1].T, aax[0:j+1,index])   
                fy[index, j] = 1.0*np.dot(smy_diag[0:j+1].T, aay[0:j+1,index])    
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
    model = ModelInfo(ndof)

    # acc = pd.DataFrame(np.hstack((time.T, gx.T)))

    # acc.to_csv("results\\Result.csv", mode='w', sep=',')

    return result, model
