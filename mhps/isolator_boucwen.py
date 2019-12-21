
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


class IsoBoucWenModel:
    def __init__(self, rmbm, tbx, zetabx, rtytxb, rzyzxb, bt, g2, a, nt, f0, g1, kg, cg, dg):
        self.rmbm = rmbm
        self.tbx = tbx
        self.zetabx = zetabx
        self.rtytxb = rtytxb
        self.rzyzxb = rzyzxb
        self.bt = bt
        self.g2 = g2
        self.a = a
        self.nt = nt
        self.f0 = f0
        self.g1 = g1
        self.kg = kg
        self.cg = cg
        self.dg = dg

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.rmbm == other.rmbm and self.tbx == other.tbx and \
                self.zetabx == other.zetabx and self.rtytxb == other.rtytxb and  self.rzyzxb == other.rzyzxb and \
                    self.bt == other.bt and self.g2 == other.g2 and self.a == other.a and self.nt == other.nt and \
                        self.f0 == other.f0 and self.g1 == other.g1 and self.kg == other.kg and self.cg == other.cg and self.dg == other.dg
        return False



def read_iso_boucwen_var_param(var_param_file):

    """
    This function reads all the set of parameters for the parametric
    studies and stores it in an array. It will be a one-time allocation
    to increase speed.
    """
    print(var_param_file)
    dtype1 = np.dtype([('IJK', 'i4'), ('RMBM', 'd'), ('TBX', 'd'), ('ZETABX','d'), ('RTYTXB','d'), ('RZXZYB', 'd'), ('BT', 'd'), ('G2', 'd'), ('A', 'd'), ('NT', 'd'), ('F0', 'd'), ('G1', 'd'), ('KG', 'd'), ('CG', 'd'), ('DG', 'd')])
    
    ijk, rmbm, tbx, zetabx, rtytxb, rzyzxb, bt, g2, a, nt, f0, g1, kg, cg, dg = np.loadtxt \
    (var_param_file, delimiter=',', usecols=range(15), skiprows=1, unpack = True, dtype=dtype1)

    try:
        total_param = len(ijk)
    except:
        total_param = 1


    for i in range(0, total_param):
        try:
            yield IsoBoucWenModel(rmbm[i], tbx[i], zetabx[i], rtytxb[i], rzyzxb[i], bt[i], g2[i], a[i], nt[i], f0[i], g1[i], kg[i], cg[i], dg[i])
        except:
            yield IsoBoucWenModel(rmbm, tbx, zetabx, rtytxb, rzyzxb, bt, g2, a, nt, f0, g1, kg, cg, dg)


def addboucwen_iso(smx, skx, cdx, smy, sky, cdy, nst, rmbm, tbx, zetabx, rtytxb, rzxzyb):
                 
    """
    Calculation of superstructure mass, stiffness and damping matrices
    in X- and Y- direction with linear isolator.
    """
    print("NST = %d, RMBM = %3.2f, TBX = %3.2f, ZETABX = %3.2f, RTYTXB = %3.2f, RZXZYB = %3.2f" %(nst, rmbm, tbx, zetabx, rtytxb, rzxzyb))

    ndof = nst + 1

    bm = smx[0,0]*rmbm
    tm = np.sum(np.diag(smx[0:nst,0:nst])) + bm
    smx[0:nst, nst] = np.diag(smx[0:nst,0:nst]) 
    smx[nst, nst] = bm
    # print(smx)

    smy[0:nst, nst] = np.diag(smy[0:nst,0:nst]) 
    smy[nst, nst] = bm

    cdabx = 2.0*zetabx*(2*math.pi/tbx)*tm
    cdx[nst, nst-1] = -1.0*(cdx[nst-1,nst-2] + cdx[nst-1,nst-1])
    cdx[nst, nst] = cdabx
    # print(cdx)

    cdaby = 2.0*zetabx*rzxzyb*(2*math.pi/(tbx*rtytxb))*tm
    cdy[nst, nst-1] = -1.0*(cdy[nst-1,nst-2] + cdy[nst-1,nst-1])
    cdy[nst, nst] = cdaby

    if tbx > 10:
        ckabx = 0.0
        ckaby = 0.0
    else:
        ckabx = math.pow(2.0*math.pi/tbx, 2.0)*tm
        ckaby = math.pow(2.0*math.pi/(tbx*rtytxb), 2.0)*tm

    skx[nst, nst-1] = -1.0*(skx[nst-1,nst-2] + skx[nst-1,nst-1])
    skx[nst, nst] = ckabx

    ckaby = math.pow(2.0*math.pi/(tbx*rtytxb), 2.0)*tm
    sky[nst, nst-1] = -1.0*(sky[nst-1,nst-2] + sky[nst-1,nst-1])
    sky[nst, nst] = ckaby



    return smx, skx, cdx, smy, sky, cdy


def wen(iso, vx, vy, zx, zy, dt, alpx, alpy, kbx, kby, fyx, fyy):
    
    beta = iso.g2
    gamma = iso.bt
    A = iso.a

    phix1 = dt*(A*vx - beta*abs(vx*zx)*zx - gamma*vx*math.pow(zx,2) - beta*abs(vy*zy)*zx - gamma*vy*zx*zy)
    phiy1 = dt*(A*vy - beta*abs(vy*zy)*zy - gamma*vy*math.pow(zy,2) - beta*abs(vx*zx)*zy - gamma*vx*zy*zx)

    zx = zx + 0.5*phix1
    zy = zy + 0.5*phiy1
    phix2 = dt*(A*vx - beta*abs(vx*zx)*zx - gamma*vx*math.pow(zx,2) - beta*abs(vy*zy)*zx - gamma*vy*zx*zy)
    phiy2 = dt*(A*vy - beta*abs(vy*zy)*zy - gamma*vy*math.pow(zy,2) - beta*abs(vx*zx)*zy - gamma*vx*zy*zx)

    zx = zx + 0.5*phix2
    zy = zy + 0.5*phiy2
    phix3 = dt*(A*vx - beta*abs(vx*zx)*zx - gamma*vx*math.pow(zx,2) - beta*abs(vy*zy)*zx - gamma*vy*zx*zy)
    phiy3 = dt*(A*vy - beta*abs(vy*zy)*zy - gamma*vy*math.pow(zy,2) - beta*abs(vx*zx)*zy - gamma*vx*zy*zx)

    zx = zx + phix3
    zy = zy + phiy3
    phix4 = dt*(A*vx - beta*abs(vx*zx)*zx - gamma*vx*math.pow(zx,2) - beta*abs(vy*zy)*zx - gamma*vy*zx*zy)
    phiy4 = dt*(A*vy - beta*abs(vy*zy)*zy - gamma*vy*math.pow(zy,2) - beta*abs(vx*zx)*zy - gamma*vx*zy*zx)

    dzx = (phix1 + 2.0*(phix2 + phix3) + phix4)/(6.0*iso.g1)
    dzy = (phiy1 + 2.0*(phiy2 + phiy3) + phiy4)/(6.0*iso.g1)
    
    dpzx = (1.0 - alpx)*fyx*dzx
    dpzy = (1.0 - alpy)*fyy*dzy

    return dpzx, dpzy, dzx, dzy

#@profile 
def simulator_boucwen(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, ndof, smx, skx, cdx, smy, sky, cdy, iso, screen_on):
    
    nit = 3
    nst = ndof - 1
    gamma = 0.5
    beta = 1/6

    smx_inv = np.linalg.inv(smx)
    smx_inv_fixed = np.linalg.inv(smx[0:nst, 0:nst])
    smx_diag = np.diag(-1.0*smx).reshape(ndof,1)
    smy_inv = np.linalg.inv(smy)
    smy_inv_fixed = np.linalg.inv(smy[0:nst, 0:nst])
    smy_diag = np.diag(-1.0*smy).reshape(ndof,1)

    tm = np.sum(np.diag(smx))
    fyx = iso.f0*tm*9.81
    fyy = iso.f0*tm*9.81

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
    vx1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0   # Check how we can tackle sign value when it's zero
    ax1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    px1 = smx_diag*xg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ax1[0:ndof, 0] = np.dot(smx_inv, px1[0:ndof, 0] - np.dot(cdx[0:ndof, 0:ndof], vx1[0:ndof, 0]) - np.dot(skx[0:ndof, 0:ndof], dx1[0:ndof, 0]))
    ax1[ndof-1, 0] = 0.0
    dy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    vy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0   # Check how we can tackle sign value when it's zero
    ay1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    py1 = smy_diag*yg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ay1[0:ndof, 0] = np.dot(smy_inv, py1[0:ndof, 0] - np.dot(cdy[0:ndof, 0:ndof], vy1[0:ndof,0]) - np.dot(sky[0:ndof, 0:ndof], dy1[0:ndof,0]))
    ay1[ndof-1, 0] = 0.0
    

    # I = np.ones((ndof,1), dtype=np.dtype('d'), order ='F')
    dx2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    vx2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0   # Check how we can tackle sign value when it's zero
    px2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    ax2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    dy2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    vy2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0   # Check how we can tackle sign value when it's zero
    py2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    ay2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    
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

    knx_fixed = skx[0:nst, 0:nst] + na1x[0:nst, 0:nst]
    knx_inv_fixed = np.linalg.inv(knx_fixed)
    kny_fixed = sky[0:nst, 0:nst] + na1y[0:nst, 0:nst]
    kny_inv_fixed = np.linalg.inv(kny_fixed)

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

    alpx = (iso.g1*skx[ndof-1, ndof-1])/(iso.f0*tm*9.81)
    alpy = (iso.g1*sky[ndof-1, ndof-1])/(iso.f0*tm*9.81)

    # print('alpha')
    # print(alpx)

    r = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')
    pzx = 0.0
    pzy = 0.0
    zx = 0.0
    zy = 0.0
    dzx = 0.0
    dzy = 0.0

    dpx = 0.0
    dpy = 0.0

    for i in range(1,len(xg)):

        t += dt

        px2 = smx_diag*xg[i]
        py2 = smy_diag*yg[i]
        dpx = (px2 - px1)
        dpy = (py2 - py1)

        dpzx = 0.0
        dpzy = 0.0

        for i2 in range(nit):
            pcx1 = dpx + np.dot(na2x, vx1) + np.dot(na3x, ax1)
            pcx1[ndof-1,0] = pcx1[ndof-1,0] - dpzx
            ddx = np.dot(knx_inv, pcx1)
            dx2 = dx1 + ddx
            dvx = (gamma/beta/dt)*ddx - gamma/beta*vx1 + dt*(1.0 - gamma/2.0/beta)*ax1
            vx2 = vx1 + dvx

            pcy1 = dpy + np.dot(na2y, vy1) + np.dot(na3y, ay1)
            pcy1[ndof-1,0] = pcy1[ndof-1,0] - dpzy
            ddy = np.dot(kny_inv, pcy1)
            dy2 = dy1 + ddy
            dvy = (gamma/beta/dt)*ddy - gamma/beta*vy1 + dt*(1.0 - gamma/2.0/beta)*ay1
            vy2 = vy1 + dvy
                
            dpzx, dpzy, dzx, dzy = wen(iso, vx2[ndof-1,0], vy2[ndof-1,0], zx, zy, dt, alpx, alpy, skx[ndof-1, ndof-1], sky[ndof-1, ndof-1], fyx, fyy)

        zx = zx + dzx
        zy = zy + dzy    
        pzx = pzx + dpzx
        pzy = pzy + dpzy
        fabx = skx[ndof-1, ndof-1]*dx2[ndof-1,0] + pzx
        faby = sky[ndof-1, ndof-1]*dy2[ndof-1,0] + pzy

        epx = px2 - np.dot(cdx, vx2) - np.dot(skx, dx2)
        epy = py2 - np.dot(cdy, vy2) - np.dot(sky, dy2)
        epx[ndof-1,0] = epx[ndof-1,0] - pzx
        epy[ndof-1,0] = epy[ndof-1,0] - pzy
        ax2 = np.dot(smx_inv, epx)
        ay2 = np.dot(smy_inv, epy)

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

            for j in range(nst):
                fx[index, j] = 1.0*np.dot(smx_diag[0:j+1].T, aax[0:j+1,index])   
                fy[index, j] = 1.0*np.dot(smy_diag[0:j+1].T, aay[0:j+1,index])    
            fx[index, ndof-1] = fabx
            fy[index, ndof-1] = faby
            ek[index, 0] = eki
            ed[index, 0] = edi
            es[index, 0] = esi
            ei[index, 0] = eii

            error[index, 0] = (edi + esi + eki - eii)/(abs(edi + esi) + eki + abs(eii))
    
    peakerror = max(abs(error))
    sumerror = sum(abs(error))
    peaktopaccx = max(abs(aax[0,:]))
    peaktopaccy = max(abs(aay[0,:]))
    peaktopdispx = max(abs(dx[0,:]))
    peaktopdispy = max(abs(dy[0,:]))
    peakbasedispx = max(abs(dx[ndof-1,:]))
    peakbasedispy = max(abs(dy[ndof-1,:]))
    residualdispx = abs(dx[ndof-1,-1])
    residualdispy = abs(dy[ndof-1,-1])

    if screen_on == True:
        print(" ")
        print("Simulation" + "\033[91m" + " SET%d-%d" %(ref, ijk) + "\033[0m" + ": Earthquake #: %d, Parameter #: %d" %(ref, ijk))
        print("Peak Error: % 8.6f" %(peakerror))
        print("Absolute Sum of Errors: % 8.6f" %(sumerror))
        print("Peak Top Floor Absolute Acceleration in X-Direction: % 8.6f m/s2" %(peaktopaccx))
        print("Peak Top Floor Absolute Acceleration in Y-Direction: % 8.6f m/s2" %(peaktopaccy))
        print("Peak Top Floor Relative Displacement in X-Direction: % 8.6f cm" %(peaktopdispx*100.0))
        print("Peak Top Floor Relative Displacement in Y-Direction: % 8.6f cm" %(peaktopdispy*100.0))
        print("Peak Isolator Displacement in X-Direction: % 8.6f cm" %(peakbasedispx*100.0))
        print("Peak Isolator Displacement in Y-Direction: % 8.6f cm" %(peakbasedispy*100.0))
        print("Isolator Residual Displacement in X-Direction: % 8.6f cm" %(residualdispx*100.0))
        print("Isolator Residual Displacement in Y-Direction: % 8.6f cm" %(residualdispy*100.0))
    
    result = ResultFixedXY(ref, ijk, time.T, gx.T, dx.T, vx.T, ax.T, aax.T, gy.T, dy.T, vy.T, ay.T, aay.T, fx, fy, ek, ed, es, ei, error, smx, skx, cdx, smy, sky, cdy)
    model = ModelInfo(ndof)

    return result, model
