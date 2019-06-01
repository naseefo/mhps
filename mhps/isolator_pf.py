
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


class IsoPFModel:
    def __init__(self, rmbm, tbx, zetabx, rtytxb, rzxzyb, typevf, mu0, alpha0, alpha1, nu, umax, kg, cg, dg):
        self.rmbm = rmbm
        self.tbx = tbx
        self.zetabx = zetabx
        self.rtytxb = rtytxb
        self.rzxzyb = rzxzyb
        self.typevf = typevf
        self.mu0 = mu0
        self.alpha0 = alpha0
        self.alpha1 = alpha1
        self.nu = nu
        self.umax = umax
        self.kg = kg
        self.cg = cg
        self.dg = dg

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.rmbm == other.rmbm and self.tbx == other.tbx and \
                self.zetabx == other.zetabx and self.rtytxb == other.rtytxb and  self.rzxzyb == other.rzxzyb and \
                    self.typevf == other.typevf and self.mu0 == other.mu0 and self.alpha0 == other.alpha0 and self.alpha1 == other.alpha1 and \
                        self.nu == other.nu and self.umax == other.umax and self.kg == other.kg and self.cg == other.cg and self.dg == other.dg
        return False



def read_iso_pf_var_param(var_param_file):

    """
    This function reads all the set of parameters for the parametric
    studies and stores it in an array. It will be a one-time allocation
    to increase speed.
    """
    print(var_param_file)
    dtype1 = np.dtype([('IJK', 'i4'), ('RMBM', 'd'), ('TBX', 'd'), ('ZETABX','d'), ('RTYTXB','d'), ('RZXZYB', 'd'), ('TYPEVF', 'i4'), ('MU0', 'd'), ('ALPHA0', 'd'), ('ALPHA1', 'd'), ('NU', 'd'), ('UMAX', 'd'), ('KG', 'd'), ('CG', 'd'), ('DG', 'd')])
    print('I have reached here')
    ijk, rmbm, tbx, zetabx, rtytxb, rzxzyb, typevf, mu0, alpha0, alpha1, nu, umax, kg, cg, dg = np.loadtxt \
    (var_param_file, delimiter=',', usecols=range(15), skiprows=1, unpack = True, dtype=dtype1)

    try:
        total_param = len(ijk)
    except:
        total_param = 1


    for i in range(0, total_param):
        try:
            yield IsoPFModel(rmbm[i], tbx[i], zetabx[i], rtytxb[i], rzxzyb[i], typevf[i], mu0[i], alpha0[i], alpha1[i], nu[i], umax[i], kg[i], cg[i], dg[i])
        except:
            yield IsoPFModel(rmbm, tbx, zetabx, rtytxb, rzxzyb, typevf, mu0, alpha0, alpha1, nu, umax, kg, cg, dg)


def stat1(fabx, faby, qx, qy):
    ratio = math.sqrt(math.pow(fabx/qx, 2) + math.pow(faby/qy, 2))
    if (ratio - 1.0) < 0:
        sliding_state = False
    else:
        fabx = fabx/ratio
        faby = faby/ratio
        sliding_state = True

    return sliding_state, fabx, faby

def stat0(ddxb, ddyb, fabx, faby):
    wd = fabx*ddxb + faby*ddyb
    if wd >= 0:
        sliding_state = True
    else:
        sliding_state = False
    return sliding_state

def mu_val(iso, ub):
    L0 = iso.alpha0*iso.mu0
    L1 = (iso.alpha1-iso.alpha0)*(iso.mu0/iso.umax)*abs(ub)
    L2 = (iso.alpha1-iso.alpha0)*(iso.mu0/math.pow(iso.umax,2))*math.pow(ub,2)
    if iso.typevf == 1:
        mu = L0 + (1-iso.nu)*L1 + iso.nu*L2
    elif iso.typevf == 2:
        mu = L0 - (1-iso.nu)*L1 - iso.nu*(L2 - 2*L1)
    else:
        mu = iso.mu0

    if abs(ub) > iso.umax and abs(ub) < iso.dg:
        mu = iso.alpha1*iso.mu0

    if abs(ub) > iso.dg and iso.typevf != 0:
        mu = 100
    return mu

#@profile 
def simulator_pf(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, ndof, smx, skx, cdx, smy, sky, cdy, iso, screen_on):
    
    nit = 5
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
    vx1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001
    ax1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    px1 = smx_diag*xg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ax1[0:nst, 0] = np.dot(smx_inv_fixed, px1[0:nst, 0] - np.dot(cdx[0:nst, 0:nst], vx1[0:nst, 0]) - np.dot(skx[0:nst, 0:nst], dx1[0:nst, 0]))
    dy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    vy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001
    ay1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    py1 = smy_diag*yg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ay1[0:nst, 0] = np.dot(smy_inv_fixed, py1[0:nst, 0] - np.dot(cdy[0:nst, 0:nst], vy1[0:nst,0]) - np.dot(sky[0:nst, 0:nst], dy1[0:nst,0]))

    # I = np.ones((ndof,1), dtype=np.dtype('d'), order ='F')
    dx2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    vx2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001
    px2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    ax2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    dy2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0
    vy2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001
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

    r = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')
    dpx = 0.0
    dpy = 0.0
    sliding_state = False

    qx = 9.81*tm*iso.mu0
    qy = 9.81*tm*iso.mu0

    p_index = 1000
    for i in range(1, len(xg)):
        

        t += dt

        px2 = smx_diag*xg[i]
        py2 = smy_diag*yg[i]
        dpx = (px2 - px1)
        dpy = (py2 - py1)

        if sliding_state == False:
                
            pcx1 = dpx[0:nst] + np.dot(na2x[0:nst, 0:nst], vx1[0:nst]) + np.dot(na3x[0:nst, 0:nst], ax1[0:nst])
            ddx[0:nst] = np.dot(knx_inv_fixed, pcx1)
            dx2[0:nst, 0] = dx1[0:nst, 0] + ddx[0:nst, 0]
            dvx[0:nst, 0] = (gamma/beta/dt)*ddx[0:nst, 0] - gamma/beta*vx1[0:nst, 0] + dt*(1.0 - gamma/2.0/beta)*ax1[0:nst, 0]
            vx2[0:nst, 0] = vx1[0:nst, 0] + dvx[0:nst, 0]
        
            pcy1 = dpy[0:nst] + np.dot(na2y[0:nst, 0:nst], vy1[0:nst]) + np.dot(na3y[0:nst, 0:nst], ay1[0:nst])
            ddy[0:nst] = np.dot(kny_inv_fixed, pcy1)
            dy2[0:nst, 0] = dy1[0:nst, 0] + ddy[0:nst, 0]
            dvy[0:nst, 0] = (gamma/beta/dt)*ddy[0:nst, 0] - gamma/beta*vy1[0:nst, 0] + dt*(1.0 - gamma/2.0/beta)*ay1[0:nst, 0]
            vy2[0:nst, 0] = vy1[0:nst, 0] + dvy[0:nst, 0]

            fabx = (-1.0*cdx[nst, nst-1])*vx2[nst-1,0] + (-1.0*skx[nst, nst-1])*dx2[nst-1,0] - cdx[ndof-1, ndof-1]*vx2[ndof-1,0] - skx[ndof-1, ndof-1]*dx2[ndof-1,0] + px2[ndof-1]
            

            faby = (-1.0*cdy[nst, nst-1])*vy2[nst-1,0] + (-1.0*sky[nst, nst-1])*dy2[nst-1,0] - cdy[ndof-1, ndof-1]*vy2[ndof-1,0] - sky[ndof-1, ndof-1]*dy2[ndof-1,0] + py2[ndof-1]
        
            if i == p_index:
                print(i, sliding_state)
                print(fabx)
                print(faby)

            ub = math.sqrt(math.pow(dx2[ndof-1,0], 2) + math.pow(dy2[ndof-1,0], 2))
            mu = mu_val(iso, ub)
            
            qx = 9.81*tm*mu
            qy = 9.81*tm*mu

            sliding_state, fabx, faby = stat1(fabx, faby, qx, qy)

            if sliding_state == False:
                ax2[0:nst, 0] = np.dot(smx_inv_fixed, px2[0:nst, 0] - np.dot(cdx[0:nst, 0:nst], vx2[0:nst, 0]) - np.dot(skx[0:nst, 0:nst], dx2[0:nst, 0]))
                ay2[0:nst, 0] = np.dot(smy_inv_fixed, py2[0:nst, 0] - np.dot(cdy[0:nst, 0:nst], vy2[0:nst, 0]) - np.dot(sky[0:nst, 0:nst], dy2[0:nst, 0]))
            else:
                epx = px2 - np.dot(cdx, vx2) - np.dot(skx, dx2)
                epy = py2 - np.dot(cdy, vy2) - np.dot(sky, dy2)

                epx[ndof-1,0] = epx[ndof-1,0] - fabx
                epy[ndof-1,0] = epy[ndof-1,0] - faby

                ax2 = np.dot(smx_inv, epx)
                ay2 = np.dot(smy_inv, epy)
            
            dx1, vx1, px1, ax1 = dx2, vx2, px2, ax2 
            dy1, vy1, py1, ay1 = dy2, vy2, py2, ay2 

        else:
            dfabx = 0.0
            dfaby = 0.0

            for ii in range (nit):
                  
                pcx1 = dpx + np.dot(na2x, vx1) + np.dot(na3x, ax1)
                pcx1[ndof-1,0] = pcx1[ndof-1,0] - dfabx
                ddx = np.dot(knx_inv, pcx1)
                dx2 = dx1 + ddx
                dvx = (gamma/beta/dt)*ddx - gamma/beta*vx1 + dt*(1.0 - gamma/2.0/beta)*ax1
                vx2 = vx1 + dvx

                pcy1 = dpy + np.dot(na2y, vy1) + np.dot(na3y, ay1)
                pcy1[ndof-1,0] = pcy1[ndof-1,0] - dfaby
                ddy = np.dot(kny_inv, pcy1)
                dy2 = dy1 + ddy
                dvy = (gamma/beta/dt)*ddy - gamma/beta*vy1 + dt*(1.0 - gamma/2.0/beta)*ay1
                vy2 = vy1 + dvy

                rvl = math.sqrt(math.pow(vx2[ndof-1, 0], 2) + math.pow(vy2[ndof-1, 0], 2))

                if rvl > 1.0e-5:
                    ub = math.sqrt(math.pow(dx2[ndof-1,0], 2) + math.pow(dy2[ndof-1,0], 2))
                    mu = mu_val(iso, ub)
                    qx = 9.81*tm*mu
                    qy = 9.81*tm*mu
                    dfabx = qx*vx2[ndof-1,0]/rvl - fabx
                    dfaby = qy*vy2[ndof-1,0]/rvl - faby
                
            if nit > 1:
                fabx = fabx + dfabx
                faby = faby + dfaby
            
            sliding_state = stat0(ddx[ndof-1], ddy[ndof-1], fabx, faby)

            if sliding_state == False:
                vx2[ndof-1, 0] = 1e-10
                ax2[ndof-1, 0] = 0.0
                vy2[ndof-1, 0] = 1e-10
                ay2[ndof-1, 0] = 0.0

                ax2[0:nst, 0] = np.dot(smx_inv_fixed, px2[0:nst, 0] - np.dot(cdx[0:nst, 0:nst], vx2[0:nst, 0]) - np.dot(skx[0:nst, 0:nst], dx2[0:nst, 0]))
                ay2[0:nst, 0] = np.dot(smy_inv_fixed, py2[0:nst, 0] - np.dot(cdy[0:nst, 0:nst], vy2[0:nst, 0]) - np.dot(sky[0:nst, 0:nst], dy2[0:nst, 0]))
            else:
                epx = px2 - np.dot(cdx, vx2) - np.dot(skx, dx2)
                epy = py2 - np.dot(cdy, vy2) - np.dot(sky, dy2)

                epx[ndof-1,0] = epx[ndof-1,0] - fabx
                epy[ndof-1,0] = epy[ndof-1,0] - faby

                ax2 = np.dot(smx_inv, epx)
                ay2 = np.dot(smy_inv, epy)
            
            dx1, vx1, px1, ax1 = dx2, vx2, px2, ax2 
            dy1, vy1, py1, ay1 = dy2, vy2, py2, ay2
        

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
    
            # index += 1
            # time[0,index] += t
            # gx[0,index] = xg[i]
            # dx[:,index] = dx2[:,0]
            # vx[:,index] = vx2[:,0]
            # ax[:,index] = ax2[:,0]
            # aax[0:ndof,index] = ax2[:,0] + ax2[ndof-1,0] + xg[i]
            # aax[ndof-1,index] = ax2[ndof-1,0] + xg[i]

            # gy[0,index] = yg[i]
            # dy[:,index] = dy2[:,0]
            # vy[:,index] = vy2[:,0]
            # ay[:,index] = ay2[:,0]
            # aay[:,index] = ay2[:,0] + ay2[ndof-1,0] + yg[i]
            # aay[ndof-1,index] = ay2[ndof-1,0] + yg[i]

            # for j in range(ndof):
            #     fx[index, j] = 1.0*np.dot(smx_diag[0:j+1].T, aax[0:j+1,index])   
            #     fy[index, j] = 1.0*np.dot(smy_diag[0:j+1].T, aay[0:j+1,index])
                    
            # ek[index, 0] = eki
            # ed[index, 0] = edi
            # es[index, 0] = esi
            # ei[index, 0] = eii

            # error[index, 0] = (edi + esi + eki - eii)/(abs(edi + esi) + eki + abs(eii))
    
   
    peakerror = max(abs(error))
    sumerror = sum(abs(error))
    peaktopaccx = max(abs(aax[0,:]))
    peaktopaccy = max(abs(aay[0,:]))
    peaktopdispx = max(abs(dx[0,:]))
    peaktopdispy = max(abs(dy[0,:]))
    peakbasedispx = max(abs(dx[ndof-1,:]))
    peakbasedispy = max(abs(dy[ndof-1,:]))

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
    
    result = ResultFixedXY(ref, ijk, time.T, gx.T, dx.T, vx.T, ax.T, aax.T, gy.T, dy.T, vy.T, ay.T, aay.T, fx, fy, ek, ed, es, ei, error, smx, skx, cdx, smy, sky, cdy)
    model = ModelInfo(ndof)

    # acc = pd.DataFrame(np.hstack((time.T, gx.T)))

    # acc.to_csv("results\\Result.csv", mode='w', sep=',')

    return result, model
