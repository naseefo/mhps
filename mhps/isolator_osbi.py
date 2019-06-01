
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
import scipy.integrate as integrate
import warnings
from mhps.colorfn import *

from data.defaults.param_manager import *

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

class IsoOSBIModel:    
    def osbt2x(self): 
        default_values = get_default_param_values()
        tdiv = default_values['TDIV'] 
        fun1 = lambda x: sqrt(1-pow(self.ecc*sin(x), 2.0))
        xbtab = np.zeros((tdiv, ), dtype=np.dtype('d'), order ='F')
        i = 0
        for td0 in np.arange(0.0,89.9999*pi/180,89.9999*pi/180/tdiv):
            trd0 = atan((self.b0/self.a0)*tan(td0))
            rd0 = sqrt(pow(self.a0*sin(trd0), 2.0) + pow(self.b0*cos(trd0), 2.0))
            cd0 = self.a0*sin(td0)*cos(trd0) - self.b0*cos(td0)*sin(trd0)
            itheta = integrate.quad(fun1, 0.0, td0)
            xbtab[i] = 2.0*self.a0*itheta[0] - 2.0*cd0
            i += 1
        ttab = np.arange(0.0,89.9999*pi/180,89.9999*pi/180/tdiv, dtype=np.dtype('d'))
        return xbtab, ttab
    def u_max(self):
        umax = 2.0*pi/4.0*(3.0*(self.a0 + self.b0) - sqrt((3.0*self.a0 + self.b0)*(self.a0 + 3*self.b0)))
        if umax < self.umax:
            print(prRed('ATTENTION: ') + prCyan('Permissible max. displacement = %8.4f m | Provided max. displacement = %8.4f m'%(umax, self.umax)))
            print('Changing OSBI_VP input file with umax = %8.4f m...'%(umax))
            return umax
        return self.umax

    def __init__(self, rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, kg, cg, dg, am, niso):
        self.rmbm = rmbm
        self.tbx = tbx
        self.zetabx = zetabx
        self.rtytxb = rtytxb
        self.rzyzxb = rzyzxb
        self.typevf = typevf
        self.mu0 = mu0
        self.alpha0 = alpha0
        self.alpha1 = alpha1
        self.nu = nu
        self.umax = umax
        self.D = D
        self.ecc = ecc
        self.rmrm = rmrm
        self.kg = kg
        self.cg = cg
        self.fos_ud = dg

        self.a0 = D/2.0
        self.b0 = sqrt((D/2.0)**2.0 - (ecc*D/2.0)**2.0)
        self.Mr = am*rmrm
        self.mr = self.Mr/niso
        self.Jr = (self.Mr/4.0)*(pow(self.a0, 2.0) + pow(self.b0, 2.0))
        # self.Jr = (self.Mr/5.0)*(pow(self.a0, 2.0) + pow(self.b0, 2.0))
        self.jr = self.Jr/niso
        self.umax = self.u_max()
        self.ud = self.fos_ud*self.umax
        self.xbtab, self.ttab = self.osbt2x()

        theta_d0 = np.interp(self.ud, self.xbtab, self.ttab)
        self.theta_r_max_lower = atan((self.b0/self.a0)*tan(theta_d0))
        p_d0 = self.a0*sin(theta_d0)*sin(self.theta_r_max_lower) + self.b0*cos(theta_d0)*cos(self.theta_r_max_lower)
        self.Lc0 = sqrt(pow(self.ud, 2.0) + pow(2.0*p_d0, 2.0))

        theta_d0 = np.interp(self.umax, self.xbtab, self.ttab)
        self.theta_r_max_upper = atan((self.b0/self.a0)*tan(theta_d0))
        p_d0 = self.a0*sin(theta_d0)*sin(self.theta_r_max_upper) + self.b0*cos(theta_d0)*cos(self.theta_r_max_upper)
        self.Lcf = sqrt(pow(self.umax, 2.0) + pow(2.0*p_d0, 2.0))
        self.maxstrain = (self.Lcf - self.Lc0)/self.Lcf
    
    def __str__(self):
        line = '\nPlan diameter of OSBI (D) = %8.4f cm'%(self.D*100.0)
        line += '\nMajor radius of OSBI (a0) = %8.4f cm'%(self.a0*100.0)
        line += '\nMinor radius of OSBI (b0) = %8.4f cm'%(self.b0*100.0)
        line += '\nEccentricity of OSBI (ecc) = %8.4f'%(self.ecc)
        line += '\nRatio of mass of all OSBI to floor mass (rmrm) = %8.4f'%(self.rmrm)
        line += '\nTotal mass of all OSBI (Mr) = %8.4f kg'%(self.Mr)
        line += '\nMass of single OSBI (mr) = %8.4f kg'%(self.mr)
        line += '\nMass Moment of inertia of all OSBI (Jr) = %8.4f kg-m2'%(self.Jr)
        line += '\nMass Moment of inertia inertia of single OSBI (jr) = %8.4f kg-m2'%(self.jr)
        line += '\nMaximum displacement of OSBI (umax) = %8.4f cm'%(self.umax*100.0)
        line += '\nType of friction profile (typevf) = %d | alpha0 = %8.4f | alpha1 = %8.4f | nu = %8.4f '%(self.typevf, self.alpha0, self.alpha1, self.nu)
        line += '\nRolling friction coefficient (mu0) = %8.4f'%(self.mu0)
        # line += '\nInitial time period of OSBI (Tb0) = %8.4f'%(self.Tb0)
        line += '\n\nTime period of linear spring in X-direction (tbx) = %8.4f sec'%(self.tbx)
        line += '\nTime period of linear spring in Y-direction (tby) = %8.4f sec'%(self.tbx*self.rtytxb)
        line += '\nCritical viscous damping ratio in X-direction (zetabx) = %8.4f p.c.'%(self.zetabx)
        line += '\nCritical viscous damping ratio in Y-direction (zetaby) = %8.4f p.c.'%(self.zetabx*self.rzyzxb)
        line += '\n\nFactor of safety for maximum displacement (fos_ud)  = %8.4f '%(self.fos_ud)
        line += '\nLower maximum design displacement of OSBI (ud)  = %8.4f cm'%(self.ud*100)
        line += '\nLower maximum design rotation of OSBI (theta_r_max_lower)  = %8.4f deg'%(self.theta_r_max_lower*180.0/pi)
        line += '\nUpper maximum design displacement of OSBI (umax)  = %8.4f cm'%(self.umax*100)
        line += '\nUpper maximum design rotation of OSBI (theta_r_max_upper)  = %8.4f deg'%(self.theta_r_max_upper*180.0/pi)
        line += '\nLength of cable at lower maximum design displacement (Lc0) = %8.4f cm'%(self.Lc0*100.0)
        line += '\nLength of cable at upper maximum design displacement (Lcf)  = %8.4f p.c.'%(self.Lcf*100.0)
        line += '\nStrain in cable at maximum design displacement (maxstrain) = %8.4f'%(self.maxstrain)
        line += '\nRatio of all cable stiffness to floor  = %8.4f p.c.'%(self.zetabx*self.rzyzxb)
        line += '\n'
        return line

    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.rmbm == other.rmbm and self.tbx == other.tbx and \
                self.zetabx == other.zetabx and self.rtytxb == other.rtytxb and  self.rzyzxb == other.rzyzxb and \
                    self.typevf == other.typevf and self.mu0 == other.mu0 and self.alpha0 == other.alpha0 and self.alpha1 == other.alpha1 and \
                        self.nu == other.nu and self.umax == other.umax and self.kg == other.kg and self.cg == other.cg and self.dg == other.dg and \
                            self.D == other.D and self.ecc == ecc and self.rmrm == other.rmrm
        return False

def read_iso_osbi_var_param(var_param_file, am, niso):

    """
    This function reads all the set of parameters for the parametric
    studies and stores it in an array. It will be a one-time allocation
    to increase speed.
    """
    # print(var_param_file)
    dtype1 = np.dtype([('IJK', 'i4'), ('RMBM', 'd'), ('TBX', 'd'), ('ZETABX','d'), ('RTYTXB','d'), ('RZYZXB', 'd'), ('TYPEVF', 'i4'), ('MU0', 'd'), ('ALPHA0', 'd'), ('ALPHA1', 'd'), ('NU', 'd'), ('UMAX', 'd'), ('D', 'd'), ('ECC', 'd'), ('RMRM', 'd'), ('KG', 'd'), ('CG', 'd'), ('DG', 'd')])
    
    ijk, rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, kg, cg, dg = np.loadtxt \
    (var_param_file, delimiter=',', usecols=range(18), skiprows=1, unpack = True, dtype=dtype1)

    try:
        total_param = len(ijk)
    except:
        total_param = 1


    for i in range(0, total_param):
        try:
            yield IsoOSBIModel(rmbm[i], tbx[i], zetabx[i], rtytxb[i], rzyzxb[i], typevf[i], mu0[i], alpha0[i], alpha1[i], nu[i], umax[i], D[i], ecc[i], rmrm[i], kg[i], cg[i], dg[i], am, niso)
        except:
            yield IsoOSBIModel(rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, kg, cg, dg, am, niso)


def fs1(M, mr, y_b_d2, c_d0, p_d0, drb, dyb, dxb):

    g = 9.81 # in %m/s2

    fs1r = ((2*M + mr)*g + (2*M + 0.5*mr)*y_b_d2)*(c_d0/2.0/p_d0)

    if drb == 0:
        fs1x = 0.0
        fs1y = 0.0
    else:
        fs1x = fs1r*dxb/drb
        fs1y = fs1r*dyb/drb

    return fs1x, fs1y

def fs1fixed(M, mr, c_d0, p_d0, drb, dyb, dxb):


    g= 9.91
    fs1r = (2*M + mr)*g*(c_d0/2.0/p_d0)
    # fs1r = (2*M + 0.0)*9.81*(c_d0/2.0/p_d0)
    
    if drb == 0:
        fs1x = 0.0
        fs1y = 0.0
    else:
        fs1x = fs1r*dxb/drb
        fs1y = fs1r*dyb/drb

    return fs1x, fs1y

def fb(i, p_index, J, mr, drb, dyb, dxb, vrb, vyb, vxb, arg, ayg, axg, arb, ayb, axb, p_d0, theta_r_dot2):

    if arb == 0:     # Don't forget to change here!!!
        fbr1x = 0.0
        fbr1y = 0.0
    else:
        if i == p_index:
            print('Resultant base acceleration %8.4f | theta_r_dot_2 = %8.4f'%(arb, theta_r_dot2))
            print('I am calculating fbr1')
        fbr1 = (1/2.0/p_d0)*(J*theta_r_dot2)
        # fbr1x = fbr1*vxb/vrb 
        # fbr1y = fbr1*vyb/vrb

        # fbr1x = fbr1*dxb/drb 
        # fbr1y = fbr1*dyb/drb

        fbr1x = fbr1*axb/arb 
        fbr1y = fbr1*ayb/arb

    if arg == 0:
        fbr2x = 0.0
        fbr2y = 0.0
    else:
        if i == p_index:
            print('Resultant ground acceleration %8.4f'%(arg))
            print('I am calculating fbr1')
        fbr2 = 0.5*mr*arg
        fbr2x = fbr2*axg/arg 
        fbr2y = fbr2*ayg/arg

    if arb == 0:
        fbr3x = 0.0
        fbr3y = 0.0
    else:
        if i == p_index:
            print('I am calculating fbr1')
        fbr3 =  0.5*mr*0.5*arb
        fbr3x = fbr3*axb/arb 
        fbr3y = fbr3*ayb/arb

    fbx = fbr1x + fbr2x + fbr3x

    fby = fbr1y + fbr2y + fbr3y

    return fbx, fby

def fbfixed(mr, arg, ayg, axg):

    fbr = 0.5*mr*arg

    if fbr == 0:
        fbx = 0.0
        fby = 0.0
    else:
        fbx = fbr*axg/arg
        fby = fbr*ayg/arg

    return fbx, fby

def fs2(M, mr, mu, y_b_d2):

    g = 9.81

    fs2r = 0.5*mu*((2*M + mr)*g + (2*M + 0.5*mr)*y_b_d2)
    
    return fs2r

def fs2fixed(mu, M, mr):

    g = 9.81

    fs2r = 0.5*mu*(2*M + mr)*g

    return fs2r


def stat1(fs2x, fs2y, qx, qy):

    # Calculating ratio
    ratio = sqrt(pow(fs2x/qx, 2.0) + pow(fs2y/qy, 2.0))
    
    if (ratio - 1.0) < 0:
        # If ratio is less than then set condition for analysis as for fixed-base structure
        rolling_state = False
    else:
        # If ratio is equal to or more than one set condition for analysis as for base-isolated structure
        # fs2x and fs2y are calculated based on the circular interaction

        fs2x = fs2x/ratio  # Component of rolling resistance imparted in the x-direction
    
        
        fs2y = fs2y/ratio  # Component of rolling resistance imparted in the y-direction
      
        rolling_state = True

    return rolling_state, fs2x, fs2y

def stat0(ddxb, ddyb, fs2x, fs2y):
    wd = fs2x*ddxb + fs2y*ddyb
    if wd >= 0:
        rolling_state = True
    else:
        rolling_state = False
    return rolling_state

def mu_val(iso, ub):
    L0 = iso.alpha0*iso.mu0
    L1 = (iso.alpha1-iso.alpha0)*(iso.mu0/iso.umax)*abs(ub)
    L2 = (iso.alpha1-iso.alpha0)*(iso.mu0/pow(iso.umax,2))*pow(ub,2)
    if iso.typevf == 1:
        mu = L0 + (1-iso.nu)*L1 + iso.nu*L2
    elif iso.typevf == 2:
        mu = L0 - (1-iso.nu)*L1 - iso.nu*(L2 - 2*L1)
    else:
        mu = iso.mu0

    if abs(ub) > iso.umax and abs(ub) < iso.dos_ud*iso.umax:
        mu = iso.alpha1*iso.mu0

    if abs(ub) > iso.fos_ud*iso.umax and iso.typevf != 0:
        mu = 100
    return mu

#@profile 
def simulator_osbi(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, ndof, smx, skx, cdx, smy, sky, cdy, iso, screen_on):
    
    # C-Order: smx, smy, 

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

    time = np.zeros((1, ndt), dtype=np.dtype('d'), order ='F') # C
    
    dx = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F') #C
    vx = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F') #C
    ax = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F') #C
    aax = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F') #C
    gx = np.zeros((1, ndt), dtype=np.dtype('d'), order ='F') #C
    ddx = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F') #C
    dvx = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F') #C

    dy = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F') #C
    vy = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F') #C
    ay = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F') #C
    aay = np.zeros((ndof, ndt), dtype=np.dtype('d'), order ='F') #C
    gy = np.zeros((1, ndt), dtype=np.dtype('d'), order ='F') #C
    ddy = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F') #C
    dvy = np.zeros((ndof, 1), dtype=np.dtype('d'), order ='F') #C

    gx[0,0] = xg[0]
    gy[0,0] = yg[0]

    fx = np.zeros((ndt, ndof), dtype=np.dtype('d'), order ='F') #C
    fy = np.zeros((ndt, ndof), dtype=np.dtype('d'), order ='F') #C
    ek = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    ed = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    es = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    ei = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    error = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C

    dx1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    vx1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001 #C
    ax1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    px1 = smx_diag*xg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ax1[0:nst, 0] = np.dot(smx_inv_fixed, px1[0:nst, 0] - np.dot(cdx[0:nst, 0:nst], vx1[0:nst, 0]) - np.dot(skx[0:nst, 0:nst], dx1[0:nst, 0]))
    dy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    vy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001 #C
    ay1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    py1 = smy_diag*yg[0] # px1 = smx_diag*xg[0] # Initial earthquake acceleration is considered zero
    ay1[0:nst, 0] = np.dot(smy_inv_fixed, py1[0:nst, 0] - np.dot(cdy[0:nst, 0:nst], vy1[0:nst,0]) - np.dot(sky[0:nst, 0:nst], dy1[0:nst,0]))

    # I = np.ones((ndof,1), dtype=np.dtype('d'), order ='F')
    dx2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    vx2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001 #C
    px2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    ax2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    dy2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    vy2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001 #C
    py2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    ay2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    
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
    rolling_state = False

    mr = iso.rmrm*smx[0,0]      # Here mr is sum of masses of all OSBI
    J = (mr/4.0)*(pow(iso.a0, 2.0) + pow(iso.b0, 2.0))

    print(mr)

    fs1x = fs1y = 0.0
    fs2x = fs2y = 0.0
    fbx = fbx = 0.0
    ub = 0.0

    p_index = 29465
    for i in range(1, len(xg)):
        
        # for i in range(1,20):

        t += dt

        px2 = smx_diag*xg[i]
        py2 = smy_diag*yg[i]
        dpx = (px2 - px1)
        dpy = (py2 - py1)

        if rolling_state == False:
                
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
            
            drb = sqrt(pow(dy2[ndof-1, 0], 2.0) + pow(dx2[ndof-1, 0], 2.0))
            # if drb > iso.umax:
            #     print(prRed('WARNING: ') + prCyan('Isolator has crossed maximum displacement.'))
            # print(dy2[ndof-1, 0], dx2[ndof-1, 0])
            theta_d0 = np.interp(drb, iso.xbtab, iso.ttab)
            theta_r_d0 = atan((iso.b0/iso.a0)*tan(theta_d0))
            c_d0 = iso.a0*sin(theta_d0)*cos(theta_r_d0) - iso.b0*cos(theta_d0)*sin(theta_r_d0)
            p_d0 = iso.a0*sin(theta_d0)*sin(theta_r_d0) + iso.b0*cos(theta_d0)*cos(theta_r_d0) 
            fs1x, fs1y = fs1fixed(tm, mr, c_d0, p_d0, drb, dy2[ndof-1, 0], dx2[ndof-1, 0])

            arg = sqrt(pow(xg[i], 2.0) + pow(yg[i], 2.0))

            fbx, fby = fbfixed(mr, arg, yg[i], xg[i])

            fs2x = (-1.0*cdx[nst, nst-1])*vx2[nst-1,0] + (-1.0*skx[nst, nst-1])*dx2[nst-1,0] - cdx[ndof-1, ndof-1]*vx2[ndof-1,0] - skx[ndof-1, ndof-1]*dx2[ndof-1,0] + px2[ndof-1] - fs1x - fbx
           
            fs2y = (-1.0*cdy[nst, nst-1])*vy2[nst-1,0] + (-1.0*sky[nst, nst-1])*dy2[nst-1,0] - cdy[ndof-1, ndof-1]*vy2[ndof-1,0] - sky[ndof-1, ndof-1]*dy2[ndof-1,0] + py2[ndof-1] - fs1y - fby
           

            if i == p_index:
                print(i, rolling_state, fbx, fby, fs1x, fs1y, fs2x, fs2y)

            mu = mu_val(iso, drb)
            qr = fs2fixed(mu, tm, mr)
            rolling_state, fs2x, fs2y = stat1(fs2x, fs2y, qr, qr)

            if rolling_state == False:
                ax2[0:nst, 0] = np.dot(smx_inv_fixed, px2[0:nst, 0] - np.dot(cdx[0:nst, 0:nst], vx2[0:nst, 0]) - np.dot(skx[0:nst, 0:nst], dx2[0:nst, 0]))
                ay2[0:nst, 0] = np.dot(smy_inv_fixed, py2[0:nst, 0] - np.dot(cdy[0:nst, 0:nst], vy2[0:nst, 0]) - np.dot(sky[0:nst, 0:nst], dy2[0:nst, 0]))
            else:
                epx = px2 - np.dot(cdx, vx2) - np.dot(skx, dx2)
                epy = py2 - np.dot(cdy, vy2) - np.dot(sky, dy2)

                epx[ndof-1,0] = epx[ndof-1,0] - fs2x - fs1x - fbx    
                epy[ndof-1,0] = epy[ndof-1,0] - fs2y - fs1y - fby  

                ax2 = np.dot(smx_inv, epx)
                ay2 = np.dot(smy_inv, epy)
            
            dx1, vx1, px1, ax1 = dx2, vx2, px2, ax2 
            dy1, vy1, py1, ay1 = dy2, vy2, py2, ay2 

        else:
            
            dfbx = 0.0
            dfby = 0.0

            dfs1x = 0.0
            dfs1y = 0.0

            dfs2x = 0.0
            dfs2y = 0.0

            for ii in range (nit):
                  
                pcx1 = dpx + np.dot(na2x, vx1) + np.dot(na3x, ax1)
                pcx1[ndof-1,0] = pcx1[ndof-1,0] - dfs2x - dfs1x - dfbx
                ddx = np.dot(knx_inv, pcx1)
                dx2 = dx1 + ddx
                dvx = (gamma/beta/dt)*ddx - gamma/beta*vx1 + dt*(1.0 - gamma/2.0/beta)*ax1
                vx2 = vx1 + dvx

                pcy1 = dpy + np.dot(na2y, vy1) + np.dot(na3y, ay1)
                pcy1[ndof-1,0] = pcy1[ndof-1,0] - dfs2y - dfs1y - dfby
                ddy = np.dot(kny_inv, pcy1)
                dy2 = dy1 + ddy
                dvy = (gamma/beta/dt)*ddy - gamma/beta*vy1 + dt*(1.0 - gamma/2.0/beta)*ay1
                vy2 = vy1 + dvy

                ##@ Maybe we may need to calculate the acceleration here as its used in various calculation of isolator forces
                epx = px2 - np.dot(cdx, vx2) - np.dot(skx, dx2)
                epy = py2 - np.dot(cdy, vy2) - np.dot(sky, dy2)
                epx[ndof-1,0] = epx[ndof-1,0] - fs2x - fs1x - fbx
                epy[ndof-1,0] = epy[ndof-1,0] - fs2y - fs1y - fby
                ax2 = np.dot(smx_inv, epx)
                ay2 = np.dot(smy_inv, epy)

                drb = sqrt(pow(dx2[ndof-1, 0], 2) + pow(dy2[ndof-1, 0], 2))
                # if drb > iso.umax:
                #     print(prRed('WARNING: ') + prCyan('Isolator has crossed maximum displacement.'))
                # phi_drb = atan(dy2[ndof-1, 0]/dx2[ndof-1, 0])
                vrb = sqrt(pow(vx2[ndof-1, 0], 2) + pow(vy2[ndof-1, 0], 2))
                # phi_vrb = atan(vy2[ndof-1, 0]/vx2[ndof-1, 0])

                ##? Check if "absolute" resultant acceleration needs to be passed ????
                arb = sqrt(pow(ax2[ndof-1, 0], 2) + pow(ay2[ndof-1, 0], 2))
                # phi_arb = atan(ay2[ndof-1, 0]/ax2[ndof-1, 0])

                arg = sqrt(pow(xg[i], 2.0) + pow(yg[i], 2.0))

                if vrb > 1.0e-5:
                    
                    theta_d0 = np.interp(drb, iso.xbtab, iso.ttab)
                    theta_r_d0 = atan((iso.b0/iso.a0)*tan(theta_d0))
                    theta_r_d1 = sin(2*theta_r_d0)/sin(2*theta_d0)
                    theta_r_d2 = 2*theta_r_d1*(tan(theta_d0) - theta_r_d1*tan(theta_r_d0))
                    F_d0 = iso.a0*sqrt(1 - pow(iso.ecc*sin(theta_d0), 2.0))
                    F_d1 = -pow(iso.a0*iso.ecc, 2.0)*sin(theta_d0)*cos(theta_d0)/F_d0
                    c_d0 = iso.a0*sin(theta_d0)*cos(theta_r_d0) - iso.b0*cos(theta_d0)*sin(theta_r_d0)
                    p_d0 = iso.a0*sin(theta_d0)*sin(theta_r_d0) + iso.b0*cos(theta_d0)*cos(theta_r_d0)
                    c_d1 = iso.a0*cos(theta_d0)*cos(theta_r_d0) + iso.b0*sin(theta_d0)*sin(theta_r_d0) - p_d0*theta_r_d1
                    p_d1 = iso.a0*cos(theta_d0)*sin(theta_r_d0) - iso.b0*cos(theta_d0)*sin(theta_r_d0) + c_d0*theta_r_d1
                    c_d2 = -p_d0*theta_r_d2 - 2*p_d1*theta_r_d1 - c_d0 + c_d0*pow(theta_r_d1, 2.0)
                    p_d2 = c_d0*theta_r_d2 + 2*c_d1*theta_r_d1 - p_d0 + p_d0*pow(theta_r_d1, 2.0)
                    theta_dot1 = (0.5*vrb)/(F_d0 - c_d1)                                               #vbr
                    theta_dot2 = ((F_d1 - c_d2)*pow(theta_dot1, 2.0) - 0.5*arb)/(c_d1 - F_d0)     #abr
                    theta_r_dot1 = theta_r_d1*theta_dot1
                    theta_r_dot2 = theta_r_d2*pow(theta_dot1, 2.0) + theta_r_d1*theta_dot2
                    y_r_d2 = pow(theta_dot1, 2.0)*p_d2 + p_d1*theta_dot2
                    y_b_d2 = 2*y_r_d2

                    fbx2, fby2 = fb(i, p_index, J, mr, drb, dy2[ndof-1, 0], dx2[ndof-1, 0], vrb, vy2[ndof-1, 0], vx2[ndof-1, 0], arg, yg[i], xg[i], arb, ay2[ndof-1, 0], ax2[ndof-1, 0], p_d0, theta_r_dot2)  ##?
                    dfbx = fbx2 - fbx
                    dfby = fby2 - fby
                    
                    fs1x2, fs1y2 = fs1(tm, mr, y_b_d2, c_d0, p_d0, drb, dy2[ndof-1, 0], dx2[ndof-1, 0]) ##?
                    dfs1x = fs1x2 - fs1x
                    dfs1y = fs1y2 - fs1y
                    
                    mu = mu_val(iso, drb)    ##?

                    qr = fs2(tm, mr, mu, y_b_d2)   ##?

                    if vrb == 0:
                        fs2x2 = 0.0
                        fs2y2 = 0.0
                    else:
                        fs2x2 = qr*vx2[ndof-1, 0]/vrb
                        fs2y2 = qr*vy2[ndof-1, 0]/vrb

                    dfs2x = fs2x2 - fs2x
                    dfs2y = fs2y2 - fs2y
                
            if nit > 1:

                fbx = fbx + dfbx
                fby = fby + dfby

                fs1x = fs1x + dfs1x
                fs1y = fs1y + dfs1y

                fs2x = fs2x + dfs2x
                fs2y = fs2y + dfs2y

                if i == p_index:
                    print(i, rolling_state, fbx, fby, fs1x, fs1y, fs2x, fs2y)
            
            rolling_state = stat0(ddx[ndof-1], ddy[ndof-1], fs2x, fs2y)

            if rolling_state == False:
                vx2[ndof-1, 0] = 0.0 # 1e-10
                ax2[ndof-1, 0] = 0.0
                vy2[ndof-1, 0] = 0.0 # 1e-10
                ay2[ndof-1, 0] = 0.0

                ax2[0:nst, 0] = np.dot(smx_inv_fixed, px2[0:nst, 0] - np.dot(cdx[0:nst, 0:nst], vx2[0:nst, 0]) - np.dot(skx[0:nst, 0:nst], dx2[0:nst, 0]))
                ay2[0:nst, 0] = np.dot(smy_inv_fixed, py2[0:nst, 0] - np.dot(cdy[0:nst, 0:nst], vy2[0:nst, 0]) - np.dot(sky[0:nst, 0:nst], dy2[0:nst, 0]))
            else:
                epx = px2 - np.dot(cdx, vx2) - np.dot(skx, dx2)
                epy = py2 - np.dot(cdy, vy2) - np.dot(sky, dy2)

                epx[ndof-1,0] = epx[ndof-1,0] - fs2x - fs1x - fbx
                epy[ndof-1,0] = epy[ndof-1,0] - fs2y - fs1y - fby

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
            fx[index, ndof-1] = fs2x
            fy[index, ndof-1] = fs2y
            ek[index, 0] = eki
            ed[index, 0] = edi
            es[index, 0] = esi
            ei[index, 0] = eii

            error[index, 0] = (edi + esi + eki - eii)/(abs(edi + esi) + eki + abs(eii))
            # print(rolling_state)
   
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

    return result, model
