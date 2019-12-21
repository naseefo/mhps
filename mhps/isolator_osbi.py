import sys
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
        
        ld0 = 0.0
        ud0 = pi/2
        delta = ud0/tdiv
        sizeTable = int((ud0 - ld0)/delta + 1)
        # print(sizeTable)
        xbtab = np.zeros((sizeTable, ), dtype=np.dtype('d'), order ='F')
        ttab = np.zeros((sizeTable, ), dtype=np.dtype('d'), order ='F')
        trtab = np.zeros((sizeTable, ), dtype=np.dtype('d'), order ='F')
        cd0tab = np.zeros((sizeTable, ), dtype=np.dtype('d'), order ='F')
        pd0tab = np.zeros((sizeTable, ), dtype=np.dtype('d'), order ='F')
        
        fun1 = lambda x: sqrt(1-pow(self.ecc*sin(x), 2.0))
        i = 0
        for td0 in np.arange(ld0, ud0 + 0*delta, delta):
            # print(i)
            ttab[i] = td0
            trd0 = atan((self.b0/self.a0)*tan(td0))
            trtab[i] = trd0
            cd0 = self.a0*sin(td0)*cos(trd0) - self.b0*cos(td0)*sin(trd0)
            pd0 = self.a0*sin(td0)*sin(trd0) + self.b0*cos(td0)*cos(trd0)
            cd0tab[i] = cd0
            pd0tab[i] = pd0
            itheta = integrate.quad(fun1, 0.0, td0)
            xbtab[i] = 2.0*self.a0*itheta[0] - 2.0*cd0
            i += 1
        
        activate_1 = 0 # To print the values of the various terms in restoring force to a file
        if activate_1 == 1:
            pd.DataFrame(xbtab).to_csv("xbtab.csv")
            pd.DataFrame(ttab).to_csv("ttab.csv")
            pd.DataFrame(trtab).to_csv("trtab.csv")
            pd.DataFrame(cd0tab).to_csv("cd0tab.csv")
            pd.DataFrame(pd0tab).to_csv("pd0tab.csv")
        
        return xbtab, ttab

    def u_max(self):
        fun1 = lambda x: sqrt(1-pow(self.ecc*sin(x), 2.0))
        trd0 = pi/2
        td0 = atan((self.a0/self.b0)*tan(trd0))
        cd0 = self.a0*sin(td0)*cos(trd0) - self.b0*cos(td0)*sin(trd0)
        itheta = integrate.quad(fun1, 0.0, td0)
        umax = 2.0*self.a0*itheta[0] - 2.0*cd0
        if umax < self.umax:
            print(prRed('ATTENTION: ') + prCyan('Permissible max. displacement = %8.4f m | Provided max. displacement = %8.4f m'%(umax, self.umax)))
            print('Changing OSBI_VP input file with umax = %8.4f m...'%(umax))
            return umax
        
        return umax

    def __init__(self, rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, tc, ecrc, fos_ud, am, niso):
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
        self.tc = tc
        self.ecrc = ecrc
        self.fos_ud = fos_ud
        self.niso = niso
        self.am = am

        self.a0 = D/2.0

        self.b0 = sqrt((D/2.0)**2.0 - (ecc*D/2.0)**2.0)

        self.Mr = am*rmrm

        self.mr = self.Mr/niso

        self.Jr = (self.Mr/5.0)*(pow(self.a0, 2.0) + pow(self.b0, 2.0))

        self.jr = self.Jr/niso

        self.umax = self.u_max()

        self.ud = self.fos_ud*self.umax

        self.xbtab, self.ttab = self.osbt2x()

        self.V = 4.0/3.0*pi*pow(self.a0, 2.0)*self.b0

        self.ro = self.mr/self.V

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
        line = '\nNumber of isolators (niso) = ' + prCyan('%d ')%(self.niso)
        line += '\nPlan diameter of OSBI (D) = ' + prCyan('%8.4f cm')%(self.D*100.0)
        line += '\nMajor radius of OSBI (a0) = ' + prCyan('%8.4f cm')%(self.a0*100.0)
        line += '\nMinor radius of OSBI (b0) = ' + prCyan('%8.4f cm')%(self.b0*100.0)
        line += '\nEccentricity of OSBI (ecc) = ' + prCyan('%8.4f')%(self.ecc)
        line += '\nRatio of mass of all OSBI to floor mass (rmrm) = ' + prCyan('%8.4f')%(self.rmrm)
        line += '\nTotal mass of all OSBI (Mr) = ' + prCyan('%8.4f kg')%(self.Mr)
        line += '\nMass of single OSBI (mr) = ' + prCyan('%8.4f kg')%(self.mr)
        line += '\nVolume of single OSBI (V) = ' + prCyan('%8.4f m3'%(self.V))
        line += '\nDensity of material required for OSBI (ro) = ' + prCyan('%8.4f kg/m3')%(self.ro)
        line += '\nMass Moment of inertia of all OSBI (Jr) = ' + prCyan('%8.4f kg-m2')%(self.Jr)
        line += '\nMass Moment of inertia inertia of single OSBI (jr) = ' + prCyan('%8.4f kg-m2')%(self.jr)
        line += '\nMaximum displacement of OSBI (umax) = ' + prCyan('%8.4f cm')%(self.umax*100.0)
        line += '\nType of friction profile (typevf) = ' + prCyan('%d ')%(self.typevf) + '|alpha0 = ' + prCyan('%8.4f ')%(self.alpha0) + '| alpha1 = ' + prCyan('%8.4f ')%(self.alpha1) + '| nu = ' + prCyan('%8.4f ')%(self.nu)
        line += '\nRolling friction coefficient (mu0) = ' + prCyan('%8.4f')%(self.mu0)
        # line += '\nInitial time period of OSBI (Tb0) = %8.4f'%(self.Tb0)
        line += '\n\nTime period of linear spring in X-direction (tbx) = ' + prCyan('%8.4f sec')%(self.tbx)
        line += '\nTime period of linear spring in Y-direction (tby) = ' + prCyan('%8.4f sec')%(self.tbx*self.rtytxb)
        line += '\nCritical viscous damping ratio in X-direction (zetabx) = ' + prCyan('%8.4f p.c.')%(self.zetabx)
        line += '\nCritical viscous damping ratio in Y-direction (zetaby) = ' + prCyan('%8.4f p.c.')%(self.zetabx*self.rzyzxb)
        line += '\n\nFactor of safety for maximum displacement (fos_ud)  = ' + prCyan('%8.4f ')%(self.fos_ud)
        line += '\nLower maximum design displacement of OSBI (ud)  = ' + prCyan('%8.4f cm')%(self.ud*100)
        line += '\nLower maximum design rotation of OSBI (theta_r_max_lower)  = ' + prCyan('%8.4f deg')%(self.theta_r_max_lower*180.0/pi)
        line += '\nUpper maximum design displacement of OSBI (umax)  = ' + prCyan('%8.4f cm')%(self.umax*100)
        line += '\nUpper maximum design rotation of OSBI (theta_r_max_upper)  = ' + prCyan('%8.4f deg')%(self.theta_r_max_upper*180.0/pi)
        line += '\nLength of cable at lower maximum design displacement (Lc0) = ' + prCyan('%8.4f cm')%(self.Lc0*100.0)
        line += '\nLength of cable at upper maximum design displacement (Lcf)  = ' + prCyan('%8.4f p.c.')%(self.Lcf*100.0)
        line += '\nStrain in cable at maximum design displacement (maxstrain) = ' + prCyan('%8.4f')%(self.maxstrain)
        line += '\n'
        return line

    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.rmbm == other.rmbm and self.tbx == other.tbx and \
                self.zetabx == other.zetabx and self.rtytxb == other.rtytxb and  self.rzyzxb == other.rzyzxb and \
                    self.typevf == other.typevf and self.mu0 == other.mu0 and self.alpha0 == other.alpha0 and self.alpha1 == other.alpha1 and \
                        self.nu == other.nu and self.umax == other.umax and self.tc == other.tc and self.ecrc == other.ecrc and self.fos_ud == other.fos_ud and \
                            self.D == other.D and self.ecc == ecc and self.rmrm == other.rmrm
        return False

def read_iso_osbi_var_param(var_param_file, am, niso):

    """
    This function reads all the set of parameters for the parametric
    studies and stores it in an array. It will be a one-time allocation
    to increase speed.
    """
    # print(var_param_file)
    dtype1 = np.dtype([('IJK', 'i4'), ('RMBM', 'd'), ('TBX', 'd'), ('ZETABX','d'), ('RTYTXB','d'), ('RZYZXB', 'd'), ('TYPEVF', 'i4'), ('MU0', 'd'), ('ALPHA0', 'd'), ('ALPHA1', 'd'), ('NU', 'd'), ('UMAX', 'd'), ('D', 'd'), ('ECC', 'd'), ('RMRM', 'd'), ('TC', 'd'), ('ECRC', 'd'), ('FOS_UD', 'd')])
    
    ijk, rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, tc, ecrc, fos_ud = np.loadtxt \
    (var_param_file, delimiter=',', usecols=range(18), skiprows=1, unpack = True, dtype=dtype1)

    try:
        total_param = len(ijk)
    except:
        total_param = 1


    for i in range(0, total_param):
        try:
            yield IsoOSBIModel(rmbm[i], tbx[i], zetabx[i], rtytxb[i], rzyzxb[i], typevf[i], mu0[i], alpha0[i], alpha1[i], nu[i], umax[i], D[i], ecc[i], rmrm[i], tc[i], ecrc[i], fos_ud[i], am, niso)
        except:
            yield IsoOSBIModel(rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, tc, ecrc, fos_ud, am, niso)


def fs1(M, fcz, rMrM, y_b_d2, c_d0, p_d0, drb, dyb, dxb):

    g = 9.81 # in %m/s2

    fs1r = ((2 + rMrM)*M*g + (2 + 0)*M*y_b_d2)*(c_d0/2.0/p_d0) + fcz*c_d0/p_d0

    if drb == 0:
        fs1x = 0.0
        fs1y = 0.0
    else:
        fs1x = fs1r*dxb/drb
        fs1y = fs1r*dyb/drb

    return fs1x, fs1y

def fs1fixed(M, fcz, rMrM, c_d0, p_d0, drb, dyb, dxb):

    g= 9.81
    fs1r = (2 + rMrM)*M*g*(c_d0/2.0/p_d0) + fcz*c_d0/p_d0
    
    if drb == 0:
        fs1x = 0.0
        fs1y = 0.0
    else:
        fs1x = fs1r*dxb/drb
        fs1y = fs1r*dyb/drb

    return fs1x, fs1y

def fc(iso, dt, dyb, dxb, drb, dc1, p_d0, tm):

    # MECHANICAL PROPERTIES OF THE CABLE
    kc = pow(2.0*pi/iso.tc, 2.0)*tm
    cc = 2.0*iso.ecrc*tm*(2.0*pi/iso.tc)

    # GEOMETRICAL PROPERTY OF THE CABLE
    lc = sqrt(pow(drb, 2.0) + pow((2*p_d0), 2.0))

    # CALCULATION OF CABLE AXIAL DEFORMATION AND VELOCITY
    dc2 = lc - iso.Lc0 # Axial deformation in cable
    vc2 = (dc2 - dc1)/dt # Axial velocity in cable


    if drb >= iso.ud:
        fc_axial = kc*dc2 + cc*vc2
        fcr = fc_axial*drb/lc
        fcz = fc_axial*2*p_d0/lc
        fcx = fcr*dxb/drb
        fcy = fcr*dyb/drb
        e_axial = (lc - iso.Lc0)/iso.Lc0
    else:
        fc_axial = 0.0
        fcz = 0.0
        fcx = 0.0
        fcy = 0.0
        dc2 = 0.0
        vc2 = 0.0
        e_axial = 0.0

        fc_axial, fcx, fcy, fcz, dc2, vc2, e_axial = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    return fc_axial, fcx, fcy, fcz, dc2, vc2, e_axial



def fb(J, mr, drb, dyb, dxb, vrb, vyb, vxb, arg, ayg, axg, arb, ayb, axb, p_d0, theta_r_dot2, c_d0, y_b_d2):
    shut = 1 # 1 for on and 0 for off
    if arb == 0:
        fbr1x = 0.0
        fbr1y = 0.0
    else:
        fbr1 = (1/2.0/p_d0)*(J*theta_r_dot2)*shut
        fbr1x = fbr1*axb/arb 
        fbr1y = fbr1*ayb/arb

    if arg == 0:
        fbr2x = 0.0
        fbr2y = 0.0
    else:
        fbr2 = 0.5*mr*arg*shut
        fbr2x = fbr2*axg/arg 
        fbr2y = fbr2*ayg/arg

    if arb == 0:
        fbr3x = 0.0
        fbr3y = 0.0
    else:
        fbr3 =  0.5*mr*0.5*arb*shut
        fbr3x = fbr3*axb/arb 
        fbr3y = fbr3*ayb/arb
    
    if drb == 0:
        fbr4x = 0.0
        fbr4y = 0.0
    else:
        fbr4 = mr*y_b_d2*c_d0/4.0*p_d0*shut
        fbr4x = fbr4*dxb/drb
        fbr4y = fbr4*dyb/drb

    fbx = fbr1x + fbr2x + fbr3x + fbr4x

    fby = fbr1y + fbr2y + fbr3y + fbr4y

    fbx = 0.0
    fby = 0.0

    return fbx, fby

def fbfixed(mr, ayg, axg):
    shut = 1 # 1 for on and 0 for off

    arg = sqrt(pow(axg, 2.0) + pow(ayg, 2.0)) #CB
    fbr = 0.5*mr*arg*shut

    if arg == 0:
        fbx = 0.0
        fby = 0.0
    else:
        fbx = fbr*axg/arg
        fby = fbr*ayg/arg

    return fbx, fby

def fs2(M, fcz, rMrM, mu, y_b_d2):

    g = 9.81

    fs2r = 0.5*mu*M*((2.0 + rMrM)*g + (2.0 + 0.5*rMrM)*y_b_d2) + 0*mu*fcz
    
    return fs2r, fs2r

def fs2fixed(mu, M, fcz, rMrM):

    g = 9.81

    fs2r = 0.5*mu*M*(2.0 + rMrM)*g + 0*mu*fcz

    return fs2r, fs2r


def stat1(fabx, faby, fs2x, fs2y):

    # Calculating ratio
    ratio = sqrt(pow(fabx/fs2x, 2.0) + pow(faby/fs2y, 2.0))
    
    if (ratio - 1.0) < 0:
        # If ratio is less than then set condition for analysis as for fixed-base structure
        rolling_state = False
    else:
        # If ratio is equal to or more than one set condition for analysis as for base-isolated structure
        # fs2x and fs2y are calculated based on the circular interaction

        fabx = fabx/ratio  # Component of rolling resistance imparted in the x-direction
        faby = faby/ratio  # Component of rolling resistance imparted in the y-direction
      
        rolling_state = True
    return rolling_state, fabx, faby


def stat0(ddxb, ddyb, fs2x, fs2y):
    wd = fs2x*ddxb + fs2y*ddyb
    if wd >= 0:
        rolling_state = True
    else:
        rolling_state = False
    return rolling_state


def mu_val(iso, ub):

    ud = iso.fos_ud*iso.umax

    L0 = iso.alpha0*iso.mu0
    L1 = (iso.alpha1-iso.alpha0)*iso.mu0*abs(ub)/ud
    L2 = (iso.alpha1-iso.alpha0)*iso.mu0*pow(ub, 2.0)/pow(ud, 2.0)

    if abs(ub) <= ud:
        if iso.typevf == 1:
            mu = L0 + (1-iso.nu)*L1 + iso.nu*L2
        elif iso.typevf == 2:
            mu = L0 - (1-iso.nu)*L1 - iso.nu*(L2 - 2.0*L1)
        else:
            mu = iso.mu0
    else:
        mu = iso.alpha1*iso.mu0

    return mu

#@profile 
def simulator_osbi(ref, xg, yg, dt, ndiv, ndt, lxy, ijk, ndof, smx, skx, cdx, smy, sky, cdy, iso, screen_on):
  
    # INITIALIZATION OF PRIMARY INPUT VARIABLES
    nit = 2
    nst = ndof - 1
    gamma = 0.5
    beta = 1/6

    # INITIALIZATIONS OF PRIMARY RESPONSE VARIABLES
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

    # INITIALIZATIOS OF SECONDARY RESPONSE VARIABLES
    fx = np.zeros((ndt, ndof), dtype=np.dtype('d'), order ='F') #C # BASE SHEAR RESPONSE AT ISOLATOR LEVEL
    fy = np.zeros((ndt, ndof), dtype=np.dtype('d'), order ='F') #C

    eki = 0.0                                                      # VARIABLES FOR ENERGY CALCULATIONS 
    edi = 0.0
    esi = 0.0
    eii = 0.0
    ek = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    ed = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    es = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    ei = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C

    roll = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    theta_0 = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    theta_r = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    theta_r_dot2 = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    zbd = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    zbd_dot2 = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    Fs1x = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    Fs1y = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    Fs2x = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    Fs2y = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    Fbx = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    Fby = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    F_axial = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    Dc_axial = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C
    Strain_axial = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C

    error = np.zeros((ndt, 1), dtype=np.dtype('d'), order ='F') #C

    # INITIALIZATIONS OF TEMPORARY VARIABLES
    dx1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    vx1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001 #C
    ax1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    px1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')

    dy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    vy1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001 #C
    ay1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    py1 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')
    
    dx2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    vx2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001 #C
    px2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    ax2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C

    dy2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    vy2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0000000001 #C
    py2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C
    ay2 = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')*0.0 #C

    dpx = 0.0 #INCREMENT IN EXCITATION FORCE
    dpy = 0.0

    r = np.ones((ndof, 1), dtype=np.dtype('d'), order ='F')  # INFLUENCE COEFFICIENT USED IN ENERGY CALCULATIONS

    fs1x = fs1y = 0.0 # RESTORING FORCE
    fs2x = fs2y = 0.0 # FRICTIONAL FORCE
    fbx = fbx = 0.0 # ISOLATOR INTERTIAL FORCE
    fcx = fcy = fcz = 0.0 # ISOLATOR CABLE FORCE

    yd1, yv1, ya1 = 0.0, 0.0, 0.0 # CALCULATION OF VERTICAL ACCELERATION IN BASE SLAB
    yd2, yv2, ya2 = 0.0, 0.0, 0.0 

    theta_r_d01, theta_r_dot11, theta_r_dot21 = 0.0, 0.0, 0.0 # CALCULATION OF ROTATIONAL ACCELERATION IN ISOLATOR
    theta_r_d02, theta_r_dot12, theta_r_dot22 = 0.0, 0.0, 0.0

    dc1, vc1 = 0.0, 0.0 # CALCULATION OF AXIAL DEFORMATION AND VELOCITY IN CABLE
    dc2, vc2 = 0.0, 0.0

    
    # PRE-CALCULATIONS: MASS PARAMETERS
    smx_inv = np.linalg.inv(smx)
    smx_inv_fixed = np.linalg.inv(smx[0:nst, 0:nst])
    smx_diag = np.diag(-1.0*smx).reshape(ndof,1)

    smy_inv = np.linalg.inv(smy)
    smy_inv_fixed = np.linalg.inv(smy[0:nst, 0:nst])
    smy_diag = np.diag(-1.0*smy).reshape(ndof,1)

    tm = np.sum(np.diag(smx))    # TOTAL MASS (INCLUDING SUPERSTRUCTURE AND BASE MASS)

    smxi = np.diag(np.diag(smx)) # FOR ENERGY CALCULATIONS (MASS MATRIX ASSOCIATED WITH EXCITATION)
    smyi = np.diag(np.diag(smy)) # FOR ENERGY CALCULATIONS (MASS MATRIX ASSOCIATED WITH EXCITATION)

    rMrM = iso.Mr/tm # RATIO OF TOTAL ISOLATOR MASS TO TOTAL MASS (INCLUDING SUPERSTRUCTURE AND BASE SLAB)


    # PRE-CALCULATIONS: NEWMARK-BETA METHOD PARAMETER CALCULATIONS
    na1x = (1.0/beta/np.power(dt,2))*smx + (gamma/beta/dt)*cdx
    na2x = (1.0/beta/dt)*smx + gamma/beta*cdx
    na3x = 1.0/2.0/beta*smx + (gamma*dt/2.0/beta - dt)*cdx
    na1y = (1.0/beta/np.power(dt,2))*smy + (gamma/beta/dt)*cdy
    na2y = (1.0/beta/dt)*smy + gamma/beta*cdy
    na3y = 1.0/2.0/beta*smy + (gamma*dt/2.0/beta - dt)*cdy


    # PRE-CALCULATIONS: EFFECTIVE STIFFNESS PARAMETERS FOR NEWMARK-BETA METHOD
    knx = skx + na1x
    knx_inv = np.linalg.inv(knx)
    kny = sky + na1y
    kny_inv = np.linalg.inv(kny)

    knx_fixed = skx[0:nst, 0:nst] + na1x[0:nst, 0:nst]
    knx_inv_fixed = np.linalg.inv(knx_fixed)
    kny_fixed = sky[0:nst, 0:nst] + na1y[0:nst, 0:nst]
    kny_inv_fixed = np.linalg.inv(kny_fixed)

    # INITIAL TIME STEP CALCULATION FOR ACCELERATION
    px1 = smx_diag*xg[0]
    ax1[0:nst, 0] = np.dot(smx_inv_fixed, px1[0:nst, 0] - np.dot(cdx[0:nst, 0:nst], vx1[0:nst, 0]) - np.dot(skx[0:nst, 0:nst], dx1[0:nst, 0]))
    py1 = smy_diag*yg[0]
    ay1[0:nst, 0] = np.dot(smy_inv_fixed, py1[0:nst, 0] - np.dot(cdy[0:nst, 0:nst], vy1[0:nst,0]) - np.dot(sky[0:nst, 0:nst], dy1[0:nst,0]))
    
    
    # INITIAL PRIMARY RESPONSE VARIABLE ASSIGNMENTS

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

    rolling_state = False # INITIAL ROLLING STATE

    
    for i in range(1, len(xg)):

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
            if drb > iso.umax:
                print(prRed('WARNING: ') + prCyan('Isolator has crossed maximum displacement.'))
                print(drb, dy2[ndof-1, 0], dx2[ndof-1, 0])
                print("System terminating simulation")
                sys.exit()
            
            theta_d02 = np.interp(drb, iso.xbtab, iso.ttab)
            theta_r_d02 = atan((iso.b0/iso.a0)*tan(theta_d02))
            c_d02 = iso.a0*sin(theta_d02)*cos(theta_r_d02) - iso.b0*cos(theta_d02)*sin(theta_r_d02)
            p_d02 = iso.a0*sin(theta_d02)*sin(theta_r_d02) + iso.b0*cos(theta_d02)*cos(theta_r_d02) 
            
            
            f_axial, fcx, fcy, fcz, dc2, vc2, e_axial = fc(iso, dt, dy2[ndof-1, 0], dx2[ndof-1, 0], drb, dc1, p_d02, tm)
           
            fs1x, fs1y = fs1fixed(tm, fcz, rMrM, c_d02, p_d02, drb, dy2[ndof-1, 0], dx2[ndof-1, 0]) #CB

            fbx, fby = fbfixed(iso.Mr, yg[i], xg[i]) #CB

            fabx = (-1.0*cdx[nst, nst-1])*vx2[nst-1,0] + (-1.0*skx[nst, nst-1])*dx2[nst-1,0] - skx[ndof-1, ndof-1]*dx2[ndof-1,0] + px2[ndof-1] - fs1x - fbx - fcx
            faby = (-1.0*cdy[nst, nst-1])*vy2[nst-1,0] + (-1.0*sky[nst, nst-1])*dy2[nst-1,0] - sky[ndof-1, ndof-1]*dy2[ndof-1,0] + py2[ndof-1] - fs1y - fby - fcy

            mu = mu_val(iso, drb)  #CB

            fs2x, fs2y = fs2fixed(mu, tm, fcz, rMrM) #CB

            rolling_state, fabx, faby = stat1(fabx, faby, fs2x, fs2y)

            if rolling_state == False:
                idx = 11 # del
                ax2[0:nst, 0] = np.dot(smx_inv_fixed, px2[0:nst, 0] - np.dot(cdx[0:nst, 0:nst], vx2[0:nst, 0]) - np.dot(skx[0:nst, 0:nst], dx2[0:nst, 0]))
                ay2[0:nst, 0] = np.dot(smy_inv_fixed, py2[0:nst, 0] - np.dot(cdy[0:nst, 0:nst], vy2[0:nst, 0]) - np.dot(sky[0:nst, 0:nst], dy2[0:nst, 0]))
            else:
                idx = 12 # del
                epx = px2 - np.dot(cdx, vx2) - np.dot(skx, dx2)
                epy = py2 - np.dot(cdy, vy2) - np.dot(sky, dy2)

                epx[ndof-1,0] = epx[ndof-1,0] - fabx - fs1x - fbx - fcx    
                epy[ndof-1,0] = epy[ndof-1,0] - faby - fs1y - fby - fcy

                ax2 = np.dot(smx_inv, epx)
                ay2 = np.dot(smy_inv, epy)
            
            dx1, vx1, px1, ax1 = dx2, vx2, px2, ax2 
            dy1, vy1, py1, ay1 = dy2, vy2, py2, ay2
            yd1, yv1, ya1 = yd2, yv2, ya2
            theta_r_d01, theta_r_dot11, theta_r_dot21 = theta_r_d02, theta_r_dot12, theta_r_dot22
            dc1, vc1 = dc2, vc2

        else:
            
            dfbx = 0.0
            dfby = 0.0

            dfs1x = 0.0
            dfs1y = 0.0

            dfcx = 0.0
            dfcy = 0.0

            dfabx = 0.0
            dfaby = 0.0

            for ii in range (nit):
                  
                pcx1 = dpx + np.dot(na2x, vx1) + np.dot(na3x, ax1)
                pcx1[ndof-1,0] = pcx1[ndof-1,0] - dfabx - dfs1x - dfbx - dfcx
                ddx = np.dot(knx_inv, pcx1)
                dx2 = dx1 + ddx
                dvx = (gamma/beta/dt)*ddx - gamma/beta*vx1 + dt*(1.0 - gamma/2.0/beta)*ax1
                vx2 = vx1 + dvx

                pcy1 = dpy + np.dot(na2y, vy1) + np.dot(na3y, ay1)
                pcy1[ndof-1,0] = pcy1[ndof-1,0] - dfaby - dfs1y - dfby - dfcy
                ddy = np.dot(kny_inv, pcy1)
                dy2 = dy1 + ddy
                dvy = (gamma/beta/dt)*ddy - gamma/beta*vy1 + dt*(1.0 - gamma/2.0/beta)*ay1
                vy2 = vy1 + dvy

                vrb = sqrt(pow(vx2[ndof-1, 0], 2.0) + pow(vy2[ndof-1, 0], 2.0))

                if vrb > 1.0e-5:

                    daxb = (1/smx[ndof-1, ndof-1])*(dpx[ndof-1, 0] - cdx[ndof-1, ndof-1]*dvx[ndof-1, 0] - skx[ndof-1, ndof-1]*ddx[ndof-1, 0] + (-1.0*cdx[nst, nst-1])*dvx[nst-1, 0] + (-1.0*skx[nst, nst-1]*ddx[nst-1, 0]) - dfs1x - dfabx - dfbx - dfcx)
                    dayb = (1/smy[ndof-1, ndof-1])*(dpy[ndof-1, 0] - cdy[ndof-1, ndof-1]*dvy[ndof-1, 0] - sky[ndof-1, ndof-1]*ddy[ndof-1, 0] + (-1.0*cdy[nst, nst-1])*dvy[nst-1, 0] + (-1.0*sky[nst, nst-1]*ddy[nst-1, 0]) - dfs1y - dfaby - dfby - dfcy)

                    ax2[ndof-1, 0] = ax1[ndof-1, 0] + daxb
                    ay2[ndof-1, 0] = ay1[ndof-1, 0] + dayb

                    drb = sqrt(pow(dx2[ndof-1, 0], 2.0) + pow(dy2[ndof-1, 0], 2.0))
                    arb = sqrt(pow(ax2[ndof-1, 0], 2.0) + pow(ay2[ndof-1, 0], 2.0))
                    arg = sqrt(pow(xg[i], 2.0) + pow(yg[i], 2.0))

                    theta_d02 = np.interp(drb, iso.xbtab, iso.ttab)
                    theta_r_d02 = atan((iso.b0/iso.a0)*tan(theta_d02))
                    c_d02 = iso.a0*sin(theta_d02)*cos(theta_r_d02) - iso.b0*cos(theta_d02)*sin(theta_r_d02)
                    p_d02 = iso.a0*sin(theta_d02)*sin(theta_r_d02) + iso.b0*cos(theta_d02)*cos(theta_r_d02)
                    
                    # dyb_dt0 = 2.0*(p_d02 - iso.b0) - yd1
                    # dyb_dt1 = (gamma/beta/dt)*dyb_dt0 - gamma/beta*yv1 + dt*(1.0 - gamma/2.0/beta)*ya1
                    # dyb_dt2 = (1.0/beta/pow(dt, 2.0))*dyb_dt0 - (1.0/beta/dt)*yv1 - (1.0/2.0/beta)*ya1
                    # yd2 = yd1 + dyb_dt0
                    # yv2 = yv1 + dyb_dt1
                    # ya2 = ya1 + dyb_dt2

                    dyb_dt0 = 2.0*(p_d02 - iso.b0) - yd1
                    yd2 = yd1 + dyb_dt0
                    yv2 = (yd2 - yd1)/dt
                    ya2 = (yv2 - yv1)/dt

                    # dtheta_r_d0_dt0 = theta_r_d02 - theta_r_d01
                    # dtheta_r_d0_dt1 = (gamma/beta/dt)*dtheta_r_d0_dt0 - gamma/beta*theta_r_dot11 + dt*(1.0 - gamma/2.0/beta)*theta_r_dot21
                    # dtheta_r_d0_dt2 = (1.0/beta/pow(dt, 2.0))*dtheta_r_d0_dt0 - (1.0/beta/dt)*theta_r_dot11 - (1.0/2.0/beta)*theta_r_dot21
                    # theta_r_d02 = theta_r_d01 + dtheta_r_d0_dt0
                    # theta_r_dot12 = theta_r_dot11 + dtheta_r_d0_dt1
                    # theta_r_dot22 = theta_r_dot21 + dtheta_r_d0_dt2

                    theta_r_dot12 = (theta_r_d02 - theta_r_d01)/dt
                    theta_r_dot22 = (theta_r_dot12 - theta_r_dot11)/dt
                    
                    f_axial, fcx2, fcy2, fcz2, dc2, vc2, e_axial = fc(iso, dt, dy2[ndof-1, 0], dx2[ndof-1, 0], drb, dc1, p_d02, tm)

                    fs1x2, fs1y2 = fs1(tm, fcz2, rMrM, ya2, c_d02, p_d02, drb, dy2[ndof-1, 0], dx2[ndof-1, 0])
                    
                    mu = mu_val(iso, drb)
                    fs2x2, fs2y2 = fs2(tm, fcz2, rMrM, mu, ya2)
                    fbx2, fby2 = fb(iso.Jr, iso.Mr, drb, dy2[ndof-1, 0], dx2[ndof-1, 0], vrb, vy2[ndof-1, 0], vx2[ndof-1, 0], arg, yg[i], xg[i], arb, ay2[ndof-1, 0], ax2[ndof-1, 0], p_d02, theta_r_dot22, c_d02, ya2)

                    dfs1x = fs1x2 - fs1x
                    dfabx = fs2x2*vx2[ndof-1, 0]/vrb - fabx
                    dfbx = fbx2 - fbx
                    dfcx = fcx2 - fcx

                    dfs1y = fs1y2 - fs1y
                    dfaby = fs2y2*vy2[ndof-1, 0]/vrb - faby
                    dfby = fby2 - fby
                    dfcy = fcy2 - fcy

                    # print(fs1x2, fs2x2*vx2[ndof-1, 0]/vrb, fbx2, fcx2)
            
            # print('\n')
            if nit > 1:

                fbx = fbx + dfbx
                fby = fby + dfby

                fs1x = fs1x + dfs1x
                fs1y = fs1y + dfs1y

                fabx = fabx + dfabx
                faby = faby + dfaby

                fcx = fcx + dfcx
                fcy = fcy + dfcy
            
            rolling_state = stat0(ddx[ndof-1], ddy[ndof-1], fabx, faby)

            if rolling_state == False:
                drb = sqrt(pow(dy2[ndof-1, 0], 2.0) + pow(dx2[ndof-1, 0], 2.0))
                if drb > iso.umax:
                    print(prRed('WARNING: ') + prCyan('Isolator has crossed maximum displacement.'))
                    print(drb, dy2[ndof-1, 0], dx2[ndof-1, 0])
                    print("System terminating simulation")
                    sys.exit()

                idx = 21 # del
                vx2[ndof-1, 0] = 1e-10
                ax2[ndof-1, 0] = 0.0
                vy2[ndof-1, 0] = 1e-10
                ay2[ndof-1, 0] = 0.0
                yv2 = 0.0
                ya2 = 0.0
                theta_r_dot12 = 0.0
                theta_r_dot22 = 0.0
                vc2 = 0.0

                ax2[0:nst, 0] = np.dot(smx_inv_fixed, px2[0:nst, 0] - np.dot(cdx[0:nst, 0:nst], vx2[0:nst, 0]) - np.dot(skx[0:nst, 0:nst], dx2[0:nst, 0]))
                ay2[0:nst, 0] = np.dot(smy_inv_fixed, py2[0:nst, 0] - np.dot(cdy[0:nst, 0:nst], vy2[0:nst, 0]) - np.dot(sky[0:nst, 0:nst], dy2[0:nst, 0]))
            else:

                drb = sqrt(pow(dy2[ndof-1, 0], 2.0) + pow(dx2[ndof-1, 0], 2.0))
                if drb > iso.umax:
                    print(prRed('WARNING: ') + prCyan('Isolator has crossed maximum displacement.'))
                    print(drb, dy2[ndof-1, 0], dx2[ndof-1, 0])
                    print("System terminating simulation")
                    sys.exit()

                idx = 22 # del
                epx = px2 - np.dot(cdx, vx2) - np.dot(skx, dx2)
                epy = py2 - np.dot(cdy, vy2) - np.dot(sky, dy2)
                epx[ndof-1,0] = epx[ndof-1,0] - fabx - fs1x - fbx - fcx
                epy[ndof-1,0] = epy[ndof-1,0] - faby - fs1y - fby - fcy

                ax2 = np.dot(smx_inv, epx)
                ay2 = np.dot(smy_inv, epy)
            
            dx1, vx1, px1, ax1 = dx2, vx2, px2, ax2
            dy1, vy1, py1, ay1 = dy2, vy2, py2, ay2
            yd1, yv1, ya1 = yd2, yv2, ya2
            theta_r_d01, theta_r_dot11, theta_r_dot21 = theta_r_d02, theta_r_dot12, theta_r_dot22
            dc1, vc1 = dc2, vc2
        
        eki = eki + 0.5*dt*(np.dot(np.dot(vx2.T, smx),ax2) + np.dot(np.dot(vx1.T, smx),ax1)) + 0.5*dt*(np.dot(np.dot(vy2.T, smy),ay2) + np.dot(np.dot(vy1.T, smy),ay1))
        edi = edi + 0.5*dt*(np.dot(np.dot(vx2.T, cdx),vx2) + np.dot(np.dot(vx1.T, cdx),vx1)) + 0.5*dt*(np.dot(np.dot(vy2.T, cdy),vy2) + np.dot(np.dot(vy1.T, cdy),vy1))
        esi = esi + 0.5*dt*(np.dot(np.dot(vx2.T, skx),dx2) + np.dot(np.dot(vx1.T, skx),dx1)) + 0.5*dt*(np.dot(np.dot(vy2.T, sky),dy2) + np.dot(np.dot(vy1.T, sky),dy1))
        eii = eii - 0.5*dt*(np.dot(np.dot(vx2.T, smxi),np.dot(r, xg[i])) + np.dot(np.dot(vx1.T, smxi),np.dot(r, xg[i-1]))) - 0.5*dt*(np.dot(np.dot(vy2.T, smyi),np.dot(r, yg[i])) + np.dot(np.dot(vy1.T, smyi),np.dot(r, yg[i-1])))
        
        dx1, vx1, px1, ax1 = dx2, vx2, px2, ax2 
        dy1, vy1, py1, ay1 = dy2, vy2, py2, ay2
        yd1, yv1, ya1 = yd2, yv2, ya2
        theta_r_d01, theta_r_dot11, theta_r_dot21 = theta_r_d02, theta_r_dot12, theta_r_dot22 
        dc1, vc1 = dc2, vc2
         
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
            fx[index, ndof-1] = fabx + fs1x + fbx
            fy[index, ndof-1] = faby + fs1y + fby
            ek[index, 0] = eki
            ed[index, 0] = edi
            es[index, 0] = esi
            ei[index, 0] = eii

            roll[index, 0] = idx # del
            theta_0[index, 0] = theta_d02
            theta_r[index, 0] = theta_r_d02
            theta_r_dot2[index, 0] = theta_r_dot22
            zbd[index, 0] = yd2
            zbd_dot2[index, 0] = ya2
            Fs1x[index, 0] = fs1x
            Fs1y[index, 0] = fs1y
            Fs2x[index, 0] = fabx
            Fs2y[index, 0] = faby 
            Fbx[index, 0] = fbx 
            Fby[index, 0] = fby
            F_axial[index, 0] = f_axial
            Dc_axial[index, 0] = dc2
            Strain_axial[index, 0] = e_axial 

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
    
    result = ResultFixedXY(ref, ijk, time.T, gx.T, dx.T, vx.T, ax.T, aax.T, gy.T, dy.T, vy.T, ay.T, aay.T, fx, fy, ek, ed, es, ei, error, smx, skx, cdx, smy, sky, cdy, roll, theta_0, theta_r, theta_r_dot2, zbd, zbd_dot2, Fs1x, Fs1y, Fs2x, Fs2y, Fbx, Fby, F_axial, Dc_axial, Strain_axial)
    model = ModelInfo(ndof)

    return result, model