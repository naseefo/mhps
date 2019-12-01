
import numpy as np
from scipy.optimize import fsolve
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

class IsoOSBITorsionModel: 

    def k0fn(e, *data):
        D, ro2, k0i = data
        fn_value = k0i - (1+0.5*ro2)*tm*9.81/D*(pow(e, 2.0)/pow((1-pow(e, 2.0)), 3.0/2.0))
        return fn_value

    def t0fn(e, *data):
        # print(data)
        i, tbe0x, D, ro2 = data
        # print(tbe0x, D, ro2)
        # print(tbe0x - sqrt(D)*sqrt(4.0*pow(pi, 2.0)/((1 + ro2)*9.81*(pow(e, 2.0)/pow((1-pow(e, 2.0)), 3.0/2.0))))  )
        fn_value = tbe0x - sqrt(D)*sqrt(4.0*pow(pi, 2.0)/((1 + ro2)*9.81*(pow(e, 2.0)/pow((1-pow(e, 2.0)), 3.0/2.0))))  
        return fn_value
        
    def calc_inital_stiffness(self):
        ecc0 = fsolve(self.t0fn, 0.5, args=(self.tbe0x, self.D, self.ro2))
        tm = (self.am + self.am*self.rmbm)/4.0
        k0 = (1 + 0.5*self.ro2)*tm*9.81/self.D*(pow(e, 2.0)/pow((1-pow(e, 2.0)), 3.0/2.0))
        return k0

    def calc_ecc(self, k0i):
        data = (self.D, self.ro2, k0i)
        e = fsolve(self.k0fn, 0.5, args=data)
        return e


    def osbt2i(self, ecc, a0, b0): 
        default_values = get_default_param_values()
        tdiv = default_values['TDIV']
        
        ld0 = 0.0
        ud0 = pi/2
        delta = ud0/tdiv
        sizeTable = int((ud0 - ld0)/delta + 1)
        xbtab = np.zeros((sizeTable, 1), dtype=np.dtype('d'), order ='F')
        ttab = np.zeros((sizeTable, 1), dtype=np.dtype('d'), order ='F')
        
        fun1 = lambda x: sqrt(1-pow(ecc*sin(x), 2.0))
        i = 0
        for td0 in np.arange(ld0, ud0 + delta, delta):
            ttab[i, 0] = td0
            trd0 = atan((b0/a0)*tan(td0))
            cd0 = a0*sin(td0)*cos(trd0) - b0*cos(td0)*sin(trd0)
            itheta = integrate.quad(fun1, 0.0, td0)
            xbtab[i, 0] = 2.0*a0*itheta[0] - 2.0*cd0
            i += 1
        return xbtab, ttab


    def u_max(self, ecc, a0, b0):
        fun1 = lambda x: sqrt(1-pow(ecc*sin(x), 2.0))
        trd0 = pi/2
        td0 = atan((a0/b0)*tan(trd0))
        cd0 = a0*sin(td0)*cos(trd0) - b0*cos(td0)*sin(trd0)
        itheta = integrate.quad(fun1, 0.0, td0)
        umax = 2.0*a0*itheta[0] - 2.0*cd0
        return umax

    def __init__(self, rmbm, tbx, zetabx, ebxd, wrwxb, efxd, f0, g1, bt, g2, a, nt, D, ecc, rmrm, am, niso, tbe0x, e0exd, mu0exd):
        self.rmbm = rmbm
        self.tbx = tbx
        self.tbe0x = tbe0x
        self.zetabx = zetabx
        self.ebxd = ebxd
        self.wrwxb = wrwxb
        self.efxd = efxd

        self.f0 = f0
        self.g1 = g1
        self.bt = bt
        self.g2 = g2
        self.a = a
        self.nt = nt
        self.D = D
        self.ecc = ecc
        self.rmrm = rmrm
        self.am = am
        self.niso = niso
        self.tbe0x = tbe0x
        self.e0exd = e0exd
        self.mu0exd = mu0exd

        self.ro2 = self.rmrm/2.0
        self.tm = (self.am + self.am*self.rmbm)
        
        S2 = lambda e: pow(e, 2.0)/pow((1 - pow(e, 2.0)), 1.5)
        
        # calc_ecc0 = lambda e: self.tbe0x - sqrt(self.D)*sqrt(4.0*pow(pi, 2.0)/((1 + self.ro2)*9.81*(pow(e, 2.0)/pow((1-pow(e, 2.0)), 3.0/2.0)))) 

        calc_ecc0_tb0 = lambda e: self.tbe0x - 2.0*pi*sqrt(self.D/(1+0.5*self.ro2)/9.81/S2(e))

        


        self.ecc0 = fsolve(calc_ecc0_tb0, 0.2)
        self.val = 2.0*pi*sqrt(self.D/(1+0.5*self.ro2)/9.81/S2(self.ecc0))

        self.k0 = (1 + 0.5*self.ro2)*self.tm*9.81/self.D*(pow(self.ecc0, 2.0)/pow((1-pow(self.ecc0, 2.0)), 3.0/2.0))

        self.k0i = np.zeros((4,),dtype=np.dtype('d'), order='F')
        self.ecc_i = np.zeros((4,),dtype=np.dtype('d'), order='F')
        self.mr = np.zeros((4,),dtype=np.dtype('d'), order='F')
        self.a0 = np.zeros((4,),dtype=np.dtype('d'), order='F')
        self.b0 = np.zeros((4,),dtype=np.dtype('d'), order='F')
        self.umax = np.zeros((4,),dtype=np.dtype('d'), order='F')
        default_values = get_default_param_values()
        tdiv = default_values['TDIV']
        ld0 = 0.0
        ud0 = pi/2
        delta = ud0/tdiv
        sizeTable = int((ud0 - ld0)/delta + 1)
        self.xbtab = np.zeros((sizeTable, 4), dtype=np.dtype('d'), order ='F')
        self.ttab = np.zeros((sizeTable, 4), dtype=np.dtype('d'), order ='F')

        self.k0i[0] = self.k0*0.25*(1 + 2.0*self.e0exd)
        self.k0i[1] = self.k0*0.25*(1 - 2.0*self.e0exd)
        self.k0i[2] = self.k0*0.25*(1 - 2.0*self.e0exd)
        self.k0i[3] = self.k0*0.25*(1 + 2.0*self.e0exd)

        calc_ecc0_k0i = lambda e, k0val: k0val - self.tm/4.0/self.D*(1+0.5*self.ro2)*9.81*S2(e)
        
        for i in range(4):
            self.mr[i] = self.rmrm*am/4.0
            self.ecc_i[i] = fsolve(calc_ecc0_k0i, 0.5, args=self.k0i[i])
            self.a0[i] = self.D/2.0
            self.b0[i] = self.a0[i]*sqrt(1 - pow(self.ecc_i[i], 2.0))
            # print(self.ecc_i[i], self.a0[i], self.b0[i])
            temp1, temp2 = self.osbt2i(self.ecc_i[i], self.a0[i], self.b0[i])
            # print(temp1.shape, temp2.shape, self.xbtab[:,i].shape)
            self.xbtab[:,i] = temp1.reshape(sizeTable,)
            self.ttab[:, i] = temp2.reshape(sizeTable,)
            
            self.umax[i] = self.u_max(self.ecc_i[i], self.a0[i], self.b0[i])

        
        # self.ecctb0 = self.ecc_calc_tbe0x()

        # self.a0 = D/2.0
        # self.b0 = sqrt((D/2.0)**2.0 - (ecc*D/2.0)**2.0)
        # self.Mr = am*rmrm
        
        # self.Jr = (self.Mr/5.0)*(pow(self.a0, 2.0) + pow(self.b0, 2.0))
        # self.jr = self.Jr/niso
        # self.umax = self.u_max()
        # self.xbtab, self.ttab = self.osbt2x()

        # self.V = 4.0/3.0*pi*pow(self.a0, 2.0)*self.b0

        # self.ro = self.mr/self.V
    
    def __str__(self):
        print("Isolator Initial Stiffness")
        print(self.k0i)
        print("Isolator oblateness ratio")
        print(self.ecc_i)
        print("Maximum Displacement")
        print(self.umax)
        print("a0")
        print(self.a0)
        print("b0")
        print(self.b0)
        print("mr")
        print(self.mr)
        print("xbtab")
        print(self.xbtab)
        print("ttab")
        print(self.ttab)
        line = '\nNumber of isolators (niso) = ' + prCyan('%d ')%(self.niso)
        line += '\nPlan diameter of OSBI (D) = ' + prCyan('%8.4f cm')%(self.D*100.0)
        # line += '\nMajor radius of OSBI (a0) = ' + prCyan('%8.4f cm')%(self.a0*100.0)
        # line += '\nMinor radius of OSBI (b0) = ' + prCyan('%8.4f cm')%(self.b0*100.0)
        # line += '\nEccentricity of OSBI (ecc) = ' + prCyan('%8.4f')%(self.ecc)
        line += '\nRatio of mass of all OSBI to floor mass (rmrm) = ' + prCyan('%8.4f')%(self.rmrm)
        # line += '\nTotal mass of all OSBI (Mr) = ' + prCyan('%8.4f kg')%(self.Mr)
        # line += '\nMass of single OSBI (mr) = ' + prCyan('%8.4f kg')%(self.mr)
        # line += '\nVolume of single OSBI (V) = ' + prCyan('%8.4f m3'%(self.V))
        # line += '\nDensity of material required for OSBI (ro) = ' + prCyan('%8.4f kg/m3')%(self.ro)
        # line += '\nMass Moment of inertia of all OSBI (Jr) = ' + prCyan('%8.4f kg-m2')%(self.Jr)
        # line += '\nMass Moment of inertia inertia of single OSBI (jr) = ' + prCyan('%8.4f kg-m2')%(self.jr)
        # line += '\nMaximum displacement of OSBI (umax) = ' + prCyan('%8.4f cm')%(self.umax*100.0)
        # line += '\nRolling friction coefficient (mu0) = ' + prCyan('%8.4f')%(self.f0)
        # line += '\n\nTime period of linear spring in X-direction (tbx) = ' + prCyan('%8.4f sec')%(self.tbx)
        # line += '\nCritical viscous damping ratio in X-direction (zetabx) = ' + prCyan('%8.4f p.c.')%(self.zetabx)
        line += '\n'
        return line

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.rmbm == other.rmbm and self.tbx == other.tbx and \
                self.zetabx == other.zetabx and self.ebxd == other.ebxd and  self.wrwxb == other.wrwxb and \
                    self.efxd == other.efxd and self.f0 == other.f0 and self.g1 == other.g1 and self.bt == other.bt and \
                        self.g2 == other.g2 and self.a == other.a and self.nt == other.nt and self.D == other.D and self.ecc == other.ecc and \
                            self.rmrm == other.rmrm and self.am == am and self.niso == other.niso and self.tbe0x == tbe0x and self.e0exd == other.e0exd and self.mu0exd == mu0exd

                            
        return False

def read_iso_osbi_torsion_var_param(var_param_file, fm, niso):

    """
    This function reads all the set of parameters for the parametric
    studies and stores it in an array. It will be a one-time allocation
    to increase speed.
    """

    dtype1 = np.dtype([('IJK', 'i4'), ('RMBM', 'd'), ('TBX', 'd'), ('ZETABX','d'), ('EBXD','d'), ('WRWXB', 'd'), ('efxd', 'd'), ('F0', 'd'), ('G1', 'd'), ('BT', 'd'), ('G2', 'd'), ('A', 'd'), ('NT', 'd'), ('D', 'd'), ('ECC', 'd'), ('RMRM', 'd'), ('TBE0XD', 'd'), ('E0EXD', 'd'), ('MU0EXD', 'd') ])
    
    ijk, rmbm, tbx, zetabx, ebxd, wrwxb, efxd, f0, g1, bt, g2, a, nt, D, ecc, rmrm, tbe0x, e0exd, mu0exd  = np.loadtxt \
    (var_param_file, delimiter=',', usecols=range(19), skiprows=1, unpack = True, dtype=dtype1)

    try:
        total_param = len(ijk)
    except:
        total_param = 1


    for i in range(0, total_param):
        try:
            yield IsoOSBITorsionModel(rmbm[i], tbx[i], zetabx[i], ebxd[i], wrwxb[i], efxd[i], f0[i], g1[i], bt[i], g2[i], a[i], nt[i], D[i], ecc[i], rmrm[i], fm, niso, tbe0x[i], e0exd[i], mu0exd[i])
        except:
            yield IsoOSBITorsionModel(rmbm, tbx, zetabx, ebxd, wrwxb, efxd, f0, g1, bt, g2, a, nt, D, ecc, rmrm, fm, niso, tbe0x, e0exd, mu0exd)


# superstructure_propxy_t
# IJK, TX1, ZETA, EXD, WRWX
# ijk, tx1, zeta, exd, wrwx, fm, nb, x, y, xb, yb





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

def fs1(dt, iso, idb, tm, dbcx2, dbcy2, yd1, yv1, fs1x1, fs1y1):
    # idb stands for the index of the isolator being passed

    g = 9.81 # in %m/s2

    drb = sqrt(pow(dbcx2, 2.0) + pow(dbcy2, 2.0))

    theta_d02 = np.interp(drb, iso.xbtab[:, idb], iso.ttab[:, idb])
    theta_r_d02 = atan((iso.b0[idb]/iso.a0[idb])*tan(theta_d02))

    c_d02 = iso.a0[idb]*sin(theta_d02)*cos(theta_r_d02) - iso.b0[idb]*cos(theta_d02)*sin(theta_r_d02)
    p_d02 = iso.a0[idb]*sin(theta_d02)*sin(theta_r_d02) + iso.b0[idb]*cos(theta_d02)*cos(theta_r_d02)

    # dyb_dt0 = 2.0*(p_d02 - iso.b0[idb]) - yd1
    # yd2 = yd1 + dyb_dt0
    # yv2 = (yd2 - yd1)/dt
    # ya2 = (yv2 - yv1)/dt
    yd2=0
    yv2=0
    ya2=0
    fs1r = 0.5*tm*((1 + 2.0*iso.mr[idb]/tm)*g + ya2)*(c_d02/2.0/p_d02)
    # fs1r = ((2 + rMrM)*M*g + (2 + 0)*M*y_b_d2)*(c_d0/2.0/p_d0) + fcz*c_d0/p_d0

    if drb == 0:
        fs1x2 = 0.0
        fs1y2 = 0.0
    else:
        fs1x2 = fs1r*dbcx2/drb
        fs1y2 = fs1r*dbcy2/drb
    
    dfs1x = fs1x2 - fs1x1
    dfs1y = fs1y2 - fs1y1

    # yd2, yv2 are passed to store it in the main subrouting to be used here again as yd1
    # ya2 is passed as it will be required in the wen's subroutine in the next step in the main subroutine
    
    return dfs1x, dfs1y, yd2, yv2, ya2

def wen(iso, vx, vy, zx, zy, dt, alpx, alpy, fyx, fyy):
    
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
def simulator_osbi_tor(ref, xg, yg, dt, ndiv, ndt, ijk, sm, sk, cd, x, y, xb, yb, nb, iso):

    gamma = 0.5
    beta = 1/6
    nit = 1

    xb = np.asarray(xb, dtype=np.dtype('d'), order='F')
    yb = np.asarray(yb, dtype=np.dtype('d'), order='F')
    
    wbx = 2.0*pi/iso.tbx
    fm = sm[0,0]
    bm = sm[0,0]*iso.rmbm
    tm = fm + bm

    ckabx = (fm + bm)*wbx*wbx
    
    bkx = np.zeros(shape=(4, ), dtype=np.dtype('d'), order='F')
    bky = np.zeros(shape=(4, ), dtype=np.dtype('d'), order='F')

    bkx[0] = 0.25*ckabx*(1.0 + 2.0*iso.ebxd)
    bkx[1] = 0.25*ckabx*(1.0 - 2.0*iso.ebxd)
    bkx[2] = 0.25*ckabx*(1.0 - 2.0*iso.ebxd)
    bkx[3] = 0.25*ckabx*(1.0 + 2.0*iso.ebxd)

    bky = bkx

    fyx = np.zeros(shape=(4, ), dtype=np.dtype('d'), order='F')
    fyy = np.zeros(shape=(4, ), dtype=np.dtype('d'), order='F')

    qyf = iso.f0*(fm + bm)*9.81
    alp = iso.g1*ckabx/qyf
    print(iso.rmbm, fm, bm, wbx, ckabx, iso.g1, qyf, alp)

    fyx[0] = qyf*0.25*(1 + 2.0*iso.mu0exd)
    fyx[1] = qyf*0.25*(1 - 2.0*iso.mu0exd)
    fyx[2] = qyf*0.25*(1 - 2.0*iso.mu0exd)
    fyx[3] = qyf*0.25*(1 + 2.0*iso.mu0exd)


    print("f0")
    print(iso.f0)
    print("g1")
    print(iso.g1)
    print("qyf")
    print(qyf)
    print("alp")
    print(alp)
    print(fyx)
    
    print(iso)

    


    sm_inv = np.linalg.inv(sm)
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
    a1 = np.dot(sm_inv, p1 - np.dot(cd, v1) - np.dot(sk, d1))

    d2 = np.zeros((6, 1), dtype=np.dtype('d'), order='F')
    v2 = np.zeros((6, 1), dtype=np.dtype('d'), order='F')
    p2 = np.zeros((6, 1), dtype=np.dtype('d'), order='F')
    a2 = np.zeros((6, 1), dtype=np.dtype('d'), order='F')

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
    iwbx = 2.0*pi/iso.tbe0x
    iwrb = iwbx*iso.wrwxb
    ikttb = np.sum(iso.k0i*np.square(yb)) + np.sum(iso.k0i*np.square(xb))
    bmr = ikttb/pow(iwrb, 2.0)
    sm[5,5] = bmr
    print(sm)
    print("Linear-Stiffness Matrix")
    print(sk)
    print("Linear-Damping Matrix")
    print(cd)

    pzx = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    pzy = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    zx = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    zy = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    dzx = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    dzy = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    dpz = np.zeros((3, 1), dtype=np.dtype('d'), order='F')
    pz = np.zeros((3, 1), dtype=np.dtype('d'), order='F')
    dpzx = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    dpzy = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    
    dbcx2 = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    dbcy2 = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    vbcx2 = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    vbcy2 = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    abcx2 = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    abcy2 = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    aabcx2 = np.zeros((4, ), dtype=np.dtype('d'), order='F')
    aabcy2 = np.zeros((4, ), dtype=np.dtype('d'), order='F')

    fabx = np.zeros((4, 1), dtype=np.dtype('d'), order='F')
    faby = np.zeros((4, 1), dtype=np.dtype('d'), order='F')

    

    yd1 = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added
    yv1 = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added
    ya1 = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added

    yd2 = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added
    yv2 = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added
    ya2 = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added

    fs1x1 = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added
    fs1y1 = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added
    pfs1 = np.zeros((3, 1), dtype=np.dtype('d'), order='F') # added
    dpfs1x = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added
    dpfs1y = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added
    dpfs1 = np.zeros((3, 1), dtype=np.dtype('d'), order='F') #added
    pfs1x = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added
    pfs1y = np.zeros((4, ), dtype=np.dtype('d'), order='F') #added

    
    for i in range(1,len(xg)):

        t += dt
        
        
        ug2 = np.array([[xg[i]],[yg[i]], [0.0]])
        p2 = -1.0*np.dot(np.dot(sm, r), ug2)
        dp = p2 - p1

        for wc in range(4):
            dpzx[wc] = 0.0
            dpzy[wc] = 0.0

            dpfs1x[wc] = 0.0 #added
            dpfs1y[wc] = 0.0 #added

        
        for bc in range(3):
            dpz[bc] = 0.0
            dpfs1[bc] = 0.0 # added

        pzt = pz
        pfs1t = pfs1 # added
        for i2 in range(nit):
            pcx1 = dp + np.dot(na2x, v1) + np.dot(na3x, a1)
            pcx1[3:6,0] = pcx1[3:6,0] - dpz.T - dpfs1.T  # added
            dd = np.dot(knx_inv, pcx1)
            d2 = d1 + dd
            dv = (gamma/beta/dt)*dd - gamma/beta*v1 + dt*(1.0 - gamma/2.0/beta)*a1
            v2 = v1 + dv

            ep = p2 - np.dot(cd, v2) - np.dot(sk, d2)
            ep[3:6,0] = ep[3:6,0] - pzt.T - pfs1t.T  # added
            a22 = np.dot(sm_inv, ep)
            
            # CALCULATION OF FORCES FROM BEARING
            dpz1, dpz2, dpz3 = 0.0, 0.0, 0.0
            dpfs11, dpfs12, dpfs13 = 0.0, 0.0, 0.0 # added
            for bc in range(4):
                dbcx2[bc] = d2[3,0] - yb[bc]*d2[5,0]
                dbcy2[bc] = d2[4,0] + xb[bc]*d2[5,0]

                vbcx2[bc] = v2[3,0] - yb[bc]*v2[5,0]
                vbcy2[bc] = v2[4,0] + xb[bc]*v2[5,0]

                abcx2[bc] = a22[3,0] - yb[bc]*a22[5,0]
                abcy2[bc] = a22[4,0] + xb[bc]*a22[5,0]

                aabcx2[bc] = a22[3,0] + xg[i] - yb[bc]*(a22[5,0] + 0.0)
                aabcy2[bc] = a22[4,0] + yg[i] + xb[bc]*(a22[5,0] + 0.0)

                # added
                dpfs1x[bc], dpfs1y[bc], yd2[bc], yv2[bc], ya2[bc] = fs1(dt, iso, bc, tm, dbcx2[bc], dbcy2[bc], yd1[bc], yv1[bc], fs1x1[bc], fs1y1[bc])
                dpfs11 = dpfs11 + dpfs1x[bc]
                dpfs12 = dpfs12 + dpfs1y[bc]
                dpfs13 = dpfs13 + (-1.0*dpfs1x[bc]*yb[bc] + dpfs1y[bc]*xb[bc])
                # added

                # added
                # print("velocity")
                # print(ya2[bc])
                fyx[bc] = fyx[bc]*((1 + 2.0*iso.mr[bc]/tm) + (1 + iso.mr[bc]/tm)*ya2[bc]/9.81)
                fyy[bc] = fyx[bc]
                # added
                # print(fyx[bc])
                # print(vbcx2[bc], vbcy2[bc], zx[bc], zy[bc], dt, alp, alp, fyx[bc], fyy[bc])
                dpzx[bc], dpzy[bc], dzx[bc], dzy[bc] = wen(iso, vbcx2[bc], vbcy2[bc], zx[bc], zy[bc], dt, alp, alp, fyx[bc], fyy[bc])
                dpz1 = dpz1 + dpzx[bc]
                dpz2 = dpz2 + dpzy[bc]
                dpz3 = dpz3 + (-1.0*dpzx[bc]*yb[bc] + dpzy[bc]*xb[bc])
                
                

            dpz[0,0] = dpz1
            dpz[1,0] = dpz2
            dpz[2,0] = dpz3

            # added
            dpfs1[0,0] = dpfs11
            dpfs1[1,0] = dpfs12
            dpfs1[2,0] = dpfs13
            # added

            pzt[0,0] += dpz[0,0]
            pzt[1,0] += dpz[1,0]
            pzt[2,0] += dpz[2,0]

            # added
            pfs1t[0,0] += dpfs1[0,0]
            pfs1t[1,0] += dpfs1[1,0]
            pfs1t[2,0] += dpfs1[2,0]
            # added
             
        
        pz[0,0] += dpz[0,0]
        pz[1,0] += dpz[1,0]
        pz[2,0] += dpz[2,0]

        # added
        pfs1[0,0] += dpfs1[0,0]
        pfs1[1,0] += dpfs1[1,0]
        pfs1[2,0] += dpfs1[2,0]
        # added

        for wc in range(4):
            zx[wc] = zx[wc] + dzx[wc]
            zy[wc] = zy[wc] + dzy[wc]

            pzx[wc] = pzx[wc] + dpzx[wc]
            pzy[wc] = pzy[wc] + dpzy[wc]

            pfs1x[wc] = pfs1x[wc] + dpfs1x[wc] # added
            pfs1y[wc] = pfs1y[wc] + dpfs1y[wc] # added

            fabx[wc, 0] = bkx[wc]*dbcx2[wc] + pzx[wc] + pfs1x[wc] # added
            faby[wc, 0] = bky[wc]*dbcy2[wc] + pzy[wc] + pfs1y[wc] # added
        
        ep = p2 - np.dot(cd, v2) - np.dot(sk, d2)
        ep[3:6,0] = ep[3:6,0] - pz.T - pfs1.T  # added
        a2 = np.dot(sm_inv, ep)

        eki = eki + 0.5*dt*(np.dot(np.dot(v2.T, sm),a2) + np.dot(np.dot(v1.T, sm),a1))
        edi = edi + 0.5*dt*(np.dot(np.dot(v2.T, cd),v2) + np.dot(np.dot(v1.T, cd),v1))
        esi = esi + 0.5*dt*(np.dot(np.dot(v2.T, sk),d2) + np.dot(np.dot(v1.T, sk),d1))
        eii = eii - 0.5*dt*(np.dot(np.dot(v2.T, sm), np.array([[xg[i]], [yg[i]], [0.0], [xg[i]], [yg[i]], [0.0]])) + np.dot(np.dot(v1.T, sm), np.array([[xg[i-1]], [yg[i-1]], [0.0], [xg[i-1]], [yg[i-1]], [0.0]])))


        d1, v1, p1, a1 = d2, v2, p2, a2 # added
        yd1, yv1, ya1 =  yd2, yv2, ya2 # added
        fs1x1, fs1y1 = pfs1x, pfs1y # added

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
            # print(time[0,index], db[2, index])
            vb[:, index] = v2[3:6,0]
            ab[:, index] = a2[3:6,0]
            aab[0,index] = a2[3,0] + xg[i]
            aab[1,index] = a2[4,0] + yg[i]
            aab[2,index] = a2[5,0] + 0.0
            
            # DO WORK HERE
            
            # baseshear = np.dot(sk, d2) + np.dot(r,np.dot(sm, a2) + np.dot(np.dot(sm, r), ug2) # np.dot(sk, d2) + np.dot(cd, v2) + np.dot(r, pz) # + 
           
            # for wc in range(3):
            #     # f[index, wc] = sm[wc,wc]*aa[wc,index] + sm[wc+3, wc+3]*aab[wc, index]
            #     f[index, wc] = baseshear[wc + 3, 0] 
                
      
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
                
                # DO WORK HERE

                # fcx[index, nc] = -1.0*(f[index, 0] - yb[nc]*f[index, 2]) # fabx[nc] # 
                # fcy[index, nc] = f[index, 1] + xb[nc]*f[index, 2] # faby[nc] #

                fcx[index, nc] = fabx[nc] # 
                fcy[index, nc] = faby[nc] #

                f[index, 0] += fabx[nc]
                f[index, 1] += faby[nc]
                f[index, 2] += -fabx[nc]*yb[nc] + faby[nc]*xb[nc]


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
    peakbaseshearcenterX = max(abs(f[:, 0]))
    peakbaseshearcenterY = max(abs(f[:, 1]))
    peakbaseshearcornerX1 = max(abs(fcx[:, 0]))
    peakbaseshearcornerY1 = max(abs(fcy[:, 0]))
    peakbaseshearcornerX2 = max(abs(fcx[:, 1]))
    peakbaseshearcornerY2 = max(abs(fcy[:, 1]))
    


    print("Velocity in Structure")
    print(v.shape)
    print(v)

    print("Velocity in Isolator")
    print(vb.shape)
    print(vb)

    print("Velocity Corner-X in Structure")
    print(vcx.shape)
    print(vcx)
    print("Velocity Corner-Y in Structure")
    print(vcy.shape)
    print(vcy)

    print("Velocity Corner-X in Isolator")
    print(vbcx.shape)
    print(vbcx)
    print("Velocity Corner-Y in Isolator")
    print(vbcy.shape)
    print(vbcy)

    print("Force Corner-X in Isolator")
    print(fcx.shape)
    print(fcx)
    print("Force Corner-Y in Isolator")
    print(fcy.shape)
    print(fcy)

    print(" ")
    print("Simulation" + "\033[91m" + " SET%d-%d" %(ref, ijk) + "\033[0m" + ": Earthquake #: %d, Parameter #: %d" %(ref, ijk))
    print("Peak Error: % 8.6f" %(peakerror))
    print("Absolute Sum of Errors: % 8.6f" %(sumerror))
    print("Peak Top Floor Absolute Acceleration in X-Direction: % 8.6f m/s2 | % 8.6f g" %(peaktopaccx, peaktopaccx/9.81))
    print("Peak Top Floor Absolute Acceleration in Y-Direction: % 8.6f m/s2 | % 8.6f g" %(peaktopaccy, peaktopaccy/9.81))
    print("Peak Top Floor Absolute Acceleration in Theta-Direction: % 8.6f rad/s2" %(peaktopacctheta))
    print("Peak Top Floor Relative Displacement in X-Direction: % 8.6f cm" %(peaktopdispx*100.0))
    print("Peak Top Floor Relative Displacement in Y-Direction: % 8.6f cm" %(peaktopdispy*100.0))
    print("Peak Top Floor Relative Displacement in Theta-Direction: % 8.6f rad" %(peaktopdisptheta))
    print("Peak Isolator Displacement in X-Direction: % 8.6f cm" %(peakbasedispx*100.0))
    print("Peak Isolator Displacement in Y-Direction: % 8.6f cm" %(peakbasedispy*100.0))
    print("Peak Isolator Displacement in Theta-Direction: % 8.6f rad" %(peakbasedisptheta))
    print("Peak Resultant Base Shear in X-Direction: % 8.6f N" %(peakbaseshearcenterX))
    print("Peak Resultant Base Shear in Y-Direction: % 8.6f N" %(peakbaseshearcenterY))
    print("Peak Corner 1 Base Shear in X-Direction: % 8.6f N" %(peakbaseshearcornerX1))
    print("Peak Corner 1 Base Shear in Y-Direction: % 8.6f N" %(peakbaseshearcornerY1))
    print("Peak Corner 2 Base Shear in X-Direction: % 8.6f N" %(peakbaseshearcornerX2))
    print("Peak Corner 2 Base Shear in Y-Direction: % 8.6f N" %(peakbaseshearcornerY2))
    
    
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