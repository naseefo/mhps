from scipy.optimize import fsolve
from math import pow
import numpy as np
from mhps.isolator_tor_osbi import IsoOSBITorsionModel
from mhps.isolator_osbi import IsoOSBIModel
from math import sin, cos, tan, atan, pow, exp, sqrt, pi
import scipy.integrate as integrate



def fs1(dt, iso, idb, tm, dbcx2, dbcy2, yd1, yv1, fs1x1, fs1y1):
    # idb stands for the index of the isolator being passed

    g = 9.81 # in %m/s2

    drb = sqrt(pow(dbcx2, 2.0) + pow(dbcy2, 2.0))

    print('drb %8.4f'%drb)

    theta_d02 = np.interp(drb, iso.xbtab[:, idb], iso.ttab[:, idb])

    print('theta_d02 %8.4f'%theta_d02)

    theta_r_d02 = atan((iso.b0[idb]/iso.a0[idb])*tan(theta_d02))

    print('theta_r_d02 %8.4f'%theta_r_d02)

    c_d02 = iso.a0[idb]*sin(theta_d02)*cos(theta_r_d02) - iso.b0[idb]*cos(theta_d02)*sin(theta_r_d02)

    print(c_d02)

    p_d02 = iso.a0[idb]*sin(theta_d02)*sin(theta_r_d02) + iso.b0[idb]*cos(theta_d02)*cos(theta_r_d02)

    print(p_d02)

    dyb_dt0 = 2.0*(p_d02 - iso.b0[idb]) - yd1

    print(dyb_dt0)

    yd2 = yd1 + dyb_dt0
    
    print(yd2)

    yv2 = (yd2 - yd1)/dt
    
    print(yv2)

    ya2 = (yv2 - yv1)/dt

    print(ya2)

    ya2 = 0

    fs1r = 0.5*tm*((1 + 2.0*iso.mr[idb]/tm)*g + ya2)*(c_d02/2.0/p_d02)

    print(iso.mr[idb])

    print("------")
    print(fs1r)
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
    return dfs1x, dfs1y#, yd2, yv2, ya2

def fs11(M, fcz, rMrM, y_b_d2, c_d0, p_d0, drb, dyb, dxb):

    g = 9.81 # in %m/s2

    fs1r = ((2 + rMrM)*M*g + (2 + 0)*M*y_b_d2)*(c_d0/2.0/p_d0) + fcz*c_d0/p_d0

    if drb == 0:
        fs1x = 0.0
        fs1y = 0.0
    else:
        fs1x = fs1r*dxb/drb
        fs1y = fs1r*dyb/drb

    return fs1x, fs1y



iso = IsoOSBITorsionModel(rmbm=1.0, tbx=50.0, zetabx=0.0, ebxd=0.001, wrwxb=1.0, efxd=0.001, f0=0.05, g1=0.001, bt=0.5, g2=0.5, a=1.0, nt=15, D=0.8, ecc=43.0, rmrm=0.05, am=1.0, niso=4, tbe0x=2.5, e0exd=0.3, mu0exd=0.0)

print("max dis %8.4f"%iso.umax[1])

iso1 = IsoOSBIModel(rmbm=1.0, tbx=50.0, zetabx=0.0, rtytxb=1, rzyzxb=0, typevf=0, mu0=0.05, alpha0=1, alpha1=1, nu=1, umax=1, D=0.8, ecc=0.43, rmrm=0.05, tc=50, ecrc=0, fos_ud=1, am=1, niso=4)


a, b = fs1(dt=0.005/20, iso=iso, idb=1, tm=2, dbcx2=0.3, dbcy2=0.3, yd1=0, yv1=0, fs1x1=0, fs1y1=0)

print(a,b)

dbcx2=0.3
dbcy2=0.3
drb = sqrt(pow(dbcx2, 2.0) + pow(dbcy2, 2.0))
theta_d02 = np.interp(drb, iso1.xbtab, iso1.ttab)
theta_r_d02 = atan((iso1.b0/iso1.a0)*tan(theta_d02))
c_d0 = iso1.a0*sin(theta_d02)*cos(theta_r_d02) - iso1.b0*cos(theta_d02)*sin(theta_r_d02)
p_d0 = iso1.a0*sin(theta_d02)*sin(theta_r_d02) + iso1.b0*cos(theta_d02)*cos(theta_r_d02)
a, b = fs11(M= 2, fcz=0, rMrM=0.05, y_b_d2=0.0, c_d0=c_d0, p_d0=p_d0, drb=drb, dyb=dbcx2, dxb=dbcx2)

print(a,b)