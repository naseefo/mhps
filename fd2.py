

from mhps.isolator_osbi import *
from math import sin, cos, tan, atan, pow, exp, sqrt, pi
import scipy.integrate as integrate
import numpy as np
import math
import matplotlib.pyplot as plt

rmbm = 1.0
tbx = 50.0
zetabx = 0.0
rtytxb = 1.0
rzyzxb = 1.0
typevf = 1
mu0 = 0.01
alpha0 = 1
alpha1 = 1
nu = 1
umax = 0.85
D = 0.4
ecc = 0.5
rmrm = 0.05
kg = 100
cg = 100
dg = 0.8
am = 4140/4
niso = 1

iso = IsoOSBIModel(rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, kg, cg, dg, am, niso)

print(iso)


xb = 0.3 # 0.41500097
yb = 0.3 # 0.41500097

M = 6000
rMrM = iso.Mr/M
ya2 = 0.1

vyb = 0.6
vxb = 0.6
vrb = sqrt(pow(vxb, 2.0) + pow(vyb, 2.0))

axg = 5.0
ayg = 5.0
arg = sqrt(pow(axg, 2.0) + pow(ayg, 2.0))

axb = 2.0
ayb = 2.0
arb = sqrt(pow(axb, 2.0) + pow(ayb, 2.0))

theta_r_dot2 = 2.0


drb = sqrt(xb**2.0 + yb**2.0)
print('Resultant displacement = %8.4f m'%(drb))

print('Major radius of osbi (a0) = %8.4f'%(iso.a0))
print('Minor radius of osbi (b0) = %8.4f'%(iso.b0))

theta_d0 = np.interp(drb, iso.xbtab, iso.ttab)
print('Rotation angle (theta_d0) = %8.4f rad'%theta_d0)

theta_r_d0 = atan((iso.b0/iso.a0)*tan(theta_d0))
print('Eccentric angle (theta_r_d0) = %8.4f rad'%(theta_r_d0))

c_d0 = iso.a0*sin(theta_d0)*cos(theta_r_d0) - iso.b0*cos(theta_d0)*sin(theta_r_d0)
print('Horizontal offset (c_d0) = %8.4f m'%(c_d0))

p_d0 = iso.a0*sin(theta_d0)*sin(theta_r_d0) + iso.b0*cos(theta_d0)*cos(theta_r_d0)
print('Vertical offset (p_d0) = %8.4f m'%(p_d0))

print('\n')

fs1x, fs1y = fs1fixed(M, 0, rMrM, c_d0, p_d0, drb, yb, xb)

print("Mass of structure (M) = %8.4f kg"%(M))
print("Mass of all osbi balls (mr) = %8.4f kg"%(iso.Mr))
print("Restoring force in X-direction (fs1x) = %8.4f N"%(fs1x))
print("Restoring force in Y-direction (fs1y) = %8.4f N"%(fs1y))

fs1x, fs1y = fs1(M, 0, rMrM, ya2, c_d0, p_d0, drb, yb, xb)
print("Mass of structure (M) = %8.4f kg"%(M))
print("Mass of all osbi balls (mr) = %8.4f kg"%(iso.Mr))
print("Restoring force in X-direction (fs1x) = %8.4f N"%(fs1x))
print("Restoring force in Y-direction (fs1y) = %8.4f N"%(fs1y))

print('\n')

mu = mu_val(iso, drb)

fs2x, fs2y = fs2fixed(mu, M, 0, rMrM)
print("Mass of structure (M) = %8.4f kg"%(M))
print("Mass of all osbi balls (mr) = %8.4f kg"%(iso.Mr))
print("Frictional force in X-direction (fs1x) = %8.4f N"%(fs2x))
print("Frictional force in Y-direction (fs1y) = %8.4f N"%(fs2y))

fs2x, fs2y = fs2(M, 0, rMrM, mu, ya2)
print("Mass of structure (M) = %8.4f kg"%(M))
print("Mass of all osbi balls (mr) = %8.4f kg"%(iso.Mr))
print("Frictional force in X-direction (fs1x) = %8.4f N"%(fs2x))
print("Frictional force in Y-direction (fs1y) = %8.4f N"%(fs2y))

print('\n')

fbx, fby = fbfixed(iso.Mr, 0.35*9.81, 0.35*9.81)
print("Mass of structure (M) = %8.4f kg"%(M))
print("Mass of all osbi balls (mr) = %8.4f kg"%(iso.Mr))
print("Inertial force in X-direction (fs1x) = %8.4f N"%(fbx))
print("Inertial force in Y-direction (fs1y) = %8.4f N"%(fby))

fbx, fby = fb(iso.Jr, iso.Mr, drb, yb, xb, vrb, vyb, vxb, arg, ayg, axg, arb, ayb, axb, p_d0, theta_r_dot2, c_d0, ya2)
print("Mass of structure (M) = %8.4f kg"%(M))
print("Mass of all osbi balls (mr) = %8.4f kg"%(iso.Mr))
print("Inertial force in X-direction (fs1x) = %8.4f N"%(fbx))
print("Inertial force in Y-direction (fs1y) = %8.4f N"%(fby))

print('\n')
print('P of Fs1')

print(np.interp(0.4, iso.xbtab, iso.ttab)*180/pi)
plt.plot(iso.xbtab, iso.ttab*180/pi)
plt.xlabel('xbtab')
plt.ylabel('theta_0_or_r')
plt.show()