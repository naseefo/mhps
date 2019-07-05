
import numpy as np
from math import *
from mhps.isolator_osbi import IsoOSBIModel
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import os


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
ecc = 0.6
rmrm = 0.1
kg = 100
cg = 100
dg = 0.8
am = 4140/4
niso = 1

iso = IsoOSBIModel(rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, kg, cg, dg, am, niso)

M_dash = 4*iso.am + iso.rmbm*am
ro_1 = iso.Mr/M_dash
g = 9.81

fun1 = lambda x: sqrt(1-pow(iso.ecc*sin(x), 2.0))
# theta_r = np.arange(-pi/2 + pi/10000, pi/2 - pi/1000, pi/100)
theta_r = np.arange(-pi/2, pi/2, pi/100)
xb = np.zeros((len(theta_r), 1), dtype=np.dtype('d'), order='F')
kt = np.zeros((len(theta_r), 1), dtype=np.dtype('d'), order='F')
Tt = np.zeros((len(theta_r), 1), dtype=np.dtype('d'), order='F')
idx = 0
for val in theta_r:
    theta_d0 = atan(iso.a0/iso.b0*tan(val))
    
    fv0 = sin(theta_d0)*sin(val) + sqrt(1 - pow(iso.ecc, 2.0))*cos(theta_d0)*cos(val)
    fh0 = sin(theta_d0)*cos(val) - sqrt(1 - pow(iso.ecc, 2.0))*cos(theta_d0)*sin(val)

    fv1 = cos(theta_d0)*sin(val) - sqrt(1 - pow(iso.ecc, 2.0))*sin(theta_d0)*cos(val) + fh0*sqrt(1 - pow(iso.ecc, 2.0))*pow(cos(val)/cos(theta_d0), 2.0)
    fh1 = cos(theta_d0)*cos(val) + sqrt(1 - pow(iso.ecc, 2.0))*sin(theta_d0)*sin(val) - fv0*sqrt(1 - pow(iso.ecc, 2.0))*pow(cos(val)/cos(theta_d0), 2.0)

    S_1 = (fv0*fh1 - fv1*fh0)/(pow(fv0, 2.0)*(sqrt(1 - pow(iso.ecc*sin(theta_d0), 2.0)) - fh1))
    
    
    itheta = integrate.quad(fun1, 0.0, theta_d0)
    xb[idx, 0] = (itheta[0] - fh0) # iso.D*(itheta[0] - fh0)
    kt[idx, 0] = S_1# ((1 + 0.5*ro_1)*M_dash*g/iso.D)*S_1
    if kt[idx, 0] > 0:
        Tt[idx, 0] = (2*pi*sqrt(1/S_1)) # sqrt(iso.D/(1 + 0.5*ro_1)/g)*(2*pi*sqrt(1/S_1))
    else:
        Tt[idx, 0] = 0
    file = os.path.join('studies', 'Kt','kt.csv')
    data = np.hstack((xb, kt, Tt))
    np.savetxt(file, data, delimiter=",")
    idx = idx + 1

plt.plot(xb, kt)
plt.xlabel('Base Displacement')
plt.ylabel('Tangent Stiffness')
plt.show()

# print(iso.mr, iso.V, iso.ro)