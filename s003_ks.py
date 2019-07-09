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
D = 0.3
ecc = 0.5
rmrm = 0.
kg = 100
cg = 100
dg = 0.8
am = 4140/4
niso = 1

iso = IsoOSBIModel(rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, kg, cg, dg, am, niso)
print(iso.mr, iso.Mr, iso.V, iso.ro)