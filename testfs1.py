from scipy.optimize import fsolve
from math import pow
import numpy as np
from mhps.isolator_tor_osbi import IsoOSBITorsionModel, fs1

iso = IsoOSBITorsionModel(rmbm=1.0, tbx=50.0, zetabx=0.0, ebxd=0.001, wrwxb=1.0, efxd=0.001, f0=0.05, g1=0.001, bt=0.5, g2=0.5, a=1.0, nt=15, D=0.4, ecc=0.0, rmrm=0.05, am=1.0, niso=4, tbe0x=2.5, e0exd=0.3, mu0exd=0.0)


dpfs1x[bc], dpfs1y[bc], yd2[bc], yv2[bc], ya2[bc] = fs1(dt, iso, bc, tm, dbcx2[bc], dbcy2[bc], yd1[bc], yv1[bc], fs1x1[bc], fs1y1[bc])