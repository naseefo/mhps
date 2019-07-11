from scipy.optimize import fsolve
from math import pow
import numpy as np
from mhps.isolator_tor_osbi import IsoOSBITorsionModel

iso = IsoOSBITorsionModel(rmbm=1.0, tbx=50.0, zetabx=0.0, ebxd=0.001, wrwxb=1.0, efxd=0.001, f0=0.05, g1=0.001, bt=0.5, g2=0.5, a=1.0, nt=15, D=0.4, ecc=0.0, rmrm=0.05, am=1.0, niso=4, tbe0x=2.5, e0exd=0.3, mu0exd=0.0)






print("Isolator Initial Stiffness")
print(iso.k0i)
print("Isolator oblateness ratio")
print(iso.ecc_i)
print("Maximum Displacement")
print(iso.umax)
print("a0")
print(iso.a0)
print("b0")
print(iso.b0)
print("mr")
print(iso.mr)
print("xbtab")
print(iso.xbtab)
print("ttab")
print(iso.ttab)