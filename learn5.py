


import scipy.integrate as integrate
import math
import numpy as np
import matplotlib.pyplot as plt

a0 = 10.0
b0 = 2.0

ecc = math.sqrt((a0**2.0 - b0**2.0)/a0**2.0)
# print('ecc = %3.2f' %(ecc))

# fun1 = lambda x: math.sqrt(1-math.pow(ecc*math.sin(x), 2.0))
# fun1 = lambda x: math.sqrt((a0*math.cos(x))**2.0 + (b0*math.sin(x))**2.0)
fun1 = lambda x: a0*math.sqrt(1 - (ecc*math.sin(x))**2.0)
result = integrate.quad(fun1, 0, math.pi/6)
# result = integrate.quad(fun1, math.pi/2, math.pi/2+math.pi/6)

print(result)


# class pingpong:

#     def prop(self):
#         area = self.a*self.b
#         perimeter = 2.0*(self.a+self.b)
#         return area, perimeter

#     def __init__(self, a, b):
#         self.a = a
#         self.b = b
#         self.area, self.perimeter = self.prop()

# shape1 = pingpong(5.0, 3.0)
# shape2 = pingpong(4.0, 6.0)

# print('SHAPE 1: area = %3.2f | perimeter = %3.2f' %(shape1.area, shape1.perimeter))
# print('SHAPE 2: area = %3.2f | perimeter = %3.2f' %(shape2.area, shape2.perimeter))

# from data.defaults.param_manager import *

# class IsoOSBIModel:
    
#     def osbt2x(self): 
#         default_values = get_default_param_values()
#         tdiv = default_values['TDIV'] 
#         fun1 = lambda x: math.sqrt(1-math.pow(self.ecc*math.cos(x), 2.0))
#         xbtab = np.zeros((tdiv, ), dtype=float)
#         i = 0
#         for ttab in np.arange(0.0,89.9999*math.pi/180,89.9999*math.pi/180/tdiv):
#             trd0 = math.atan((self.b0/self.a0)*math.tan(ttab))
#             rd0 = math.sqrt(math.pow(self.a0*math.sin(ttab), 2.0) + math.pow(self.b0*math.cos(ttab), 2.0))
#             cd0 = rd0*(math.sin(ttab)*math.cos(trd0)-math.cos(ttab)*math.sin(trd0))
#             itheta = integrate.quad(fun1, 0.0, ttab)
#             xbtab[i] = 2.0*self.a0*itheta[0] - 2.0*cd0
#             i += 1
#         ttab = np.arange(0.0,89.9999*math.pi/180,89.9999*math.pi/180/tdiv)
#         return xbtab, ttab

#     def __init__(self, rmbm, tbx, zetabx, rtytxb, rzxzyb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, kg, cg, dg):
#         self.rmbm = rmbm
#         self.tbx = tbx
#         self.zetabx = zetabx
#         self.rtytxb = rtytxb
#         self.rzxzyb = rzxzyb
#         self.typevf = typevf
#         self.mu0 = mu0
#         self.alpha0 = alpha0
#         self.alpha1 = alpha1
#         self.nu = nu
#         self.umax = umax
#         self.a0 = D/2.0
#         self.ecc = ecc
#         self.rmrm = rmrm
#         self.kg = kg
#         self.cg = cg
#         self.dg = dg
#         self.b0 = math.sqrt(math.pow(D/2.0, 2.0) - math.pow(ecc*D/2.0, 2.0))
#         self.xbtab, self.ttab = self.osbt2x()

    
#     def __eq__(self, other):
#         if isinstance(other, self.__class__):
#             return self.rmbm == other.rmbm and self.tbx == other.tbx and \
#                 self.zetabx == other.zetabx and self.rtytxb == other.rtytxb and  self.rzxzyb == other.rzxzyb and \
#                     self.typevf == other.typevf and self.mu0 == other.mu0 and self.alpha0 == other.alpha0 and self.alpha1 == other.alpha1 and \
#                         self.nu == other.nu and self.umax == other.umax and self.kg == other.kg and self.cg == other.cg and self.dg == other.dg and \
#                             self.D == other.D and self.ecc == ecc and self.rmrm == other.rmrm
#         return False


# osbi1 = IsoOSBIModel(1.0, 2.5, 0.0, 1.0, 1.0, 0, 0.1, 1.0, 1.0, 1.0, 0.8155, 0.4, 0.9, 1.0, 100.0, 0.05, 1000)

# print(osbi1.a0, osbi1.b0)
# print(osbi1.xbtab)
# print(osbi1.ttab)

# plt.plot(osbi1.ttab, osbi1.xbtab, label='linear')
# plt.legend()
# plt.show()