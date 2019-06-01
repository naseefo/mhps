

from mhps.isolator_osbi import fb, fbfixed, fs1, fs1fixed, fs2, fs2fixed, IsoOSBIModel
from math import sin, cos, tan, atan, pow, exp, sqrt, pi
import scipy.integrate as integrate
import numpy as np
import math

iso = IsoOSBIModel(rmbm = 1.0, tbx = 50.0, zetabx = 0.0, rtytxb = 1.0, rzyzxb = 0.0, typevf = 1, mu0 = 0.1, alpha0 = 1, alpha1 = 1, nu = 1, umax = 0.85, D = 0.4, ecc = 0.5, rmrm = 0.05, kg = 100, cg = 100, dg = 1000)

drb = 10  # -0.5869 # -0.6283

xb = 0.2 # 0.41500097
yb = 0.2 # 0.41500097

M = 1000
mr = iso.rmrm*M


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

fs1x, fs1y = fs1fixed(M, mr, c_d0, p_d0, drb, yb, xb)
print("Mass of structure (M) = %8.4f kg"%(M))
print("Mass of all osbi balls (mr) = %8.4f kg"%(mr))
print("Restoring force in X-direction (fs1x) = %8.4f N"%(fs1x))
print("Restoring force in Y-direction (fs1y) = %8.4f N"%(fs1y))



# theta_r_d0 = atan((iso.b0/iso.a0)*tan(theta_d0))
# theta_r_d1 = sin(2*theta_r_d0)/sin(2*theta_d0)
# theta_r_d2 = 2*theta_r_d1*(tan(theta_d0) - theta_r_d1*tan(theta_r_d0))
# F_d0 = iso.a0*sqrt(1 - pow(iso.ecc*sin(theta_d0), 2.0))
# F_d1 = -pow(iso.a0*iso.ecc, 2.0)*sin(theta_d0)*cos(theta_d0)/F_d0
# c_d0 = iso.a0*sin(theta_d0)*cos(theta_r_d0) - iso.b0*cos(theta_d0)*sin(theta_r_d0)
# p_d0 = iso.a0*sin(theta_d0)*sin(theta_r_d0) + iso.b0*cos(theta_d0)*cos(theta_r_d0)
# c_d1 = iso.a0*cos(theta_d0)*cos(theta_r_d0) + iso.b0*sin(theta_d0)*sin(theta_r_d0) - p_d0*theta_r_d1
# p_d1 = iso.a0*cos(theta_d0)*sin(theta_r_d0) - iso.b0*cos(theta_d0)*sin(theta_r_d0) + c_d0*theta_r_d1
# c_d2 = -p_d0*theta_r_d2 - 2*p_d1*theta_r_d1 - c_d0 + c_d0*pow(theta_r_d1, 2.0)
# p_d2 = c_d0*theta_r_d2 + 2*c_d1*theta_r_d1 - p_d0 + p_d0*pow(theta_r_d1, 2.0)
# theta_dot1 = (0.5*vrb)/(F_d0 - c_d1)                                               #vbr
# theta_dot2 = ((F_d1 - c_d2)*pow(theta_dot1, 2.0) - 0.5*arb)/(c_d1 - F_d0)     #abr
# theta_r_dot1 = theta_r_d1*theta_dot1
# theta_r_dot2 = theta_r_d2*pow(theta_dot1, 2.0) + theta_r_d1*theta_dot2
# y_r_d2 = pow(theta_dot1, 2.0)*p_d2 + p_d1*theta_dot2
# y_b_d2 = 2*y_r_d2

# fs1(M, mr, y_b_d2, c_d0, p_d0, phi_drb)
# fs1fixed(M, mr, c_d0, p_d0, phi_drb)
# fb(J, mr, arg, arb, phi_thetadot2, phi_arg, phi_arb, p_d0, theta_r_dot2)
# fbfixed(mr, arg, phi_arg)
# fs2(M, mr, mu, y_b_d2)
# fs2fixed(mu, M, mr)