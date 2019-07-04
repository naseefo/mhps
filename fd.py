



from mhps.isolator_osbi import fb, fbfixed, fs1, fs1fixed, fs2, fs2fixed, IsoOSBIModel
from math import sin, cos, tan, atan, pow, exp, sqrt, pi
import scipy.integrate as integrate
import numpy as np
import math
import matplotlib.pyplot as plt
                 
# iso = IsoOSBIModel(1.0, 50.0, 0.0, 1.0, 0.0, 1, 0.1, 1, 1, 1, 0.85, 0.4, 0.5, 0.05, 100, 100, 0.8, 4140, 4)
# print(iso)

# drb = 10  # -0.5869 # -0.6283

# xb = 0.2 # 0.41500097
# yb = 0.2 # 0.41500097

# M = 4140 + 4140



# drb = sqrt(xb**2.0 + yb**2.0)
# print('Resultant displacement = %8.4f m'%(drb))

# print('Major radius of osbi (a0) = %8.4f'%(iso.a0))
# print('Minor radius of osbi (b0) = %8.4f'%(iso.b0))

# theta_d0 = np.interp(drb, iso.xbtab, iso.ttab)
# print('Rotation angle (theta_d0) = %8.4f rad'%theta_d0)

# theta_r_d0 = atan((iso.b0/iso.a0)*tan(theta_d0))
# print('Eccentric angle (theta_r_d0) = %8.4f rad'%(theta_r_d0))

# c_d0 = iso.a0*sin(theta_d0)*cos(theta_r_d0) - iso.b0*cos(theta_d0)*sin(theta_r_d0)
# print('Horizontal offset (c_d0) = %8.4f m'%(c_d0))

# p_d0 = iso.a0*sin(theta_d0)*sin(theta_r_d0) + iso.b0*cos(theta_d0)*cos(theta_r_d0)
# print('Vertical offset (p_d0) = %8.4f m'%(p_d0))

# fs1x, fs1y = fs1fixed(M, iso.Mr, c_d0, p_d0, drb, yb, xb)
# print("Mass of structure (M) = %8.4f kg"%(M))
# print("Mass of all osbi balls (mr) = %8.4f kg"%(iso.Mr))
# print("Restoring force in X-direction (fs1x) = %8.4f N"%(fs1x))
# print("Restoring force in Y-direction (fs1y) = %8.4f N"%(fs1y))


# print(np.interp(0.628319, iso.xbtab, iso.ttab)*180/pi)
# plt.plot(iso.xbtab, iso.ttab*180/pi)
# plt.xlabel('xbtab')
# plt.ylabel('theta_0_or_r')
# plt.show()



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


# drb = 10  # -0.5869 # -0.6283

# xb = 0.2 # 0.41500097
# yb = 0.2 # 0.41500097

# M = 4140 + 4140



# drb = sqrt(xb**2.0 + yb**2.0)
# print('Resultant displacement = %8.4f m'%(drb))

# print('Major radius of osbi (a0) = %8.4f'%(iso.a0))
# print('Minor radius of osbi (b0) = %8.4f'%(iso.b0))


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
ecc = 0.4
rmrm = 0.05
kg = 100
cg = 100
dg = 0.8
am = 4140/4
niso = 1

iso = IsoOSBIModel(rmbm, tbx, zetabx, rtytxb, rzyzxb, typevf, mu0, alpha0, alpha1, nu, umax, D, ecc, rmrm, kg, cg, dg, am, niso)
# print(iso)

# trd0 = 90
# print('Theta_r = %6.3f deg'%(trd0))
# trd0 = trd0*pi/180.0

# td0 = atan(iso.a0/iso.b0*tan(trd0))
# print('Theta = %6.3f deg'%(td0*180.0/pi))

# c = iso.a0*sin(td0)*cos(trd0) - iso.b0*cos(td0)*sin(trd0) # with modification
# print('c = %6.3f m'%(c))

# p = iso.a0*sin(td0)*sin(trd0) + iso.b0*cos(td0)*cos(trd0)
# print('p = %6.3f'%(p))

# fun1 = lambda x: iso.a0*sqrt(1-pow(iso.ecc*sin(x), 2.0)) # with modification

# I = integrate.quad(fun1, 0, td0)
# print('I = %6.3f'%(I[0]))

# xb = 2.0*(I[0] - c)
# print('xb = %6.3f'%(xb))

# If = integrate.quad(fun1, 0, pi/2)

# print('Half-way mark = %6.3f'%(I[0]/If[0]))

Dt = 2.0*iso.a0
ro = 3400
mr = pi/6.0*ro*pow(Dt, 3.0)*sqrt(1 - pow(iso.ecc, 2.0))
M = 1*4140.0/4*5
rmrM = mr/M
g = 9.81


t_d0 = lambda tr: atan((iso.a0/iso.b0)*tan(tr))
T_d0 = lambda tr: atan((1.0/sqrt(1 - pow(iso.ecc, 2.0)))*tan(tr))

F = lambda tr : iso.a0*sqrt(1-pow(iso.ecc*sin(tr), 2.0))
FF = lambda tr: sqrt(1-pow(iso.ecc*sin(T_d0(tr)), 2.0))

c_d0 = lambda tr : iso.a0*sin(atan(iso.a0/iso.b0*tan(tr)))*cos(tr) - iso.b0*cos(atan(iso.a0/iso.b0*tan(tr)))*sin(tr)
C_d0 = lambda tr: sin(T_d0(tr))*cos(tr) - sqrt(1 - pow(iso.ecc, 2.0))*cos(T_d0)*sin(tr)

p_d0 = lambda tr : iso.a0*sin(atan(iso.a0/iso.b0*tan(tr)))*sin(tr) + iso.b0*cos(atan(iso.a0/iso.b0*tan(tr)))*cos(tr)
P_d0 = lambda tr: sin(T_d0(tr))*sin(tr) + sqrt(1 - pow(iso.ecc, 2.0))*cos(T_d0)*cos(tr)

tr_d1 = lambda tr : (iso.b0/iso.a0)*pow(cos(tr)/cos(t_d0(tr)),2.0) #sin(2*tr)/sin(2*atan(iso.a0/iso.b0*tan(tr)))
Tr_d1 = lambda tr : sqrt(1 - pow(iso.ecc, 2.0))*pow(cos(tr)/cos(t_d0(tr)),2.0)

c_d1 = lambda tr : (iso.a0*cos(atan(iso.a0/iso.b0*tan(tr)))*cos(tr) + iso.b0*sin(atan(iso.a0/iso.b0*tan(tr)))*sin(tr)) - p_d0(tr)*tr_d1(tr)
C_d1 = lambda tr : cos(T_d0(tr))*cos(tr) + sqrt(1 - pow(iso.ecc, 2.0))*sin(T_d0(tr))*sin(tr) - P_d0(tr)*Tr_d1(tr)

p_d1 = lambda tr : (iso.a0*cos(atan(iso.a0/iso.b0*tan(tr)))*sin(tr) - iso.b0*sin(atan(iso.a0/iso.b0*tan(tr)))*cos(tr)) + c_d0(tr)*tr_d1(tr)
P_d1 = lambda tr : cos(T_d0(tr))*sin(tr) - sqrt(1 - pow(iso.ecc, 2.0))*sin(T_d0(tr))*cos(tr) + C_d0(tr)*Tr_d1(tr)




ld0 = 0.0 # Lower limit of td0
ud0 = pi/2 # Upper limit of td0
tdiv = (ud0)/10000
sizeTable = int((ud0 - ld0)/tdiv + 1)
print(sizeTable)

xbtab = np.zeros((sizeTable, ), dtype=np.dtype('d'), order ='F')
kt = np.zeros((sizeTable, ), dtype=np.dtype('d'), order ='F')
tr = np.zeros((sizeTable, ), dtype=np.dtype('d'))
t0 = np.zeros((sizeTable, ), dtype=np.dtype('d'))


# i = 0
# for tr in np.arange(0.0,89.9999*pi/180,89.9999*pi/180/tdiv):
    
#     trs[i] = tr
#     kt[i] = (M*g*(p_d0(tr)*c_d1(tr) - c_d0(tr)*p_d1(tr)))/(pow(p_d0(tr),2)*(2.0*tr_d1(tr)*F(tr) - 2*c_d1(tr)))
#     itheta = integrate.quad(F, 0.0, t_d0(tr))
#     xbtab[i] = 2.0*itheta[0] - 2.0*c_d0(tr)
#     if i >= 1:
#         if kt[i-1] >= 0 and kt[i] < 0: 
#             Fb = (2*M + iso.Mr)*g*c_d0(tr)/p_d0(tr)
#             print('kt = %8.3f N/m | xb = %8.3f m | tr = %8.3f deg'%(kt[i-1], xbtab[i-1], trs[i-1]*180/pi))
#             print('kt = %8.3f N/m | xb = %8.3f m | tr = %8.3f deg'%(kt[i], xbtab[i], trs[i]*180/pi))
#             print('Fb = %8.3f N'%Fb)
#             break

#     i = i+1

i = 0

for td0 in np.arange(ld0, ud0 + tdiv, tdiv):
    t0[i] = td0
    tr[i] = atan(iso.b0/iso.a0*tan(td0))
    kt[i] = (M*g*(p_d0(tr[i])*c_d1(tr[i]) - c_d0(tr[i])*p_d1(tr[i])))/(pow(p_d0(tr[i]),2)*(2.0*tr_d1(tr[i])*F(tr[i]) - 2*c_d1(tr[i])))
    itheta = integrate.quad(F, 0.0, td0)
    xbtab[i] = 2.0*itheta[0] - 2.0*c_d0(tr[i])
    if i >= 1:
        if kt[i-1] >= 0 and kt[i] < 0: 
            Fb = (2*M + iso.Mr)*g*c_d0(tr[i])/p_d0(tr[i])
            print('kt = %8.3f N/m | xb = %8.3f m | tr = %8.3f deg'%(kt[i-1], xbtab[i-1], tr[i-1]*180/pi))
            print('kt = %8.3f N/m | xb = %8.3f m | tr = %8.3f deg'%(kt[i], xbtab[i], tr[i]*180/pi))
            print('Fb = %8.3f N'%Fb)
    i = i+1
print(i)
print(tr, xbtab)

plt.plot(t0*180/pi, xbtab, iso.ttab*180/pi, iso.xbtab)
plt.xlabel('trs (deg)')
plt.ylabel('xbtab (m)')
plt.show()

# itheta = integrate.quad(F, 0.0, t_d0(20*pi/180))
# xbtab = 2.0*itheta[0] - 2.0*c_d0(20*pi/180)
# print(itheta)
# print(c_d0(20*pi/180))
# print(xbtab)


# fun1 = lambda x: sqrt(1 + pow((iso.b0/iso.a0)*((x*sqrt(pow(iso.a0, 2.0) - pow(x, 2.0)))/2.0 + (pow(iso.a0, 2.0)*np.sinh(x/iso.a0))/2.0), 2.0))
# itheta = integrate.quad(fun1, 0.0, )