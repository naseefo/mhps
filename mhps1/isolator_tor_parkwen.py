
import numpy as np
import math
from data.defaults.param_manager import default
from data.defaults.param_manager import get_default_param_values
import cProfile
import io
import pstats
import pandas as pd
from mhps.postprocessor import ResultFixedXY, ModelInfo
from math import sin, cos, tan, atan, pow, exp, sqrt, pi
import matplotlib.pyplot as plt
import sys, os


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
def simulator_parkwen_tor(ref, xg, yg, dt, ndiv, ndt, ijk, sm, sk, cd, x, y, xb, yb, nb, iso):

    gamma = 0.5
    beta = 1/6
    nit = 1

    

    xb = np.asarray(xb, dtype=np.dtype('d'), order='F')
    yb = np.asarray(yb, dtype=np.dtype('d'), order='F')
    
    wbx = 2.0*pi/iso.tbx
    fm = sm[0,0]
    bm = sm[0,0]*iso.rmbm
    tm = fm + bm

    if iso.tbx > 10:
        ckabx = 0.0
    else:
        ckabx = (fm + bm)*wbx*wbx
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
    print("ALPX")
    print("%20.16f %20.16f %20.16f %20.16f %20.16f %20.16f %20.16f"%(iso.tbx, wbx, tm, iso.g1, ckabx, iso.f0, fm+bm))

    # print(iso.rmbm, fm, bm, wbx, ckabx, iso.g1, qyf, alp)

    fyx[0] = qyf*0.25*(1 + 2.0*iso.efxd)
    fyx[1] = qyf*0.25*(1 - 2.0*iso.efxd)
    fyx[2] = qyf*0.25*(1 - 2.0*iso.efxd)
    fyx[3] = qyf*0.25*(1 + 2.0*iso.efxd)

    fyy = fyx

    # print(fyx)

    # print("f0")
    # print(iso.f0)
    # print("g1")
    # print(iso.g1)
    # print("qyf")
    # print(qyf)
    # print("alp")
    # print(alp)
    # print("fyx")
    # print(fyx)

    


    sm_inv = np.linalg.inv(sm)
    r = np.zeros((6,3), dtype=np.dtype('d'), order='F')
    r[3:6,0:3] = np.diag([1.0, 1.0, 1.0])
    print(r)

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
    a1[3,0] = 0.0
    a1[4,0] = 0.0
    a1[5,0] = 0.0


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
    print(sm)
    print("Linear-Stiffness Matrix")
    print(sk)
    # print("Linear-Damping Matrix")
    # print(cd)

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

   
    for i in range(1,len(xg)):

        t += dt
        
        ug2 = np.array([[xg[i]],[yg[i]], [0.0]])
        p2 = -1.0*np.dot(np.dot(sm, r), ug2)
        # print(p2)
        dp = p2 - p1

        for wc in range(4):
            dpzx[wc] = 0.0
            dpzy[wc] = 0.0
        
        for bc in range(3):
            dpz[bc] = 0.0

        pzt = pz
        dpz1, dpz2, dpz3 = 0.0, 0.0, 0.0
        for i2 in range(nit):
            # print(dp)
            # sys.exit()
            pcx1 = dp + np.dot(na2x, v1) + np.dot(na3x, a1)
            pcx1[3:6,0] = pcx1[3:6,0] - dpz.T
            dd = np.dot(knx_inv, pcx1)
            # print("Debugger")
            # print(dd)
            # sys.exit()
            d2 = d1 + dd
            dv = (gamma/beta/dt)*dd - gamma/beta*v1 + dt*(1.0 - gamma/2.0/beta)*a1
            v2 = v1 + dv

            ep = p2 - np.dot(cd, v2) - np.dot(sk, d2)
            ep[3:6,0] = ep[3:6,0] - pzt.T
            a22 = np.dot(sm_inv, ep)
            
            # CALCULATION OF FORCES FROM BEARING
           
            for bc in range(4):
                dbcx2[bc] = d2[3,0] - yb[bc]*d2[5,0]
                dbcy2[bc] = d2[4,0] + xb[bc]*d2[5,0]

                vbcx2[bc] = v2[3,0] - yb[bc]*v2[5,0]
                vbcy2[bc] = v2[4,0] + xb[bc]*v2[5,0]

                abcx2[bc] = a22[3,0] - yb[bc]*a22[5,0]
                abcy2[bc] = a22[4,0] + xb[bc]*a22[5,0]

                aabcx2[bc] = a22[3,0] + xg[i] - yb[bc]*(a22[5,0] + 0.0)
                aabcy2[bc] = a22[4,0] + yg[i] + xb[bc]*(a22[5,0] + 0.0)
                
                dpzx[bc], dpzy[bc], dzx[bc], dzy[bc] = wen(iso, vbcx2[bc], vbcy2[bc], zx[bc], zy[bc], dt, alp, alp, fyx[bc], fyy[bc])
                dpz1 = dpz1 + dpzx[bc]
                dpz2 = dpz2 + dpzy[bc]
                dpz3 = dpz3 + (-1.0*dpzx[bc]*yb[bc] + dpzy[bc]*xb[bc])
            
            # if i%ndiv == 0:
            #     fffyx = fyx[0]+fyx[1]+fyx[2]+fyx[3]
            #     fffyy = fyy[0]+fyy[1]+fyy[2]+fyy[3]
            #     print("time: %8.6f s | fyx in X : %18.12f | fyy: %18.12f | alp: %18.12f"%(t, fffyx, fffyy, alp))
            #     print("time: %8.6f s | dpz in X : %18.12f | Y: %18.12f"%(t, dpz1, dpz2))

            dpz[0,0] = dpz1
            dpz[1,0] = dpz2
            dpz[2,0] = dpz3

            pzt[0,0] += dpz[0,0]
            pzt[1,0] += dpz[1,0]
            pzt[2,0] += dpz[2,0]
             
        
        pz[0,0] += dpz[0,0]
        pz[1,0] += dpz[1,0]
        pz[2,0] += dpz[2,0]

        for wc in range(4):
            zx[wc] = zx[wc] + dzx[wc]
            zy[wc] = zy[wc] + dzy[wc]    
            pzx[wc] = pzx[wc] + dpzx[wc]
            pzy[wc] = pzy[wc] + dpzy[wc]
            fabx[wc, 0] = bkx[wc]*dbcx2[wc] + pzx[wc]
            faby[wc, 0] = bky[wc]*dbcy2[wc] + pzy[wc]
        
        ep = p2 - np.dot(cd, v2) - np.dot(sk, d2)
        ep[3:6] = ep[3:6] - pz
        a2 = np.dot(sm_inv, ep)

        eki = eki + 0.5*dt*(np.dot(np.dot(v2.T, sm),a2) + np.dot(np.dot(v1.T, sm),a1))
        edi = edi + 0.5*dt*(np.dot(np.dot(v2.T, cd),v2) + np.dot(np.dot(v1.T, cd),v1))
        esi = esi + 0.5*dt*(np.dot(np.dot(v2.T, sk),d2) + np.dot(np.dot(v1.T, sk),d1))
        eii = eii - 0.5*dt*(np.dot(np.dot(v2.T, sm), np.array([[xg[i]], [yg[i]], [0.0], [xg[i]], [yg[i]], [0.0]])) + np.dot(np.dot(v1.T, sm), np.array([[xg[i-1]], [yg[i-1]], [0.0], [xg[i-1]], [yg[i-1]], [0.0]])))

        d1, v1, p1, a1 = d2, v2, p2, a2
        

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
    


    # print("Velocity in Structure")
    # print(v.shape)
    # print(v)

    # print("Velocity in Isolator")
    # print(vb.shape)
    # print(vb)

    # print("Velocity Corner-X in Structure")
    # print(vcx.shape)
    # print(vcx)
    # print("Velocity Corner-Y in Structure")
    # print(vcy.shape)
    # print(vcy)

    # print("Velocity Corner-X in Isolator")
    # print(vbcx.shape)
    # print(vbcx)
    # print("Velocity Corner-Y in Isolator")
    # print(vbcy.shape)
    # print(vbcy)

    # print("Force Corner-X in Isolator")
    # print(fcx.shape)
    # print(fcx)
    # print("Force Corner-Y in Isolator")
    # print(fcy.shape)
    # print(fcy)

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

