
import numpy as np
import math
from data.defaults.param_manager import default
from data.defaults.param_manager import get_default_param_values
import cProfile
import io
import pstats
from scipy import linalg

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

def read_const_param(const_param_file):

    """
    This function reads all the constant set of parameters for the
    program which does not require to be read in iteration of the
    parametric study.
    """
    f1  = open(const_param_file, 'r')
    
    maxnst = int(f1.readline().split(':')[1])
    am = np.array([float(x.strip()) for x in f1.readline().split(':')[1:][0].strip().split(',')])
    ak = np.array([float(x.strip()) for x in f1.readline().split(':')[1:][0].strip().split(',')])
    
    f1.close()
    
    return maxnst, am, ak

def read_ss_var_param(var_param_file):

    """
    This function reads all the set of parameters for the parametric
    studies and stores it in an array. It will be a one-time allocation
    to increase speed.
    """
    
    dtype1 = np.dtype([('IJK', 'i4'), ('NST', 'i4'), ('TX1', 'f4'), ('ZETA','f4'), ('RTYTX','f4')])

    ijk, nst, tx1, zeta, rtytx = np.loadtxt \
    (var_param_file, delimiter=',', usecols=range(5),\
     skiprows=1, unpack = True, dtype=dtype1)

    for i in range(0, len(ijk)):
        yield ijk[i], nst[i], tx1[i], zeta[i], rtytx[i]


def get_total_variable_parameter_set(var_param_file):
 
    dtype1 = np.dtype([('IJK', 'i4'), ('NST', 'i4'), ('TX1', 'f4'), ('ZETA','f4')])

    ijk, nst, tx1, zeta = np.loadtxt(var_param_file, delimiter=',', usecols=range(4), skiprows=1, unpack = True, dtype=dtype1)

    return len(ijk)

def superstructure_propxy(nst, tx1, rtytx, am, ak, zeta, knor):
    
    """
    Calculation of superstructure mass, stiffness and damping matrices
    in X- and Y- direction.
    """
    
    zeta = np.ones(nst)*zeta

    # Calculation of msss matrix in x-direction
    smx = np.zeros(shape=(nst+1,nst+1), dtype=np.dtype('d'), order='C')
    smx[0:nst,0:nst] = np.diag(am[0:nst])   # ------> Mass matrix in x-direction

    # Calculation of msss matrix in y-direction
    smy = np.zeros(shape=(nst+1,nst+1), dtype=np.dtype('d'), order='C')
    smy = smx                  # ------> Mass matrix in y-direction

    # Calculation of stiffness matrix in x- and y-direction
    temp = np.array([[1.,-1.],[-1.,1.]], dtype=np.dtype('d'))
    skx = np.zeros(shape=(nst+1,nst+1), dtype=np.dtype('d'), order='C')
    for i in range(nst-1):
        skx[i:i+2,i:i+2] += temp*ak[i]
    skx[-2,-2] += ak[-1]
    sky = np.zeros(shape=(nst+1,nst+1), dtype=np.dtype('d'), order='C')
    sky = skx
    del temp
    
    D = np.dot(skx[0:nst,0:nst],np.linalg.inv(smx[0:nst,0:nst]))
    e,v = np.linalg.eig(D)

    wx1 = 2.0*np.pi/tx1
    RT = (wx1**2.0)/min(e)
    skx = RT*skx                # -------> Stiffness matrix in x-direction
    
    wy1 = wx1/rtytx
    RT = (wy1**2.0)/min(e)
    sky = RT*sky                # -------> Stiffness matrix in y-direction
    del e, RT, v

    # Calculation of damping matrix in x-direction
    D = np.dot(skx[0:nst,0:nst],np.linalg.inv(smx[0:nst,0:nst]))
    e,v = np.linalg.eig(D)
    idx = e.argsort()[::1]
    e = e[idx]
    #v = v[:,idx]

    T = 2.0*math.pi/np.sqrt(e)
    del e, idx, v

    W = 2.0*math.pi/T
    c_w = (np.ones((nst,nst))*W[None,:]).T
    c_w = np.power(c_w,np.arange(-1,2*nst-2,2))
    a_i = np.dot(np.linalg.inv(c_w),2.0*zeta[0:nst].reshape(nst,1))
    cdx = np.zeros(shape=(nst+1,nst+1),dtype=np.dtype('d'))
    for i in range(nst):
        cdx[0:nst,0:nst] += np.dot(np.linalg.matrix_power(D,i),smx[0:nst,0:nst])*a_i[i]  # Damping matrix in x-direction
    del W, c_w, a_i, D, T

    # Calculation of damping matrix in y-direction
    D = np.dot(sky[0:nst,0:nst],np.linalg.inv(smx[0:nst,0:nst]))
    e,v = np.linalg.eig(D)
    idx = e.argsort()[::1]
    e = e[idx]
    #v = v[:,idx]

    T = 2.0*math.pi/np.sqrt(e)
    del e, idx, v

    W = 2.0*math.pi/T
    c_w = (np.ones((nst,nst))*W[None,:]).T
    c_w = np.power(c_w,np.arange(-1,2*nst-2,2))
    a_i = np.dot(np.linalg.inv(c_w),2.0*zeta[0:nst].reshape(nst,1))
    cdy = np.zeros(shape=(nst+1,nst+1),dtype=np.dtype('d'))
    for i in range(nst):
        cdy[0:nst,0:nst] += np.dot(np.linalg.matrix_power(D,i),smx[0:nst,0:nst])*a_i[i]  # Damping matrix in y-direction
    del W, c_w, a_i, D, T

    return smx, skx, cdx, smy, sky, cdy

@profile
def fixed_simulator(ref, xg, yg, dt, ndiv, ijk, nst, smx, skx, cdx, smy, sky, cdy):
    
    gamma = default_values['NEWMARK_LINEAR'][0]
    beta = default_values['NEWMARK_LINEAR'][1]

    knx = skx + (gamma/beta/dt)*cdx + (1.0/beta/np.power(dt,2))*smx

    dx1 = np.ones((nst,1), dtype=float)*0.0
    vx1 = np.ones((nst,1), dtype=float)*0.0
    px1 = np.diag(-1.0*smx).reshape(nst,1)*xg[0]
    ax1 = np.matmul(linalg.inv(smx), px1 - np.matmul(cdx, vx1) - np.matmul(skx, dx1))

    # I = np.ones((nst,1), dtype=float)

    dx2 = np.zeros((nst,1), dtype=float)
    vx2 = np.zeros((nst,1), dtype=float)
    px2 = np.zeros((nst,1), dtype=float)
    ax2 = np.zeros((nst,1), dtype=float)

    na1x = (1.0/beta/np.power(dt,2))*smx + (gamma/beta/dt)*cdx
    na2x = (1.0/beta/dt)*smx + (gamma/beta - 1.0)*cdx
    na3x = (1.0/2.0/beta - 1.0)*smx + (gamma*dt/2.0/beta - dt)*cdx

    for i in range(1,len(xg)):
        
        px2 = np.diag(-1.0*smx).reshape(nst,1)*xg[i]
        
        pcx1 = px2 + np.matmul(na1x, dx1) + np.matmul(na2x, vx1) + np.matmul(na3x, ax1)
        
        dx2 =  np.matmul(np.linalg.inv(smx), pcx1)

        vx2 = (gamma/beta/dt)*(dx2 - dx1) + (1.0 - gamma/beta)*vx1 + dt*(1.0 - gamma/2.0/beta)*ax1

        ax2 = np.matmul(np.linalg.inv(smx), px2 - np.matmul(cdx, vx2) - np.matmul(skx, dx2))

        dx1, vx1, px1, ax1 = dx2, vx2, px2, ax2 
        

    pass

