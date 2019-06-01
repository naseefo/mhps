
from __future__ import print_function
import numba
import numpy as np
import math
import llvmlite
import ctypes
print("numba version: %s \nNumPy version: %s\nllvm version: %s" % (numba.__version__,np.__version__, llvmlite.__version__))

import progressbar
import time
import logging

from mhps.earthquake import *


# def genit():
#     for i in range(0,5):
#         x = np.array([1,2,3,4])*i
#         yield i, x


# def raf_txt(file):
#     path = 'data\\earthquakes\\not_classified\\'
#     xg = pd.read_csv(path + file, header=None).values
#     return xg

# xg = raf_txt("AX.txt")
# print(xg)



# Link: https://www.youtube.com/watch?v=zCxElmuOD8s

# import multiprocessing as mp

# p = mp.pool(2)
# isPrimes = p.map(isPrime, range(2,m+1))
# nPrimes = sum(isPrimes)
# p.close()
# p.join()


#Link: https://player.oreilly.com/videos/9781787280496    SUPER!!!!