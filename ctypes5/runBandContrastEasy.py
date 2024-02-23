import ctypes
from ctypes import *
from scipy.optimize import direct


lib = ctypes.CDLL(r"C:\Users\mhaja\OneDrive\Desktop\Github\Honours\ctypes5\lib.so")


amoeba_chisq = lib.amoeba_chisq
amoeba_chisq.argtypes = [
    ctypes.c_double,               # mStdDev
    ctypes.c_int,               # simStdDev
    POINTER(c_int),             # asbs
    ctypes.c_int                   # ndim
]

amoeba_chisq.restype = ctypes.c_double

python_list = [640, 3.220126, 0.000781, 0.000002, 0.000002, 0.000002, 640.000000, 0.000781, 1.199413, 0.000002, 0.000002, 0.000002]
num = ctypes.c_double(29.5)
asbs = (c_int * len(python_list))(*python_list)
amoeba_chisq(num, 24, asbs, 12)
#def foo(arr, args):
    #a0, a1, a2, a3, a4, a5, b0, b1, b2, b3, b4, b5 = arr
    #mStd, simStd, simplexCorner, ndim = args
    #return amoeba_chisq(mStd, simStd, simplexCorner, ndim)
# put mSTDdev and simSTDDev and simplexCorner and ndim into the args thing which will never change

