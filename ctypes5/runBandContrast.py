import ctypes
from ctypes import *
from scipy.optimize import direct



# creating the AFMData structure here
class AFMData(Structure):
    _fields_ = [
        ("xResolution", c_int),
        ("yResolution", c_int),
        ("zValues", POINTER(POINTER(c_double)))  # 2D array of doubles
    ]

# creating the BandContrast structure here
class BandContrast(Structure):
    _fields_ = [
        ("nrow", c_int),
        ("ncol", c_int),
        ("greyScale", POINTER(POINTER(c_double))),
        ("EBSDred", POINTER(POINTER(c_int))),
        ("EBSDgreen", POINTER(POINTER(c_int))),
        ("EBSDblue", POINTER(POINTER(c_int))),
    ]
    
class BandContrastAFMMapper(Structure):
    _fields_ = [
        ("nrow", c_int),
        ("ncol", c_int),
        # Assuming NUMBER_OF_LAYERS_IN_BCAFMM is a constant defining how many layers there are
        # "map" is a pointer to a pointer to a pointer of c_double, representing a 3D array
        ("map", POINTER(POINTER(POINTER(c_double))))
    ]

lib = ctypes.CDLL(r"C:\Users\mhaja\OneDrive\Desktop\Github\Honours\ctypes5\lib.so")


amoeba_chisq = lib.amoeba_chisq
amoeba_chisq.argtypes = [
    ctypes.POINTER(BandContrast),  # bcMeasured
    AFMData,                       # afm (assuming by value, adjust if it's a pointer)
    ctypes.POINTER(BandContrast),  # bcTilted
    ctypes.POINTER(BandContrastAFMMapper),  # bcAFMmOut
    ctypes.c_double,               # mStdDev
    ctypes.c_double,               # simStdDev
    ctypes.POINTER(ctypes.c_double),  # simplexCorner
    ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),  # asbs
    ctypes.c_int                   # ndim
]



