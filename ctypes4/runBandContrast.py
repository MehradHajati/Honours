import ctypes
from ctypes import *
from scipy.optimize import direct

lib = ctypes.CDLL(r"C:\Users\mhaja\OneDrive\Desktop\Github\Honours\ctypes4\lib.so")


main = lib.main

argc = c_int(7) # name of program + 6 arguments
argv = (c_char_p * 7)(b".\BandContrastSim", b"EC2_000", b"70", b"20", b"8", b"0.7", b"1")

main(argc, argv)



