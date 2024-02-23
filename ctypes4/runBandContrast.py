import ctypes
from ctypes import *

# getting the library file
lib = ctypes.CDLL(r"C:\Users\mhaja\OneDrive\Desktop\Github\Honours\ctypes4\lib.so")

# reading the function from the library and defining it
main = lib.main

# name of program + 6 arguments
argc = c_int(7) 
# the arguments to be passed into the program: sample name, gamma, beta, widthInUM, alpha, readFacets (0 to calculate facets, 1 to read from files)
argv = (c_char_p * 7)(b".\BandContrastSim", b"EC2_000", b"70", b"20", b"8", b"0.7", b"1")

#calling the main function
main(argc, argv)



