import ctypes
from ctypes import *

# need to put r before the path of the library so that it takes it as a raw path
clibrary = ctypes.CDLL(r"C:\Users\mhaja\OneDrive\Desktop\Github\Honours\ctypes1\clibrary.so")

# defining the function in python by reading it from the library
func = clibrary.display

# assigning the arg types that are passed into the function
# need to use the types from the ctypes library which converts automatically for us
func.argtypes = [ctypes.c_char_p, ctypes.c_int]
# assiging the return type so that we can get it from the function
func.restype = ctypes.c_char_p

#need to put the b before the string so that it takes it in as a binary string
# in python string are immutable but in c strings are mutable which is why we need to add the b before string
msg = func(b"mehrad", 22)
print(msg)




