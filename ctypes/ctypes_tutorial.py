import ctypes
from ctypes import *
import numpy as np
from numpy import ctypeslib
from scipy.optimize import minimize
from scipy.optimize import direct, Bounds, basinhopping

# need to put before the path of the library so that it takes it as a raw path
clibrary = ctypes.CDLL(r"C:\Users\mhaja\OneDrive\Desktop\Github\Honours\ctypes\clibrary.so")

#need to put the b before the string so that it takes it in as a binary string
# in python string are immutable but in c strings are mutable which is why we need to add the b before string
#clibrary.display(b"mehrad", 22)

# 
func = clibrary.display

# assigning the arg types that are passed into the function
# need to use the types from the ctypes library which converts automatically for us
func.argtypes = [ctypes.c_char_p, ctypes.c_int]
# assiging the return type so that we can get it from the function
func.restype = ctypes.c_char_p

func(b"mehrad", 22)


# trying to use scipy here 

styblinski_tang = clibrary.styblinski_tang

# we have a an array here to we have to pass a pointer to it
styblinski_tang.argtypes = [ctypes.c_float, ctypes.c_float]
styblinski_tang.restype = ctypes.c_float

# we have to create a c int array to be passed to the function
#arr = [1,1]
#arr_c = (ctypes.c_int * 2)(*arr)
#print(styblinski_tang(arr_c))

print(styblinski_tang(1, 1))
#bounds = ([-4, -4], [4, 4])
#bounds_c = np.ctypeslib.as_ctypes(bounds)

# Using a wrapper function 
def foo(xy):
    x,y = xy
    # need to add float to conver it cause they are np.float objects but for it to work with the compiled library it needs to a float object
    return styblinski_tang(x, y)

print(foo([1,1]))

bounds = [[-4, 4], [-4, 4]]

result = direct(foo, bounds, len_tol=1e-3)

print(result.x, result.fun, result.nfev)



