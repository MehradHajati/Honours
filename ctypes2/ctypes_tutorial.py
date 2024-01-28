import ctypes
from ctypes import *
from scipy.optimize import direct

clibrary = ctypes.CDLL(r"C:\Users\mhaja\OneDrive\Desktop\Github\Honours\ctypes2\clibrary.so")

# Importing function from library
styblinski_tang = clibrary.styblinski_tang

# converting the parameter and return types
styblinski_tang.argtypes = [ctypes.c_float, ctypes.c_float]
styblinski_tang.restype = ctypes.c_float


# Using a wrapper function 
def foo(xy):
    x,y = xy
    return styblinski_tang(x, y)

# creating the bounds for the minimization function
bounds = [[-4, 4], [-4, 4]]
# calling the direct function with given bounds and small tolerance
result = direct(foo, bounds, len_tol=1e-3)
# printing results
print(result.x, result.fun, result.nfev)



