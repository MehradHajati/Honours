import numpy as np
from scipy.optimize import minimize
from scipy.optimize import direct, Bounds

def rosen(x):
    """The Rosenbrock function"""
    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)


def styblinski_tang(pos):
    x, y = pos
    return 0.5 * (x**4 - 16*x**2 + 5*x + y**4 - 16*y**2 + 5*y)


#x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
#res = minimize(rosen, x0, method='nelder-mead', options={'xatol': 1e-8, 'disp': True})


bounds = [[-4, 4], [-4, 4]]
#result = direct(rosen, bounds)


result = direct(styblinski_tang, bounds, len_tol=1e-3)
print(result.x, result.fun, result.nfev)

print(styblinski_tang([1,1]))