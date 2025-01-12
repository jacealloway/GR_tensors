import sympy as sp
from functions import calcTensor as t
import numpy as np 

# xp, xm, y, z = sp.symbols('xplus xminus y z')

# W, L = sp.symbols('W L', cls=sp.Function)

# W = L(xp, y, z)

# A = t.retrieveAll(np.array([[W, 1, 0, 0], [1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]), np.array([xp, xm, y, z]), 'ALL')





time, r, theta, phi, m = sp.symbols('t r theta phi mass')

F = sp.symbols('F', cls=sp.Function)


t.retrieveAll([[1-2*m/r, 0, 0, 0],
               [0, -1/(1-2*m/r), 0, 0],
               [0, 0, -r**2, 0],
               [0, 0, 0, -r**2 * sp.sin(theta)**2]], 
               [time, r, theta, phi], 'SCI')








# time, r, theta, phi, m, a = sp.symbols('t r theta phi m a')

# P, D, X = sp.symbols('Rho Delta Xi', cls=sp.Function)


# gdndn, gupup = t.metricGenerator([  [1 - 2*m*r/P(r, theta),   0,                            0,                      X(r, theta)   ],
#                                     [0,                       -P(r, theta)/D(r),            0,                      0             ],
#                                     [0,                       0,                            -P(r, theta),           0             ],
#                                     [X(r, theta),             0,                            0,                      -(sp.sin(theta)**2)*(r**2 + a**2 + a*X(r, theta))]])

# sp.init_printing()

# print(t.chrisUDD(gupup, gdndn, [time, r, theta, phi])[1])




# print(t.classInfo())


















