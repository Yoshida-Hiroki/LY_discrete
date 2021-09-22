# Focus on static 2level system
import numpy as np
from sympy import *

z = Symbol('z')

# initial distrib
p1 = 0.5
p2 = 0.5

# time
dt = 1

# transition matrix elements
a_R = 0.4
a_L = 0.3
b_R = 0.1
b_L = 0.2

a = a_R+a_L
b = b_R+b_L

# continuum zeros
z1_cont = (-((b-a)**2/4+b_L*a_L+b_R*a_R)+np.sqrt(((b-a)**2/4+b_L*a_L+b_R*a_R)**2-4*a_R*a_L*b_R*b_L))/(2*a_R*b_L)
z2_cont = (-((b-a)**2/4+b_L*a_L+b_R*a_R)-np.sqrt(((b-a)**2/4+b_L*a_L+b_R*a_R)**2-4*a_R*a_L*b_R*b_L))/(2*a_R*b_L)

print(z1_cont)
print(z2_cont)

# discrete zeros
a_z = a_L + a_R * z
b_z = b_L + b_R / z

w_z = [[1-b*dt, a_z],[b_z, 1-a*dt]]
print(w_z)

def w(w_z,k):
    w_n = [[1,0],[0,1]]
    for i in range(k):
        w_n = np.dot(w_n,w_z)
    return simplify(w_n)


print(solveset(np.dot(np.dot([1,1],w(w_z,3)),[[p1],[p2]])[0],z))
print(solveset(np.dot(np.dot([1,1],w(w_z,5)),[[p1],[p2]])[0],z))
print(solveset(np.dot(np.dot([1,1],w(w_z,10)),[[p1],[p2]])[0],z))
print(solveset(np.dot(np.dot([1,1],w(w_z,20)),[[p1],[p2]])[0],z))
