# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import csv
import cmath
import time


type = r"\exact_zeros"
date = "211105"
ver = "7"

z = Symbol('z')

# transition matrix elements
r = lambda x: np.pi/100*np.sin(x)
phi_a = lambda x: np.pi*3/4+np.pi/50*np.cos(x)
phi_b = lambda x: np.pi/4+np.pi/50*np.sin(x)

a_L = lambda x: 1/2*(1+r(x))*np.sin(phi_a(x)/2)**2
a_R = lambda x: 1/2*(1+r(x))*np.cos(phi_a(x)/2)**2
b_L = lambda x: 1/2*(1-r(x))*np.sin(phi_b(x)/2)**2
b_R = lambda x: 1/2*(1-r(x))*np.cos(phi_b(x)/2)**2
# a_L = lambda x: 0.3
# a_R = lambda x: 0.3
# b_L = lambda x: 0.2
# b_R = lambda x: 0.2

a = lambda x : a_R(x)+a_L(x)
b = lambda x : b_R(x)+b_L(x)

# initial distrib
# p1 = 0.8
# p2 = 0.2
p1 = a(0)/(a(0)+b(0))
p2 = b(0)/(a(0)+b(0))
print(p1,p2)

N = 5
M = 2
iter = M*N
dt = 1
W = [np.eye(2)]
V1 = np.array([[0,1],[0,0]])
V2 = np.array([[0,0],[1,0]])


for s in range(iter):
    theta = (2*np.pi/N*s)%(2*np.pi)

    X = np.array([[1-b(theta)*dt,a_L(theta)*dt],[b_L(theta)*dt,1-a(theta)*dt]])
    w_z_1 = X +z*a_R(theta)*dt*V1+1/z*b_R(theta)*dt*V2    # normal way of definition is not good.

    W.append(np.dot(W[-1],w_z_1))


Z=[]
for i in range(iter):
    Z.append(np.dot(np.dot([1,1],W[i+1]),[[p1],[p2]])[0])

# start = time.time()
z_disc=[]
for i in range(iter):
    z_disc.append(list(solveset(Z[i],z)))
# print(time.time()-start)

f_zero = open(r"C:\Users\hyoshida\Desktop\timedep"+str(type)+"_"+str(date)+"_"+str(ver)+r"_zeros.dat",'w')
for i in range(N*M):
    for j in range(2*(int(i/2)+1)):
        f_zero.write(str(z_disc[i][j]))
        f_zero.write(' ')
    f_zero.write('\n')
f_zero.close()


########## current density #############

Z = np.exp(np.linspace(-3,2,100))
J = [[0]*100 for i in range(N*M)]
for j in range(N*M):
    n = 0
    for z in Z:
        sum = 0
        for i in range(2*(int(j/2)+1)):
            sum += 2/(2*(int(j/2)+1))*z/(z-z_disc[j][i])
        J[j][n] = 0.5*(sum-1)
        n += 1

Phi = [[0]*100 for i in range(N*M)]
for j in range(N*M):
    n = 0
    for z in Z:
        sum = 0
        for i in range(2*(int(j/2)+1)):
            sum += 2/(2*(int(j/2)+1))*log((z-z_disc[j][i])/(1-z_disc[j][i]))
        Phi[j][n] = 0.5*(sum-2*J[j][n]*np.log(z)-np.log(z))
        n += 1

# fig,ax = plt.subplots()
# for i in range(N*M):
#     ax.plot(J[i],Phi[i],label=i)
# plt.legend()
# plt.show()

f_zero = open(r"C:\Users\hyoshida\Desktop\timedep"+str(type)+"_"+str(date)+"_"+str(ver)+r"_current.dat",'w')
for j in range(100):
    for i in range(N*M):
        f_zero.write(str(J[i][j])+' '+str(Phi[i][j])+' ')
    f_zero.write('\n')
f_zero.close()
