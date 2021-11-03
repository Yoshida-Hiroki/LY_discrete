# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import csv
import cmath
import time


type = r"\exact_zeros"
date = "211104"
ver = "1"

z = Symbol('z')

# initial distrib
p1 = 0.5
p2 = 0.5

# transition matrix elements
r = lambda x: np.pi/100*np.sin(x)
phi_a = lambda x: np.pi*3/4+np.pi/10*np.cos(x)
phi_b = lambda x: np.pi/4+np.pi/10*np.sin(x)

a_L = lambda x: 1/2*(1+r(x))*np.sin(phi_a(x)/2)**2
a_R = lambda x: 1/2*(1+r(x))*np.cos(phi_a(x)/2)**2
b_L = lambda x: 1/2*(1-r(x))*np.sin(phi_b(x)/2)**2
b_R = lambda x: 1/2*(1-r(x))*np.cos(phi_b(x)/2)**2

a = lambda x : a_R(x)+a_L(x)
b = lambda x : b_R(x)+b_L(x)

N=5
cycle_num = 2
iter = cycle_num*N
dt = 1
W = [np.eye(2)]
V1 = np.array([[0,1],[0,0]])
V2 = np.array([[0,0],[1,0]])

start = time.time()
for s in range(iter):
    theta = (2*np.pi/N*s)%(2*np.pi)

    X = np.array([[1-b(theta)*dt,a_L(theta)*dt],[b_L(theta)*dt,1-a(theta)*dt]])
    w_z_1 = X +z*a_R(theta)*dt*V1+1/z*b_R(theta)*dt*V2    # normal way of definition is not good.

    W.append(np.dot(W[-1],w_z_1))
print(time.time()-start)

Z=[]
for i in range(iter):
    Z.append(np.dot(np.dot([1,1],W[i+1]),[[p1],[p2]])[0])

start = time.time()
z_disc=[]
for i in range(iter):
    z_disc.append(list(solveset(Z[i],z)))
print(time.time()-start)

outfile = open(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".csv",'w', newline='')
writer = csv.writer(outfile)
writer.writerows(z_disc)
outfile.close()

# f = open(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".txt",'w')
# f.write('p1='+str(p1)+'\n'+'p2='+str(p2)+'\n\n'+'a_R='+str(a_R)+'\n'+'a_L='+str(a_L)+'\n'+'b_R='+str(b_R)+'\n'+'b_L='+str(b_L)+'\n\n'+'T=0.5s'+'\n'+'N='+str(N)+'\n\n'+'s=1 to '+str(iter))
# f.close()
