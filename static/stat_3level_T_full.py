# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import csv
import cmath


type = r"\discrete3_full"
date = "210928"
ver = "1"

z = Symbol('z')

# initial distrib
p1 = 0.5
p2 = 0.3
p3 = 0.2

# transition matrix elements
a1_R = 0.4
a1_L = 0.3
a2_R = 0.2
a2_L = 0.3
a3_R = 0.4
a3_L = 0.1
b1_R = 0.7
b1_L = 0.2
b2_R = 0.3
b2_L = 0.2
b3_R = 0.6
b3_L = 0.2

a1 = a1_R+a1_L
b1 = b1_R+b1_L
a2 = a2_R+a2_L
b2 = b2_R+b2_L
a3 = a3_R+a3_L
b3 = b3_R+b3_L

# continuum zeros does not exist

# discrete zeros

def w(s):
    # params
    T = 100
    M = s
    dt = T/M
    X = np.array([[1-(b1+b2)*dt,a1_L*dt,a2_L*dt],[b1_L*dt,1-(a1+b3)*dt,a3_L*dt],[b2_L*dt,b3_L*dt,1-(a2+a3)*dt]])
    V1 = np.array([[0,1,0],[0,0,0],[0,0,0]])
    V2 = np.array([[0,0,1],[0,0,0],[0,0,0]])
    V3 = np.array([[0,0,0],[0,0,1],[0,0,0]])
    V4 = np.array([[0,0,0],[1,0,0],[0,0,0]])
    V5 = np.array([[0,0,0],[0,0,0],[1,0,0]])
    V6 = np.array([[0,0,0],[0,0,0],[0,1,0]])
    w_z_1 = X +z*a1_R*dt*V1+z*a2_R*dt*V2+z*a3_R*dt*V3+1/z*b1_R*dt*V4+1/z*b2_R*dt*V5+1/z*b3_R*dt*V6    # normal way of definition is not good.
    w_n = np.eye(3)

    for i in range(M):
        w_n = np.dot(w_n,w_z_1)
    return simplify(w_n)

def Z(M):
    return np.dot(np.dot([1,1,1],w(M)),[[p1],[p2],[p3]])[0]


iter = 4
z_disc=[]
for i in range(1,iter+1):
    z_disc.append(list(solveset(Z(i),z)))

def num_root(i):
    if i%3==0:
        return int(2+4*i/3)
    else:
        return 4+4*int(i/3)

for i in range(iter):
    for j in range(num_root(i)):
        z_disc[i][j] = complex(z_disc[i][j])

fig = plt.figure()
ax1 = plt.subplot2grid((1,1),(0,0))

for i in range(iter):
    ax1.plot(np.real(z_disc[i]),np.imag(z_disc[i]),linestyle="None",marker="o",color=[0.9-i/10,0.9-i/10,0.9-i/10],label="Discrete k="+str(i+1))
# for i in range(iter):
#     ax1.plot([np.real(z_disc[i][int(i/2)]),np.real(z_disc[i][int(i/2)+1])],[np.imag(z_disc[i][int(i/2)]),np.imag(z_disc[i][int(i/2)+1])],linestyle="None",marker="o",color=[0.9-i/20,0.9-i/20,0.9-i/20],label="Discrete k="+str(i+1))
ax1.legend()
# ax1.set_xlim([-10,1])
# plt.show()
plt.savefig(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".png")

outfile = open(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".csv",'w', newline='')
writer = csv.writer(outfile)
writer.writerows(z_disc)
outfile.close()
