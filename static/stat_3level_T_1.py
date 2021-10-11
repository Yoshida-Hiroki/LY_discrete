# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import csv
import cmath


type = r"\discrete3_1"
date = "211011"
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
    w_z_1 = X +z*a1_R*dt*V1+a2_R*dt*V2+a3_R*dt*V3+1/z*b1_R*dt*V4+b2_R*dt*V5+b3_R*dt*V6    # normal way of definition is not good.
    w_n = np.eye(3)

    for i in range(M):
        w_n = np.dot(w_n,w_z_1)
    return simplify(w_n)

def Z(M):
    return np.dot(np.dot([1,1,1],w(M)),[[p1],[p2],[p3]])[0]

c2 = a1+a2+a3+b1+b2+b3
c1 = -((a1_L+a1_R*z)*(b1_L+b1_R/z)+a2*b2+a3*b3-(a1+b3)*(a2+a3)-(b1+b2)*(a1+b3)+(b1+b2)*(a1+b3))
c0 = -((a1_L+a1_R*z)*b2*a3+(b1_L+b1_R/z)*a2*b3-(b1+b2)*(a1+b3)*(a2+a3)+(a1_L+a1_R*z)*(b1_L+b1_R/z)*(a2+a3)+a2*b2*(a1+b3)+a3*b3*(b1+b2))

X = (9*c1*c1-27*c0-2*c2**3)/54
Y = (3*c1-c2**2)/9

Zp = (X+(X**2+Y**3))**(1/3)
Zn = (X-(X**2+Y**3))**(1/3)

omega = np.exp(1.j*(2*np.pi)/3)

Lambda_1 = Zp+Zn
Lambda_2 = omega*Zp+omega**(-1)*Zn
Lambda_3 = omega**(-1)*Zp+omega*Zn

solve(Lambda_1-Lambda_2,z)
solve(Lambda_1-Lambda_3,z)
# solve(Lambda_2-Lambda_3,z)

simplify(Lambda_1.subs(z,min_z[0])-Lambda_2.subs(z,min_z[0]))

np.real(np.complex(Lambda_1.subs(z,-9)))

hoge = np.linspace(-10,0,100)
Lambda = [[0]*100 for i in range(3)]
i=0
for x in hoge:
    Lambda[0][i] = np.abs(np.complex(Lambda_1.subs(z,x)))
    Lambda[1][i] = np.abs(np.complex(Lambda_2.subs(z,x)))
    Lambda[2][i] = np.abs(np.complex(Lambda_3.subs(z,x)))
    i+=1

for i in range(3):
    plt.plot(hoge,Lambda[i])


def diff(x):
    W = np.matrix([[-(b1+b2),a1_L+a1_R*x,a2],[b1_L+b1_R/x,-(a1+b3),a3],[b2,b3,-(a2+a3)]])
    return np.linalg.eig(W)[0][0]-np.linalg.eig(W)[0][1],np.linalg.eig(W)[0][0]-np.linalg.eig(W)[0][2],np.linalg.eig(W)[0][1]-np.linalg.eig(W)[0][2]

X = np.linspace(-10,0,100)
Y = np.linspace(-5,5,100)

min_z = [0]*3
min1 =10000
i=0
for x in X:
    j=0
    for y in Y:
        if(np.abs(diff(x+1.j*y)[0])<min1):
            min1 = np.abs(diff(x+1.j*y)[0])
            min_z[0] = x+1.j*y
        j+=1
    i+=1

min2 =10000
i=0
for x in X:
    j=0
    for y in Y:
        if(np.abs(diff(x+1.j*y)[1])<min2):
            min2 = np.abs(diff(x+1.j*y)[1])
            min_z[1] = x+1.j*y
        j+=1
    i+=1

min3 =10000
i=0
for x in X:
    j=0
    for y in Y:
        if(np.abs(diff(x+1.j*y)[2])<min3):
            min3 = np.abs(diff(x+1.j*y)[2])
            min_z[2] = x+1.j*y
        j+=1
    i+=1

print(min_z,min1,min2,min3)

# list(W.eigenvals().keys())
#
# def dif(x):
#     return list(W.eigenvals().keys())[0].subs(z,x)-list(W.eigenvals().keys())[1].subs(z,x)
#
# x1=np.linspace(-10,0,101)
# for x in x1:
#     plt.plot(x,dif(x))
# list(W.eigenvals().keys())[0].subs(z,9)

iter = 6
z_disc=[]
for i in range(1,iter+1):
    z_disc.append(list(solveset(Z(i),z)))

for i in range(iter):
    for j in range(int(i/2+1)*2):
        z_disc[i][j] = complex(z_disc[i][j])

fig = plt.figure()
ax1 = plt.subplot2grid((1,1),(0,0))

ax1.plot(np.real(min_z[0]),np.imag(min_z[0]),linestyle="None",marker="s",color="blue",label="min$|\Lambda_1-\Lambda_2|$")
ax1.plot(np.real(min_z[1]),np.imag(min_z[1]),linestyle="None",marker="s",color="blue",label="min$|\Lambda_1-\Lambda_3|$")
ax1.plot(np.real(min_z[2]),np.imag(min_z[2]),linestyle="None",marker="s",color="blue",label="min$|\Lambda_2-\Lambda_3|$")
for i in range(iter):
    ax1.plot(np.real(np.array(z_disc[i])),np.imag(np.array(z_disc[i])),linestyle="None",marker="o",color=[0.9-i/10,0.9-i/10,0.9-i/10],label="Discrete k="+str(i+1))
# for i in range(iter):
#     ax1.plot([np.real(z_disc[i][int(i/2)]),np.real(z_disc[i][int(i/2)+1])],[np.imag(z_disc[i][int(i/2)]),np.imag(z_disc[i][int(i/2)+1])],linestyle="None",marker="o",color=[0.9-i/20,0.9-i/20,0.9-i/20],label="Discrete k="+str(i+1))
ax1.legend()
ax1.set_xlim([-10,1])
# plt.show()
plt.savefig(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".png")

outfile = open(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".csv",'w', newline='')
writer = csv.writer(outfile)
writer.writerows(z_disc)
outfile.close()
