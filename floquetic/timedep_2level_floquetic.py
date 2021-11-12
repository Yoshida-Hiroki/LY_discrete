# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import csv
import cmath
import time


type = r"\floquetic_zeros"
date = "211112"
ver = "1"

z = Symbol('z')

# # transition matrix elements
# r = lambda x: np.pi/100*np.sin(x)
# phi_a = lambda x: np.pi*3/4+np.pi/10*np.cos(x)
# phi_b = lambda x: np.pi/4+np.pi/10*np.sin(x)
#
a1_L = lambda x: 0.3
a1_R = lambda x: 0.3
b1_L = lambda x: 0.2
b1_R = lambda x: 0.2

a1 = lambda x : a1_R(x)+a1_L(x)
b1 = lambda x : b1_R(x)+b1_L(x)

a2_L = lambda x: 0.31
a2_R = lambda x: 0.31
b2_L = lambda x: 0.19
b2_R = lambda x: 0.19

a2 = lambda x : a2_R(x)+a2_L(x)
b2 = lambda x : b2_R(x)+b2_L(x)

# N = 1
# M = 10
# iter = M*N
dt = 1
W = [np.eye(2)]
V1 = np.array([[0,1],[0,0]])
V2 = np.array([[0,0],[1,0]])

theta = np.pi

W_1 = lambda z : np.array([[1-b1(theta)*dt,a1_L(theta)*dt+a1_R(theta)*dt*z],[b1_L(theta)*dt+b1_R(theta)*dt/z,1-a1(theta)*dt]])
W_2 = lambda z : np.array([[1-b2(theta)*dt,a2_L(theta)*dt+a2_R(theta)*dt*z],[b2_L(theta)*dt+b2_R(theta)*dt/z,1-a2(theta)*dt]])

U_2 = lambda z : np.dot(W_2(z),W_1(z))

Trace = lambda z : np.trace(U_2(z))
Det = lambda z : U_2(z)[0][0]*U_2(z)[1][1]-U_2(z)[0][1]*U_2(z)[1][0]

z_disc = list(solveset(Trace(z)**2-4*Det(z),z))

z_disc = np.array(z_disc)
z_disc = z_disc.astype(np.float64)
z_disc
def R(x):
    # return (a1(theta)*(1-b2(theta)*dt)+a2(theta)*(1-a1(theta)*dt)+b1(theta)*(1-a2(theta)*dt)+b2(theta)*(1-b1(theta)*dt))*2/Trace(x)
    return np.abs(Trace(1)-2)/Trace(x)
R(1)
def root(x):
    return np.complex(sqrt(-(x-z_disc[0])*(x-z_disc[1])*(x-z_disc[2])*(x-z_disc[3])/(x**2*(1-z_disc[0])*(1-z_disc[1])*(1-z_disc[2])*(1-z_disc[3]))))
root(-0.1)
def rho(x):
    temp = 1/np.pi*(R(x)*root(x))/(1+R(x)**2*root(x)**2)*(1/(x-z_disc[0])+1/(x-z_disc[1])+1/(x-z_disc[2])+1/(x-z_disc[3])-2/x-(a1_R(theta)*b2_L(theta)-(a1_L(theta)*b2_R(theta)+a2_L(theta)*b1_R(theta))/x**2)*2/Trace(x))
    if (np.imag(temp)!=0):
        # return str(nan)
        return 0
    else:
        return np.abs(temp)
rho(-1.5)
N=10000
Z = np.linspace(-2,-0.001,N)
Rho = []
for z in Z:
    Rho.append(rho(z))
Rho[9000]
plt.ylim([0,5])
plt.plot(Z,Rho)
plt.show()

def J(z):
    sum = 0
    i = 0
    for x in Z:
        sum += 1/N*Rho[i]*z/(z-x)
        i += 1
    return 0.5*(sum-1)

Chi = np.linspace(-4,5,100)
J_dat = []
for chi in Chi:
    J_dat.append(J(np.exp(chi)))

def integ(z):
    sum = 0
    i = 0
    for x in Z:
        sum += 1/N*Rho[i]*np.log((z-x)/(1-x))
        i += 1
    return 0.5*sum

Phi_dat = []
i = 0
for chi in Chi:
    Phi_dat.append(integ(np.exp(chi))-(J_dat[i]+0.5)*chi)
    i += 1

f_zero = open(r"C:\Users\hyoshida\Desktop\floquetic\phi_"+str(date)+"_"+str(ver)+r".dat",'w')
for j in range(100):
    f_zero.write(str(J_dat[j])+' '+str(Phi_dat[j]))
    f_zero.write('\n')
f_zero.close()
