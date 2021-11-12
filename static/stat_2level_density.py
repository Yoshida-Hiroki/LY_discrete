# Focus on static 2level system
import numpy as np
from sympy import *
from scipy import integrate
import matplotlib.pyplot as plt
import csv
import cmath


type = r"\static_density"
date = "211014"
ver = "1"

# initial distrib
p1 = 0.5
p2 = 0.5

# transition matrix elements
a_R = 0.3
a_L = 0.3
b_R = 0.2
b_L = 0.2

a = a_R+a_L
b = b_R+b_L

z1 = (-((b-a)**2/4+b_L*a_L+b_R*a_R)+np.sqrt(((b-a)**2/4+b_L*a_L+b_R*a_R)**2-4*a_R*a_L*b_R*b_L))/(2*a_R*b_L)
z2 = (-((b-a)**2/4+b_L*a_L+b_R*a_R)-np.sqrt(((b-a)**2/4+b_L*a_L+b_R*a_R)**2-4*a_R*a_L*b_R*b_L))/(2*a_R*b_L)

dt = 1
R = (a+b)/2*dt/(1-(a+b)/2*dt)
R
def root(x):
    return np.real(np.complex(sqrt(-(x-z1)*(x-z2)/(x*(1-z1)*(1-z2)))))
root(-0.1)
def rho1(x):
    return R*root(x)/(np.pi*(1+R**2*root(x)**2))*(1/(x-z1)+1/(x-z2)-1/x)
rho1(-0.5)
def rho2(x):
    return -R*root(x)/(np.pi*(1+R**2*root(x)**2))*(1/(x-z1)+1/(x-z2)-1/x)
N=10000
Z = np.linspace(-2,-0.001,N)
Rho = []
for z in Z:
    Rho.append(np.abs(rho1(z)))
Rho[9000]
plt.ylim([0,5])
plt.xlim([-2,0])
plt.plot(Z,Rho)
plt.show()
#################################################

dx = 0.00001

def integ1(x,z):
    return rho1(x)*np.log((z-x)/(1-x))

def integ2(x,z):
    return rho2(x)*np.log((z-x)/(1-x))

def integ3(x,z):
    return rho1(x)*z/(z-x)

def integ4(x,z):
    return rho2(x)*z/(z-x)

def J(z):
    return (integrate.quad(integ3,z1+dx,0-dx,args=z)[0]+integrate.quad(integ4,-1000,z2-dx,args=z)[0]-1)*0.5
J(1)
def phi(z):
    return (integrate.quad(integ1,z1+dx,0-dx,args=z)[0]+integrate.quad(integ2,-1000,z2-dx,args=z)[0]-np.log(z))*0.5-J(z)*np.log(z)

Chi = np.linspace(-4,5,100)
J_dat=[]
phi_dat=[]
for chi in Chi:
    J_dat.append(J(np.exp(chi)))
    phi_dat.append(phi(np.exp(chi)))
J_dat
######################## simulation ####################################
iter = 100000
M = 100
J_sim = [0]*(4*M+1)
for i in range(iter):
    n = [1]+[0]*(2*M-1)
    j = [0]*(2*M)
    for k in range(1,2*M):
        r = np.random.rand()

        n[k] = (1-n[k-1])*0.5*(1+np.sign(b*dt-r))+n[k-1]*0.5*(1-np.sign(a*dt-r))
        j[k] = -(1-n[k-1])*0.5*(1+np.sign(b_R*dt-r))+n[k-1]*0.5*(1+np.sign(a_R*dt-r))

    J_sim[int(np.sum(j)+2*M)] += 1

# S = iter*2/(2*M)
# MAX=np.max(J_sim)

phi_sim = [0]*(4*M+1)
for i in range(4*M+1):
    phi_sim[i] = np.log(J_sim[i]/iter)/(2*M)

######################## plot ################################

fig = plt.figure()
ax1 = plt.subplot2grid((1,1),(0,0))

ax1.plot(J_dat,phi_dat)
ax1.set_xlim([-0.2,0.2])
ax1.set_ylim([-0.3,0.05])
ax1.plot(np.linspace(-1,1,4*M+1),phi_sim,linestyle="None",marker="s")
ax1.hlines(0,-0.2,0.2,color="gray")
ax1.vlines(0,-0.3,0.05,color="gray")
plt.show()
