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

def root(x):
    return np.sqrt(-(x-z1)*(x-z2)/(x*(1-z1)*(1-z2)))

def rho1(x):
    return R*root(x)/(np.pi*(1+R**2*root(x)**2))*(1/(x-z1)+1/(x-z2)-1/x)

def rho2(x):
    return -R*root(x)/(np.pi*(1+R**2*root(x)**2))*(1/(x-z1)+1/(x-z2)-1/x)

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
    return (integrate.quad(integ3,z1+dx,0-dx,args=z)[0]+integrate.quad(integ4,-100,z2-dx,args=z)[0]-1)*0.5

def phi(z):
    return (integrate.quad(integ1,z1+dx,0-dx,args=z)[0]+integrate.quad(integ2,-100,z2-dx,args=z)[0]-np.log(z))*0.5-J(z)*np.log(z)

Chi = np.linspace(-4,5,100)
J_dat=[]
phi_dat=[]
for chi in Chi:
    J_dat.append(J(np.exp(chi)))
    phi_dat.append(phi(np.exp(chi)))

plt.plot(J_dat,phi_dat)
plt.xlim([-0.2,0.2])
plt.ylim([-0.3,0.05])
plt.show()
J_dat

# def g(x):
#     return 0.5*integrate()
#
#
# x1 = np.linspace(z1+0.000001,-0.001,1000)
# x2 = np.linspace(-2,z2-0.000001,1000)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# fig = plt.figure()
# ax = plt.subplot2grid((1,1),(0,0))
#
# ax.plot(x1,rho1(x1))
# ax.plot(x2,rho2(x2))
# ax.set_ylim([0,10])
# plt.show()
