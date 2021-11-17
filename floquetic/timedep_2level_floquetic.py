# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy import integrate
import csv
import cmath
import time


type = r"\floquetic_zeros"
date = "211117"
ver = "3"

z = Symbol('z')

# # transition matrix elements
# r = lambda x: np.pi/100*np.sin(x)
# phi_a = lambda x: np.pi*3/4+np.pi/10*np.cos(x)
# phi_b = lambda x: np.pi/4+np.pi/10*np.sin(x)
#
a1_L = lambda x: 0.4
a1_R = lambda x: 0.1
b1_L = lambda x: 0.2
b1_R = lambda x: 0.3

a1 = lambda x : a1_R(x)+a1_L(x)
b1 = lambda x : b1_R(x)+b1_L(x)

a2_L = lambda x: 0.3
a2_R = lambda x: 0.4
b2_L = lambda x: 0.2
b2_R = lambda x: 0.1

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
z_1 = list(solveset(np.trace(W_1(z))**2-4*(W_1(z)[0][0]*W_1(z)[1][1]-W_1(z)[0][1]*W_1(z)[1][0]),z))
z_2 = list(solveset(np.trace(W_2(z))**2-4*(W_2(z)[0][0]*W_2(z)[1][1]-W_2(z)[0][1]*W_2(z)[1][0]),z))

z_disc = list(solveset(Trace(z)**2-4*Det(z),z))

z_1
z_2
z_disc = np.array(z_disc)
z_disc = z_disc.astype(np.float64)
z_disc
Trace(0.1)
f_zero = open(r"C:\Users\hyoshida\Desktop\floquetic\zero_"+str(date)+"_"+str(ver)+r".dat",'w')
f_zero.write(str(z_disc[0])+" "+str(z_disc[1])+" "+str(z_disc[2])+" "+str(z_disc[3]))
f_zero.close()
def R(x):
    # return (a1(theta)*(1-b2(theta)*dt)+a2(theta)*(1-a1(theta)*dt)+b1(theta)*(1-a2(theta)*dt)+b2(theta)*(1-b1(theta)*dt))*dt/Trace(x)
    return np.abs(Trace(1)-2)/Trace(x)
R(1)
def root(x):
    return np.complex(sqrt(-(x-z_disc[0])*(x-z_disc[1])*(x-z_disc[2])*(x-z_disc[3])/(x**2*(1-z_disc[0])*(1-z_disc[1])*(1-z_disc[2])*(1-z_disc[3]))))
root(-0.1)
def rho(x):
    temp = 1/np.pi*(R(x)*root(x))/(1+R(x)**2*root(x)**2)*(1/(x-z_disc[0])+1/(x-z_disc[1])+1/(x-z_disc[2])+1/(x-z_disc[3])-2/x-2*(a1_R(theta)*b2_L(theta)+a2_R(theta)*b1_L(theta)-(a1_L(theta)*b2_R(theta)+a2_L(theta)*b1_R(theta))/x**2)*dt**2/Trace(x))
    if (np.imag(temp)!=0):
        # return str(nan)
        return 0
    else:
        return np.abs(temp)
rho(-0.5)
N=100000
dx = 0.00001
Z1 = np.linspace(z_disc[0]-dx,z_disc[1]+dx,N)
Z2 = np.linspace(z_disc[2]-dx,z_disc[3]+dx,N)
Rho = []
for z in Z1:
    Rho.append(rho(z))
# sum(Rho)*(Z1[1]-Z1[0])
for z in Z2:
    Rho.append(rho(z))
# sum(Rho)*(Z2[1]-Z2[0])

Z = np.append(Z1,Z2)
plt.ylim([0,5])
plt.xlim([-25,0])
plt.plot(Z,Rho)
plt.show()
#
# plt.figure(figsize=(4,3))
# plt.xlim([-4.1,0])
# plt.yticks([0,1,2],["0","1","2"])
# plt.xlabel(r"$\Re (z)$")
# plt.plot(z_1,[2,2],label=r"$W_1$ only", linestyle="None",marker="+",color="blue")
# plt.plot(z_2,[1,1],label=r"$W_2$ only", linestyle="None",marker="*",color = "blue")
# plt.plot(z_disc,[0,0,0,0],label=r"$U_2$", linestyle="None",marker=".",color="red")
# plt.legend()
# plt.savefig(r"C:\Users\hyoshida\Desktop\fig.png")
# plt.show()
#
# def integ1(x,z):
#     return rho(x)*np.log((z-x)/(1-x))
#
# def integ2(x,z):
#     return rho(x)*z/(z-x)
#
# def J(z):
#     return (integrate.quad(integ2,z_disc[2]+dx,z_disc[3]-dx,args=z)[0]+integrate.quad(integ2,z_disc[0]+dx,z_disc[1]-dx,args=z)[0])*0.5-1
#
# def phi(z):
#     return (integrate.quad(integ1,z_disc[2]+dx,z_disc[3]-dx,args=z)[0]+integrate.quad(integ1,z_disc[0]+dx,z_disc[1]-dx,args=z)[0])*0.5-np.log(z)-J(z)*np.log(z)
#
#
#
# # def J(z):
# #     sum = 0
# #     i = 0
# #     for x in Z1:
# #         sum += (z_disc[1]-z_disc[0])/N*Rho[i]*z/(z-x)
# #         i += 1
# #     for x in Z2:
# #         sum += (z_disc[3]-z_disc[2])/N*Rho[i]*z/(z-x)
# #         i += 1
# #     return 0.5*sum-1
#
# Chi = np.linspace(-4,5,200)
# J_dat = []
# phi_dat=[]
# for chi in Chi:
#     J_dat.append(J(np.exp(chi)))
#     phi_dat.append(phi(np.exp(chi)))
#
# # def integ(z):
# #     sum = 0
# #     i = 0
# #     for x in Z1:
# #         sum += (z_disc[1]-z_disc[0])/N*Rho[i]*np.log((z-x)/(1-x))
# #         i += 1
# #     for x in Z2:
# #         sum += (z_disc[3]-z_disc[2])/N*Rho[i]*np.log((z-x)/(1-x))
# #         i += 1
# #     return sum
# #
# # Phi_dat = []
# # i = 0
# # for chi in Chi:
# #     Phi_dat.append(0.5*integ(np.exp(chi))-(J_dat[i]+1)*chi)
# #     i += 1
#
# # plt.xlim([-0.6,0.2])
# # plt.ylim([-0.005,0.0001])
# # plt.plot(J_dat,Phi_dat)
# # plt.show()

# f_zero = open(r"C:\Users\hyoshida\Desktop\floquetic\phi_"+str(date)+"_"+str(ver)+r".dat",'w')
# for j in range(100):
#     f_zero.write(str(J_dat[j])+' '+str(phi_dat[j]))
#     f_zero.write('\n')
# f_zero.close()
