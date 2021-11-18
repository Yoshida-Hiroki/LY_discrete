# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy import integrate
import csv
import cmath
import time


type = r"\floquetic_N4_zeros"
date = "211118"
ver = "1"

z = Symbol('z')

# transition matrix elements
r = lambda x: 0.5+1/3*np.sin(x)
phi_a = lambda x: np.pi*0.5+np.pi/3*np.cos(x)
phi_b = lambda x: phi_a(x)

a_L = lambda x: 1/2*(1+r(x))*np.sin(phi_a(x)/2)**2
a_R = lambda x: 1/2*(1+r(x))*np.cos(phi_a(x)/2)**2
b_L = lambda x: 1/2*(1-r(x))*np.sin(phi_b(x)/2)**2
b_R = lambda x: 1/2*(1-r(x))*np.cos(phi_b(x)/2)**2

a = lambda x : a_R(x)+a_L(x)
b = lambda x : b_R(x)+b_L(x)

dt = 1
W = [np.eye(2)]
V1 = np.array([[0,1],[0,0]])
V2 = np.array([[0,0],[1,0]])

W_1 = lambda z : np.array([[1-b(0)*dt,a_L(0)*dt+a_R(0)*dt*z],[b_L(0)*dt+b_R(0)*dt/z,1-a(0)*dt]])
W_2 = lambda z : np.array([[1-b(np.pi*0.5)*dt,a_L(np.pi*0.5)*dt+a_R(np.pi*0.5)*dt*z],[b_L(np.pi*0.5)*dt+b_R(np.pi*0.5)*dt/z,1-a(np.pi*0.5)*dt]])
W_3 = lambda z : np.array([[1-b(np.pi)*dt,a_L(np.pi)*dt+a_R(np.pi)*dt*z],[b_L(np.pi)*dt+b_R(np.pi)*dt/z,1-a(np.pi)*dt]])
W_4 = lambda z : np.array([[1-b(np.pi*1.5)*dt,a_L(np.pi*1.5)*dt+a_R(np.pi*1.5)*dt*z],[b_L(np.pi*1.5)*dt+b_R(np.pi*1.5)*dt/z,1-a(np.pi*1.5)*dt]])

U_4 = lambda z : np.dot(W_4(z),np.dot(W_3(z),np.dot(W_2(z),W_1(z))))

Trace = lambda z : np.trace(U_4(z))
Det = lambda z : U_4(z)[0][0]*U_4(z)[1][1]-U_4(z)[0][1]*U_4(z)[1][0]
# z_1 = list(solveset(np.trace(W_1(z))**2-4*(W_1(z)[0][0]*W_1(z)[1][1]-W_1(z)[0][1]*W_1(z)[1][0]),z))
# z_2 = list(solveset(np.trace(W_2(z))**2-4*(W_2(z)[0][0]*W_2(z)[1][1]-W_2(z)[0][1]*W_2(z)[1][0]),z))

z_disc = list(solveset(Trace(z)**2-4*Det(z),z))

z_disc = np.array(z_disc)
z_disc = z_disc.astype(np.float64)
z_disc

simplify(diff(log(Trace(z))))
simplify(diff(log(Trace(z)))*Trace(z)*z**3).coeff(z,0)
simplify(Trace(z))
R(1)
nume = []
denom = []
for i in range(5):
    nume.append(simplify(diff(log(Trace(z)))*Trace(z)*z**3).coeff(z,4-i))
    denom.append(simplify(Trace(z)*z**2).coeff(z,4-i))
sum(denom)
f_zero = open(r"C:\Users\hyoshida\Desktop\floquetic\zero_"+str(date)+"_"+str(ver)+r".dat",'w')
f_zero.write(str(z_disc[0])+" "+str(z_disc[1])+" "+str(z_disc[2])+" "+str(z_disc[3])+" "+str(z_disc[4])+" "+str(z_disc[5])+" "+str(z_disc[6])+" "+str(z_disc[7])+" ")
f_zero.write(str(nume[0])+" "+str(nume[1])+" "+str(0)+" "+str(nume[3])+" "+str(nume[4])+" ")
f_zero.write(str(denom[0])+" "+str(denom[1])+" "+str(denom[2])+" "+str(denom[3])+" "+str(denom[4]))
f_zero.close()


# def R(x):
#     # return (a1(theta)*(1-b2(theta)*dt)+a2(theta)*(1-a1(theta)*dt)+b1(theta)*(1-a2(theta)*dt)+b2(theta)*(1-b1(theta)*dt))*dt/Trace(x)
#     return np.abs(Trace(1)-2)/Trace(x)
#
# def root(x):
#     return np.complex(sqrt(-(x-z_disc[0])*(x-z_disc[1])*(x-z_disc[2])*(x-z_disc[3])/(x**2*(1-z_disc[0])*(1-z_disc[1])*(1-z_disc[2])*(1-z_disc[3]))))
#
# def rho(x):
#     temp = 1/np.pi*(R(x)*root(x))/(1+R(x)**2*root(x)**2)*(1/(x-z_disc[0])+1/(x-z_disc[1])+1/(x-z_disc[2])+1/(x-z_disc[3])-2/x-2*diff(log(Trace(x))))
#     if (np.imag(temp)!=0):
#         # return str(nan)
#         return 0
#     else:
#         return np.abs(temp)
#
# N=1000
# dx = 0.00001
# Z1 = np.linspace(z_disc[0]-dx,z_disc[1]+dx,N)
# Z2 = np.linspace(z_disc[2]-dx,z_disc[3]+dx,N)
# Z3 = np.linspace(z_disc[4]-dx,z_disc[5]+dx,N)
# Z4 = np.linspace(z_disc[6]-dx,z_disc[7]+dx,N)
# Rho = []
# for z in Z1:
#     Rho.append(rho(z))
# for z in Z2:
#     Rho.append(rho(z))
# for z in Z3:
#     Rho.append(rho(z))
# for z in Z4:
#     Rho.append(rho(z))
#
#
# Z = np.append(np.append(np.append(Z1,Z2),Z3),Z4)
# plt.ylim([0,0.5])
# plt.xlim([-0.004,0])
# plt.plot(Z,Rho,linestyle="None",marker="s")
# plt.show()
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
