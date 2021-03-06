# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy import integrate
import csv
import cmath
import time


type = r"\floquetic_zeros"
date = "211228"
ver = "N2_2"

z = Symbol('z')

r_base = 0.5
# r_coef = 0
r_coef = 0.4
phi_coef = 5
r = lambda x: r_base+r_coef*np.sin(x)
# phi_a = lambda x: 0.75*np.pi+np.pi/phi_coef*np.cos(x)
phi_a = lambda x: 0.5*np.pi+np.pi/phi_coef*np.cos(x)
# phi_b = lambda x: 0.5*np.pi+2*np.pi/phi_coef*np.sin(x)
phi_b = lambda x: phi_a(x)

r_prime = lambda x: r_coef*np.cos(x)

a_L = lambda x: 1/2*(1+r(x))*np.sin(phi_a(x)/2)**2
a_R = lambda x: 1/2*(1+r(x))*np.cos(phi_a(x)/2)**2
b_L = lambda x: 1/2*(1-r(x))*np.sin(phi_b(x)/2)**2
b_R = lambda x: 1/2*(1-r(x))*np.cos(phi_b(x)/2)**2

a = lambda x : a_R(x)+a_L(x)
b = lambda x : b_R(x)+b_L(x)

z1 = lambda x : (-((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))+np.sqrt(((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))**2-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x))
z2 = lambda x : (-((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))-np.sqrt(((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))**2-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x))

dt = 1

theta = 0
W_1 = lambda z : np.array([[1-b(theta)*dt,a_L(theta)*dt+a_R(theta)*dt*z],[b_L(theta)*dt+b_R(theta)*dt/z,1-a(theta)*dt]])
W_2 = lambda z : np.array([[1-b(theta+np.pi)*dt,a_L(theta+np.pi)*dt+a_R(theta+np.pi)*dt*z],[b_L(theta+np.pi)*dt+b_R(theta+np.pi)*dt/z,1-a(theta+np.pi)*dt]])

U_2 = lambda z : np.dot(W_2(z),W_1(z))

Trace = lambda z : np.trace(U_2(z))
Det = lambda z : U_2(z)[0][0]*U_2(z)[1][1]-U_2(z)[0][1]*U_2(z)[1][0]
z_1 = list(solveset(np.trace(W_1(z))**2-4*(W_1(z)[0][0]*W_1(z)[1][1]-W_1(z)[0][1]*W_1(z)[1][0]),z))
z_2 = list(solveset(np.trace(W_2(z))**2-4*(W_2(z)[0][0]*W_2(z)[1][1]-W_2(z)[0][1]*W_2(z)[1][0]),z))

coeff = [0]*5
for i in range(5):
    coeff[i] = simplify((Trace(z)**2-4*Det(z))*z**2).coeff(z,4-i)
z_disc = np.roots(coeff)

z_disc


j_d = lambda x : (1-r(x)**2)/(16*np.pi)*(np.cos(phi_a(x))-np.cos(phi_b(x)))
J_d = integrate.quad(j_d,0,2*np.pi)[0]

j_ad = lambda x : r_prime(x)*np.cos(phi_a(x))/(8*np.pi)
J_ad = integrate.quad(j_ad,0,2*np.pi)[0]

nume = []
denom = []
for i in range(3):
    nume.append(simplify(diff(log(Trace(z)))*Trace(z)*z**2).coeff(z,2-i))
    denom.append(simplify(Trace(z)*z).coeff(z,2-i))

# ################## file make ##########################
# f_zero = open(r"C:\Users\hyoshida\Desktop\floquetic\zero_"+str(date)+"_"+str(ver)+r".dat",'w')
# f_zero.write(str(z_disc[0])+" "+str(z_disc[1])+" "+str(z_disc[2])+" "+str(z_disc[3])+" ")
# f_zero.write(str(nume[0])+" "+str(nume[1])+" "+str(nume[2])+" ")
# f_zero.write(str(denom[0])+" "+str(denom[1])+" "+str(denom[2])+" ")
# f_zero.write(str(z_1[0])+" "+str(z_1[1])+" "+str(z_2[0])+" "+str(z_2[1]))
# f_zero.close()

############### graph plot #######################
x = np.linspace(0,2*np.pi,1000)
fig = plt.figure()
ax1 = plt.subplot2grid((2,2),(0,0),rowspan=2)
ax2 = plt.subplot2grid((2,2),(0,1))
ax3 = plt.subplot2grid((2,2),(1,1))

# ###### J_ad = 0 ############
# ax1.set_title(f"$r = {r_base:.2f}$"+"\n"+r"$\phi_a=3\pi/4+$"+f"$\pi/{phi_coef:.3g}\cos$"+"\n"+r"$\phi_b=\pi/2+$"+f"2$\pi/{phi_coef:.3g}\sin$")

###### J_d = 0 ############
ax1.set_title(f"$r = {r_base:.2f}+{r_coef:.4f}\sin$"+"\n"+r"$\phi=\pi/2+$"+f"$\pi/{phi_coef:.3g}\cos$")

ax1.plot(x,np.real(z1(x)),color="red")
ax1.plot(x,np.real(z2(x)),color="blue")
ax1.plot(theta,z1(theta),marker="s",color="red",label=f"{z_1[1]:.3g}",linestyle="None")
ax1.plot(theta,z2(theta),marker="s",color="blue",label=f"{z_1[0]:.3g}",linestyle="None")
ax1.plot(theta+np.pi,z1(theta+np.pi),marker="o",color="red",label=f"{z_2[1]:.3g}",linestyle="None")
ax1.plot(theta+np.pi,z2(theta+np.pi),marker="o",color="blue",label=f"{z_2[0]:.3g}",linestyle="None")
ax1.legend()
# ax1.set_ylim([-10,0.1])
ax1.set_xlabel(r"$\theta$")
ax1.set_ylabel("$z$")
ax1.set_xticks([0,np.pi/2,np.pi,np.pi*3/2,np.pi*2])
ax1.set_xticklabels([r"$0$",r"$\frac{\pi}{2}$",r"$\pi$",r"$\frac{3\pi}{2}$",r"$2\pi$"])
ax1.hlines(0,0,2*np.pi,color="gray")

ax2.set_title(f"$J_d=${J_d:.2g}"+"\n"+f"$J_{{ad}}=${J_ad:.2g}")

# ########### phi_a-phi_b ################
# ax2.plot(phi_a(x),phi_b(x))
# ax2.plot(phi_a(theta),phi_b(theta),marker="s",color="orange")
# ax2.plot(phi_a(theta+np.pi),phi_b(theta+np.pi),marker="s",color="orange")
# ax2.set_xlabel(r"$\phi_a(\theta)$")
# ax2.set_ylabel(r"$\phi_b(\theta)$")
# ax2.set_xlim([0,np.pi])
# ax2.set_xticks([0,np.pi/2,np.pi])
# ax2.set_ylim([0,np.pi])
# ax2.set_yticks([0,np.pi/2,np.pi])
# ax2.set_xticklabels([r"$0$",r"$\frac{\pi}{2}$",r"$\pi$"])
# ax2.set_yticklabels([r"$0$",r"$\frac{\pi}{2}$",r"$\pi$"])


########### r-phi ################
ax2.plot(phi_a(x),r(x))
ax2.plot(phi_a(theta),r(theta),marker="s",color="orange")
ax2.plot(phi_a(theta+np.pi),r(theta+np.pi),marker="s",color="orange")
ax2.set_xlabel(r"$\phi_a(\theta)$")
ax2.set_ylabel(r"$r(\theta)$")
ax2.set_xlim([0,np.pi])
ax2.set_xticks([0,np.pi/2,np.pi])
ax2.set_ylim([0,1])
ax2.set_yticks([0,1/2,1])
ax2.set_xticklabels([r"$0$",r"$\frac{\pi}{2}$",r"$\pi$"])
ax2.set_yticklabels([r"$0$",r"$\frac{1}{2}$",r"$1$"])

############ affinity #############
# ax2.plot(x,z1(x)*z2(x),label="$z_1z_2$")
# ax2.plot(x,-np.log(z1(x)*z2(x)),label="A")
# ax2.set_ylim([-5,5])
# ax2.set_xlabel(r"$\theta$")
# ax2.legend()
# ax2.set_xticks([0,np.pi/2,np.pi,np.pi*3/2,np.pi*2])
# ax2.set_xticklabels([r"$0$",r"$\frac{\pi}{2}$",r"$\pi$",r"$\frac{3\pi}{2}$",r"$2\pi$"])


ax3.text(0.2, 0.66, f"z0={z_disc[0]:.2e}" , va="center", ha="center",size=10)
ax3.text(0.7, 0.66, f"z1={z_disc[1]:.2e}" , va="center", ha="center",size=10)
ax3.text(0.2, 0.33, f"z2={z_disc[2]:.3f}" , va="center", ha="center",size=10)
ax3.text(0.7, 0.33, f"z3={z_disc[3]:.3f}" , va="center", ha="center",size=10)
ax3.tick_params(labelbottom=False, labelleft=False)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)
ax3.tick_params('x', length=0, which='major')
ax3.tick_params('y', length=0, which='major')

plt.tight_layout()
plt.savefig(r"C:/Users/hyoshida/Desktop/floquetic/z_"+str(date)+"_"+str(ver)+".png")
plt.clf()
plt.close()
# def R(x):
#     # return (a1(theta)*(1-b2(theta)*dt)+a2(theta)*(1-a1(theta)*dt)+b1(theta)*(1-a2(theta)*dt)+b2(theta)*(1-b1(theta)*dt))*dt/Trace(x)
#     return np.abs(Trace(1)-2)/Trace(x)
# R(1)
# def root(x):
#     return np.complex(sqrt(-(x-z_disc[0])*(x-z_disc[1])*(x-z_disc[2])*(x-z_disc[3])/(x**2*(1-z_disc[0])*(1-z_disc[1])*(1-z_disc[2])*(1-z_disc[3]))))
# root(-0.1)
# def rho(x):
#     temp = 1/np.pi*(R(x)*root(x))/(1+R(x)**2*root(x)**2)*(1/(x-z_disc[0])+1/(x-z_disc[1])+1/(x-z_disc[2])+1/(x-z_disc[3])-2/x-2*(a1_R(theta)*b2_L(theta)+a2_R(theta)*b1_L(theta)-(a1_L(theta)*b2_R(theta)+a2_L(theta)*b1_R(theta))/x**2)*dt**2/Trace(x))
#     if (np.imag(temp)!=0):
#         # return str(nan)
#         return 0
#     else:
#         return np.abs(temp)
# rho(-0.5)
# N=100000
# dx = 0.00001
# Z1 = np.linspace(z_disc[0]-dx,z_disc[1]+dx,N)
# Z2 = np.linspace(z_disc[2]-dx,z_disc[3]+dx,N)
# Rho = []
# for z in Z1:
#     Rho.append(rho(z))
# # sum(Rho)*(Z1[1]-Z1[0])
# for z in Z2:
#     Rho.append(rho(z))
# # sum(Rho)*(Z2[1]-Z2[0])
#
# Z = np.append(Z1,Z2)
# plt.ylim([0,5])
# plt.xlim([-5,0])
# plt.plot(Z,Rho)
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
