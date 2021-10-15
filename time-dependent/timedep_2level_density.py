# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import csv
import cmath


type = r"\time_density"
date = "211015"
ver = "1"

# definitions of functions
R = lambda x: np.ones_like(x)*1
r = lambda x: np.pi/4*np.sin(x)/50
phi_a = lambda x: np.pi*3/4+np.pi/10*np.cos(x)
phi_b = lambda x: np.pi/4+np.pi/10*np.cos(x)

# r_prime = lambda x: np.pi/4*np.cos(x)
# phi_a_prime = lambda x: -np.pi/4*np.sin(x)
# phi_b_prime = lambda x: np.pi/4*np.cos(x)

a_L = lambda x: R(x)/2*(1+r(x))*np.sin(phi_a(x)/2)**2
a_R = lambda x: R(x)/2*(1+r(x))*np.cos(phi_a(x)/2)**2
b_L = lambda x: R(x)/2*(1-r(x))*np.sin(phi_b(x)/2)**2
b_R = lambda x: R(x)/2*(1-r(x))*np.cos(phi_b(x)/2)**2

a = lambda x : a_R(x)+a_L(x)
b = lambda x : b_R(x)+b_L(x)

# continuum zeros
z1 = lambda x : (-((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))+np.sqrt(((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))**2-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x))
z2 = lambda x : (-((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))-np.sqrt(((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))**2-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x))


######################## function checking (moving zeros plot) ################################

# fig = plt.figure()
# ax1 = plt.subplot2grid((1,2),(0,0))
# ax2 = plt.subplot2grid((1,2),(0,1))
#
# x = np.linspace(0,2*np.pi,1000)
# ax1.set_title("moving zeros")
# ax1.plot(x,np.real(-z1(x)),color="black",label="$-z_1$")
# ax1.plot(x,np.real(-z2(x)),color="blue",label="$-z_2$")
# ax1.legend()
# ax1.set_xlabel(r"$\theta$")
# ax1.set_xticks([0,np.pi/2,np.pi,np.pi*3/2,np.pi*2])
# ax1.set_xticklabels([r"$0$",r"$\frac{\pi}{2}$",r"$\pi$",r"$\frac{3\pi}{2}$",r"$2\pi$"])
# # ax1.set_xlim([-10,0])
# ax1.hlines(0,0,2*np.pi,color="gray")
# # ax1.set_yscale("log")
#
# # ax2.set_title("affinity $z_1z_2$")
# ax2.set_yscale("log")
# ax2.plot(x,z1(x)*z2(x),label="$z_1z_2$")
# ax2.plot(x,1/(a_R(x)*b_L(x)),label="$1/(a_Rb_L)$")
# ax2.plot(x,(a(x)+b(x))**2/(4*a_R(x)*b_L(x))-1-z1(x)*z2(x),label=r"$(\frac{a+b}{4a_Rb_L}-1-z_1z_2)$",color="black")
# ax2.legend()
# ax2.set_xticks([0,np.pi/2,np.pi,np.pi*3/2,np.pi*2])
# ax2.set_xticklabels([r"$0$",r"$\frac{\pi}{2}$",r"$\pi$",r"$\frac{3\pi}{2}$",r"$2\pi$"])
# plt.show()

# plt.plot(np.real(z1(x)),np.imag(z1(x)),linestyle="None",marker="s")
# plt.plot(np.real(z2(x)),np.imag(z2(x)),linestyle="None",marker="s")

###############################################

M=100
dt = 1
R = (0.5*dt)/(1-0.5*dt)

def root(x):
    sum = 0
    for i in range(2*M):
        sum += 1/(2*M)*(0.5*np.sqrt(-(x-z1(i))*(x-z2(i))/(x*(1-z1(i))*(1-z2(i)))))
    return sum

def rho1(x):
    sum = 0
    for i in range(2*M):
        sum += 1/(2*M)*(np.sqrt(-(x-z1(i))*(x-z2(i))/(x*(1-z1(i))*(1-z2(i)))))*(1/(x-z1(i))+1/(x-z2(i))-1/x)
    return R/(np.pi*(1+R**2*root(x)**2))*sum

def rho2(x):
    return -R*root(x)/(np.pi*(1+R**2*root(x)**2))*(1/(x-z1)+1/(x-z2)-1/x)
root(-0.5)


# # discrete zeros
# a_z = a_L + a_R * z
# b_z = b_L + b_R / z
#
# def w(s):
#     # params
#     T = 100
#     M = s
#     dt = T/M
#     X = np.array([[1-b*dt,a_L*dt],[b_L*dt,1-a*dt]])
#     V1 = np.array([[0,1],[0,0]])
#     V2 = np.array([[0,0],[1,0]])
#     w_z_1 = X +z*a_R*dt*V1+1/z*b_R*dt*V2    # normal way of definition is not good.
#     w_n = np.eye(2)
#
#     for i in range(s):
#         w_n = np.dot(w_n,w_z_1)
#     return simplify(w_n)
#
# def Z(M):
#     return np.dot(np.dot([1,1],w(M)),[[p1],[p2]])[0]
#
# iter = 16
# z_disc=[]
# for i in range(1,iter+1):
#     z_disc.append(list(solveset(Z(i),z)))
#
# for i in range(iter):
#     for j in range((int(i/2)+1)*2):
#         z_disc[i][j] = complex(z_disc[i][j])
#
# fig = plt.figure()
# ax1 = plt.subplot2grid((1,1),(0,0))
#
# ax1.plot([np.real(z1_cont),np.real(z2_cont)],np.full(np.size([np.real(z1_cont),np.real(z2_cont)]),iter),linestyle="None",marker="s",color="blue",label="Continuous",markersize=8)
# for i in range(iter):
#     ax1.plot(np.real(z_disc[i]),np.full(np.size(np.real(z_disc[i])),i),linestyle="None",marker="o",color=[0.9-i/25,0.9-i/25,0.9-i/25],label="Discrete k="+str(i+1))
# # for i in range(iter):
# #     ax1.plot([np.real(z_disc[i][int(i/2)]),np.real(z_disc[i][int(i/2)+1])],[np.imag(z_disc[i][int(i/2)]),np.imag(z_disc[i][int(i/2)+1])],linestyle="None",marker="o",color=[0.9-i/20,0.9-i/20,0.9-i/20],label="Discrete k="+str(i+1))
# # ax1.legend()
# ax1.set_xlim([-10,3])
# ax1.set_yticks([x for x in range(iter+1)])
# ax1.set_yticklabels([x for x in range(1,iter+1)]+['continuous'])
# # plt.show()
# plt.savefig(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".png")
#
# outfile = open(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".csv",'w', newline='')
# writer = csv.writer(outfile)
# writer.writerows(z_disc)
# writer.writerow([z1_cont, z2_cont])
# outfile.close()
