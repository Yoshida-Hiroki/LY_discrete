# Focus on static 2level system
import numpy as np
from sympy import *
from scipy import integrate
import matplotlib.pyplot as plt
import csv
import cmath


type = r"\time_density"
date = "211019"
ver = "3"

# definitions of functions
# R = lambda x: np.ones_like(x)*1
r = lambda x: np.pi/7*np.sin(x)
phi_a = lambda x: np.pi*3/4+np.pi/10*np.cos(x)
phi_b = lambda x: np.pi/4+np.pi/10*np.cos(x)

# r_prime = lambda x: np.pi/4*np.cos(x)
# phi_a_prime = lambda x: -np.pi/4*np.sin(x)
# phi_b_prime = lambda x: np.pi/4*np.cos(x)

a_L = lambda x: 1/2*(1+r(x))*np.sin(phi_a(x)/2)**2
a_R = lambda x: 1/2*(1+r(x))*np.cos(phi_a(x)/2)**2
b_L = lambda x: 1/2*(1-r(x))*np.sin(phi_b(x)/2)**2
b_R = lambda x: 1/2*(1-r(x))*np.cos(phi_b(x)/2)**2

a = lambda x : a_R(x)+a_L(x)
b = lambda x : b_R(x)+b_L(x)

# continuum zeros
z1 = lambda x : (-((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))+np.sqrt(((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))**2-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x))
z2 = lambda x : (-((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))-np.sqrt(((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))**2-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x))

x = np.linspace(0,2*np.pi,1000)
z_p = np.max(z1(x))
z_m = np.min(z2(x))

print(z_p)
print(z_m)

######################## function checking (moving zeros plot) ################################

fig = plt.figure()
ax1 = plt.subplot2grid((1,2),(0,0))
ax2 = plt.subplot2grid((1,2),(0,1))

ax1.set_title("moving zeros")
ax1.plot(x,np.real(z1(x)),color="black",label="$z_1$")
ax1.plot(x,np.real(z2(x)),color="blue",label="$z_2$")
ax1.legend()
ax1.set_xlabel(r"$\theta$")
ax1.set_ylabel("$z$")
ax1.set_xticks([0,np.pi/2,np.pi,np.pi*3/2,np.pi*2])
ax1.set_xticklabels([r"$0$",r"$\frac{\pi}{2}$",r"$\pi$",r"$\frac{3\pi}{2}$",r"$2\pi$"])
# ax1.set_xlim([-10,0])
ax1.hlines(0,0,2*np.pi,color="gray")
# ax1.set_yscale("log")

ax2.set_title("$z_1z_2$")
# ax2.set_yscale("log")
ax2.plot(x,z1(x)*z2(x),label="$z_1z_2$")
ax2.plot(x,-np.log(z1(x)*z2(x)),label="A")
ax2.set_xlabel(r"$\theta$")
# ax2.plot(x,1/(a_R(x)*b_L(x)),label="$1/(a_Rb_L)$")
# ax2.plot(x,(a(x)+b(x))**2/(4*a_R(x)*b_L(x))-1-z1(x)*z2(x),label=r"$(\frac{a+b}{4a_Rb_L}-1-z_1z_2)$",color="black")
ax2.legend()
ax2.set_xticks([0,np.pi/2,np.pi,np.pi*3/2,np.pi*2])
ax2.set_xticklabels([r"$0$",r"$\frac{\pi}{2}$",r"$\pi$",r"$\frac{3\pi}{2}$",r"$2\pi$"])
plt.tight_layout()
# plt.show()
#
# ax1.plot(np.real(z1(x)),np.imag(z1(x)),linestyle="None",marker="s")
# ax1.plot(np.real(z2(x)),np.imag(z2(x)),linestyle="None",marker="s")

plt.savefig(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r"_zeros.png")
plt.clf()
plt.close()
# plt.show()

#################### zero density definition #########################

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
    sum = 0
    for i in range(2*M):
        sum += 1/(2*M)*(np.sqrt(-(x-z1(i))*(x-z2(i))/(x*(1-z1(i))*(1-z2(i)))))*(1/(x-z1(i))+1/(x-z2(i))-1/x)
    return -R/(np.pi*(1+R**2*root(x)**2))*sum

x1 = np.linspace(z_p+0.0001,-0.1,10000)
x2 = np.linspace(-50,z_m-0.0001,10000)

############# zero density plot ####################

fig = plt.figure()
ax1 = plt.subplot2grid((1,1),(0,0))
ax1.plot(x1,rho1(x1),color="blue",label=r"$\rho_+$")
ax1.plot(x2,rho2(x2),color="g",label=r"$\rho_-$")
ax1.plot(np.real(z1(x)),np.imag(z1(x)),linestyle="None",marker="s",markersize=2,color="blue",label=r"$z_1$")
ax1.plot(np.real(z2(x)),np.imag(z2(x)),linestyle="None",marker="s",markersize=2,color="g",label=r"$z_2$")
plt.legend()
ax1.set_xlabel(r"$z$")
ax1.set_ylabel(r"$\rho$")
# ax1.set_ylim([-0.05,1])
plt.savefig(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r"_density.png")
plt.clf()
plt.close()
# plt.show()

# print(integrate.quad(rho1,z_p+0.0001,-0.1))
# rho1(x1).sum()*(x1[1]-x1[0])

############# current density calculation ##############

dx1 = x1[1]-x1[0]
dx2 = x2[1]-x2[0]

def integ1(x,z):
    return rho1(x)*np.log((z-x)/(1-x))

def integ2(x,z):
    return rho2(x)*np.log((z-x)/(1-x))

def integ3(x,z):
    return rho1(x)*z/(z-x)

def integ4(x,z):
    return rho2(x)*z/(z-x)

def J(z):
    return (integ3(x1,z).sum()*dx1+integ4(x2,z).sum()*dx2-1)*0.5

def phi(z):
    return (integ1(x1,z).sum()*dx1+integ2(x2,z).sum()*dx2-np.log(z))*0.5-J(z)*np.log(z)

Chi = np.linspace(-4,5,100)
J_dat=[]
phi_dat=[]
for chi in Chi:
    J_dat.append(J(np.exp(chi)))
    phi_dat.append(phi(np.exp(chi)))


################## simulation #############################

iter = 100000
J_sim = [0]*(4*M+1)
for i in range(iter):
    n_before = 0
    n_after = 0
    j = 0
    # n = [0]*(2*M)
    # j = [0]*(2*M)
    for k in range(1,2*M):
        rand = np.random.rand()
        theta = k*2*np.pi/(2*M)

        n_after = (1-n_before)*0.5*(1+np.sign(b(theta)*dt-rand))+n_before*0.5*(1-np.sign(a(theta)*dt-rand))
        j += -(1-n_before)*0.5*(1+np.sign(b_R(theta)*dt-rand))+n_before*0.5*(1+np.sign(a_R(theta)*dt-rand))
        # n[k] = (1-n[k-1])*0.5*(1+np.sign(b(theta)*dt-rand))+n[k-1]*0.5*(1-np.sign(a(theta)*dt-rand))
        # j[k] = -(1-n[k-1])*0.5*(1+np.sign(b_R(theta)*dt-rand))+n[k-1]*0.5*(1+np.sign(a_R(theta)*dt-rand))

        n_before = n_after

    J_sim[int(j+2*M)] += 1
    # J_sim[int(np.sum(j)+2*M)] += 1

# plt.plot(np.linspace(-1,1,4*M+1),J_sim)
# S = iter*2/(2*M)
# MAX=np.max(J_sim)

phi_sim = [0]*(4*M+1)
for i in range(4*M+1):
    phi_sim[i] = np.log(J_sim[i]/iter)/(2*M)


################### plot ##########################
fig = plt.figure()
ax1 = plt.subplot2grid((1,1),(0,0))

ax1.plot(J_dat,phi_dat)
ax1.set_xlim([-0.4,0.1])
ax1.set_ylim([-0.3,0.05])
ax1.plot(np.linspace(-1,1,4*M+1),phi_sim,linestyle="None",marker="+")
ax1.hlines(0,-0.4,0.1,color="gray")
ax1.vlines(0,-0.3,0.05,color="gray")
ax1.set_xlabel(r"$J$")
ax1.set_ylabel(r"$\phi$")
plt.savefig(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r"_current.png")
plt.clf()
plt.close()
# plt.show()

##################### csv ##########################

outfile = open(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".csv",'w', newline='')
writer = csv.writer(outfile)
writer.writerow(z1(x))
writer.writerow(z2(x))
writer.writerow(J_dat)
writer.writerow(phi_dat)
writer.writerow(np.linspace(-1,1,4*M+1))
writer.writerow(phi_sim)
outfile.close()

##################### txt ###########################
f = open(r"G:\マイドライブ\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".txt",'w')
f.write('R=1'+'\n'+'r=pi/7*sin(x)'+'\n'+'phi_a=pi*3/4+pi/10*cos(x)'+'\n'+'phi_b=pi/4+pi/10*cos(x)'+'\n\n'+'M='+str(M)+'\n'+'iteration='+str(iter)+'\n\n'+'-------scv data------'+"\n\n"+'z1:[0,2pi]'+"\n"+'z2:[0,2pi]'+"\n"+'J:Derived from density of zeros'+"\n"+'phi:Derived from density of zeros to be plotted with J above'+'\n'+"J_sim:4*M+1 points in [-1,1]"+"\n"+'phi_sim:Phi derived from simulation. to be plotted with J_sim above.'+"\n\n"+"------Remarks-------")
f.close()
