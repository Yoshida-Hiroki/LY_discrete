# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import csv
import cmath


type = r"\discrete3_1"
date = "211011"
ver = "1"

# transition matrix elements
a1_R = 0.4
a1_L = 0.3
a2_R = 0.2
a2_L = 0.3
a3_R = 0.4
a3_L = 0.1
b1_R = 0.7
b1_L = 0.2
b2_R = 0.2
b2_L = 0.4
b3_R = 0.5
b3_L = 0.1


a1 = a1_R+a1_L
b1 = b1_R+b1_L
a2 = a2_R+a2_L
b2 = b2_R+b2_L
a3 = a3_R+a3_L
b3 = b3_R+b3_L

def c2(z):
    return a1+a2+a3+b1+b2+b3

def c1(z):
    return -((a1_L+a1_R*z)*(b1_L+b1_R/z)+a2*b2+a3*b3-(a1+b3)*(a2+a3)-(b1+b2)*(a1+b3)-(b1+b2)*(a1+b3))

def c0(z):
    return -((a1_L+a1_R*z)*b2*a3+(b1_L+b1_R/z)*a2*b3-(b1+b2)*(a1+b3)*(a2+a3)+(a1_L+a1_R*z)*(b1_L+b1_R/z)*(a2+a3)+a2*b2*(a1+b3)+a3*b3*(b1+b2))

def X(z):
    return (9*c2(z)*c1(z)-27*c0(z)-2*c2(z)**3)/54

def Y(z):
    return (3*c1(z)-c2(z)**2)/9

omega = np.exp(1.j*(2*np.pi)/3)

def Zp(z):
    return (X(z)+(X(z)**2+Y(z)**3)**(1/2))**(1/3)

def Zn(z):
    return (X(z)-(X(z)**2+Y(z)**3)**(1/2))**(1/3)

def Lambda1(z):
    return Zp(z)+Zn(z)

def Lambda2(z):
    return omega*Zp(z)+Zn(z)/omega

def Lambda3(z):
    return Zp(z)/omega+omega*Zn(z)

def Delta12(z):
    return np.abs(Lambda1(z)-Lambda2(z))

def Delta13(z):
    return np.abs(Lambda1(z)-Lambda3(z))

def Delta23(z):
    return np.abs(Lambda2(z)-Lambda3(z))
Lambda1(-2.468+0.00001*1.j)
Lambda2(-2.468+0.00001*1.j)
Lambda3(-2.468+0.00001*1.j)
x = np.linspace(-5,-2,1000)

plt.plot(x,(X(x)**2+Y(x)**3))
plt.xlabel("$\Re(z)$")
plt.ylabel("$X^2+Y^3$")
plt.axhline(0,0,1,color="gray")
plt.savefig(r"C:\Users\hyoshida\Desktop\XY.png")
plt.show()



fig = plt.figure()
ax = plt.subplot2grid((1,1),(0,0))
# for i in range(1,100):
ax.plot(x,Delta12(x+0.0001*1.j),color="red",label="$|\Lambda_1-\Lambda_2|$")
ax.plot(x,Delta13(x+0.0001*1.j),color="blue",label="$|\Lambda_1-\Lambda_2|$")
ax.plot(x,Delta23(x+0.0001*1.j),color="black",label="$|\Lambda_2-\Lambda_3|$")
ax.set_title("$|\Lambda_n^{\chi}-\Lambda_m^{\chi}|$ at $\Im(z)=0.0001$")
ax.set_xlabel("$\Re(z)$")
plt.axhline(0,0,1,color="gray")
plt.legend()
plt.savefig(r"C:\Users\hyoshida\Desktop\Lambda.png")
plt.show()
