# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from scipy import integrate
import csv
import cmath
import time


type = r"\floquetic_zeros"
date = "211223"
ver = "N1_2"

z = Symbol('z')

a_L = lambda x: 0.3
# a_R = lambda x: 0.3
a_R = lambda x: 0.4
b_L = lambda x: 0.2
# b_R = lambda x: 0.2
b_R = lambda x: 0.1

a = lambda x : a_R(x)+a_L(x)
b = lambda x : b_R(x)+b_L(x)

z1 = lambda x : (-((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))+np.sqrt(((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))**2-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x))
z2 = lambda x : (-((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))-np.sqrt(((b(x)-a(x))**2/4+b_L(x)*a_L(x)+b_R(x)*a_R(x))**2-4*a_R(x)*a_L(x)*b_R(x)*b_L(x)))/(2*a_R(x)*b_L(x))

dt = 1

theta = 0
W_1 = lambda z : np.array([[1-b(theta)*dt,a_L(theta)*dt+a_R(theta)*dt*z],[b_L(theta)*dt+b_R(theta)*dt/z,1-a(theta)*dt]])
W_2 = lambda z : np.eye(2)

U_2 = lambda z : np.dot(W_2(z),W_1(z))

Trace = lambda z : np.trace(U_2(z))
Det = lambda z : U_2(z)[0][0]*U_2(z)[1][1]-U_2(z)[0][1]*U_2(z)[1][0]
z_1 = list(solveset(np.trace(W_1(z))**2-4*(W_1(z)[0][0]*W_1(z)[1][1]-W_1(z)[0][1]*W_1(z)[1][0]),z))
# z_2 = list(solveset(np.trace(W_2(z))**2-4*(W_2(z)[0][0]*W_2(z)[1][1]-W_2(z)[0][1]*W_2(z)[1][0]),z))

z_disc = z_1
z_disc

################## file make ##########################
f_zero = open(r"C:\Users\hyoshida\Desktop\floquetic\zero_"+str(date)+"_"+str(ver)+r".dat",'w')
f_zero.write(str(z_disc[0])+" "+str(z_disc[1]))
f_zero.close()
