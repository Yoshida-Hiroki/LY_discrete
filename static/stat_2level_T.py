# Focus on static 2level system
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import csv
import cmath


type = r"\discrete2T"
date = "211003"
ver = "4"

z = Symbol('z')

# initial distrib
p1 = 0.5
p2 = 0.5

# transition matrix elements
a_R = 0.4
a_L = 0.3
b_R = 0.4
b_L = 0.2
# 1/(a+b)=0.77

a = a_R+a_L
b = b_R+b_L

# continuum zeros
z1_cont = (-((b-a)**2/4+b_L*a_L+b_R*a_R)+np.sqrt(((b-a)**2/4+b_L*a_L+b_R*a_R)**2-4*a_R*a_L*b_R*b_L))/(2*a_R*b_L)
z2_cont = (-((b-a)**2/4+b_L*a_L+b_R*a_R)-np.sqrt(((b-a)**2/4+b_L*a_L+b_R*a_R)**2-4*a_R*a_L*b_R*b_L))/(2*a_R*b_L)

print(z1_cont)
print(z2_cont)

# discrete zeros
a_z = a_L + a_R * z
b_z = b_L + b_R / z

T = 100
def w(s):
    # params
    # T = 100
    M = s
    dt = T/M
    X = np.array([[1-b*dt,a_L*dt],[b_L*dt,1-a*dt]])
    V1 = np.array([[0,1],[0,0]])
    V2 = np.array([[0,0],[1,0]])
    w_z_1 = X +z*a_R*dt*V1+1/z*b_R*dt*V2    # normal way of definition is not good.
    w_n = np.eye(2)

    for i in range(M):
        w_n = np.dot(w_n,w_z_1)
    return simplify(w_n)

def Z(M):
    return np.dot(np.dot([1,1],w(M)),[[p1],[p2]])[0]

iter = 10
z_disc=[]
for i in range(1,iter+1):
    z_disc.append(list(solveset(Z(i),z)))

for i in range(iter):
    for j in range((int(i/2)+1)*2):
        z_disc[i][j] = complex(z_disc[i][j])

# z_prod = [[]*(iter-2)]
# for i in range(2,iter):
#     for j in range((int((i-2)/2)+1)):
#         z_prod[i-2][j].append(z_disc[i][j]*z_disc[i][i-i%2-1-j])

################### current ########################
# J_cont = (a+b)/2*(1-z1_cont*z2_cont)/(2*(1-z1_cont)*(1-z2_cont))
# def J_disc(T,s):
#     temp = 0
#     for i in range((int(s/2)+1)*2):
#         temp = temp + 1/(1-z_disc[s][i])
#     return (temp-(int(s/2)+1))/T
#
# J_cont
# J_disc(T,2)
#
# J = [0]*(iter+1)
# J[0]=J_cont
# for i in range(1,iter+1):
#     J[i] = J_disc(T,i-1)
#
# print(J)

fig = plt.figure()
ax1 = plt.subplot2grid((1,1),(0,0))

ax1.set_title('M=s, T='+str(T))
ax1.plot([np.real(z1_cont),np.real(z2_cont)],np.full(np.size([np.real(z1_cont),np.real(z2_cont)]),iter),linestyle="None",marker="s",color="blue",label="Continuous")
for i in range(iter):
    ax1.plot(np.real(z_disc[i]),np.full(np.size(np.real(z_disc[i])),i),linestyle="None",marker="o",markersize=2,color=[0,0,0],label="Discrete k="+str(i+1))
# for i in range(iter):
#     ax1.plot([np.real(z_disc[i][int(i/2)]),np.real(z_disc[i][int(i/2)+1])],[np.imag(z_disc[i][int(i/2)]),np.imag(z_disc[i][int(i/2)+1])],linestyle="None",marker="o",color=[0.9-i/20,0.9-i/20,0.9-i/20],label="Discrete k="+str(i+1))
# ax1.legend()
ax1.set_xlim([-10,3])
ax1.set_xlabel('$\Re (z)$')
ax1.set_ylabel('s')
ax1.set_yticks([x for x in range(iter+1)])
ax1.set_yticklabels([x for x in range(1,iter+1)]+['continuous'])
# plt.show()
plt.savefig(r"G:\??????????????????\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".png")

outfile = open(r"G:\??????????????????\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".csv",'w', newline='')
writer = csv.writer(outfile)
writer.writerows(z_disc)
writer.writerow([z1_cont, z2_cont])
outfile.close()

f = open(r"G:\??????????????????\research"+str(type)+"_"+str(date)+"_"+str(ver)+r".txt",'w')
f.write('p1='+str(p1)+'\n'+'p2='+str(p2)+'\n\n'+'a_R='+str(a_R)+'\n'+'a_L='+str(a_L)+'\n'+'b_R='+str(b_R)+'\n'+'b_L='+str(b_L)+'\n\n'+'T='+str(T)+'\n'+'M=s'+'\n\n'+'s=1 to '+str(iter))
f.close()
