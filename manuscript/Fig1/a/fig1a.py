import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt


os.chdir(os.path.realpath(''))

################## Parameter definitions
c = 10
k = 5
T0 = 1.33*10**3
mu = 0.001

dt = 0.0001
t = np.arange(0,40+dt,dt)

################# Initial condition
V0 = 1/30
E0 = 0
I0 = 0
W0 = 0

################# Import data + plot data (only loading patient ID, time, viral load)
my_data = genfromtxt('data.txt', delimiter=',', skip_header = 1,usecols=(0,2,4))

############### Retrieve patient IDs
p_id = np.unique(my_data[:,0])

###############
############### first patient
###############
tdata = []
vdata = []
ind = np.where((my_data == p_id[0]))[0]
for j in range(len(ind)):
     tdata.append(my_data[ind[j],1])
     vdata.append(10**(my_data[ind[j],2]))

plt.semilogy(tdata,vdata,marker='o',linestyle='dotted',markersize=10)

############### Individual estimate for the patients
Tfull, Efull, Ifull, Vfull, Wfull = [], [], [], [], []

Tfull.append(T0)
Efull.append(E0)
Ifull.append(I0)
Vfull.append(V0)
Wfull.append(W0)

R0 = 9.77
delta = 0.71
p = 11016
beta = c *delta*R0/((mu*p-delta*R0)*T0)

################# ODE simulation (Euler scheme) - full model
for i in range(len(t)-1):
    Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
    Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] *dt))
    Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta *Ifull[-1] * dt))
    Vfull.append(max(0,Vfull[-1] + mu * p * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
    Wfull.append(max(0,Wfull[-1] + (1-mu) * p * Ifull[-1] * dt - c * Wfull[-1] * dt))

Vtot = np.asarray(Vfull) + np.asarray(Wfull)
plt.plot(t,Vtot,linewidth=2,color='C0')

###################
################### second patient
###################
tdata = []                                                  
vdata = []
ind = np.where((my_data == p_id[2]))[0]
for j in range(len(ind)):
     tdata.append(my_data[ind[j],1])
     vdata.append(10**(my_data[ind[j],2]))

plt.plot(tdata,vdata,marker='o',linestyle='dotted',markersize=10,color='C1')

############### Individual estimate for the patients
Tfull, Efull, Ifull, Vfull, Wfull = [], [], [], [], []

Tfull.append(T0)
Efull.append(E0)
Ifull.append(I0)
Vfull.append(V0)
Wfull.append(W0)

R0 = 8.73
delta = 0.66
p = 11104
beta = c *delta*R0/((mu*p-delta*R0)*T0)

################# ODE simulation (Euler scheme) - full model
for i in range(len(t)-1):
    Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
    Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] *dt))
    Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta *Ifull[-1] * dt))
    Vfull.append(max(0,Vfull[-1] + p * mu * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
    Wfull.append(max(0,Wfull[-1] + (1-mu) * p * Ifull[-1] * dt - c * Wfull[-1] * dt))

Vtot = np.asarray(Vfull) + np.asarray(Wfull)
plt.plot(t,Vtot,linewidth=2,color='C1')

###################
################### third patient
###################
tdata = []
vdata = []
ind = np.where((my_data == p_id[12]))[0]
for j in range(len(ind)):
     tdata.append(my_data[ind[j],1])
     vdata.append(10**(my_data[ind[j],2]))

plt.plot(tdata,vdata,marker='o',linestyle='dotted',markersize=10,color='C2')

############### Individual estimate for the patients
Tfull, Efull, Ifull, Vfull, Wfull = [], [], [], [], []

Tfull.append(T0)
Efull.append(E0)
Ifull.append(I0)
Vfull.append(V0)
Wfull.append(W0)

R0 = 7.12
delta = 0.5
p = 11031
beta = c *delta*R0/((mu*p-delta*R0)*T0)

################# ODE simulation (Euler scheme) - full model
for i in range(len(t)-1):
    Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
    Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] *dt))
    Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta *Ifull[-1] * dt))
    Vfull.append(max(0,Vfull[-1] + p * mu * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
    Wfull.append(max(0,Wfull[-1] + (1-mu) * p * Ifull[-1] * dt - c * Wfull[-1] * dt))

Vtot = np.asarray(Vfull) + np.asarray(Wfull)
plt.plot(t,Vtot,linewidth=2,color='C2')

###################
################### fourth patient
###################
tdata = []
vdata = []
ind = np.where((my_data == p_id[7]))[0]
for j in range(len(ind)):
     tdata.append(my_data[ind[j],1])
     vdata.append(10**(my_data[ind[j],2]))

plt.plot(tdata,vdata,marker='o',linestyle='dotted',markersize=10,color='C3')

############### Individual estimate for the patients
Tfull, Efull, Ifull, Vfull, Wfull = [], [], [], [], []

Tfull.append(T0)
Efull.append(E0)
Ifull.append(I0)
Vfull.append(V0)
Wfull.append(W0)

R0 = 9.2
delta = 0.86
p = 11060
beta = c *delta*R0/((mu*p-delta*R0)*T0)

################# ODE simulation (Euler scheme) - full model
for i in range(len(t)-1):
    Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
    Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] *dt))
    Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta *Ifull[-1] * dt))
    Vfull.append(max(0,Vfull[-1] + p * mu * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
    Wfull.append(max(0,Wfull[-1] + (1-mu) * p * Ifull[-1] * dt - c * Wfull[-1] * dt))

Vtot = np.asarray(Vfull) + np.asarray(Wfull)
plt.plot(t,Vtot,linewidth=2,color='C3')

################# Load plot
#plt.ylim((0,5*10**8))
plt.xlim((0,30))
plt.ylim((10**(-2),10**8))
plt.plot(t,[10**(1.84)]*len(t), linewidth=2,linestyle='dotted',color='black')
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()


