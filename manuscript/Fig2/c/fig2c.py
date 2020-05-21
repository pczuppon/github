import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt
from scipy.special import lambertw
from scipy.optimize import fsolve

################################
################################ data
################################

os.chdir(os.path.realpath(''))

################################
################################ Parameter definitions
################################

################################ within-host dynamics
R0 = 7. 
k = 5
delta = 0.58
c = 10.
p = 13.2
B = p/delta
T = 4.9*10**5

beta = c*delta*R0/((p-delta*R0)*T)

################################ drug efficacy
eps1 = 0.5
eps2 = 0.25
eps3 = 0.75
V0 = np.arange(1,1000,1)


################################
################################ Import Data
################################

################################ neutral
my_data = genfromtxt('surv_LN_eps_0_sc_0_model_1.txt', delimiter=',')

data1 = []
eps_data1 = []
for i in range(len(my_data)):
    data1.append(float(my_data[i][1]))
    eps_data1.append(float(my_data[i][0]))

################################ Reducing burst size, eps = 0.2
my_data = genfromtxt('surv_LN_eps_025_sc_0_model_1.txt', delimiter=',')

data2 = []
eps_data2 = []
for i in range(len(my_data)):
    data2.append(float(my_data[i][1]))
    eps_data2.append(float(my_data[i][0]))


################################ Reducing infectivity, eps = 0.2
my_data = genfromtxt('surv_LN_eps_025_sc_1_model_1.txt', delimiter=',')

data3 = []
eps_data3 = []
for i in range(len(my_data)):
    data3.append(float(my_data[i][1]))
    eps_data3.append(float(my_data[i][0]))

################################ Reducing burst size, eps = 0.5
my_data = genfromtxt('surv_LN_eps_05_sc_0_model_1.txt', delimiter=',')

data4 = []
eps_data4 = []
for i in range(len(my_data)):
    data4.append(float(my_data[i][1]))
    eps_data4.append(float(my_data[i][0]))

################################ Reducing infectivity, eps = 0.5
my_data = genfromtxt('surv_LN_eps_05_sc_1_model_1.txt', delimiter=',')

data5 = []
eps_data5 = []
for i in range(len(my_data)):
    data5.append(float(my_data[i][1]))
    eps_data5.append(float(my_data[i][0]))


################################ Reducing burst size, eps = 0.8
my_data = genfromtxt('surv_LN_eps_075_sc_0_model_1.txt', delimiter=',')

data6 = []
eps_data6 = []
for i in range(len(my_data)):
    data6.append(float(my_data[i][1]))
    eps_data6.append(float(my_data[i][0]))

################################ Reducing infectivity, eps = 0.8
my_data = genfromtxt('surv_LN_eps_075_sc_1_model_1.txt', delimiter=',')

data7 = []
eps_data7 = []
for i in range(len(my_data)):
    data7.append(float(my_data[i][1]))
    eps_data7.append(float(my_data[i][0]))


################################
################################ Surv prob of virus
################################

################################ Reducing burst size
surv1 = []
surv12 = []
surv13 = []
survneu = []

for i in range(len(V0)):
    p2 = (1-eps1)*p
    p21 = (1-eps2)*p
    p22 = (1-eps3)*p
    
    surv1.append(max(0,1-(c/(c+beta*T)+delta/p2)**V0[i]))
    surv12.append(max(0,1-(c/(c+beta*T)+delta/p21)**V0[i]))
    surv13.append(max(0,1-(c/(c+beta*T)+delta/p22)**V0[i]))
    survneu.append(max(0,1-(c/(c+beta*T)+delta/p)**V0[i]))
    

################################ Reducing infectivity
surv2 = []
surv22 = []
surv23 = []

for i in range(len(V0)):
    beta2 = (1-eps1)*beta
    beta21 = (1-eps2)*beta
    beta22 = (1-eps3)*beta
    
    surv2.append(max(0,1-(c/(c+beta2*T)+delta/p)**V0[i]))
    surv22.append(max(0,1-(c/(c+beta21*T)+delta/p)**V0[i]))
    surv23.append(max(0,1-(c/(c+beta22*T)+delta/p)**V0[i]))

################################
################################ Plot surv prob
################################
plt.semilogx(V0,np.array(surv1)/np.array(survneu),linewidth=3,color='C0',linestyle='dashed')
plt.plot(eps_data2,np.array(data2)/np.array(data1),'o',color='C0',markersize=10)
plt.plot(V0,np.array(surv2)/np.array(survneu),linewidth=3,color='C1',linestyle='dashed')
plt.plot(eps_data3,np.array(data3)/np.array(data1),'o',color='C1',markersize=10)
plt.plot(V0,np.array(surv12)/np.array(survneu),linewidth=3,color='C0',linestyle='solid')
plt.plot(eps_data4,np.array(data4)/np.array(data1),'o',color='C0',markersize=10)
plt.plot(V0,np.array(surv22)/np.array(survneu),linewidth=3,color='C1',linestyle='solid')
plt.plot(eps_data5,np.array(data5)/np.array(data1),'o',color='C1',markersize=10)
plt.plot(V0,np.array(surv13)/np.array(survneu),linewidth=3,color='C0',linestyle='dotted')
plt.plot(eps_data6,np.array(data6)/np.array(data1),'o',color='C0',markersize=10)
plt.plot(V0,np.array(surv23)/np.array(survneu),linewidth=3,color='C1',linestyle='dotted')
plt.plot(eps_data7,np.array(data7)/np.array(data1),'o',color='C1',markersize=10)
#plt.plot(V0,survneu,linewidth=2,color='black')
#plt.plot(eps_data1,data1,'o',color='black',markersize=10)

plt.ylim((-0.05,1.05))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()

