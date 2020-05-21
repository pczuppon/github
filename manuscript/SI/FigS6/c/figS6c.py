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
R0 = 4.3 
k = 5
delta = 0.58
c = 10.
p = 11.9
B = p/delta
T = 6000000

beta = R0*c*delta/(T*(p-R0*delta))

################################ drug efficacy
eps = np.arange(0,1,0.0001)


################################
################################ Import Data
################################

################################ Short
################################

################################ Reducing burst size
my_data = genfromtxt('surv_LR_LN_V0_1_sc_0_model_1.txt', delimiter=',')

data1 = []
eps_data1 = []
for i in range(len(my_data)):
    data1.append(float(my_data[i][1]))
    eps_data1.append(float(my_data[i][0]))

################################ Reducing infectivity
my_data = genfromtxt('surv_LR_LN_V0_1_sc_1_model_1.txt', delimiter=',')

data2 = []
eps_data2 = []
for i in range(len(my_data)):
    data2.append(float(my_data[i][1]))
    eps_data2.append(float(my_data[i][0]))

################################ Reducing burst size, I0=10
my_data = genfromtxt('surv_LR_LN_V0_10_sc_0_model_1.txt', delimiter=',')

data3 = []
eps_data3 = []
for i in range(len(my_data)):
    data3.append(float(my_data[i][1]))
    eps_data3.append(float(my_data[i][0]))

################################ Reducing infectivity, I0=10
my_data = genfromtxt('surv_LR_LN_V0_10_sc_1_model_1.txt', delimiter=',')

data4 = []
eps_data4 = []
for i in range(len(my_data)):
    data4.append(float(my_data[i][1]))
    eps_data4.append(float(my_data[i][0]))

################################ Reducing burst size, I0=100
my_data = genfromtxt('surv_LR_LN_V0_100_sc_0_model_1.txt', delimiter=',')

data5 = []
eps_data5 = []
for i in range(len(my_data)):
    data5.append(float(my_data[i][1]))
    eps_data5.append(float(my_data[i][0]))

################################ Reducing infectivity, I0=10
my_data = genfromtxt('surv_LR_LN_V0_100_sc_1_model_1.txt', delimiter=',')

data6 = []
eps_data6 = []
for i in range(len(my_data)):
    data6.append(float(my_data[i][1]))
    eps_data6.append(float(my_data[i][0]))

################################
################################ Surv prob of virus
################################

################################ Reducing burst size
surv1 = []
surv11 = []
surv12 = []
l1 =[]
for i in range(len(eps)):
    p2 = (1-eps[i])*p
    
    surv1.append(max(0,1-(c/(c+beta*T)+delta/p2)))
    surv11.append(max(0,1-(c/(c+beta*T)+delta/p2)**10))
    surv12.append(max(0,1-(c/(c+beta*T)+delta/p2)**100))

################################ Reducing infectivity
surv2 = []
surv21 = []
surv22 = []
l2 = []
for i in range(len(eps)):
    beta2 = (1-eps[i])*beta
    
    surv2.append(max(0,1-(c/(c+beta2*T)+delta/p)))
    surv21.append(max(0,1-(c/(c+beta2*T)+delta/p)**10))
    surv22.append(max(0,1-(c/(c+beta2*T)+delta/p)**100))


################################
################################ Plot surv prob
################################
plt.plot(eps,surv1,linewidth=2,color='C0')
plt.plot(eps_data1,data1,'o',color='C0',markersize=10)
plt.plot(eps,surv2,linewidth=2,color='C1')
plt.plot(eps_data2,data2,'o',color='C1',markersize=10)
plt.plot(eps,surv11,linewidth=2,color='C0',linestyle='dashed')
plt.plot(eps_data3,data3,'o',color='C0',markersize=10)
plt.plot(eps,surv21,linewidth=2,color='C1',linestyle='dashed')
plt.plot(eps_data4,data4,'o',color='C1',markersize=10)
plt.plot(eps,surv12,linewidth=2,color='C0',linestyle='dotted')
plt.plot(eps_data5,data5,'o',color='C0',markersize=10)
plt.plot(eps,surv22,linewidth=2,color='C1',linestyle='dotted')
plt.plot(eps_data6,data6,'o',color='C1',markersize=10)

plt.ylim((-0.05,1.05))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()

