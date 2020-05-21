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
R0 = 7 
k = 5
delta = 0.58
c = 10.
p = 116
B = p/delta
T = 490000

beta = c*delta*R0/((p-delta*R0)*T)

################################ drug efficacy
eps = np.arange(0,1,0.0001)
V0 = 100

################################
################################ Import Data
################################

################################ Reducing burst size (2x)
my_data = genfromtxt('comb_HN_V0_%s_sc_0.txt' % V0, delimiter=',')

data1 = []
eps_data1 = []
for i in range(len(my_data)):
    data1.append(float(my_data[i][1]))
    eps_data1.append(float(my_data[i][0]))

################################ Reducing burst size + infectivity
my_data = genfromtxt('comb_HN_V0_%s_sc_1.txt' % V0, delimiter=',')

data2 = []
eps_data2 = []
for i in range(len(my_data)):
    data2.append(float(my_data[i][1]))
    eps_data2.append(float(my_data[i][0]))

################################ Reducing infectivity (2x)
my_data = genfromtxt('comb_HN_V0_%s_sc_2.txt' % V0, delimiter=',')

data3 = []
eps_data3 = []
for i in range(len(my_data)):
    data3.append(float(my_data[i][1]))
    eps_data3.append(float(my_data[i][0]))


################################
################################ Surv prob of virus
################################

################################ Combination therapy
surv1 = []
for i in range(len(eps)):
    p2 = (1-eps[i])*p
    delta2 = delta/(1-eps[i])
    surv1.append(max(0,1-(c/(c+beta*T)+delta2/p2)**V0))

surv2 = []
for i in range(len(eps)):
    p2 = (1-eps[i])*p
    beta2 = beta*(1-eps[i])
    surv2.append(max(0,1-(c/(c+beta2*T)+delta/p2)**V0))

surv3 = []
for i in range(len(eps)):
    c2 = c/(1-eps[i])
    beta2 = beta*(1-eps[i])
    surv3.append(max(0,1-(c2/(c2+beta2*T)+delta/p)**V0))

surv4 = []
surv5 = []
for i in range(len(eps)):
    p2 = p*(1-eps[i])
    beta2 = beta*(1-eps[i])
    surv4.append(max(0,1-(c/(c+beta*T)+delta/p2)**V0))
    surv5.append(max(0,1-(c/(c+beta2*T)+delta/p)**V0))

################################
################################ Plot surv prob
################################
plt.plot(eps,surv4,linewidth=3,color='C0')
plt.plot(eps,surv5,linewidth=3,color='C1')
plt.plot(eps,surv1,linewidth=3,color='black')
plt.plot(eps_data1,data1,'o',color='black',markersize=10)
plt.plot(eps,surv2,linewidth=3,color='black',linestyle='dashed')
plt.plot(eps_data2,data2,'o',color='black',markersize=10)
plt.plot(eps,surv3,linewidth=3,color='black',linestyle='dotted')
plt.plot(eps_data3,data3,'o',color='black',markersize=10)

plt.ylim((-0.05,1.05))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()

