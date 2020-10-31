import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt

################################
################################ data
################################

os.chdir(os.path.realpath(''))

################################
################################ Parameter definitions
################################

################################ within-host dynamics
R0 = 7.69 
k = 5
delta = 0.595
c = 10.
<<<<<<< HEAD
p = 11200
B = p/delta
mu = 0.001
T = 4*10**4
=======
p = 112000
mu = 0.001
B = mu*p/delta
T = 4000
>>>>>>> master

beta = c*delta*R0/((mu*p-delta*R0)*T)

################################ drug efficacy
eps = np.arange(0,1,0.0001)
<<<<<<< HEAD

=======
V0 = 1
>>>>>>> master

################################
################################ Import Data
################################

<<<<<<< HEAD
################################ Reducing burst size
my_data = genfromtxt('surv_LN_I0_1_sc_0_model_1.txt', delimiter=',')
=======
################################ Reducing burst size (2x)
my_data = genfromtxt('comb_HN_V0_%s_sc_0.txt' % V0, delimiter=',')
>>>>>>> master

data1 = []
eps_data1 = []
for i in range(len(my_data)):
    data1.append(float(my_data[i][1]))
    eps_data1.append(float(my_data[i][0]))

<<<<<<< HEAD
################################ Reducing infectivity
my_data = genfromtxt('surv_LN_I0_1_sc_1_model_1.txt', delimiter=',')
=======
################################ Reducing burst size + infectivity
my_data = genfromtxt('comb_HN_V0_%s_sc_1.txt' % V0, delimiter=',')
>>>>>>> master

data2 = []
eps_data2 = []
for i in range(len(my_data)):
    data2.append(float(my_data[i][1]))
    eps_data2.append(float(my_data[i][0]))

<<<<<<< HEAD
=======
################################ Reducing infectivity (2x)
my_data = genfromtxt('comb_HN_V0_%s_sc_2.txt' % V0, delimiter=',')

data3 = []
eps_data3 = []
for i in range(len(my_data)):
    data3.append(float(my_data[i][1]))
    eps_data3.append(float(my_data[i][0]))
>>>>>>> master


################################
################################ Surv prob of virus
################################

################################ Reducing burst size
surv1 = []

for i in range(len(eps)):
    p2 = (1-eps[i])*p*mu
<<<<<<< HEAD
    
    surv1.append(max(0,p2*(1-(c/(c+beta*T)+delta/p2))/(p2+delta-p2*(c/(c+beta*T)+delta/p2))))

=======
    delta2 = delta/(1-eps[i])
    surv1.append(max(0,1-(c/(c+beta*T)+delta2/p2)**V0))
>>>>>>> master

################################ Reducing infectivity
surv2 = []
<<<<<<< HEAD
=======
for i in range(len(eps)):
    p2 = (1-eps[i])*p*mu
    beta2 = beta*(1-eps[i])
    surv2.append(max(0,1-(c/(c+beta2*T)+delta/p2)**V0))
>>>>>>> master

for i in range(len(eps)):
<<<<<<< HEAD
    beta2 = (1-eps[i])*beta
    
    surv2.append(max(0,p*mu*(1-(c/(c+beta2*T)+delta/(p*mu)))/(p*mu+delta-p*mu*(c/(c+beta2*T)+delta/(p*mu)))))


=======
    c2 = c/(1-eps[i])
    beta2 = beta*(1-eps[i])
    surv3.append(max(0,1-(c2/(c2+beta2*T)+delta/(mu*p))**V0))

surv4 = []
surv5 = []
for i in range(len(eps)):
    p2 = p*(1-eps[i])*mu
    beta2 = beta*(1-eps[i])
    surv4.append(max(0,1-(c/(c+beta*T)+delta/p2)**V0))
    surv5.append(max(0,1-(c/(c+beta2*T)+delta/(mu*p))**V0))
>>>>>>> master

################################
################################ Plot surv prob
################################
plt.plot(eps,surv1,linewidth=3,color='C0')
plt.plot(eps_data1,data1,'o',color='C0',markersize=10)
plt.plot(eps,surv2,linewidth=3,color='C1')
plt.plot(eps_data2,data2,'o',color='C1',markersize=10)
#plt.plot(eps,surv11,linewidth=3,color='C0',linestyle='dashed')
#plt.plot(eps_data3,data3,'o',color='C0',markersize=10)
#plt.plot(eps,surv21,linewidth=3,color='C1',linestyle='dashed')
#plt.plot(eps_data4,data4,'o',color='C1',markersize=10)
#plt.plot(eps,surv12,linewidth=3,color='C0',linestyle='dotted')
#plt.plot(eps_data5,data5,'o',color='C0',markersize=10)
#plt.plot(eps,surv22,linewidth=3,color='C1',linestyle='dotted')
#plt.plot(eps_data6,data6,'o',color='C1',markersize=10)

plt.ylim((-0.05,0.4))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()
