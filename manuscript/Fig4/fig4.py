import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt
from scipy.integrate import solve_ivp

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
p = 11200
N = p/delta
mu = 0.001
T0 = 40000
I0 = 0
E0 = 0
W0 = 0

beta = c*delta*R0/((mu*p-delta*R0)*T0)

eps = np.arange(0,1,0.01)

################################
################################ Theory
################################
########################### Defining function to compute the deterministic time to reach 1000 individuals
def system(t,z,betas,ps,deltas,cs,ks,mus):
    T, E, I, V, W = z
    return [-betas*T*V,betas*T*V-ks*E,ks*E-deltas*I,ps*mus*I-cs*V-betas*T*V,ps*(1-mus)*I-cs*W]

V0 = 1        #10,
t1, t2, t3, t4 = [], [], [], []
for j in range(len(eps)):
    print(j)
    beta2 = beta
    p2 = (1-eps[j])*p
    delta2 = delta
    c2 = c
    
    phi = 1-(c2/(c2+beta2*T0)+delta2/(mu*p2))**V0
    thresh = 2000*phi
    if(thresh >= 100):
        def hit_thresh(t,z,betas,ps,deltas,cs,ks,mus): return(thresh - z[3] - z[4])
        hit_thresh.terminal = True
        z0 = [T0,E0,I0,V0,W0]
        sol = solve_ivp(system,[0,50],z0,args=(beta2,p2,delta2,c2,k,mu),events=hit_thresh)
        t1.append(sol.t_events[0][0])

for j in range(len(eps)):
    print(j)
    beta2 = (1-eps[j])*beta
    p2 = p
    delta2 = delta
    c2 = c
    
    phi = 1-(c2/(c2+beta2*T0)+delta2/(mu*p2))**V0
    thresh = 2000*phi
    if(thresh >= 100):
        def hit_thresh(t,z,betas,ps,deltas,cs,ks,mus): return(thresh - z[3] - z[4])
        hit_thresh.terminal = True
        z0 = [T0,E0,I0,V0,W0]
        sol = solve_ivp(system,[0,50],z0,args=(beta2,p2,delta2,c2,k,mu),events=hit_thresh)
        t2.append(sol.t_events[0][0])    

for j in range(len(eps)):
    print(j)
    beta2 = beta
    p2 = p
    delta2 = delta/(1-eps[j])
    c2 = c
    
    phi = 1-(c2/(c2+beta2*T0)+delta2/(mu*p2))**V0
    thresh = 2000*phi
    if(thresh >= 100):
        def hit_thresh(t,z,betas,ps,deltas,cs,ks,mus): return(thresh - z[3] - z[4])
        hit_thresh.terminal = True
        z0 = [T0,E0,I0,V0,W0]
        sol = solve_ivp(system,[0,50],z0,args=(beta2,p2,delta2,c2,k,mu),events=hit_thresh)
        t3.append(sol.t_events[0][0])    

for j in range(len(eps)):
    print(j)
    beta2 = beta
    p2 = p
    delta2 = delta
    c2 = c/(1-eps[j])
    
    phi = 1-(c2/(c2+beta2*T0)+delta2/(mu*p2))**V0
    thresh = 2000*phi
    if(thresh >= 100):
        def hit_thresh(t,z,betas,ps,deltas,cs,ks,mus): return(thresh - z[3] - z[4])
        hit_thresh.terminal = True
        z0 = [T0,E0,I0,V0,W0]
        sol = solve_ivp(system,[0,50],z0,args=(beta2,p2,delta2,c2,k,mu),events=hit_thresh)
        t4.append(sol.t_events[0][0])   


plt.plot(eps[0:80],t1[0:80],linewidth=3,color='C0')
plt.plot(eps[0:80],t2[0:80],linewidth=3,color='C1')
plt.plot(eps[0:80],t3[0:80],linewidth=3,color='C2')
plt.plot(eps[0:80],t4[0:80],linewidth=3,color='C3')

#plt.axvspan(eps[len(t2)-1],1.,facecolor = 'grey',alpha=0.75)
#plt.axvspan(eps[len(t1)-1],1.,facecolor = 'grey',alpha=0.25)
#plt.ylim((0,15))
#plt.show()

################################
################################ Import Data
################################
################################ Reducing burst size
eps_label = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
V0 = [1]
eps_plot = np.arange(0,0.9,0.1)
for j in range(len(eps_plot)):
    for l in range(len(V0)):
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_0_model_1.txt' % (V0[l], eps_label[j]), delimiter=',')
        #plt.plot(eps_data[j],np.mean(data),'o',color='C0',markersize=10)
        error = np.transpose([[np.mean(data)-np.percentile(data,5),np.percentile(data,95)-np.mean(data)]])
        plt.errorbar(eps_plot[j]-0.015,np.mean(data),yerr=error,fmt='o',color='C0', markersize=10, elinewidth=2,capsize=4)
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_1_model_1.txt' % (V0[l], eps_label[j]), delimiter=',')
        #plt.plot(eps_data[j],np.mean(data),'o',color='C1',markersize=10)
        error = np.transpose([[np.mean(data)-np.percentile(data,5),np.percentile(data,95)-np.mean(data)]])
        plt.errorbar(eps_plot[j]-0.005,np.mean(data),yerr=error,fmt='o',color='C1', markersize=10, elinewidth=2,capsize=4)
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_2_model_1.txt' % (V0[l], eps_label[j]), delimiter=',')
        #plt.plot(eps_data[j],np.mean(data),'o',color='C2',markersize=10)
        error = np.transpose([[np.mean(data)-np.percentile(data,5),np.percentile(data,95)-np.mean(data)]])
        plt.errorbar(eps_plot[j]+0.005,np.mean(data),yerr=error,fmt='o',color='C3', markersize=10, elinewidth=2,capsize=4)
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_3_model_1.txt' % (V0[l], eps_label[j]), delimiter=',')
        #plt.plot(eps_data[j],np.mean(data),'o',color='C2',markersize=10)
        error = np.transpose([[np.mean(data)-np.percentile(data,5),np.percentile(data,95)-np.mean(data)]])
        plt.errorbar(eps_plot[j]+0.015,np.mean(data),yerr=error,fmt='o',color='C2', markersize=10, elinewidth=2,capsize=4)
        

################################ Plot
#plt.ylim((0,10))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()

