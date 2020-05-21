import numpy as np
import matplotlib.pyplot as plt
import os
from numpy import genfromtxt
from scipy.special import binom 

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
p = 13.2
N = p/delta
T0 = 6*10**6
I0 = 0
E0 = 0

beta = c*delta*R0/((p-delta*R0)*T0)

eps = np.arange(0,1,0.01)

################################
################################ Theory
################################
########################### Defining function to compute the deterministic time to reach 1000 individuals
def euler(sc,V0,det):
    tpred = []
    for i in range(len(eps)):
        T = T0
        E = E0
        I = I0
        V = V0
        
        dt = 0.01
        
        if (sc == 0):
            beta2 = beta
            p2 = (1-eps[i])*p
            delta2 = delta
        
        if (sc == 1):
            beta2 = (1-eps[i])*beta
            p2 = p
            delta2 = delta
        
        if (sc == 2):
            beta2 = beta
            p2 = p
            delta2 = delta/(1-eps[i])
           
        phi = 1-(c/(c+beta2*T)+delta2/p2)**V0
        
        t = 0
        
        if (det == 0):
            thresh = 2000*phi
        
        else:
            thresh = 2000
        
        ################# ODE simulation (Euler scheme) - full model
        if (phi>0):
            while (V < thresh):
                Tnew = max(0,T - beta2 * T * V * dt)
                Enew = max(0,E + beta2 * T * V * dt - k * E *dt)
                Inew = max(0,I + k * E *dt - delta2 *I * dt)
                Vnew = max(0,V + p2 * I * dt - c * V * dt - beta2 * T * V * dt)
                t += dt
                
                T = Tnew
                E = Enew
                I = Inew
                V = Vnew
            
            tpred.append(t)
        
        else:
            break
        
    return(tpred)


V0 = [10]        #10,
styles = ['-','dotted'] #'--',
for j in range(len(V0)):
    t1 = euler(0,V0[j],0)
    t2 = euler(1,V0[j],0)
    t3 = euler(2,V0[j],0)
    #det1 = euler(0,V0[j],1)
    #det2 = euler(1,V0[j],1)
    plt.plot(eps[0:len(t1)],t1,linewidth=3,linestyle=styles[j],color='C0')
    plt.plot(eps[0:len(t2)],t2,linewidth=3,linestyle=styles[j],color='C1')
    plt.plot(eps[0:len(t3)],t3,linewidth=3,linestyle=styles[j],color='C2')
    #plt.plot(eps[0:len(t1)],[t1[0]]*(len(t1)),linewidth=2,linestyle='dashed',color='grey')
    #plt.plot(eps[0:len(t2)],[t2[0]]*(len(t2)),linewidth=2,linestyle='dashed',color='grey')
    #plt.plot(eps[0:len(t1)],det1,linewidth=3,linestyle=styles[j],color='black')
    #plt.plot(eps[0:len(t2)],det2,linewidth=3,linestyle=styles[j],color='grey')

plt.axvspan(eps[len(t2)-1],1.,facecolor = 'grey',alpha=0.75)
plt.axvspan(eps[len(t1)-1],1.,facecolor = 'grey',alpha=0.25)
#plt.ylim((0,15))
#plt.show()

################################
################################ Import Data
################################
################################ Reducing burst size
eps_data = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
for j in range(len(eps_data)):
    for k in range(len(V0)):
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_0_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        #plt.plot(eps_data[j],np.mean(data),'o',color='C0',markersize=10)
        plt.errorbar(eps_data[j],np.mean(data),yerr=2*np.std(data),fmt='o',color='C0', markersize=10, elinewidth=2,capsize=4)
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_1_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        #plt.plot(eps_data[j],np.mean(data),'o',color='C1',markersize=10)
        plt.errorbar(eps_data[j],np.mean(data),yerr=2*np.std(data),fmt='o',color='C1', markersize=10, elinewidth=2,capsize=4)
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_3_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        #plt.plot(eps_data[j],np.mean(data),'o',color='C2',markersize=10)
        plt.errorbar(eps_data[j],np.mean(data),yerr=2*np.std(data),fmt='o',color='C2', markersize=10, elinewidth=2,capsize=4)

################################ Plot
plt.ylim((0,30))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()

