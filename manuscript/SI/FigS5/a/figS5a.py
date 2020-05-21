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
        
        if (sc == 1):
            beta2 = (1-eps[i])*beta
            p2 = p
           
        phi = 1-(c/(c+beta2*T)+delta/p2)**V0
        
        t = 0
        
        if (det == 0):
            thresh = 2000*phi
        
        else:
            thresh = 2000
        
        print(eps[i])
        ################# ODE simulation (Euler scheme) - full model
        ######## phi >0.02 to speed up simulations!
        if (det == 1):
            if (phi>0.02):
                while (V < thresh):
                    T = max(0,T - beta2 * T * V * dt)
                    E = max(0,E + beta2 * T * V * dt - k * E *dt)
                    I = max(0,I + k * E *dt - delta *I * dt)
                    V = max(0,V + p2 * I * dt - c * V * dt - beta2 * T * V * dt)
                    t += dt
            
                tpred.append(t)
            
            else:
                break
        
        else:
            if (phi > 0):
                while (V < thresh):
                    T = max(0,T - beta2 * T * V * dt)
                    E = max(0,E + beta2 * T * V * dt - k * E *dt)
                    I = max(0,I + k * E *dt - delta *I * dt)
                    V = max(0,V + p2 * I * dt - c * V * dt - beta2 * T * V * dt)
                    t += dt
            
                tpred.append(t)
            
            else:
                break
        
    return(tpred)


V0 = [1]        #10,
styles = ['-','dotted'] #'--',
for j in range(len(V0)):
    t1 = euler(0,V0[j],0)
    t2 = euler(1,V0[j],0)
    det1 = euler(0,V0[j],1)
    det2 = euler(1,V0[j],1)
    plt.plot(eps[0:len(t1)],t1,linewidth=3,linestyle=styles[j],color='C0')
    plt.plot(eps[0:len(t2)],t2,linewidth=3,linestyle=styles[j],color='C1')
    plt.plot(eps[0:len(det1)],det1,linewidth=3,linestyle=styles[j],color='C0',alpha=0.5)
    plt.plot(eps[0:len(det2)],det2,linewidth=3,linestyle=styles[j],color='C1',alpha=0.5)

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
        plt.plot(eps_data[j],np.mean(data),'o',color='C0',markersize=10)
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_1_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        plt.plot(eps_data[j],np.mean(data),'o',color='C1',markersize=10)

################################ Plot
plt.ylim((0,15))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()

