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
p = 116
N = p/delta
T = 4.9*10**5

beta = c*delta*R0/((p-delta*R0)*T)

################################
################################ Theory
################################
def exttime(sc,V0,eps):
    
    if (sc==0):
        p2 = p*(1-eps)
        beta2 = beta
        c2 = c
    
    if (sc==1):
        p2 = p
        beta2 = (1-eps)*beta
        c2 = c
    
    if (sc==2):
        p2 = p
        beta2 = beta
        c2 = c/(1-eps)
    
    def g(z):
        return( c2/(c2+ beta2*T) + beta2*T/(c2+beta2*T) * delta/(p2+delta) / (1-p2/(p2+delta) * z) )
    
    gens = np.arange(1,50,1)
    proba = []
    proba.append(g(0)**V0)
    gn = g(0)
    
    indvec = np.arange(1,V0+1,1)
    
    for i in range(len(gens)-1):
        p_before = gn**(V0-indvec)
        p_at = (g(gn)-gn)**indvec
        poss = binom(V0,indvec)
        proba.append(np.sum(p_before*p_at*poss))
        gn = g(gn)
    
    phi = min(1,max(0,1-(c2/(c2+beta2*T)+delta/p2)**V0))
    print(phi)
    
    return(np.array(proba)/(1-phi))

V0 = 100
eps = 0.9
pc = exttime(2,V0,eps)
pbeta = exttime(1,V0,eps)

gens = np.arange(1,50,1)
gens_plotc = np.arange(1,50,1)

gens_plotc[0] = gens[0]/(beta*T + c/(1-eps))
gens_plotc[1] = gens[1]*(1/k + 1/delta + 1/(beta*T+c/(1-eps))) - 1/k - 1/delta
gens_plotc[2] = gens[2]*(1/k + 1/delta + 1/(beta*T+c/(1-eps))) - 3/(4*k) - 3/(4*delta)
gens_plotc[3] = gens[3]*(1/k + 1/delta + 1/(beta*T+c/(1-eps))) - 1/(2*k) - 1/(2*delta)
gens_plotc[4] = gens[4]*(1/k + 1/delta + 1/(beta*T+c/(1-eps))) - 1/(4*k) - 1/(4*delta)
#gens_plot[0:5] = gens[0:5]*(1/k + 1/delta + 1/(beta*T+c)) - 1/k - 1/delta
gens_plotc[5: ] = gens[5: ]*(1/k + 1/(beta*T + c/(1-eps)) + 1/delta)

gens_plotbeta = np.arange(1,50,1)
gens_plotbeta[0] = gens[0]/(beta*(1-eps)*T + c)
gens_plotbeta[1] = gens[1]*(1/k + 1/delta + 1/(beta*(1-eps)*T+c)) - 1/k - 1/delta
gens_plotbeta[2] = gens[2]*(1/k + 1/delta + 1/(beta*(1-eps)*T+c)) - 3/(4*k) - 3/(4*delta)
gens_plotbeta[3] = gens[3]*(1/k + 1/delta + 1/(beta*(1-eps)*T+c)) - 1/(2*k) - 1/(2*delta)
gens_plotbeta[4] = gens[4]*(1/k + 1/delta + 1/(beta*(1-eps)*T+c)) - 1/(4*k) - 1/(4*delta)
#gens_plot[0:5] = gens[0:5]*(1/k + 1/delta + 1/(beta*T+c)) - 1/k - 1/delta
gens_plotbeta[5: ] = gens[5: ]*(1/k + 1/(beta*(1-eps)*T + c) + 1/delta)

plt.plot(gens_plotc,np.array(pc)/(1/k + 1/(beta*T + c) + 1/delta),linewidth=3,color='C3')
plt.plot(gens_plotbeta,np.array(pbeta)/(1/k + 1/(beta*T + c) + 1/delta),linewidth=3,color='C1')
#plt.show()

################################
################################ Import Data
################################
################################ Reducing burst size
data1 = genfromtxt('ext_HN_V0_100_eps_0.9_sc_2_model_1.txt', delimiter=',')
data2 = genfromtxt('ext_HN_V0_100_eps_0.9_sc_1_model_1.txt', delimiter=',')

plt.hist(data1,bins=12*len(gens),density=True,color='C3',alpha = 0.7)  
plt.hist(data2,bins=12*len(gens),density=True,color='C1',alpha = 0.7)
#plt.axvline(gens[14]*(1/k + 1/(beta*T + c) + 1/delta),color='C0',linewidth=2)
#plt.axvline(np.percentile(my_data,95),color='C1',linewidth=2)

plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)

plt.xlim((0,2))
plt.show()

