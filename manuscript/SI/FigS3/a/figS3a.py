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
T = 6*10**6

beta = c*delta*R0/((p-delta*R0)*T)

################################
################################ Theory
################################
def exttime(sc,V0,eps):
    
    if (sc==0):
        p2 = p*(1-eps)
        beta2 = beta
        delta2 = delta
    
    if (sc==1):
        p2 = p
        beta2 = (1-eps)*beta
        delta2 = delta
    
    if (sc==2):
        p2 = p
        beta2 = beta
        delta2 = delta/(1-eps)
    
    def g(z):
        return( c/(c+ beta2*T) + beta2*T/(c+beta2*T) * delta2/(p2+delta2) / (1-p2/(p2+delta2) * z) )
    
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
    
    phi = min(1,max(0,1-(c/(c+beta2*T)+delta2/p2)**V0))
    print(phi)
    
    return(np.array(proba)/(1-phi))

V0 = 10
eps = 0.9
pN = exttime(0,V0,eps)
pbeta = exttime(1,V0,eps)
pdelta = exttime(2,V0,eps)

gens = np.arange(1,50,1)
gens_plotN = [0.]*len(gens)

gens_plotN[0] = gens[0]/(beta*T + c)
gens_plotN[1] = gens[1]*(1/k + 1/delta + 1/(beta*T+c)) - 1/k - 1/delta
gens_plotN[2] = gens[2]*(1/k + 1/delta + 1/(beta*T+c)) - 3/(4*k) - 3/(4*delta)
gens_plotN[3] = gens[3]*(1/k + 1/delta + 1/(beta*T+c)) - 1/(2*k) - 1/(2*delta)
gens_plotN[4] = gens[4]*(1/k + 1/delta + 1/(beta*T+c)) - 1/(4*k) - 1/(4*delta)
#gens_plot[0:5] = gens[0:5]*(1/k + 1/delta + 1/(beta*T+c)) - 1/k - 1/delta
gens_plotN[5: ] = gens[5: ]*(1/k + 1/(beta*T + c) + 1/delta)

gens_plotbeta = [0.]*len(gens)
gens_plotbeta[0] = gens[0]/(beta*(1-eps)*T + c)
gens_plotbeta[1] = gens[1]*(1/k + 1/delta + 1/(beta*(1-eps)*T+c)) - 1/k - 1/delta
gens_plotbeta[2] = gens[2]*(1/k + 1/delta + 1/(beta*(1-eps)*T+c)) - 3/(4*k) - 3/(4*delta)
gens_plotbeta[3] = gens[3]*(1/k + 1/delta + 1/(beta*(1-eps)*T+c)) - 1/(2*k) - 1/(2*delta)
gens_plotbeta[4] = gens[4]*(1/k + 1/delta + 1/(beta*(1-eps)*T+c)) - 1/(4*k) - 1/(4*delta)
#gens_plot[0:5] = gens[0:5]*(1/k + 1/delta + 1/(beta*T+c)) - 1/k - 1/delta
gens_plotbeta[5: ] = gens[5: ]*(1/k + 1/(beta*(1-eps)*T + c) + 1/delta)

gens_plotdelta = [0]*len(gens)
gens_plotdelta[0] = gens[0]/(beta*T + c)
gens_plotdelta[1] = gens[1]*(1/k + (1-eps)/delta + 1/(beta*T+c)) - 1/k - (1-eps)/delta
gens_plotdelta[2] = gens[2]*(1/k + (1-eps)/delta + 1/(beta*T+c)) - 3/(4*k) - 3*(1-eps)/(4*delta)
gens_plotdelta[3] = gens[3]*(1/k + (1-eps)/delta + 1/(beta*T+c)) - 1/(2*k) - (1-eps)/(2*delta)
gens_plotdelta[4] = gens[4]*(1/k + (1-eps)/delta + 1/(beta*T+c)) - 1/(4*k) - (1-eps)/(4*delta)
#gens_plot[0:5] = gens[0:5]*(1/k + 1/delta + 1/(beta*T+c)) - 1/k - 1/delta
gens_plotdelta[5: ] = gens[5: ]*(1/k + 1/(beta*T + c) + (1-eps)/delta)

plt.plot(gens_plotN,np.array(pN)/(1/k+1/(beta*T+c)+1/delta),linewidth=3)
plt.plot(gens_plotbeta,np.array(pbeta)/(1/k+1/((1-eps)*beta*T+c)+1/delta),linewidth=3)
plt.plot(gens_plotdelta,np.array(pdelta)/(1/k + 1/(beta*T + c) + (1-eps)/delta),linewidth=3)

################################
################################ Import Data
################################
################################ Reducing burst size
data1 = genfromtxt('ext_LN_V0_10_eps_0.9_sc_0_model_1.txt', delimiter=',')
data2 = genfromtxt('ext_LN_V0_10_eps_0.9_sc_1_model_1.txt', delimiter=',')
data4 = genfromtxt('ext_LN_V0_10_eps_0.9_sc_3_model_1.txt', delimiter=',')

plt.hist(data1,bins=4*len(gens),density=True,color='C0',alpha = 0.7)  
plt.hist(data4,bins=len(gens),density=True,color='C2',alpha = 0.7)
plt.hist(data2,bins=8*len(gens),density=True,color='C1',alpha = 0.7)
#plt.axvline(gens[14]*(1/k + 1/(beta*T + c) + 1/delta),color='C0',linewidth=2)
#plt.axvline(np.percentile(my_data,95),color='C1',linewidth=2)

plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)

plt.xlim((0,20))
plt.show()

