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
################################ Import Data
################################
################################ eps = 0.75
model_label = [0,1,2,3]
reps = np.arange(0,10,1)
colors = ['C0','C1','C3','C2']
t_plot = np.arange(0,5,0.01)

for j in range(len(model_label)):
    for k in range(len(reps)):
        data = genfromtxt('traj_time_LN_eps_075_ind_%s_sc_%s.txt' % (reps[k],model_label[j]), delimiter=',')
        plt.plot(data[:,0],data[:,1],color=colors[j],linewidth=1,alpha=0.75)

################################ eps = 0.5
#for j in range(len(model_label)):
#    data = genfromtxt('traj_time_LN_eps_05_sc_%s.txt' % (model_label[j]), delimiter=',')
#    prop = []
#    total = data[:,1].tolist().count(0.0)
#    for l in range(len(t_plot)):
#        prop.append(1-sum(map(lambda x : x>t_plot[l], data[:,0].tolist()))/1000)
#    plt.plot(t_plot,prop,'x',color=colors[j],markersize=6)

################################ Plot
#plt.ylim((0,10))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10)
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()

