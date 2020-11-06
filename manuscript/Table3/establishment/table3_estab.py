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
################################ Import Data
################################
################################ R0 = 7, LowN
print("R0=7.69,LowN")
V0 = [1,10]
eps_data = [0]
for j in range(len(eps_data)):
    print(eps_data[j])
    for k in range(len(V0)):
        print(V0[k])
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_0_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')[:,0]
        print(np.percentile(data,5))
        print(np.percentile(data,50))
        print(np.percentile(data,95))
        print("\n")
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_3_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')[:,0]
        print(np.percentile(data,5))
        print(np.percentile(data,50))
        print(np.percentile(data,95))
        print("\n")
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_1_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')[:,0]
        print(np.percentile(data,5))
        print(np.percentile(data,50))
        print(np.percentile(data,95))
        print("\n")
        data = genfromtxt('esttime_LN_V0_%s_eps_%s_sc_2_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')[:,0]
        print(np.percentile(data,5))
        print(np.percentile(data,50))
        print(np.percentile(data,95))
        print("\n")

################################ R0 = 7, HighN
print("R0=7.69,HighN")
V0 = [1,10]
eps_data = [0]
for j in range(len(eps_data)):
    print(eps_data[j])
    for k in range(len(V0)):
        print(V0[k])
        data = genfromtxt('esttime_HN_V0_%s_eps_%s_sc_0_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')[:,0]
        print(np.percentile(data,5))
        print(np.percentile(data,50))
        print(np.percentile(data,95))
        print("\n")
        data = genfromtxt('esttime_HN_V0_%s_eps_%s_sc_3_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')[:,0]
        print(np.percentile(data,5))
        print(np.percentile(data,50))
        print(np.percentile(data,95))
        print("\n")
        data = genfromtxt('esttime_HN_V0_%s_eps_%s_sc_1_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')[:,0]
        print(np.percentile(data,5))
        print(np.percentile(data,50))
        print(np.percentile(data,95))
        print("\n")
        data = genfromtxt('esttime_HN_V0_%s_eps_%s_sc_2_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')[:,0]
        print(np.percentile(data,5))
        print(np.percentile(data,50))
        print(np.percentile(data,95))
        print("\n")
