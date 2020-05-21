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
################################ R0 = 7, LN
print("R0=7, LN")
V0 = [10,100]
eps_data = [0.9]
for j in range(len(eps_data)):
    print(eps_data[j])
    for k in range(len(V0)):
        print(V0[k])
        data = genfromtxt('ext_LN_V0_%s_eps_%s_sc_0_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")
        data = genfromtxt('ext_LN_V0_%s_eps_%s_sc_3_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")
        data = genfromtxt('ext_LN_V0_%s_eps_%s_sc_1_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")

################################ R0 = 7, HN
print("R0=7,HN")
V0 = [10,100]
eps_data = [0.9]
for j in range(len(eps_data)):
    print(eps_data[j])
    for k in range(len(V0)):
        print(V0[k])
        data = genfromtxt('ext_HN_V0_%s_eps_%s_sc_0_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")
        data = genfromtxt('ext_HN_V0_%s_eps_%s_sc_3_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")
        data = genfromtxt('ext_HN_V0_%s_eps_%s_sc_1_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")

################################ R0 = 4.3
print("R0=4.3")
V0 = [100]
eps_data = [0.9]
for j in range(len(eps_data)):
    print(eps_data[j])
    for k in range(len(V0)):
        print(V0[k])
        data = genfromtxt('ext_LR_LN_V0_%s_eps_%s_sc_0_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")
        data = genfromtxt('ext_LR_LN_V0_%s_eps_%s_sc_3_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")
        data = genfromtxt('ext_LR_LN_V0_%s_eps_%s_sc_1_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")

################################ R0 = 12
print("R0=12")
V0 = [10]
eps_data = [0, 0.75, 0.9]
for j in range(len(eps_data)):
    print(eps_data[j])
    for k in range(len(V0)):
        print(V0[k])
        data = genfromtxt('ext_HR_LN_V0_%s_eps_%s_sc_0_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")
        data = genfromtxt('ext_HR_LN_V0_%s_eps_%s_sc_3_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")
        data = genfromtxt('ext_HR_LN_V0_%s_eps_%s_sc_1_model_1.txt' % (V0[k], eps_data[j]), delimiter=',')
        print(np.percentile(data,10))
        print(np.percentile(data,50))
        print(np.percentile(data,90))
        print("\n")


