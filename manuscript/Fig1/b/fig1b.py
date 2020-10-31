import numpy as np
import matplotlib.pyplot as plt

################## Parameter definitions
c = 10
delta = 0.595
k = 5
p = 11200
T0 = 1.33*10**7
mu = 10**(-3)

dt = 0.001
t = np.arange(0,10+dt,dt)

frac = [0.001,0.0015,0.002,0.003,0.004,0.005,0.006,0.008,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1,1.5,2,3,4,6,8,10,15,20,30,40,60,80,100]
################# Contour lines time to viral peak
days = [3,4,5,6,7,8,9]
res_days = np.zeros((len(days),len(frac)))

l, j = 0, 0

for l in range(len(days)):
    print("days")
    print(days[l])
    R0 = 1.1
    for j in range(len(frac)):
        print("frac")
        print(frac[j])
        update = 0
        while (update == 0):
            print(R0)
            beta = c *delta*R0/((p*mu-delta*R0)*frac[j]*T0/100)
            #
            ################# ODE Simulation
            ################# Initial condition
            V0 = 1
            E0 = 0
            I0 = 0
            NV0 = 0
            T0_upd = T0*frac[j]/100
            #
            ################# Trajectories
            Tfull, Efull, Ifull, Vfull, NVfull = [], [], [], [], []
            #
            Tfull.append(T0_upd)
            Efull.append(E0)
            Ifull.append(I0)
            Vfull.append(V0)
            NVfull.append(NV0)
            #
            i = 0
            for i in range(len(t)-1):
                Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
                Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] * dt))
                Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta * Ifull[-1] * dt))
                Vfull.append(max(0,Vfull[-1] + mu * p * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
                NVfull.append(max(0,NVfull[-1] + (1-mu) * p * Ifull[-1] * dt - c * NVfull[-1] * dt))
            #
            Vwork = np.asarray(Vfull)
            NVwork = np.asarray(NVfull)
            Vtot = Vwork + NVwork
            #
            if (np.max(Vtot)==Vtot[-1] or np.max(Vtot)==Vtot[0]):
                R0 += 0.1
            else:
                tmax = np.where(Vtot == np.max(Vtot))[0][0]*dt
                if (tmax <= days[l]):
                    res_days[l,j] = R0
                    update = 1
                else:
                    R0 += 0.1
                    if (R0 >= 20):
                        res_days[l,j]=21
                        update = 1

################# Contour lines peak viral load
load = [5,6,7,8,9,10]
res_load = np.zeros((len(load),len(frac)))

t = np.arange(0,20+dt,dt)
l, j = 0, 0

for l in range(len(load)):
    print("load")
    print(load[l])
    R0 = 18.7
    for j in range(len(frac)):
        print("frac")
        print(frac[j])
        update = 0
        while (update == 0):
            print(R0)
            beta = c *delta*R0/((mu*p-delta*R0)*frac[j]*T0/100)
            #
            ################# ODE Simulation
            ################# Initial condition
            V0 = 1/30
            E0 = 0
            I0 = 0
            T0_upd = T0*frac[j]/100
            NV0 = 0
            #
            ################# Trajectories
            Tfull, Efull, Ifull, Vfull, NVfull = [], [], [], [], []
            #
            Tfull.append(T0_upd)
            Efull.append(E0)
            Ifull.append(I0)
            Vfull.append(V0)
            NVfull.append(NV0)
            #
            i = 0
            for i in range(len(t)-1):
                Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
                Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] * dt))
                Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta * Ifull[-1] * dt))
                Vfull.append(max(0,Vfull[-1] + mu * p * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
                NVfull.append(max(0,NVfull[-1] + (1-mu) * p * Ifull[-1] *dt - c * NVfull[-1] * dt ))
            #
            Vwork = np.asarray(Vfull)
            NVwork = np.asarray(NVfull)
            Vtot = Vwork + NVwork
            #
            if (np.max(Vtot)==Vtot[-1]):
                res_load[l,j] = 22
                update = 1
            elif (np.max(Vtot) <= 10**(load[l]) and np.max(Vtot)!=Vtot[0]):
                if (R0 == 20):
                    res_load[l,j] = 22
                else:
                    res_load[l,j] = R0
                update = 1
            elif (np.max(Vtot)==Vtot[0]):
                res_load[l,j] = 0
                update = 1
                R0 = 20
            else:
                R0 -= 0.1


################# Load plot
prop = np.asarray(frac)/100

for i in range(len(days)):
    if (max(res_days[i,:])>18.7):
        ind_max = min(np.where(res_days[i,:]>18.7)[0])-1
    else:
        ind_max = len(res_days[i,:])
    plt.semilogx(prop[0:ind_max],res_days[i,0:ind_max],linewidth=2,color='C0')

for i in range(len(load)):
    if (min(res_load[i,:])<18.7):
        ind_min = min(np.where(res_load[i,:]<18.7)[0])
        ind_max = max(np.where(res_load[i,:]<18.7)[0])
    else:
        ind_min = 0
        ind_max = len(res_load[i,:])
    if (i == 1 or i==0):
        plt.plot(prop[ind_min:ind_max+1],res_load[i,ind_min:ind_max+1],linewidth=2,color='C1')
    else:
    	plt.plot(prop[ind_min-1:ind_max+1],res_load[i,ind_min-1:ind_max+1],linewidth=2,color='C1')

################ add data points
######### France
t = np.arange(0,7.5+dt,dt)

R0 = 3
frac_France = 10**(-2)
update = 0
while (update == 0):
    print(R0)
    beta = c *delta*R0/((p*mu-delta*R0)*frac_France*T0/100)
    #
    ################# ODE Simulation
    ################# Initial condition
    V0 = 1
    E0 = 0
    I0 = 0
    NV0 = 0
    T0_upd = T0*frac_France/100
    #
    ################# Trajectories
    Tfull, Efull, Ifull, Vfull, NVfull = [], [], [], [], []
    #
    Tfull.append(T0_upd)
    Efull.append(E0)
    Ifull.append(I0)
    Vfull.append(V0)
    NVfull.append(NV0)
    #
    i = 0
    for i in range(len(t)-1):
        Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
        Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] * dt))
        Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta * Ifull[-1] * dt))
        Vfull.append(max(0,Vfull[-1] + mu * p * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
        NVfull.append(max(0,NVfull[-1] + (1-mu) * p * Ifull[-1] * dt - c * NVfull[-1] * dt))
    #
    Vwork = np.asarray(Vfull)
    NVwork = np.asarray(NVfull)
    Vtot = Vwork + NVwork
    #
    if (np.max(Vtot)==Vtot[-1] or np.max(Vtot)==Vtot[0]):
        R0 += 0.1
    else:
        tmax = np.where(Vtot == np.max(Vtot))[0][0]*dt
        if (np.max(Vtot) < 2*10**6):
            if (tmax > 7):
                R0 += 0.1
            else:
                frac_France += 10**(-3)
                R0 = 3
        else:
            if (tmax <= 7):
                res_France = R0
                update = 1
            else:
                R0 += 0.1

#plt.plot(frac_France/100,res_France,marker="o",color="gray",markersize=10)

################# Singapore
R0 = 3
frac_Sing = 10**(-2)
update = 0
while (update == 0):
    print(R0)
    beta = c *delta*R0/((p*mu-delta*R0)*frac_Sing*T0/100)
    #
    ################# ODE Simulation
    ################# Initial condition
    V0 = 1
    E0 = 0
    I0 = 0
    NV0 = 0
    T0_upd = T0*frac_Sing/100
    #
    ################# Trajectories
    Tfull, Efull, Ifull, Vfull, NVfull = [], [], [], [], []
    #
    Tfull.append(T0_upd)
    Efull.append(E0)
    Ifull.append(I0)
    Vfull.append(V0)
    NVfull.append(NV0)
    #
    i = 0
    for i in range(len(t)-1):
        Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
        Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] * dt))
        Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta * Ifull[-1] * dt))
        Vfull.append(max(0,Vfull[-1] + mu * p * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
        NVfull.append(max(0,NVfull[-1] + (1-mu) * p * Ifull[-1] * dt - c * NVfull[-1] * dt))
    #
    Vwork = np.asarray(Vfull)
    NVwork = np.asarray(NVfull)
    Vtot = Vwork + NVwork
    #
    if (np.max(Vtot)==Vtot[-1] or np.max(Vtot)==Vtot[0]):
        R0 += 0.1
    else:
        tmax = np.where(Vtot == np.max(Vtot))[0][0]*dt
        if (np.max(Vtot) < 10**6):
            if (tmax > 5):
                R0 += 0.1
            else:
                frac_Sing += 10**(-3)
                R0 = 3
        else:
            if (tmax <= 5):
                res_Sing = R0
                update = 1
            else:
                R0 += 0.1

plt.plot(frac_Sing/100,res_Sing,marker="o",color="black",markersize=10)

################# Germany
R0 = 3
frac_Ger = 10**(-2)
update = 0
while (update == 0):
    print(R0)
    beta = c *delta*R0/((p*mu-delta*R0)*frac_Ger*T0/100)
    #
    ################# ODE Simulation
    ################# Initial condition
    V0 = 1
    E0 = 0
    I0 = 0
    NV0 = 0
    T0_upd = T0*frac_Ger/100
    #
    ################# Trajectories
    Tfull, Efull, Ifull, Vfull, NVfull = [], [], [], [], []
    #
    Tfull.append(T0_upd)
    Efull.append(E0)
    Ifull.append(I0)
    Vfull.append(V0)
    NVfull.append(NV0)
    #
    i = 0
    for i in range(len(t)-1):
        Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
        Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] * dt))
        Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta * Ifull[-1] * dt))
        Vfull.append(max(0,Vfull[-1] + mu * p * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
        NVfull.append(max(0,NVfull[-1] + (1-mu) * p * Ifull[-1] * dt - c * NVfull[-1] * dt))
    #
    Vwork = np.asarray(Vfull)
    NVwork = np.asarray(NVfull)
    Vtot = Vwork + NVwork
    #
    if (np.max(Vtot)==Vtot[-1] or np.max(Vtot)==Vtot[0]):
        R0 += 0.1
    else:
        if (np.max(Vtot) < 7*10**5):
            frac_Ger += 10**(-3)
            R0 = 3
        else:
            tmax = np.where(Vtot == np.max(Vtot))[0][0]*dt
            if (tmax <= 7):
                res_Ger = R0
                update = 1
            else:
                R0 += 0.1

plt.plot(frac_Ger/100,res_Ger,marker="o",color="gray",markersize=10)

################ add data points
######### Hongkong
t = np.arange(0,7.5+dt,dt)

R0 = 3
frac_Hongkong = 10**(-2)
update = 0
while (update == 0):
    print(R0)
    beta = c *delta*R0/((p*mu-delta*R0)*frac_Hongkong*T0/100)
    #
    ################# ODE Simulation
    ################# Initial condition
    V0 = 1
    E0 = 0
    I0 = 0
    NV0 = 0
    T0_upd = T0*frac_Hongkong/100
    #
    ################# Trajectories
    Tfull, Efull, Ifull, Vfull, NVfull = [], [], [], [], []
    #
    Tfull.append(T0_upd)
    Efull.append(E0)
    Ifull.append(I0)
    Vfull.append(V0)
    NVfull.append(NV0)
    #
    i = 0
    for i in range(len(t)-1):
        Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
        Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] * dt))
        Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta * Ifull[-1] * dt))
        Vfull.append(max(0,Vfull[-1] + mu * p * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
        NVfull.append(max(0,NVfull[-1] + (1-mu) * p * Ifull[-1] * dt - c * NVfull[-1] * dt))
    #
    Vwork = np.asarray(Vfull)
    NVwork = np.asarray(NVfull)
    Vtot = Vwork + NVwork
    #
    if (np.max(Vtot)==Vtot[-1] or np.max(Vtot)==Vtot[0]):
        R0 += 0.1
    else:
        tmax = np.where(Vtot == np.max(Vtot))[0][0]*dt
        if (np.max(Vtot) < 10**6):
            if (tmax > 4):
                R0 += 0.1
            else:
                frac_Hongkong += 10**(-3)
                R0 = 3
        else:
            if (tmax <= 4):
                res_Hongkong = R0
                update = 1
            else:
                R0 += 0.1

plt.plot(frac_Hongkong/100,res_Hongkong,marker="o",color="darkgray",markersize=10)

################ add data points
######### NBA
t = np.arange(0,7.5+dt,dt)

R0 = 3
frac_nba = 10**(-3)
update = 0
while (update == 0):
    print(R0)
    beta = c *delta*R0/((p*mu-delta*R0)*frac_nba*T0/100)
    #
    ################# ODE Simulation
    ################# Initial condition
    V0 = 1
    E0 = 0
    I0 = 0
    NV0 = 0
    T0_upd = T0*frac_nba/100
    #
    ################# Trajectories
    Tfull, Efull, Ifull, Vfull, NVfull = [], [], [], [], []
    #
    Tfull.append(T0_upd)
    Efull.append(E0)
    Ifull.append(I0)
    Vfull.append(V0)
    NVfull.append(NV0)
    #
    i = 0
    for i in range(len(t)-1):
        Tfull.append(max(0,Tfull[-1] - beta * Tfull[-1] * Vfull[-1] * dt))
        Efull.append(max(0,Efull[-1] + beta * Tfull[-1] * Vfull[-1] * dt - k * Efull[-1] * dt))
        Ifull.append(max(0,Ifull[-1] + k * Efull[-1] *dt - delta * Ifull[-1] * dt))
        Vfull.append(max(0,Vfull[-1] + mu * p * Ifull[-1] * dt - c * Vfull[-1] * dt - beta * Tfull[-1] * Vfull[-1] * dt))
        NVfull.append(max(0,NVfull[-1] + (1-mu) * p * Ifull[-1] * dt - c * NVfull[-1] * dt))
    #
    Vwork = np.asarray(Vfull)
    NVwork = np.asarray(NVfull)
    Vtot = Vwork + NVwork
    #
    if (np.max(Vtot)==Vtot[-1] or np.max(Vtot)==Vtot[0]):
        R0 += 0.1
    else:
        tmax = np.where(Vtot == np.max(Vtot))[0][0]*dt
        if (np.max(Vtot) < 4*10**5):
            if (tmax > 3):
                R0 += 0.1
            else:
                frac_nba += 10**(-4)
                R0 = 3
        else:
            if (tmax <= 3):
                res_nba = R0
                update = 1
            else:
                R0 += 0.1

plt.plot(frac_nba/100,res_nba,marker="o",color="olive",markersize=10)


########################## Finish plot
plt.ylim((0,20))
plt.tick_params(axis='both', which='major', labelsize=15, width=1, length=10) #top='True', right='True',
plt.tick_params(axis='both', which='minor', labelsize=10, width=1, length=5)
plt.show()


