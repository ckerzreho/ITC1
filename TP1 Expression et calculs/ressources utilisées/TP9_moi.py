# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 14:54:43 2017

@author: kerzrehocorentin
"""
import matplotlib.pyplot as plt
import os
import numpy as np
import math as m
from scipy.integrate import odeint

os.chdir("/Users/kerzrehocorentin/Documents/cours/_PCSI/LE RESTE/Informatique/PCSI3/TP9")

t,uc=np.loadtxt("TP9_acquisition2.txt",skiprows=1,delimiter=';',unpack=True)
t,uc=np.loadtxt("TP9_acquisition2.txt",skiprows=1,delimiter=';',unpack=True)

##

plt.plot(t,uc,'-x')
plt.grid(True)
plt.xlabel("Temps")
plt.ylabel("Uc")

def pseudoT(uc,t): # 3 méthodes toutes les valeurs de T, valeur moyenne T, comptage de nT
    max=determination_max(uc)
    Tmesur=(t[max[-1]]-t[max[0]])/(len(max)-1)
    return Tmesur
    
def determination_max(uc):
    max=[]
    for i in range(3,len(uc)-3):
        if uc[i]>uc[i-3] and uc[i]>uc[i+3]:
            max.append(i)
    return max


def decrement_log(uc):
    max=determination_max(uc)
    for i in max:
        plt.plot(t[i],uc[i],'o')
    max1=max[0]
    max4=max[3]
    return 1/3*m.log(uc[max1]/uc[max4])
#############################################
def Q_et_w0(uc,t):
    delta = decrement_log(uc) 
    T = pseudoT(uc,t)
    Q = m.sqrt((m.pi/delta)**2+1/4)
    w0= 2*m.pi/T / m.sqrt(1 - 1/(4*Q**2))
    return Q,w0
#############################################
def R_et_L(uc,t):
    C=100e-9 #○ valeur de la capacité en farad
    Q,w0=Q_et_w0(uc,t)
    L=1/(w0**2*C)
    R=1/Q*m.sqrt(L/C)
    return R,L
############################################
Q,omega0=Q_et_w0(uc,t)
uc0=9.8
uc0point=0
def phiX(X,t):
    uc=X[0]
    ucpoint=X[1]
    return ([ucpoint,-omega0**2*uc-omega0/Q*ucpoint])
#########################################################
X0=[uc0,uc0point]
#########################################################
sol=odeint(phiX,X0,t)
ucdiff=sol[:,0]
ucdiffpoint=sol[:,1]
#########################################################
plt.figure()
plt.plot(t,uc,'+-b',t,ucdiff,'x-g')
plt.xlabel('$t$ (s)')
plt.ylabel('$u_c(t)$ (V)')
plt.title('Tension aux bornes de C')
plt.legend(['exp','odeint'],loc=0)
plt.grid()
##########################################################
##### Trajectoire de phase ###############################
plt.figure()
plt.plot(ucdiff,ucdiffpoint,'-b')
plt.xlabel('$u_c$ (V)')
#plt.ylabel('$\frac{du_c}{dt}$ (V/s)')
plt.title('Trajectoire de phase')
plt.grid()
##########################################################

