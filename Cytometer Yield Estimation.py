# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:31:02 2018

@author: Ming Lun Wu

Project:  Flow cytometer model for calculate volume and yield 
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# the equation is
occrT=4.8E-4
occrN = 2.94E-5
Purity = 0.11
#LamdaT = sep/occrT
#LamdaN = sep/occrN

a = 0.5
R = 1.6

#func = lambda tau : R - ((1.0 - np.exp(-tau))/(1.0 - np.exp(-a*tau))) 
func = lambda tau : Purity - 1/(1 + tau*np.exp(-tau/occrT)/occrN)

#func = lambda sep : Purity - 1.0/((sep/occrN)*np.exp(-1.0*(sep/occrT))

tau = np.linspace(1E-1,5E-10,2000)
y = func(tau)
max(y)
plt.plot(tau, func(tau))

#print(func)
sep_solution = fsolve(func,0.0000006)
k = sep_solution
#sep = np.linspace(-5,1.5)
#plt.plot(sep,func(sep))
#print(str(Purity))

nLamN = k/(1.06E-4)
nLamT = k/(2.35E-3)
Newpur = 1/(1+nLamN*np.exp(-nLamT))