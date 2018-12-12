# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 18:08:13 2018

@author: Ming Lun Wu
"""
import matplotlib.pyplot as plt
import pandas as pd


Data = pd.read_csv('C:/Users/Ming Lun Wu/Desktop/Table.csv')
T2LNum = [3.4E6]
T1Num = [2.34E6]
T2NLNum = [2.76E6]
for i in range(7):
    T2LNum.append(T2LNum[i]/2)
    T1Num.append(T1Num[i]/2)
    T2NLNum.append(T2NLNum[i]/2)


"""
T2Loaded1 = [18811, 16802, 15075, 13956,10279, 6840, 4811,3503]
T2Unloaded1 = [7498,8468,7517,6689,6806,4334,3157,2180]
T11 = [16986, 22670, 21303, 13479, 11298, 6891, 5066, 3564]

T12 = [6639, 10305,11048,11439,4657,3726,2010,909]
T2Unloaded2 = [14272, 39199,27595,13412,8856,3681,2070,778]
T2Loaded2 = [4634, 9353, 5168, 5142, 6402, 2695, 1679, 904]
"""

T1 = [1515,9830,20369,23768,43412,38522,140295,182690]
T2L = [535,6287,6087,10027,13681,20727,28503,26255]
T2NL = [799,3801,9541,24856,14415,27322,61980,67794]

plt.plot(T2LNum, T2L, '--mo', label = 'T2 Loaded')
plt.plot(T2NLNum, T2NL, '--ro',label = 'T2 Unloaded')
plt.plot(T1Num, T1, '--bo', label ='T1')

plt.xlabel('Cell Number')
plt.ylabel('Median Fluorescence Value')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.title('Gathered without Auto Sampler')
plt.legend()
plt.xscale('log')
plt.show()
plt.savefig('T2 Titration Result.png')
#plt.cla()
"""
plt.plot(CellNum1, T2Loaded2, '--mo', label = 'T2 Loaded')
plt.plot(CellNum1, T2Unloaded2, '--ro',label = 'T2 Unloaded')
plt.plot(CellNum1, T12, '--bo', label ='T1')

plt.xlabel('Cell Number')
plt.ylabel('Median Fluorescence Value')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.title('Gathered with Auto Sampler')
plt.legend()
plt.xscale('log')
#plt.show()
plt.savefig('Gathered With Auto Sampler.png')
"""
