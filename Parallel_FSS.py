#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 12:05:28 2019

@author: chase
"""
import time
import numpy as np
import matplotlib.pyplot as plt
import Parallel

plt.style.use('seaborn-deep')
start = time.time()

FSSLists = []  
for l in (0,1,2,3,4):
    Data = open('TriZoomedIn' + str(l) + '.csv', 'r')
    RawData=Data.readlines()
    
    SplitData = []
    Data = []
    
    for i in range(50):
        SplitData = (float(RawData[i]))
        Data.append(float(SplitData))
    FSSLists.append(Data) 
#FSSLists = np.array(FSSLists) 

#adjust these to be appropriate for your lattices (pc range should be different, but I suspect nu and zeta will not)
pcSpace = np.linspace(.32,.325, 50)
nuSpace = np.linspace(1.2, 1.4, 200)
#nuSpace = [4/3]
zetaSpace = np.linspace(-1.4, -1.2, 200)

minDifference = 100000000
rho = np.linspace(.25, .45, 50)

#FSSLists are your data it is a list of lists, FSSErrors are the apropriate errors. If you did not include errors, don't worry
#just do a regular plot instead of an errorbar ploti

results = Parallel.fssRoutine(minDifference, np.array(FSSLists), pcSpace, nuSpace, zetaSpace, rho)
pc = results[0]
nu = results[1]
zeta = results[2]
print(pc)
print(nu)
print(zeta)
#print(time.time() - start)



x5=[]
y5=[]
x40=[]
y40=[]
x10=[]
y10=[]
x20=[]
y20=[]
y60=[]
x60=[]

for i in range(len(rho)):
    x5.append((rho[i]-pc)/pc*5**(1/nu))
    x10.append((rho[i]-pc)/pc*10**(1/nu))
    x20.append((rho[i]-pc)/pc*20**(1/nu))
    x40.append((rho[i]-pc)/pc*40**(1/nu))
    x60.append((rho[i]-pc)/pc*60**(1/nu))
    y5.append(FSSLists[0][i]*5**(-zeta/nu))
    y10.append(FSSLists[1][i]*10**(-zeta/nu))
    y20.append(FSSLists[2][i]*20**(-zeta/nu))
    y40.append(FSSLists[3][i]*40**(-zeta/nu))
    y60.append(FSSLists[4][i]*60**(-zeta/nu))


plt.plot(x5,y5, marker='.',linestyle='None',label='5')
plt.plot(x10,y10, marker='.',linestyle='None',label='10')
plt.plot(x20,y20, marker='.',linestyle='None',label='20')
plt.plot(x40,y40, marker='.',linestyle='None',label='40')
plt.plot(x60,y60, marker='.',linestyle='None',label='60')
plt.xlabel("$(p-p_c)L^{1/v}/p_c$")
plt.ylabel("\u03C3 $L^{-K/v}$")
plt.legend()
plt.title("$v$ = %5.3f, $K$ = %5.3f, $P_c$ = %5.3f" % (nu, zeta, pc))
#plt.show()
print(pc)
print(nu)
print(zeta)
