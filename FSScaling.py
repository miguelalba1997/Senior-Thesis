#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 12:05:28 2019

@author: chase
"""

import numpy as np
import matplotlib.pyplot as plt


def compareCurves(low,high,yfunc1,yfunc2,xfunc1,xfunc2):
    #assume step2 is larger than step1
    Total=0
    for i in range(len(xfunc2)):

        if xfunc2[i] >=low and xfunc2[i] <= high:
            minX=xfunc2[i]
            maxX=xfunc2[i+1]
            Avg2=.5*(yfunc2[i]+yfunc2[i+1])
            Avg1=0
            n=0
            for j in range(len(xfunc1)):
                if xfunc1[j] <= maxX and xfunc1[j] >= minX:
                    n+=1
                    Avg1+=yfunc1[j]
            if n>0:
                Avg1=Avg1/n
            Total+=np.abs(Avg2-Avg1)/((Avg2+Avg1)/2)
    return Total



FSSLists=[]
FSSErrors=[]
for l in (4,16,32,64):

    Dat=open('FiniteSizeScaling' + str(l) + '.txt', 'r')
    RawData=Dat.readlines()
        
     
    
    SplitData=[];
    Data=[];
    Errors=[]
    for i in range(len(RawData)):
        SplitData=(RawData[i].split())
        Data.append(float(SplitData[1]))
        Errors.append(float(SplitData[2]))
    
    minErr=100
    for er in Errors:
        if er>0 and er<minErr:
            minErr=er
    for i in range(len(Errors)):
        if Errors[i]==0:
            #print("ERROR 0 DETECTED")
            Errors[i]=minErr/100
           # print(Errors[i])
    FSSLists.append(Data)
    FSSErrors.append(Errors)


#adjust these to be appropriate for your lattices (pc range should be different, but I suspect nu and zeta will not)
pcSpace=np.linspace(.48,.52,100)
nuSpace=np.linspace(1,3,100)
zetaSpace=np.linspace(-3,-1,100)

minDifference=100000000
rho=np.linspace(.01,1,100)

#FSSLists are your data it is a list of lists, FSSErrors are the apropriate errors. If you did not include errors, don't worry
#just do a regular plot instead of an errorbar plot
for pc in pcSpace:
    print(pc)
    for nu in nuSpace:
        for zeta in zetaSpace:
            #Below, the 4 is a stand in for the smallest linear lattice size, and the 64 is a stand in for the largest one.
            lowLim=(rho[0]-pc)/pc*4**(1/nu)
            highLim=(rho[99]-pc)/pc*4**(1/nu)
            smallxFunc=[]
            bigxFunc=[]
            smallyFunc=[]
            bigyFunc=[]
            for i in range(len(rho)):
                smallxFunc.append((rho[i]-pc)/pc*4**(1/nu))
                bigxFunc.append((rho[i]-pc)/pc*64**(1/nu))
                smallyFunc.append(FSSLists[0][i]*4**(-zeta/nu))
                bigyFunc.append(FSSLists[1][i]*64**(-zeta/nu))
            
           
            Dif=compareCurves(lowLim,highLim,smallyFunc,bigyFunc,smallxFunc,bigxFunc,pc,nu,4)
            if Dif<minDifference:
                minDifference=Dif
                minPc=pc
                minNu=nu
                minZeta=zeta
                
print(minPc)
print(minNu)
print(minZeta)
        

pc=minPc
nu=minNu
zeta=minZeta




x4=[]
y4=[]
x64=[]
y64=[]
x16=[]
y16=[]
x32=[]
y32=[]

for i in range(len(rho)):
    x4.append((rho[i]-pc)/pc*4**(1/nu))
    x16.append((rho[i]-pc)/pc*16**(1/nu))
    x32.append((rho[i]-pc)/pc*32**(1/nu))
    x64.append((rho[i]-pc)/pc*64**(1/nu))
    y4.append(FSSLists[0][i]*4**(-zeta/nu))
    y16.append(FSSLists[1][i]*16**(-zeta/nu))
    y32.append(FSSLists[2][i]*32**(-zeta/nu))
    y64.append(FSSLists[3][i]*64**(-zeta/nu))

plt.figure(figsize=(10,7))
plt.title('nu = ' +str(nu) + ' zeta = ' + str(zeta) + ' pc = ' + str(pc))
plt.errorbar(x4,y4,yerr=FSSErrors[0],marker='.',linestyle='None',label='4')
plt.errorbar(x16,y16,yerr=FSSErrors[1],marker='.',linestyle='None',label='16')
plt.errorbar(x32,y32,yerr=FSSErrors[2],marker='.',linestyle='None',label='32')
plt.errorbar(x64,y64,yerr=FSSErrors[3],marker='.',linestyle='None',label='64')
plt.legend(loc=2)
plt.xlim(-5,5)
plt.ylim(0,3)
plt.show()
print(pc)
print(nu)
print(zeta)