import numpy as np
import math
import time
import csv
import torch

Dtype = torch.float
gpu = torch.device("cuda:0")

""" We are intrested in how non square lattice geometries will 
    behave in a resistor network. The idea is to identify how 
    the percolation threshold is effected. In addition, we will
    be using laplacian matrices to map the graphs of the 
    resistor networks to a matrix that we can analyze more easily.
"""
start_time = time.time()

#top = 25 middle = 24 and 40 rows gives about 1000 nodes
top = [5, 11, 21, 41, 61] 
middle = [4, 10, 20, 40, 60]
row = [5, 11, 21, 41, 61]
rundist = [100, 1000, 10, 200, 100]

#generate the number of nodes for each run
def NumNodes(top, middle, row):
    n = math.floor(row/2)
    nodes = int((top * (n+1)) + (middle * n))
    nodes = nodes + 2
    return nodes

#make the first and last row of nodes equipotential
def equipotential(L, nodes, top, p):
    C = link(p)
    for i in range(nodes):
        for j in range(nodes):
            if (i == 0 and ( 1 <= j <= top)):
                    L[i][j] = C
                    L[j][i] = C
            if (i == nodes - 1 and ( nodes - top - 1 <= j <= nodes - 2)):
                    L[i][j] = C
                    L[j][i] = C
 
#this function randomly places down resistors
def link(p):
    C = np.random.random()
    if C > p:
        C = -1/1000 
    else:
        C = -1
    return C

#this function picks out a node and makes that connection zero so that we can get a zero eigenvector.
def zerovec(L):    
    L = np.delete(L, 0, 1)
    L = np.delete(L, 0, 0)    
    return L

def bulk(L, nodes, p, top, middle):
    x = 2
    sep = top + middle
    y = top - 1
    w = top + 1
    z = sep - 1
    #fill out the adjacency matrix
    for i in range(1, nodes - 1):
        for j in range(1, nodes - 1):
            C = link(p)
            
            if (i-1)%sep == 0:#this one is working
                if (j == i + 1 or j == i + top) and L[i][j] == 0:
                    L[i][j] = C
                    L[j][i] = C
            
            if (x <= i <= y) or (w <= i <= z):#this one is working
                if ( j == i + 1 or j == i + top or j == i + middle) and L[i][j] == 0:
                    L[i][j] = C
                    L[j][i] = C
            
            if i%sep == 0:#this one is working
                if (j == i + middle or j == i + top) and L[i][j] == 0:
                    L[i][j] = C
                    L[j][i] = C
            
            if (i + top - 1)%sep == 0:#this one is working
                if (j == i + middle) and L[i][j] == 0:
                    L[i][j] = C
                    L[j][i] = C
            
        if (i == y):
            x += sep
            y += sep
        if (i == z):
            w += sep
            z += sep                  
               
    #fill out the degree matrix
    for i in range(nodes):
        sum = 0.
        for j in range(nodes):
            sum -= L[i][j]
        L[i][i] = sum

for n in range(4,5):    
    nodes = NumNodes(top[n], middle[n], row[n]) 
    #iterate through all the percentages of conducting nodes versus total nodes to find the percolation threshold.
    avg_mob = []
    for p in range(100):
            print(p)
            r_avg = []
            runs = rundist[n]
            for i in range(runs):
                L = np.zeros(shape = (nodes, nodes))
                #create the equipotential
                equipotential(L, nodes, top[n], p/100)
                #create the bulk of the resistor ntwork matrix
                bulk(L, nodes, p/100, top[n], middle[n])
                #create the minor matrix
                minor = zerovec(L)
                minor_gpu = torch.tensor(minor, device = gpu, dtype = Dtype)
               
		for r in range(3):
                    try:
                        #find the eigenvalues and eigenvectors of this hermitian vector
                	(w,v) = torch.symeig(minor_gpu, eigenvectors = True, upper = True)
			#this returns a tuple of Tensors that is stored on the GPU. 
                	w = w.cpu().detach().numpy() #Transfer the data from the cpu to the gpu
                	v = v.cpu().detach().numpy()
                	
                    except:
			L = np.zeros(shape = (nodes, nodes))
                	#create the equipotential
                	equipotential(L, nodes, top[n], p/100)
                	#create the bulk of the resistor ntwork matrix
                	bulk(L, nodes, p/100, top[n], middle[n])
                	#create the minor matrix
                	minor = zerovec(L)
                	minor_gpu = torch.tensor(minor, device = gpu, dtype = Dtype)
			#find the eigenvalues and eigenvectors of this hermitian vector
                	(w,v) = torch.symeig(minor_gpu, eigenvectors = True, upper = True)
			#this returns a tuple of Tensors that is stored on the GPU. 
                	w = w.cpu().detach().numpy() #Transfer the data from the cpu to the gpu
                	v = v.cpu().detach().numpy()
                        
                R=0
                for j in range(nodes - 1):
                    R+=(1/w[j])*np.abs(v[nodes - 2][j])**2
                r_avg.append(1/R)
            avg_mob.append(np.mean(r_avg)) 
    
    with open("mobility1.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(avg_mob)
print(time.time() - start_time, 's')
