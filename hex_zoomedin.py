import sys
import numpy as np
import math
import time
import csv
import torch
from google.colab import files

Dtype = torch.float
gpu = torch.device("cuda:0")

top = [4, 10, 20, 40, 60]
middle = [5, 11, 21, 41, 61]
rows = [5, 10, 20, 40, 60]

rundist = [1000, 100000, 1000, 200, 3000]
start = time.time()

def numNodes(top, middle, rows):
    n = rows - 2
    return (2 * top) + (n * middle)

def link(p):
    C = np.random.random()
    if C > p:
        C = -1/1000 
    else:
        C = -1
    return C

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
def zerovec(L):
    L = np.delete(L, 0, 1)
    L = np.delete(L, 0, 0)
    return L

def bulk(L, nodes, top, middle, p):
    x = 0
    y = top - 1
    v = top 
    u = 2 * top
    z = top + middle
    w = (2 * top) + middle - 1
    a = top + (2 * middle) - 1
    b = top + (3 * middle)
    sep = top + (3 * middle) + 1
    for i in range(1, nodes-1):
        for j in range(1, nodes-1):
            C = link(p)
            if (x <= i <= y):
                if (j == i + middle or j == i + top) and L[i][j] == 0:
                    L[i][j] = C
                    L[j][i] = C
            
            if (z <= i <= w):
                if (j == i + middle + 1 or j == i + top + 1) and L[i][j] == 0:
                    L[i][j] = C
                    L[j][i] = C
            if (v <= i <= u):
                if( j == i + middle and L[i][j] == 0):
                    L[i][j] = C
                    L[j][i] = C
            if (a <= i <= b):
                if( j == i + middle and L[i][j] == 0):
                    L[i][j] = C
                    L[j][i] = C
        
        if (i == y):
            x += sep
            y += sep
        
        if ( i == w ):
            z += sep    
            w += sep
        if ( i == u ):
            v += sep
            u += sep
        if ( i == b):
            a += sep
            b += sep
        
    for i in range(nodes):
        sum = 0.
        for j in range(nodes):
            sum -= L[i][j]
        L[i][i] = sum

for n  in range(1,2):
    nodes = numNodes(top[n], middle[n], rows[n])
    avg_mob = []
    print(nodes)
    for p in range(200,400,4):
        print(p)
        r_avg = []
        runs = rundist[n]
        for i in range(runs):
            L = np.zeros(shape = (nodes, nodes))
            equipotential(L, nodes, top[n], p/1000)
            bulk(L, nodes, top[n], middle[n], p/1000)
            minor = zerovec(L)
            minor_gpu = torch.tensor(minor, device = gpu, dtype = Dtype)
            #find the eigenvalues and eigenvectors of this hermitian vector
            (w,v) = torch.symeig(minor_gpu, eigenvectors = True, upper = True)
            #this returns a tuple of Tensors that is stored on the GPU. 
            w = w.cpu().detach().numpy() #Transfer the data from the cpu to the gpu
            v = v.cpu().detach().numpy()

            R = 0
            for j in range(nodes - 1):
                R += (1/w[j])*np.abs(v[nodes -2][j])**2
            r_avg.append(1/R)
            #rint(n, i)
        avg_mob.append(np.mean(r_avg))
        
    
    with open("HexZoomedIn1.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect = "excel")
        wr.writerow(avg_mob)
 
print( time.time() - start)          
files.download('HexZoomedIn1.csv')
   