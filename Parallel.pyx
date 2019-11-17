from cython.parallel import prange
import numpy as np
cimport numpy as np

#@cython.boundscheck(False)
#@cython.wraparound(False)

cpdef double compareCurves(double low, double high, 
                            np.ndarray[double, ndim = 1] yfunc1,
                            np.ndarray[double, ndim = 1] yfunc2, 
                            np.ndarray[double, ndim = 1] xfunc1, 
                            np.ndarray[double, ndim = 1] xfunc2):
    
    #assume step2 is larger than step1
    cdef double Total = 0
    cdef double minX
    cdef double maxX
    cdef double Avg2
    cdef double Avg1
    cdef double n
    cdef int length1 = len(xfunc1)
    cdef int length2 = len(xfunc2)
    cdef int i = 0
    cdef int j = 0
    cdef double Low = low
    cdef double High = high

    for i in prange(length2, nogil = True): 
        if xfunc2[i] >= Low and xfunc2[i] <= High:
            minX = xfunc2[i]
            maxX = xfunc2[i+1]
            Avg2 = 0.5*(yfunc2[i]+yfunc2[i+1])
            Avg1 = 0
            n = 0
            for j in xrange(length1):
                if xfunc1[j] <= maxX and xfunc1[j] >= minX:
                    n = n + 1
                    Avg1 = Avg1 + yfunc1[j]
            if n>0:
                Avg1 = Avg1/n
            Total += (Avg2-Avg1)**2/((Avg2+Avg1)/2)
    return Total

cpdef np.ndarray[double, ndim = 1] fssRoutine(float mindifference,
                                                np.ndarray[double, ndim = 2]FSSLists,
                                                np.ndarray[double, ndim = 1]pcSpace, 
                                                np.ndarray[double, ndim = 1]nuSpace,
                                                np.ndarray[double, ndim = 1]zetaSpace, 
                                                np.ndarray[double, ndim = 1]rho): 
    """
    cdef int length = len(pcSpace)
    cdef int lengthNu = len(nuSpace)
    cdef int lengthZeta = len(zetaSpace)
    """
    cdef int lengthRho = len(rho)
    cdef float minPc, minZeta, minNu = 0.0
    cdef double minDifference = mindifference
    cdef float dif = 0.0
    cdef float loLim, highLim = 0.0
    cdef list smallxFunc, smallyFunc, bigxFunc, bigyFunc = [] 
    cdef np.ndarray[double, ndim = 1] exponents

    for pc in pcSpace:
        print(pc)
        for nu in nuSpace:
            for zeta in zetaSpace:
                #Below, the 4 is a stand in for the smallest linear lattice size, and the 64 is a stand in for the largest one.
                lowLim=(rho[0]-pc)/pc*5**(1/nu)
                highLim=(rho[49]-pc)/pc*5**(1/nu)
                smallxFunc=[]
                bigxFunc=[]
                smallyFunc=[]
                bigyFunc=[]
                for i in range(lengthRho):
                    smallxFunc.append((rho[i]-pc)/pc*5**(1/nu))
                    bigxFunc.append((rho[i]-pc)/pc*60**(1/nu))
                    smallyFunc.append(FSSLists[0][i]*5**(-zeta/nu))
                    bigyFunc.append(FSSLists[4][i]*60**(-zeta/nu))
                
                Dif = compareCurves(lowLim,highLim,
                                        np.array(smallyFunc),
                                        np.array(bigyFunc),
                                        np.array(smallxFunc),
                                        np.array(bigxFunc))
                if (Dif < minDifference):
                    minDifference=Dif
                    minPc=pc
                    minNu=nu
                    minZeta=zeta
                
                #print(minDifference, minPc, minNu, minZeta)
    #print(minPc, minNu, minZeta)
    exponents = np.array([minPc, minNu, minZeta], dtype = float)
    return exponents
                
