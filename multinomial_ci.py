# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
from scipy import stats
from scipy import special

def moments(c, _lambda):
    a = _lambda + c
    b = _lambda - c
    if b < 0:
        b = 0
    if b > 0:
        den = stats.poisson.cdf(a, _lambda) - stats.poisson.cdf(b-1, _lambda)
    if b == 0:
        den=stats.poisson.cdf(a, _lambda)
    mu = np.zeros(4)
    mom = np.zeros(5)
    for r in range(1,5):
        poisA = 0
        poisB = 0
        if (a-r) >= 0:
            poisA = stats.poisson.cdf(a,_lambda) - stats.poisson.cdf(a-r,_lambda)
        if (a-r) < 0:
            poisA = stats.poisson.cdf(a,_lambda)
        if (b-r-1) >= 0:
            poisB = stats.poisson.cdf(b-1,_lambda) - stats.poisson.cdf(b-r-1,_lambda)
        if (b-r-1) < 0 and (b-1) >= 0:
            poisB = stats.poisson.cdf(b-1,_lambda) 
        if (b-r-1) < 0 and (b-1) < 0: 
            poisB = 0
        mu[r-1]=(_lambda**r)*(1-(poisA-poisB)/den)
    mom[0]=mu[0]
    mom[1]=mu[1]+mu[0]-mu[0]**2
    mom[2]=mu[2]+mu[1]*(3-3*mu[0])+(mu[0]-3*mu[0]**2+2*mu[0]**3)
    mom[3]=mu[3]+mu[2]*(6-4*mu[0])+mu[1]*(7-12*mu[0]+6*mu[0]**2)+mu[0]-4*mu[0]**2+6*mu[0]**3-3*mu[0]**4
    mom[4]=den
    return mom
# end def moments(c, _lambda)

def truncpoi(c,x,n,k):
    m=np.zeros((k,5))       
    for i in range(k):      #
        _lambda=x[i]
        mom = moments( c, _lambda )
        for j in range(5):
            m[i,j] = mom[j]
        # end inner for
    # end outer for
    for i in range(k):
        m[i,3]=m[i,3]-3*m[i,1]**2
    # end for
    
    s = m.sum(axis=0)
    s1 = s[0]
    s2 = s[1]
    s3 = s[2]
    s4 = s[3]
    probn=1/(stats.poisson.cdf(n,n)-stats.poisson.cdf(n-1,n)) 
    z=(n-s1)/np.sqrt(s2) # ok
    g1=s3/(s2**(3/2))
    g2=s4/(s2**2)
    poly = 1 + g1*(z**3-3*z)/6 + g2*(z**4-6*z**2+3)/24 + g1**2*(z**6-15*z**4+45*z**2-15)/72
    f=poly*np.exp(-z**2/2)/(np.sqrt(2)*special.gamma(0.5))
    probx=1
    for i in range(k):
        probx = probx*m[i,4]
    return probn*probx*f/np.sqrt(s2)
# end def truncpoi(c,x,n,k)
    
def sison(x,alpha,verbose=False):
    n = int(sum(filter(None,x))) 
    k = len(x)             
    p = x/n
    c = 0
    pold = 0
    for cc in range(n):
        p = truncpoi(cc+1,x,n,k)
        if p > 1-alpha and pold < 1-alpha:
            c = cc + 1
            break
        pold = p
    # end for
    
    salida = np.zeros((k,2))
    delta = (1-alpha-pold)/(p-pold)
    #print(delta)
    out = np.zeros((k,5))
    num = np.zeros((k,1))
    c = c-1
    vol1 = 1
    vol2 = 1
    for i in range(k):
        num[i,0]=i 
        obsp=x[i]/n 
        out[i,0]=obsp
        out[i,1]=obsp-c/n 
        out[i,2]=obsp+c/n+2*delta/n 
        if out[i,1]<0:  
            out[i,1]=0 
        if out[i,2]>1: 
            out[i,2]=1 
        out[i,3]=obsp-c/n-1/n 
        out[i,4]=obsp+c/n+1/n 
        if out[i,1]<0: 
            out[i,1]=0 
        if out[i,2]>1: 
            out[i,2]=1 
        vol1=vol1*(out[i,2]-out[i,1]) 
        vol2=vol2*(out[i,4]-out[i,3]) 
        
        salida[i,0] = out[i,1]
        salida[i,1] = out[i,2]
    # end for
    return salida
# end def sison(x,alpha,verbose=False)