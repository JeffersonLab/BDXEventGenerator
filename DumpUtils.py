#!/usr/local/bin/python

import math
import sys
from scipy.integrate import quad
from scipy.integrate import dblquad
from scipy.integrate import quadrature
from scipy.integrate import romberg
from scipy.special import sph_jn
from scipy.special import gamma
import numpy


#The following two functions are,respectively:

## dNdEIntegrand is dn/dE(t), i.e. the number of electrons in the dump,
## differential wrt the electron energy, at depth t (in rad. units), normalized to a single impinging electron
## Parameters are: 
## 1)x:   the depth  in the dump
## 2)Ein: the energy in GeV 
## dNdE Integral

def dNdEIntegrand(x,Ein):
    par=[]
    par.append(0.)  #useless parameter
    par.append(787.)
    par.append(-108.6)
    par.append(1.564)
    par.append(2445)
    par.append(29920)
    par.append(0.1122)
    par.append(632.5)
    par.append(-583.1)
    par.append(7208)
    par.append(1.055)
    par.append(0.8282)
    par.append(1.795)
    par.append(-0.07442)
    par.append(0.0287)


    ret=0.
    termp=0.
    N1=0.
    E1=1000000.
    N2=0.
    E2=0.009
    E0=11000.
    Ecut=1000.
    Norm=1.
    E=Ein*1000.

    #N1
    N1p0 = par[1]
    N1p1 = par[2]
    N1p2 = par[3]
    u = x-0.1    
    f= u*(N1p2*u*u+N1p1*u+N1p0)*(1-math.exp(-u/0.4))*math.exp(-u/8);
    if(f>0. and u>0.):
        N1  =  f

    #E1
    E1p0 = par[4]
    if(x>0.01):
        E1 = E1p0/math.pow(x,3./4.);
        
    # N2 
    N2p0 = par[5];
    N2p1 = par[6];
    N2p2 = par[7];
    N2p3 = par[8];
    if(x < 0.5):
        N2  =  N2p0*math.exp(-x/N2p1);
    else:
        temp = N2p2 + N2p3*x;
        if(temp>0):
            N2 = temp
      
    #E2
    E2p0 = par[9];
    if(x>0.01):
        E2 = (E2p0*x*x)*(1-math.exp(-x/0.09))

    #normalization
    Normp0 = par[10]
    Normp1 = par[11]
    Normp2 = par[12]
    Normp3 = par[13]
    Normp4 = par[14]
        
    Norm = math.exp(-x*Normp0)*(1+Normp1*x+Normp2*x*x+Normp3*x*x*x+Normp4*x*x*x*x);

    deltaE = E0-Ecut;

    ret = Norm*(N1*math.exp(-(E-Ecut)/E1) + N2*math.exp((E-E0)/E2))/(E1*N1*(1 - math.exp(-deltaE/E1)) + E2*N2*(1 - math.exp(-deltaE/E2)));
    
    ret*=1000; #This is necessary, since the dN/dE was with E in MeV (there is a differential factor)
    return ret; 



def dNdEIntegrand_nopositrons(x,Ein):
    par=[]
    par.append(0.)     # Par0
    par.append(481.2)  # Par1
    par.append(-82.23) # Par2
    par.append(2.878)  # Par3
    par.append(2600.)  # Par4
    par.append(30410.) # Par5
    par.append(0.1102) # Par6
    par.append(602.2)  # Par7
    par.append(-552.3) # Par8
    par.append(8232.)  # Par9
    par.append(1.785)  # Par10
    par.append(1.573)  # Par11
    ret=0.
    termp=0.
    N1=0.
    E1=1000000.
    N2=0.
    E2=0.009
    E0=11000.
    Ecut=1000.
    Norm=1.
    E=Ein*1000.
    
    N1p0 = par[1]
    N1p1 = par[2]
    N1p2 = par[3]
    u = x-0.1
    f = u*(N1p2*u*u+N1p1*u+N1p0)*(1.-math.exp(-u/0.4))*math.exp(-u/8.)
    if f > 0. and u > 0.:
        N1=f
    E1p0 = par[4]
    if x > 0.01:
        E1 = E1p0/math.pow(x,3./4.)
    N2p0 = par[5]
    N2p1 = par[6]
    N2p2 = par[7]
    N2p3 = par[8]
    if x < 0.5:
        N2  =  N2p0*math.exp(-x/N2p1)
    else:
        temp = N2p2 + N2p3*x
        if temp > 0.:
            N2 = temp
    E2p0 = par[9]
    if x > 0.01:
        E2 = (E2p0*x*x)*(1.-math.exp(-x/0.09))
    Normp0 = par[10]
    Normp1 = par[11]
    u = x - Normp0
    if u > 0.:
        Norm = Norm*math.exp(-u/Normp1)
    deltaE = E0-Ecut
    ret = Norm*(N1*math.exp(-(E-Ecut)/E1) + N2*math.exp((E-E0)/E2))/(E1*N1*(1 - math.exp(-deltaE/E1)) + E2*N2*(1 - math.exp(-deltaE/E2)))
    
    ret*=1000; #This is necessary, since the dN/dE was with E in MeV (there is a differential factor)
    
    return ret



def dNdEIntegrandTsaiFormula(x,Ein0):
    E0=11000.0
    Ein=Ein0*1000.
    numer=math.log(E0/Ein)
    numer=math.pow(numer,4.0/3.0*x-1.0)
    denom=E0*gamma(4.0*x/3.0)
    integrand=numer/denom
    integrand*=1000;#This is necessary, since the dN/dE was with E in MeV (there is a differential factor)
    return integrand



## dNdEintegral is the dn/dE(t) integrated OVER t between tMIN and tMAX. 
## Parameters are:
## 1) Ein: the electron energy in MeV
## 2) cut: the integration range in radiation units!




def dNdEIntegral(Ein,cut=[0.03,8.0]):
    dNdE=quad(lambda x: dNdEIntegrand(x,Ein),cut[0],cut[1],epsabs=1.0e-04, epsrel=1.0e-04, limit=50000)[0]   
   # dNdE=quad(lambda x: dNdEIntegrandTsaiFormula(x,Ein),cut[0],cut[1],epsabs=1.0e-04, epsrel=1.0e-04, limit=50000)[0]
    return dNdE
