# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 10:27:23 2022

@author: edf
"""
import numpy as np
import matplotlib.pyplot as plt

def riemannLeft(n,a,b,f,fArgs=[]):
    dx = np.float64(b-a)/n
    sm = 0.
    for i in range(n):
        x = a + np.float64(i)*dx
        sm = sm + f(x,*fArgs)
    return sm*dx

def riemannRight(n,a,b,f,fArgs=[]):
    dx = np.float64(b-a)/n
    sm = 0.
    for i in range(n):
        x = a + np.float64(i+1)*dx
        sm = sm + f(x,*fArgs)
    return sm*dx

def riemannMiddle(n,a,b,f,fArgs=[]):
    dx = np.float64(b-a)/n
    sm = 0.
    for i in range(n):
        x = a + np.float64(i+0.5)*dx
        sm = sm + f(x,*fArgs)
    return sm*dx

def Trapezoidal(n,a,b,f,fArgs=[]):
    dx = np.float64(b-a)/n
    sm = 0
    for i in range(1,n):
        x = a + np.float64(i)*dx
        sm += f(x,*fArgs)
    return 0.5*dx*(f(a,*fArgs) +2*sm + f(b,*fArgs))

def simpsons(n,a,b,f,fArgs=[]):
    if n%2:
        print('Warning: adjusting n to an even number')
        n += 1
    dx = np.float64(b-a)/n
    smEven = 0.
    for i in range(1,n//2):
        x = a + np.float64(2*i)*dx
        smEven = smEven + f(x,*fArgs)
    smOdd = 0.
    for i in range(n//2):
        x = a + np.float64(2*i+1)*dx
        smOdd = smOdd + f(x,*fArgs)
    return (f(a,*fArgs)+2*smEven+
            4*smOdd+f(b,*fArgs))*dx/3.

def integrand(x):
    return np.cos(x)

def power(x,m):
    return x**m

a = 1
b = 20
n = 10

exact = np.sin(b)-np.sin(a)

integralLeft = riemannLeft(n,a,b,integrand)
print('Riemann left',n,integralLeft,np.abs(integralLeft-exact))

integralRight = riemannRight(n,a,b,integrand)
print('Riemann right',n,integralRight,np.abs(integralRight-exact))

integralMiddle = riemannMiddle(n,a,b,integrand)
print('Riemann middle',n,integralMiddle,np.abs(integralMiddle-exact))

integralTrapezoidal = Trapezoidal(n,a,b,integrand)
print('Trapezoidal',n,integralTrapezoidal,np.abs(integralTrapezoidal-exact))

integralSimpsons = simpsons(n,a,b,integrand)
print('Simpsons',n,integralSimpsons,np.abs(integralSimpsons-exact))

nvalues = 2**np.arange(1,8)
residualLeft = []
residualMiddle = []
residualTrapezoidal = []
residualSimpsons = []

for n in nvalues:
    integralLeft = riemannLeft(n,a,b,integrand)
    residualLeft.append(np.abs(integralLeft-exact))
    integralMiddle = riemannMiddle(n,a,b,integrand)
    residualMiddle.append(np.abs(integralMiddle-exact))
    integralTrapezoidal = Trapezoidal(n,a,b,integrand)
    residualTrapezoidal.append(np.abs(integralTrapezoidal-exact))
    integralSimpsons = simpsons(n,a,b,integrand)
    residualSimpsons.append(np.abs(integralSimpsons-exact))

fig, ax = plt.subplots(1,1)
ax.set_yscale('log')
ax.set_xscale('log')
ax.plot(nvalues,residualLeft,'o-',label='Left')
ax.plot(nvalues,residualMiddle,'.-',label='Middle')
ax.plot(nvalues,residualTrapezoidal,'x-',label='Trapezoidal')
ax.plot(nvalues,residualSimpsons,'v-',label='Simpsons')
ax.legend()
plt.show()
