# ASTE-404 Homework 8 Part 2
# Colin Coane
# 1D Unsteady Heat Equation Solver, using FTCS scheme

import numpy as np
import matplotlib.pyplot as plt
from random import random as rand
import math

"""
Given parameters:

Domain Length: 1
Number of spatial cells: 60
Number of time steps: 15,000
Time Step: 10^-5
Boundary Conditions: T(t,0) = T(t,1) = 1
Initial Conditions: T(0,x) = 0

"""

ni = 61     # number of nodes
nk = 35001  # number of time steps plus initial time t = 0

dx = 1/(ni - 1) # node spacing
dt = 1e-3       # time step

# Set coefficients
L = 2*math.pi
q = 1.0
alpha = 0.3
uc = 1.1
s = 3.5
k0 = 1.5
eps = 0.15
r = 5
nu1 = 0.0075
nu2 = 0.0075
up = 0.5

# Create array for u and lambda
u = np.zeros((nk,ni))   # np.zeros() sets initial condition to zero
lam = np.zeros((nk,ni)) # np.zeros() sets initial condition to zero

# Initialize pulse characteristics before wave
for i in range(ni):
    u[0,i] = (3/2)*((1/np.cosh(8*(i*dx-1)))**2)
    lam[0,i] = 0.5

#for i in range(ni):
    #print(u[0,i])
    #print(lam[0,i])

# loop over k - march forward in time
for k in range(nk-1):
    #print(k)
    # loop over spatial range
    for i in range(1,ni-1):
        # FTCS algorithm
        # Solve u
        u[k+1,i] = u[k,i] + dt*nu1/(dx*dx)*(u[k,i-1] - 2*u[k,i] + u[k,i+1]) 
        u[k+1,i] += -dt*u[k,i]*(u[k,i+1] - u[k,i-1])/(2*dx) 
        u[k+1,i] += dt*k0*q*(1-lam[k,i])*np.exp((u[k,i]-uc)/alpha)
        u[k+1,i] += -dt*eps*u[k,i]*u[k,i]
        # Solve lambda
        lam[k+1,i] = lam[k,i] + dt*nu2/(dx*dx)*(lam[k,i-1] - 2*lam[k,i] + lam[k,i+1]) 
        lam[k+1,i] += dt*k0*(1-lam[k,i])*np.exp((u[k,i]-uc)/alpha) 
        lam[k+1,i] += -dt*s*up*lam[k,i]/(1 + np.exp(r*(u[k,i] - up)))

    # Apply periodic boundary conditions

    # At i = 0
    # solve u
    u[k+1,0] = u[k,0] + dt*nu1/(dx*dx)*(u[k,ni-1] - 2*u[k,0] + u[k,1]) 
    u[k+1,0] += -dt*u[k,0]*(u[k,1] - u[k,ni-1])/(2*dx) 
    u[k+1,0] += dt*k0*q*(1-lam[k,0])*np.exp((u[k,0]-uc)/alpha) 
    u[k+1,0] += -dt*eps*u[k,0]*u[k,0]
    # Solve lambda
    lam[k+1,0] = lam[k,0] + dt*nu2/(dx*dx)*(lam[k,ni-1] - 2*lam[k,0] + lam[k,1]) 
    lam[k+1,0] += dt*k0*(1-lam[k,0])*np.exp((u[k,0]-uc)/alpha) 
    lam[k+1,0] += -dt*s*up*lam[k,0]/(1 + np.exp(r*(u[k,0] - up)))

    # At i = ni-1
    # Solve u
    u[k+1,ni-1] = u[k,ni-1] + dt*nu1/(dx*dx)*(u[k,ni-2] - 2*u[k,ni-1] + u[k,0]) 
    u[k+1,ni-1] += -dt*u[k,ni-1]*(u[k,0] - u[k,ni-2])/(2*dx) 
    u[k+1,ni-1] += dt*k0*q*(1-lam[k,ni-1])*np.exp((u[k,ni-1]-uc)/alpha) 
    u[k+1,ni-1] += -dt*eps*u[k,ni-1]*u[k,ni-1]
    # Solve lambda
    lam[k+1,ni-1] = lam[k,ni-1] + dt*nu2/(dx*dx)*(lam[k,ni-2] - 2*lam[k,ni-1] + lam[k,0]) 
    lam[k+1,ni-1] += dt*k0*(1-lam[k,ni-1])*np.exp((u[k,ni-1]-uc)/alpha) 
    lam[k+1,ni-1] += -dt*s*up*lam[k,ni-1]/(1 + np.exp(r*(u[k,ni-1] - up)))



# plot results as contour plot
#lev = np.linspace(0,2,10)
uplot = np.matrix.transpose(u)
lplot = np.matrix.transpose(lam)

plt.contourf(uplot,levels=255,cmap='jet')
#plt.clim(0.1,100)
plt.colorbar()
#plt.contour(T,levels=24,colors='black')
plt.show()
