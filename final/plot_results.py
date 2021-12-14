# ASTE-404 Final Project
# Plot simulation results for RDE simulation

# Headers needed
import numpy as np
import matplotlib.pyplot as plt
import csv
import math

# Extract time history from csv file
u = np.genfromtxt("time_history.csv",delimiter=",")

# User defined values
snap_plot = 500         # plot results from this snapshot

# Get domain and number of snapshots
num_snap = u.shape[0]   # number of snapshots
ni = u.shape[1]         # domain length

if snap_plot >= num_snap:
    # Automatically redefine which snapshot is plotted if exceeds max
    snap_plot = num_snap - 1

# Locate which index corresponds to domain angle of 2pi or pi 
pi2val = int((ni-1))        # 2pi
pival = int((ni-1)/2)       # pi

# angles from 0 to 2*pi
angles = np.linspace(0,2*math.pi,num=ni)

# plot results as contour plot

uplot = np.matrix.transpose(u)  # transpose so time is on x axis

fig,ax = plt.subplots(figsize=(8,6))

plt.contourf(uplot,levels=255,cmap='jet')

ax.set_xticks([600,1200,1800,2400,3000]) # convert timestep to time
ax.set_xticklabels([1,2,3,4,5])          # convert timestep to time

ax.set_yticks([0,pival,pi2val])          # plot over angles
ax.set_yticklabels(["0","pi","2pi"])     # angles from 0 to 2*pi

ax.set_xlabel("Nondimensional time [-]")
ax.set_ylabel("theta [-]")
#plt.clim(0.1,100)
plt.colorbar()
#plt.contour(T,levels=24,colors='black')
plt.show()

# Plot a single wave over the annulus
fig,ax = plt.subplots(figsize=(8,6))

# Plot energy at selected snapshot
plot1 = ax.plot(angles,u[snap_plot,:], color="blue", label="Internal Energy")
ax.set_xlabel("theta [-]")
ax.set_ylabel("u [-]")

plt.show()
