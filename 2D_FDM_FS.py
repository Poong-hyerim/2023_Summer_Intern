import numpy as np
import matplotlib as plt

#Parameter Initialize
nx = 200
nz = 100
dx = 0.005
dz = 0.005
tmax = 1
G = 16
velocity = 1.5
fmax = np.min(velocity)/(G*dx)

dt = 0.1/fmax
nt = int(tmax/dt)

df = 1/tmax
nf = int(fmax/df)

source = np.zeros(nt)
csource = np.zeros(nt)
u = np.zeros(nt)
green = ((nx*nz, nf))