import numpy as np
import matplotlib.pyplot as plt

def fdgaus(cutoff, dt, nt):
    w = np.zeros(nt)
    phi = 4 * np.arctan(1.0)
    """ a=phi*cutoff**2 """
    a = phi * (5.0 * cutoff / 8.0) ** 2
    amp = np.sqrt(a / phi)
    for i in range(nt):
        t = i * dt
        arg = -a * t ** 2
        if arg < -32.0:
            arg = -32.0
        w[i] = amp * np.exp(arg)

    for i in range(nt):
        if w[i] < 0.001 * w[0]:
            icut = i
            t0 = icut * dt
            break

    for i in range(nt):
        t = i * dt
        t = t - t0
        arg = -a * t ** 2
        if arg < -32.0:
            arg = -32.0
        w[i] = -2.0 * np.sqrt(a) * a * t * np.exp(arg) / np.sqrt(phi)

    smax = np.max(np.abs(w))
    
    for i in range(nt):
        w[i] = w[i] / smax

    plt.plot(w)
    plt.show()
    return w

#Parameter Initialize
v = 2
dx = 0.005
dz = 0.005

nx = 70
nz = 70

tmax = 1
G = 10  #dispersion 자료 참고
fmax = np.min(v)/(G*dx)
#fmax = 50
dt = 0.1/fmax
nt = int(tmax/dt)

df = 1/tmax
nf = int(fmax/df)

pi = np.arccos(-1)
alpha = np.log(100)/tmax
source = np.zeros(nt)
u = np.zeros((nx, nt), dtype=complex)
temp = np.zeros(nx*nz, dtype=complex)
cf = np.zeros(nx*nz, dtype=complex)
mat = np.zeros((nx*nz, nx*nz), dtype=complex)
green = np.zeros((nx, nf), dtype=complex)
print(f"famx:{fmax}, dt:{dt}, df:{df}")

source = fdgaus(fmax, dt, nt)
source = np.fft.fft(source)

cf[(nx//2)*nz+2] = 1.0

for ifreq in range(1,nf):
    print(ifreq)
    w = 2.0 * pi * (ifreq) * df - 1j * alpha
    for ix in range(nx):       
        for iz in range(nz):
            m = ix*nz+iz
            mat[m,m] = (-w**2/v**2)+(2/dx**2)+(2/dz**2)         
            try: mat[m,m-1] = -(1/dz**2)
            except IndexError: 
                continue
            try: mat[m,m+1] = -(1/dz**2)
            except IndexError:
                continue
            try: mat[m,m+nz] = -(1/dx**2)
            except IndexError:
                continue
            try: mat[m,m-nz] = -(1/dx**2)  
            except IndexError:
                continue    
    temp=np.linalg.solve(mat, cf)
    for ix in range(nx):
        green[ix, ifreq] = np.copy(temp[ix*nz+2])

for i in range(nx):
    for ifreq in range(nf):
        u[i, ifreq] = green[i, ifreq]*source[ifreq]
    for ifreq in range(1, nf):
        u[i, nt-ifreq] = np.conj(u[i, ifreq])

u = np.fft.ifft(u)/nt

for it in range(nt):
    u[:, it] = u[:,it]*np.exp(alpha*it*dt)

plt.xlabel('time')
plt.ylabel('x_dist')
plt.title("2D Wave Equation FDM Modeling in F-S")
plt.imshow(np.real(u), cmap='binary', aspect='auto', extent=[0, nx, 0, tmax])
plt.colorbar()
plt.show()