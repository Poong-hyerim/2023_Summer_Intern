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

    # plt.plot(w)
    # plt.show()
    return w

delta = 20
v = 2
dx = 0.01

d0 = np.log(1000)*((3*v)/(2*delta*dx))
print(d0)

pi = np.arccos(-1)
bnx = 400
nx = bnx+2*delta
fmax = 20
tmax = 2

dt = (0.5/fmax)
nt = int(tmax/dt)
df = 1/tmax
nf = int(fmax/df)

ndt = (0.1/fmax)
nnt = int(tmax/ndt)
alpha = np.log(100)/tmax

source = np.zeros(nt, dtype=complex)
mat = np.zeros((nx,nx), dtype=complex)
csource = np.zeros(nt, dtype=complex)
u = np.zeros((nx, nnt), dtype=complex)
f = np.zeros(nx, dtype=complex)
green = np.zeros((nx,nf), dtype=complex)

csource = fdgaus(fmax, dt, nt)
# for it in range(0, nt):
#     csource[it] = source[it]*np.exp(-alpha*it*dt) 
csource = np.fft.fft(csource)

for ifreq in range(1, nf):
    print(ifreq)
    w = 2.0 * pi * (ifreq) * df - 1j * alpha
    f[:] = 0.0
    f[nx//2] = 1.0
    mat[:,:] = 0.0
    #mat[0, 0] = -(w ** 2 / v ** 2) + (2 / dx ** 2)
    #mat[0, 1] = -1.0 / dx**2
    #mat[nx - 1, nx - 1] = -(w ** 2 / v ** 2) + (2 / dx ** 2)
    #mat[nx - 1, nx - 2] = -1.0 / dx**2
    for i in range(0,nx):
        #PML 영역일 경우 감쇠함수 d(x) sx, sxp를 설정해서 감쇠를 발생시킴!
        Sx = 1
        Sxp = 0
        if i<delta:
            Sx = (w*1j)/((w*1j)+(d0*((delta-i)/delta)**2))
            Sxp = ((w*1j)/((w*1j)+(d0*((delta-i)/delta)**2))**2)*(2*d0*((delta-i)/delta**2/dx))
        elif i>(bnx+delta-1):
            Sx = (w*1j)/((w*1j)+(d0*((i-bnx-delta+1)/delta)**2))
            Sxp = -((w*1j)/((w*1j)+(d0*((i-bnx-delta+1)/delta)**2))**2)*(2*d0*((i-bnx-delta+1)/delta**2/dx))
        #print(i,Sx)
        mat[i,i] = -(w**2/v**2)+(2*Sx**2/dx**2)
        if(i!=0):
            mat[i,i-1] = ((Sx*Sxp)/(2*dx))-((Sx**2)/dx**2)
        if(i!=nx-1):
            mat[i,i+1] = -((Sx*Sxp)/(2*dx))-((Sx**2)/dx**2)
    #print(f"shape of mat{mat.shape}, mat value : {mat}")
    f = np.linalg.solve(mat, f)
    green[:, ifreq] = f[:]   
    

for ix in range(nx):
    for ifreq in range(nf):
        u[ix, ifreq] = green[ix, ifreq] * csource[ifreq]
    for ifreq in range(1,nf):
        u[ix, (nnt-1)-ifreq+1] = np.conj(u[ix, ifreq])

u = np.fft.ifft(u) / nnt

for it in range (0, nnt):
    u[:,it]=u[:,it]*np.exp(alpha*it*ndt)
    
plt.ylabel('x_dist')
plt.xlabel('time')
plt.title("1D Wave Equation FDM Modeling in F-S applying PML Boundary Condition")
plt.imshow(np.real(u), cmap='binary', aspect='auto', extent=[0, nx, 0, tmax])
plt.colorbar()
plt.show()