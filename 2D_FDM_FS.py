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

pi = np.pi
alpha = np.log(100)/tmax

source = np.zeros(nt)
u = np.zeros((nx*nz, nt), dtype=complex)
temp = np.zeros(nx*nz, dtype=complex)
cf = np.zeros(nx*nz, dtype=complex)
mat = np.zeros((nx*nz, nx*nz), dtype=complex)
green = np.zeros((nx*nz, nf), dtype=complex)

source = fdgaus(fmax, dt, nt)
source = np.fft.fft(source)

# 소스 주기!
cf[:] = 0.0
cf[(100)*nz+2] = 1.0

# k를 기준으로 넘버링을..해줘야겠죠?
for ifreq in range(1,nf):
    print(ifreq)
    w = 2.0 * pi * (ifreq) * df - 1j * alpha
    # 경계부 조건 처리
    mat[:,:] = 0.0
    # index align 하면서 mat 구성
    for ix in range(nx):
        for iz in range(nz):
            m = ix*nz+iz
            mat[m, m] = (-w**2/velocity**2)+(2/dx**2)+(2/dz**2)
            try: mat[m, m-1] = -(2/dz**2)
            except IndexError: 
                continue
            try: mat[m, m+1] = -(2/dz**2)
            except IndexError:
                continue
            try: mat[m, m+nz] = -(2/dx**2)
            except IndexError:
                continue
            try: mat[m, m-nz] = -(2/dx**2)  
            except IndexError:
                continue
    temp=np.linalg.solve(mat, cf)
    print(f"matrix=${temp}\n matrix_shape=${temp.shape}")
    green[:, ifreq] = np.copy(temp)
    
for ix in range(nx):
    for iz in range(nz):
        m = ix*nz+iz
        for ifreq in range(nf):
            u[m, ifreq] = green[m, ifreq]*source[ifreq]
        for ifreq in range(1, nf):
            u[m, nt-ifreq] = np.conj(u[m, ifreq])


u = np.fft.ifft(u)/nt
for it in range(nt):
    u[:, it] = u[:it]*np.exp(alpha*it*dt)
    
plt.xlabel('x_dist')
plt.ylabel('time')
plt.title("1D Wave Equation FDM Modeling in F-S")
plt.imshow(np.real(u), cmap='binary', aspect='auto', extent=[0, nx, 0, tmax])
plt.colorbar()  # Optionally, add a color bar for reference 
plt.show()