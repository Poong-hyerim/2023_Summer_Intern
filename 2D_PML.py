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
delta = 20

dx = 0.005
dx0 = np.log(1000)*((3*v)/(delta*dx))
dz = 0.005
dz0 = np.log(1000)*((3*v)/(delta*dz))

nx = 70
nz = 70
nnx = nx+2*delta
nnz = nz+delta

tmax = 1
G = 10  #dispersion 자료 참고
fmax = np.min(v)/(G*dx)
#fmax = 25
dt = 0.05/fmax
print(dt/dx)
nt = int(tmax/dt)

df = 1/tmax
nf = int(fmax/df)

pi = np.arccos(-1)
alpha = np.log(100)/tmax
source = np.zeros(nt)
u = np.zeros((nnx, nt), dtype=complex)
uz = np.zeros((nnz, nt), dtype=complex)
#temp = np.zeros(nnx*nnz, dtype=complex)
cf = np.zeros(nnx*nnz, dtype=complex)
mat = np.zeros((nnx*nnz, nnx*nnz), dtype=complex)
green = np.zeros((nnx, nf), dtype=complex)
#greenz = np.zeros((nnz, nf), dtype=complex)

source = fdgaus(fmax, dt, nt)
source = np.fft.fft(source)

cf[(nnx//2)*nnz+2] = 1.0
print(f"famx:{fmax}, dt:{dt}, df:{df}")

for ifreq in range(1,nf):
    print(ifreq)
    w = 2.0 * pi * (ifreq) * df - 1j * alpha
    for ix in range(nnx):       
        for iz in range(nnz):
            Sx = 1
            Sxp = 0
            Sz = 1
            Szp = 0
            if ix<delta:
                Sx = (w*1j)/((w*1j)+(dx0*((delta-ix)/delta)**2))
                Sxp = ((w*1j)/((w*1j)+(dx0*((delta-ix)/delta)**2))**2)*(2*dx0*((delta-ix)/delta**2/dx))
            elif ix>(nx+delta-1):
                Sx = (w*1j)/((w*1j)+(dx0*((ix-nx-delta+1)/delta)**2))
                Sxp = -((w*1j)/((w*1j)+(dx0*((ix-nx-delta+1)/delta)**2))**2)*(2*dx0*((ix-nx-delta+1)/delta**2/dx))
            elif iz>(nz-1):    
                Sz = (w*1j)/((w*1j)+(dz0*((iz-nz+1)/delta)**2))
                Szp = -((w*1j)/((w*1j)+(dz0*((iz-nz+1)/delta)**2))**2)*(2*dz0*((iz-nz+1)/delta**2/dz))
            m = ix*nnz+iz
            mat[m,m] = -(w**2/v**2)+(2*Sx**2/dx**2)+(2*Sz**2/dz**2)
            #mat[m,m] = (-w**2/v**2)+(2/dx**2)+(2/dz**2)        
            if(iz!=0): 
                mat[m,m-1] = ((Sz*Szp)/(2*dz))-((Sz**2)/dz**2)
            if(iz!=nnz-1):                
                mat[m,m+1] = -((Sz*Szp)/(2*dz))-((Sz**2)/dz**2)
            if(ix!=nnx-1):
                mat[m,m+nnz] = -((Sx*Sxp)/(2*dx))-((Sx**2)/dx**2)
            if(ix!=0):
                mat[m,m-nnz] = ((Sx*Sxp)/(2*dx))-((Sx**2)/dx**2)
            

    temp=np.linalg.solve(mat, cf)
    for ix in range(nnx):
        green[ix, ifreq] = np.copy(temp[ix*nnz+2])
    # for iz in range(nnz):
    #     greenz[iz, ifreq] = np.copy(temp[nnx//2*nnz+iz])

for i in range(nnx):
    for ifreq in range(nf):
        u[i, ifreq] = green[i, ifreq]*source[ifreq]
    for ifreq in range(1, nf):
        u[i, nt-ifreq] = np.conj(u[i, ifreq])
# for i in range(nnz):
#     for ifreq in range(nf):
#         uz[i, ifreq] = greenz[i, ifreq]*source[ifreq]
#     for ifreq in range(1, nf):
#         uz[i, nt-ifreq] = np.conj(uz[i, ifreq])
u = np.fft.ifft(u)/nt
#uz = np.fft.ifft(uz)/nt
for it in range(nt):
    u[:, it] = u[:,it]*np.exp(alpha*it*dt)
    #uz[:, it] = uz[:,it]*np.exp(alpha*it*dt)

plt.xlabel('time')
plt.ylabel('x_dist')
plt.title("2D Wave Equation FDM Modeling in F-S")
# plt.subplot(121)
# plt.imshow(np.real(u),cmap='binary')
# plt.subplot(122)
# plt.imshow(np.real(uz),cmap='binary')
# plt.show()
plt.imshow(np.real(u), cmap='binary', aspect='auto', extent=[0, nx, 0, tmax])
plt.colorbar()
plt.show()