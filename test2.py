import numpy as np
import matplotlib.pyplot as plt

def gaussian_derivative(t, sigma, mu):
    """
    Gaussian derivative (Gaussian function first derivative).

    Parameters:
        t (array-like): Time array where the derivative will be computed.
        sigma (float): Standard deviation of the Gaussian distribution.
        mu (float): Mean (center) of the Gaussian distribution.

    Returns:
        ndarray: Gaussian derivative values corresponding to the given time array.
    """
    gaussian = np.exp(-0.5 * ((t - mu) / sigma) ** 2)
    derivative = -((t - mu) / sigma**2) * gaussian
    return derivative

# Example usage
sampling_rate = 1000  # Hz
duration = 1.0  # seconds
t = np.linspace(0, duration, int(sampling_rate * duration))

sigma = 0.1  # Standard deviation of the Gaussian distribution
mu = 0.5     # Mean (center) of the Gaussian distribution

gaussian_derivative_wavelet = gaussian_derivative(t, sigma, mu)

# Plotting the Gaussian derivative wavelet
plt.figure()
plt.plot(t, gaussian_derivative_wavelet)
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.title('Gaussian Derivative Wavelet')
plt.grid(True)
plt.show()

""" 
xmax = 1
dx = 0.003
nx = 333
zmax = 1
dz = 0.003
nz = 333
tmax = 0.6
dt = 0.00075
nt = 800
fmax = 25
v_wave = 2v

f=np.zeros([nx,nz])
f[nx//2][1]=1

u1=np.zeros((nx,nz))
u2=np.zeros((nx,nz))
u3=np.zeros((nx,nz))
u=np.zeros((nt,nx,nz))
result_array = fdgaus(fmax, dt, nt)
#반복문
for it in range(nt):
    for ix in range(1,nx-1):
        for iz in range(1,nz-1):
            u3[ix][iz] = v_wave ** 2 * dt ** 2 * (u2[ix + 1][iz] - 2 * u2[ix][iz] + u2[ix - 1][iz]) / dx ** 2 + v_wave ** 2 * dt ** 2 * (u2[ix][iz + 1] - 2 * u2[ix][iz] + u2[ix][iz - 1]) / dz ** 2 + v_wave ** 2 * dt ** 2 * result_array[it] * f[ix][iz] + 2 * u2[ix][iz] - u1[ix][iz]
    u1[:,:]=u2[:,:]
    u2[:,:]=u3[:,:]
    u[it,:,:] = u3[:,:]
    print(it)

u_plot = np.zeros((nt,nx))
u_plot[:,:] = u[:,:,1]

plt.imshow(u_plot,cmap='binary',aspect='auto')
plt.colorbar()
plt.xlabel('x_dist')
plt.ylabel('time')
plt.title('FDM_2d_time \n receiver on ground(z=0)')
plt.show() """