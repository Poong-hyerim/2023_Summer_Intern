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

    return w
        
def main():
    #변수 선언
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
    v_wave = 2
    
    f=np.zeros([nx,nz])
    f[nx//2][1]=1
    
    result_array=fdgaus(fmax, dt, nt)
    plt.plot(result_array)
    plt.show()   
    
    u1=np.zeros((nx,nz))
    u2=np.zeros((nx,nz))
    u3=np.zeros((nx,nz))
    u=np.zeros((nt,nx,nz))

    #반복문
    for it in range(nt):
        for ix in range(1,nx-1):
            for iz in range(1,nz-1):
                u3[ix][iz] = v_wave ** 2 * dt ** 2 * ( 
                            u2[ix + 1][iz] - 2 * u2[ix][iz] + u2[ix - 1][iz]) / dx ** 2 + v_wave ** 2 * dt ** 2 * ( \
                            u2[ix][iz + 1] - 2 * u2[ix][iz] + u2[ix][iz - 1]) / dz ** 2 + v_wave ** 2 * dt ** 2 * \
                            result_array[it] * f[ix][iz] + 2 * u2[ix][iz] - u1[ix][iz]
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
    plt.show()

if __name__ == "__main__":
    main()