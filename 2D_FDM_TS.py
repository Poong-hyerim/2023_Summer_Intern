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
def show_in_plot(imgs, perc=99):
    clip = max(np.percentile(imgs, perc), -np.percentile(imgs, 100 - perc))
    #plt.xticks([])
    #plt.yticks([])
    #plt.gca().set_title(titles[i])
    plt.xlabel('x_dist')
    plt.ylabel('time')
    plt.title('FDM_2d_time \n receiver on ground(z=0)')
    plt.imshow(imgs.squeeze(), cmap='binary', vmin=-clip, vmax=clip, aspect='auto')
    plt.show()

#main 함수
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
    vmax = 2

    #np.zero로 배열 선언시 에 복소수 활용은 type을 지정해야함 
    #배열 선언
    #0~nx-1까지 nx개, 0~nz-1까지 nz개의 배열을 선언
    u1 = np.zeros((nx, nz))
    u2 = np.zeros((nx, nz))
    u3 = np.zeros((nx, nz))
    f = np.zeros((nx, nz))
    #it에 대해서 계산할거라 nt개만큼 배열 선언
    U = np.zeros((nt, nx, nz)) 

    #f설정, 소스위치 설정
    #nx 중간 지점에 소스 1.0을 넣어
    #x위치로는 중간, z=1인 지점만 1로 소스를 넣어줌
    index = int((nx)/2)
    f[index, 1] = 1.0
    print(f"source : {f[index, 1]}")
    #wavelet 계산 적용
    w=fdgaus(fmax, dt, nt)
    
    #u3 계산 식 들어가기 > 0~ix,nz 랑 nx, 0~iz 까지 값(경계 값 계산..!?)
    # 0~n-1까지 연산을 해봅시당!
    for it in range(nt):    
        print(it)    
        for ix in range(1, nx-1):
            for iz in range(1, nz-1):
                u3[ix][iz] = vmax ** 2 * dt ** 2 * (u2[ix + 1][iz] - 2 * u2[ix][iz] + u2[ix - 1][iz]) / dx ** 2 \
                    + vmax ** 2 * dt ** 2 * (u2[ix][iz + 1] - 2 * u2[ix][iz] + u2[ix][iz - 1]) / dz ** 2 \
                    + vmax ** 2 * dt ** 2 * w[it] * f[ix][iz] + 2 * u2[ix][iz] - u1[ix][iz]      
        U[it, :, :] = u3[:,:]
        u1[:,:]=u2[:,:]
        u2[:,:]=u3[:,:]
        # it에 대한 모든 u3를 다 저장하기 위함!
        # 계산결과를 U[2it, :, :]=u3에 차례로 넣어준다
    """ show_in_plot(U[:,:,1]) """
    
    result = U[:,:,1]
    perc = 99
    boundary = max(np.percentile(result, perc), -np.percentile(result, 100-perc))
    
    plt.xlabel('x_dist')
    plt.ylabel('time')
    plt.title("2D Wave Equation FDM Modeling in T-S")
    plt.imshow(result.squeeze(), cmap='binary', aspect='auto', vmin= -boundary, vmax=boundary)
    plt.show()

    print(f"boundary condition : {-boundary, boundary}")
    print(f"U.squeeze {U[:,:,1].squeeze()} U.sq_shpae : {U[:,:,1].squeeze().shape}")
    
if __name__ == "__main__":
    main()


