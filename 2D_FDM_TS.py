import numpy as np
import matplotlib.pyplot as plt

#wavelet관련 f항 계산용 함수 선언
def ricker_wavelet(t, f0, t0):
    arg = (np.pi * f0 * (t - t0)) ** 2
    return (1 - 2 * arg) * np.exp(-arg)

def fdgaus(w, cutoff, dt, nt):
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
    
    #배열 선언
    #0~nx-1까지 nx개, 0~nz-1까지 nz개의 배열을 선언
    u1 = np.zeros((nx+2, nz+2))
    u2 = np.zeros((nx+2, nz+2))
    u3 = np.zeros((nx+2, nz+2))
    f = np.zeros((nx+2, nz+2))
    #it에 대해서 계산할거라 nt개만큼 배열 선언
    w = np.zeros(nt)
    U = np.zeros((nt, nx+2, nz+2)) 

    #f설정, 소스위치 설정
    #nx 중간 지점에 소스 1.0을 넣어
    #x위치로는 중간, z=1인 지점만 1로 소스를 넣어줌
    index = int((nx+1)/2)
    f[index, 1] = 1.0
    print(f"source : {f[index, 1]}")
    #wavelet 계산 적용
    fdgaus(w, fmax, dt, nt)
    #u3 계산 식 들어가기 > 0~ix,nz 랑 nx, 0~iz 까지 값(경계 값 계산..!?)
    # 0~n-1까지 연산을 해봅시당!
    for it in range(nt):    
        print(it)    
        for ix in range(1, nx+1):
            for iz in range(1, nz+1):
                u3[ix][iz] = vmax ** 2 * dt ** 2 * (u2[ix + 1][iz] - 2 * u2[ix][iz] + u2[ix - 1][iz]) / dx ** 2 \
                    + vmax ** 2 * dt ** 2 * (u2[ix][iz + 1] - 2 * u2[ix][iz] + u2[ix][iz - 1]) / dz ** 2 \
                    + vmax ** 2 * dt ** 2 * w[it] * f[ix][iz] + 2 * u2[ix][iz] - u1[ix][iz]      
        U[it, :, :] = u3
        u1=u2
        u2=u3
        # it에 대한 모든 u3를 다 저장하기 위함!
        # 계산결과를 U[2it, :, :]=u3에 차례로 넣어준다    
    # Seismogram 시각화
    plt.imshow(U[:,:, 50], aspect='auto', extent=[0, xmax, tmax, 0], cmap='seismic')
    plt.colorbar()
    plt.xlabel('X')
    plt.ylabel('Time(s)')
    plt.title("2D_FDM Seismogram")
    plt.show()
    print(f"u3 shape: {u3.shape}, u3 value : {u3}")
    print(f"U shape: {U.shape}, U values: {U[it, :, :]}")

if __name__ == "__main__":
    main()