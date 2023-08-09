import numpy as np
import matplotlib.pyplot as plt

# Gaussian의 1차 미분형을 wavelet으로 사용하기 위한 사용자 정의 함수
def fdgaus(w, cutoff, dt, nt):
    phi = 4 * np.arctan(1.0)
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

def main():
    #변수 선언
    xmax = 1.0
    dx = 0.005
    nx = int(xmax / dx)
    vmax = 1.2
    tmax = 1.0
    dt = 0.001
    nt = int(tmax / dt)
    fmax = 25.0
    #배열 선언
    u1 = np.zeros(nx)
    u2 = np.zeros(nx)
    u3 = np.zeros(nx)
    f = np.zeros(nx)
    w = np.zeros(nt)
    #nx 중간 지점에만 source를 1로 설정(중앙부에서 source 발파의미)
    f[nx // 2] = 1.0
    #wavelet 생성
    fdgaus(w, fmax, dt, nt)
    #커브 모양 확인
    plt.plot(w)
    plt.show()
    #seismogram 추출을 위한 배열 추가 선언
    seismogram = np.zeros((nt, nx))
    #u3 계산을 위한 반복문
    for it in range(nt):
        for ix in range(1, nx - 1):
            u3[ix] = (vmax ** 2) * (dt ** 2) * (
                    (u2[ix + 1] - 2 * u2[ix] + u2[ix - 1]) / dx ** 2 + w[it] * f[ix]) + 2.0 * u2[ix] - u1[ix]
        for ix in range(nx):
            u1[ix] = u2[ix]
            u2[ix] = u3[ix]
        #매 it에 대한 u3계산 값을 seismogram에 저장
        seismogram[it, :] = u3
    # Seismogram 시각화
    plt.imshow(seismogram, aspect='auto', extent=[0, xmax, tmax, 0], cmap='binary')
    plt.colorbar(label="Amplitude")
    plt.xlabel("X")
    plt.ylabel("Time (s)")
    plt.title("1D Wave Equation FDM Modeling in T-S")
    plt.show()

if __name__ == "__main__":
    main()
