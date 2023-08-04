import numpy as np
import matplotlib.pyplot as plt

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
    xmax = 1.0
    dx = 0.005
    nx = int(xmax / dx)
    vmax = 1.2
    tmax = 1.0
    dt = 0.001
    nt = int(tmax / dt)
    fmax = 25.0

    u1 = np.zeros(nx)
    u2 = np.zeros(nx)
    u3 = np.zeros(nx)
    f = np.zeros(nx)
    w = np.zeros(nt)

    f[nx // 2] = 1.0

    fdgaus(w, fmax, dt, nt)
    plt.plot(w)
    plt.show()

    seismogram = np.zeros((nt, nx))  # 2차원 배열로 seismogram 저장

    for it in range(nt):
        for ix in range(1, nx - 1):
            u3[ix] = (vmax ** 2) * (dt ** 2) * (
                    (u2[ix + 1] - 2 * u2[ix] + u2[ix - 1]) / dx ** 2 + w[it] * f[ix]) + 2.0 * u2[ix] - u1[ix]
        for ix in range(nx):
            u1[ix] = u2[ix]
            u2[ix] = u3[ix]
            #단방향으로 푼 파동방정식 값을 대입해줘서 양 옆을 없앤다!
        u2[0] = vmax * dt / dx * (u1[1] - u1[0]) + u1[0]
        u2[nx-1] = vmax * dt / dx * (u1[nx - 2] - u1[nx-1]) + u1[nx-1]

        seismogram[it, :] = u3  # 지진 센서 위치에 해당하는 데이터를 seismogram에 저장

    # Seismogram을 2차원 그래프로 시각화합니다.
    plt.imshow(seismogram, aspect='auto', extent=[0, xmax, tmax, 0], cmap='binary')
    plt.colorbar(label="Amplitude")
    plt.xlabel("X")
    plt.ylabel("Time (s)")
    plt.title("1D Wave Equation FDM Modeling in T-S")
    plt.show()

if __name__ == "__main__":
    main()
