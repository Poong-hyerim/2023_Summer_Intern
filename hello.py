import numpy as np
import matplotlib.pyplot as plt

# 파동 모델링에 사용할 매개변수 설정
nx = 100  # x 방향 격자점 수
nz = 100  # z 방향 격자점 수
nt = 1000  # 시간 스텝 수
dx = 10.0  # x 방향 격자점 간격 (미터)
dz = 10.0  # z 방향 격자점 간격 (미터)
dt = 0.001  # 시간 간격 (초)

# 파동 속도 모델 (이 부분은 필요에 따라 변경 가능)
velocity_model = np.ones((nz, nx)) * 1500.0  # 일단 모든 점에서 파동 속도를 1500 m/s로 가정합니다.
velocity_model[50:80, 20:40] = 2500.0  # 일부 영역에서 파동 속도를 2500 m/s로 변경합니다.

# 초기 파동장 설정 (시작 조건)
u = np.zeros((nz, nx))
u_old = np.zeros((nz, nx))

# 시간-공간 영역 파동 모델링 (유한 차분법)
for n in range(1, nt):
    # 2차원 파동방정식 유한 차분법으로 풀이
    for i in range(1, nx - 1):
        for j in range(1, nz - 1):
            u[j, i] = 2 * u_old[j, i] - u[j, i] + (velocity_model[j, i] * dt / dx) ** 2 * (
                    u[j, i + 1] - 2 * u[j, i] + u[j, i - 1]) + (velocity_model[j, i] * dt / dz) ** 2 * (
                              u[j + 1, i] - 2 * u[j, i] + u[j - 1, i])

    # 경계 조건 설정 (벽으로 가정)
    u[0, :] = 0.0
    u[-1, :] = 0.0
    u[:, 0] = 0.0
    u[:, -1] = 0.0

    # 이전 파동장 저장
    u_old = u.copy()

# 시그널을 얻기 위해 특정 위치에서 파동장 기록 (z=50 위치로 가정)
seismogram = u[50, :]

# x-t 평면에 대한 seismogram 출력
x_values = np.arange(0, nx * dx, dx)
t_values = np.arange(0, nt * dt, dt)
plt.imshow(u, cmap='seismic', aspect='auto', extent=[0, nx * dx, nt * dt, 0])
plt.colorbar(label='Amplitude')
plt.xlabel('x (m)')
plt.ylabel('t (s)')
plt.title('Seismogram in x-t plane at z = 50m')
plt.show()
