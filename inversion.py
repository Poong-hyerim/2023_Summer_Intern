# 역산 코드 로직
#1.파라미터 설정(주파수 FDM+PML)
#2.초기 모델 가져오기
    #3.반복문1 진입(수렴될때까지)
    #4.행렬구성
    #5.행렬 FACTORIZE
        #6.반복문2 진입(각 샷마다 매트릭스 풀이)
        #7.source 할당
        #8.S*U=F 식 풀이
        #9.OBSERVED DATA 로딩
        #10.WAVELET 추정
        #11.Residual 계산
        #12.L2-norm목적함수 계산
        #13.Adjoint source 할당
    #14.back propagation 연산
    #15.UPDATE 방향 계산(SD)
    #16.파라미터 UPDATE > 각 frequency에 대해 계산 해 수렴될 때까지 update

import numpy as np
import matplotlib.pyplot as plt

def generate_linear_initial_model(nx, nz, delta):
    v = np.zeros((nx, nz), dtype=float)   
    gradient = 0.01
    for ix in range(nx):
        for iz in range(nz):
            # Linear gradient in both x and z directions
            v[ix, iz] = 1.5 + gradient * iz
    return v

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

#1.파라미터 설정(주파수 FDM+PML)
v = 2
delta = 30
dx = 0.005
dx0 = np.log(1000)*((3*v)/(delta*dx))
dz = 0.005
dz0 = np.log(1000)*((3*v)/(delta*dz))
nx = 200
nz = 40
nnx = nx+2*delta
nnz = nz+delta
tmax = 1
G = 10
fmax = np.min(v)/(G*dx)
dt = 0.05/fmax
nt = int(tmax/dt)
df = 1/tmax
nf = int(fmax/df)
pi = np.arccos(-1)
alpha = np.log(100)/tmax
source = np.zeros(nt)
u = np.zeros((nnx, nt), dtype=complex)
cf = np.zeros(nnx*nnz, dtype=complex)
mat = np.zeros((nnx*nnz, nnx*nnz), dtype=complex)
green = np.zeros((nnx, nf), dtype=complex)

import numpy as np
#2.초기 모델 가져오기
data = np.fromfile('cut_marmousiela_vp_400_80_20m.dat', dtype=np.dtype(np.float32))
origin_size = (400,80)
origin_model = np.transpose(np.reshape(data, origin_size))

new_size = (200,40)
new_model = np.zeros(new_size, dtype=data.dtype)
#80. 400
red_x = origin_size[0]//new_size[0]
red_z = origin_size[1]//new_size[1]
#200
for ix in range(new_size[0]):
    #40
    for iz in range(new_size[1]):
        xs = ix*red_x
        xe = (ix+1)*red_x
        zs = iz*red_z
        ze = (iz+1)*red_z
        avg = np.mean(origin_model[zs:ze, xs:xe])        
        new_model[ix,iz] = avg
new_model = np.transpose(new_model)
init = np.transpose(generate_linear_initial_model(nx, nz, delta))
# plt.subplot(121)
# plt.imshow(new_model)
# plt.colorbar()
# plt.subplot(122)
# plt.imshow(init)
# plt.colorbar()
# plt.show()


    #3.반복문1 진입(수렴될때까지)
    
    #4.행렬구성
    #5.행렬 FACTORIZE
        #6.반복문2 진입(각 샷마다 매트릭스 풀이)
        #7.source 할당
        #8.S*U=F 식 풀이
        #9.OBSERVED DATA 로딩
        #10.WAVELET 추정
        #11.Residual 계산
        #12.L2-norm목적함수 계산
        #13.Adjoint source 할당
    #14.back propagation 연산
    #15.UPDATE 방향 계산(SD)
    #16.파라미터 UPDATE > 각 frequency에 대해 계산 해 수렴될 때까지 update
