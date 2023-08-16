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

###[---주파수 영역의 1D 모델링 가즈아!---]###
#변수 설정 및 초기 설정
# 주어진 parameter로부터 필요한 값 계산해내기!(nx, dx, nt, dt, nf, df)
    #1. Nyquist frequency로 부터 dt 결정하기(fmax로부터 dt 도출)
    #2. nt 결정(tmax, dt로 nt 도출)
    #3. nf 결정(Tmax의 역수로 도출)
    #4. df 결정(fmax, df로 nf 도출)
    #5. 위의 파라미터들로 행렬식을 완성
pi = np.arccos(-1)
nx = 400
dx = 0.01
v = 2
fmax = 20
tmax = 2

dt = (0.5/fmax)
nt = int(tmax/dt)
df = 1/tmax
nf = int(fmax/df)

ndt = (0.1/fmax)
nnt = int(tmax/ndt)
alpha = np.log(100)/tmax

source = np.zeros(nt, dtype=complex)
mat = np.zeros((nx,nx), dtype=complex)
csource = np.zeros(nt, dtype=complex)
u = np.zeros((nx,nt), dtype=complex)
nu = np.zeros((nx,nnt), dtype=complex)
f = np.zeros(nx, dtype=complex)
green = np.zeros((nx,nf), dtype=complex)

#source를 선언하고 복소수 형태로 바꿔준 후 푸리에 변환
#t domain의 소스를 생성!
source = fdgaus(fmax, dt, nt)
#복소수로 source를 변환해준다
for it in range(nt):
    csource[it] = source[it]*np.exp(-alpha*it*dt) 

#csource를 주파수 영역의 소스로 바꿔준다!
csource = np.fft.fft(csource)

#얻어진 식은 행렬 matrix연산을 통해서 u에 대한 값계산이 가능하다.
#[계수 행렬]*[u] = [f]이므로 우리는 [u] = [계수행렬-1(역행렬)]*[f]를 진행할 예정
    #1. 계수 행렬을 계산하고 역행렬 취하기 tridiagonal matrix > CSY subroutine(파이썬에서 방법으로 변경)
    #2. f항에 wavelet 설정하기
    #3. 1에서 구한 역행렬*f를 곱해서 u에 대한 식 풀이하기
#주파수에따라 이를 반복적으로 계산(주파수를 변화시키며 tridiagonal matrix 생성 후 풀이)
#w는 1~df 범위내에서 df 간격으로 변화 시키는 중

for ifreq in range(1, nf):
    #w를 변화시키면서
    print(ifreq)
    w = 2.0 * pi * (ifreq) * df + 1j * alpha
    #매번 f는 초기화 후 중앙 부분에서 발파하는 것으로 설정
    f[:] = 0.0
    f[nx//2] = 1.0
    #행렬 구성하기
    mat[:] = 0.0
    mat[0, 0] = -(w ** 2 / v ** 2) + (2 / dx ** 2)
    mat[0, 1] = -1.0 / dx**2
    mat[nx - 1, nx - 1] = -(w ** 2 / v ** 2) + (2 / dx ** 2)
    mat[nx - 1, nx - 2] = -1.0 / dx**2
    
    for i in range(1, nx - 1):
        mat[i, i] = -(w ** 2 / v ** 2) + (2 / dx ** 2)
        mat[i, i - 1] = -(1.0 / dx ** 2)
        mat[i, i + 1] = -(1.0 / dx ** 2)

    #mat*u = f에서 u를 구해줌!
    f = np.linalg.solve(mat, f)
    #설정한 f(중앙발파)와 행렬(tridiagonal matrix)를 대입!
    
    green[:, ifreq] = f[:]   
    
    #green func로 주파수마다 내역을 저장해주기

#fdgaus 파형*green 배열로 fdgaus source 설정
for ix in range(nx):
    #파형적용
    for ifreq in range(nf):
        nu[ix, ifreq] = green[ix, ifreq] * csource[ifreq]
    #conjugate로 반쪽 만들어주깅!
    for ifreq in range(1,nf):
        nu[ix, (nt-1)-ifreq+1] = np.conj(nu[ix, ifreq])
        
    # inverse fourier를 통해서 다시 time domain으로 돌아간 후 u를 seismogram으로 뽑아보기
    # *** inverse 후에 u결과를 nt로 나눠주어야 처음과 진폭을 동일하게 맞출 수 있다(왜!?)
    #u = np.fft.fft(u)
nu = np.fft.ifft(nu) / nt

for it in range (nt):
    nu[:,it]=nu[:,it]*np.exp(-alpha*it*dt)

plt.xlabel('x_dist')
plt.ylabel('time')
plt.title("1D Wave Equation FDM Modeling in F-S")
plt.imshow(np.real(nu), cmap='binary', aspect='auto', extent=[0, nx, 0, tmax])
plt.colorbar()  # Optionally, add a color bar for reference 
plt.show()