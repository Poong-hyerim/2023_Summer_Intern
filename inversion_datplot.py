import numpy as np
import matplotlib.pyplot as plt

rectype = np.dtype(np.float32)
data = np.fromfile('cut_marmousiela_vp_400_80_20m.dat',dtype=rectype)
model_data = np.transpose(np.reshape(data,(400,80)))
plt.imshow(model_data)
plt.show()

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