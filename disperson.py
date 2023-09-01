import numpy as np
import matplotlib.pyplot as plt


theta = [0, 15, 30, 45]
v=1.5
fmax=1.5
pi=np.arccos(-1)

G = np.linspace(2,100,500)
G_i = np.zeros(500)
Vph_v = np.zeros(500)

for t in theta:
    for i in range(len(G)):
        G_i[i] = 1/G[i]
        Vph_v[i] =np.sqrt(np.sin(pi/G[i]*np.cos(t))**2+np.sin(pi/G[i]*np.sin(t))**2)/(pi/G[i])
    X = G_i[:]
    Y = Vph_v[:]
    plt.plot(X, Y)
    
plt.xlim(np.min(G_i),np.max(G_i)/2)
plt.ylim(0.93,1.01)

plt.yticks(np.arange(0.94,1.01,0.01))
plt.xticks(np.arange(np.min(G_i),np.max(G_i)/2,0.01))
#print([np.min(Vph_v),np.max(Vph_v)])
#print([np.min(G_i),np.max(qG_i)])
plt.xlabel('1/G')
plt.ylabel('Vph/v')
plt.legend(theta, loc='best')
plt.grid()
plt.show()

print(1/0.11)