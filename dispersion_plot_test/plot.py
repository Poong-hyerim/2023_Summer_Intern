import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('disper1d09.dat')


x=data[:,0]
y=data[:,1]
plt.plot(x,y, marker='o')
plt.legend
plt.grid(True)
plt.show()