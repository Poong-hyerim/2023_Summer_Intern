import numpy as np
import matplotlib.pyplot as plt

datas = [
    "C:/Users/owner/Desktop/2023_Summer_Intern/dispersion_plot_test/dispersion/1Dstable/disper1d01.dat",
    "C:/Users/owner/Desktop/2023_Summer_Intern/dispersion_plot_test/dispersion/1Dstable/disper1d02.dat",
    "C:/Users/owner/Desktop/2023_Summer_Intern/dispersion_plot_test/dispersion/1Dstable/disper1d03.dat",
    "C:/Users/owner/Desktop/2023_Summer_Intern/dispersion_plot_test/dispersion/1Dstable/disper1d04.dat",
    "C:/Users/owner/Desktop/2023_Summer_Intern/dispersion_plot_test/dispersion/1Dstable/disper1d05.dat",
    "C:/Users/owner/Desktop/2023_Summer_Intern/dispersion_plot_test/dispersion/1Dstable/disper1d06.dat",
    "C:/Users/owner/Desktop/2023_Summer_Intern/dispersion_plot_test/dispersion/1Dstable/disper1d07.dat",
    "C:/Users/owner/Desktop/2023_Summer_Intern/dispersion_plot_test/dispersion/1Dstable/disper1d08.dat",
    "C:/Users/owner/Desktop/2023_Summer_Intern/dispersion_plot_test/dispersion/1Dstable/disper1d09.dat"
]
legend = []
for data in datas:
    d = np.loadtxt(data)
    x = d[:,0]
    y = d[:,1]
    plt.plot(x,y)

plt.legend()
plt.grid(True)
plt.show()