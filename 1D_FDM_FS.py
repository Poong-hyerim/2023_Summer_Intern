import numpy as np
import matplotlib.pyplot as plt

def ricker_wavelet(t, f0, t0):
    arg = (np.pi * f0 * (t - t0)) ** 2
    return (1 - 2 * arg) * np.exp(-arg)

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
    