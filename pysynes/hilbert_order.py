# -*- coding: utf-8 -*-
"""
Convert from ordering in a Hilbert-like curve and a 1D array and
vice-versa.

Taken from:

   https://en.wikipedia.org/wiki/Hilbert_curve

@author: Nicolas Guarin-Zapata
"""
from __future__ import division
import matplotlib.pyplot as plt


def xy2d(n, x, y):
    d = 0
    s = n//2
    while s>0:
        rx = (x & s) > 0
        ry = (x & s) > 0
        d += s * s * ((3 * rx) ^ ry)
        x, y = rot(s, x, y, rx, ry)
        s //= 2
    return d


def d2xy(n, d):
    x = 0
    y = 0
    t = d
    s = 1
    while s < n:
        rx = 1 & (t//2)
        ry = 1 & (t ^ rx)
        x, y = rot(s, x, y, rx, ry)
        x += s * rx
        y += s * ry
        t //= 4
        s *= 2
    return x, y


def rot(n, x, y, rx, ry):
    if ry == 0:
        if rx == 1:
            x = n - 1 - x
            y = n - 1 - y

        x, y = y, x
    return x, y


if __name__ == "__main__":

    from skimage import data
    from skimage.color import rgb2gray
    import numpy as np
    from scipy.fftpack import ifft,  ifftshift
    from scipy.io import wavfile
    
    astro = rgb2gray(data.astronaut())
    n = astro.shape[0]
    astro_1D = np.zeros((2 * n**2))
    
    for row in range(n):
        for col in range(n):
            d = xy2d(n, row, col)
            astro_1D[n**2 + d + 1] = astro[row, col]
            astro_1D[n**2 - d] = astro[row, col]
    fmax = 2**14
    freq = np.linspace(-2**14, 2**14, 2*n**2)
    astro_sound = np.real(ifft(astro_1D, n=2 * n**2))
    astro_sound = ifftshift(astro_sound)
    time = 1/fmax * np.linspace(0, 2**12, 2*n**2)
    # Export the sound file
    # wavfile.write("astro_sound.wav", 2**12, astro_sound*1e6)
    
    
    plt.figure(figsize=(10,5))
    plt.subplot(1, 2, 1)
    plt.imshow(astro, cmap="gray")
    plt.subplot(2, 2, 2)
    plt.plot(freq, astro_1D)
    plt.subplot(2, 2, 4)
    plt.plot(time, astro_sound)
