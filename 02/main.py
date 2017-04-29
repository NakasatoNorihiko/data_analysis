#!/usr/bin/env python

import numpy as np
from numpy.random import *
import matplotlib.pyplot as plt

#%% 三角多項式モデル
seed(0)

w = 0.1
W = 0.001

x = np.arange(-3, 3, w)
X = np.arange(-3, 3, W)

n = len(x)
N = len(X)

pix = np.pi * x
y = np.sin(pix) / (pix) + 0.1 * x + 0.05 * randn(n)

p = np.zeros([n, 32])
P = np.zeros([N, 32])

p[:, 1] = np.ones(n)
P[:, 1] = np.ones(N)

for i in range(1, 16):
    p[:, 2 * i] = np.sin(i / 2 * x)
    p[:, 2 * i + 1] = np.cos(i / 2 * x)
    P[:, 2 * i] = np.sin(i / 2 * X)
    P[:, 2 * i + 1] = np.cos(i / 2 * X)

t = np.linalg.pinv(p).dot(y)
F = P.dot(t) 

print(F, X.shape, F.shape, P.shape, t.shape)

plt.grid()
plt.xlim(-2.8, 2.8)
plt.ylim(-0.5, 1.2)
plt.plot(x, y, 'o')
plt.plot(X, F, '-')
plt.show()

#%% 
