#!/usr/bin/env python

import numpy as np
from numpy.random import *
import matplotlib.pyplot as plt

#%% 三角多項式モデル
seed(0)

w = 0.1
W = 0.001

x = np.arange(-5, 5, w)
X = np.arange(-5, 5, W)

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

pT = p.T
p2 = pT.dot(p)
t = np.zeros(32) 
z = np.zeros(32)
u = np.zeros(32)
g = 0.01
for i in range(100):
    t0 = t
    z0 = z
    u0 = u
    t = np.linalg.inv(p2 + np.identity(len(p2))).dot(pT.dot(y) + z - u)
    for j in range(len(t)):
        z[j] = np.max([0, t[j] + u0[j] - g]) + np.min([0, t[j] + u0[j] + g])
    u = u0 + t - z

for i in range(len(t)):
    if np.abs(t[i]) < 0.001:
        t[i] = 0

F = P.dot(t)

print(t)

plt.grid()
plt.xlim(-4.8, 4.8)
plt.ylim(-0.5, 1.2)
plt.plot(x, y, 'o')
plt.plot(X, F, '-')
# plt.savefig('plot.png')
plt.show()
