#!/usr/bin/env ipython
#--------------------
import numpy as np
import matplotlib.pyplot as plt

din = np.loadtxt('out.dat')

x = din[:,0]
theta = din[:, 1]
y = din[:, 2]*-1.0

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.plot(x, theta, label='Theta')
ax.plot(x, y, label='y')
ax.legend()
plt.show()
