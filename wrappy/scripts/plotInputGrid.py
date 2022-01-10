import matplotlib.pyplot as plt
import numpy as np
import sys

X = np.loadtxt(sys.argv[1], skiprows=3)

phi0c = X[:,0]
rho0c = X[:,1]

phi0c = phi0c/np.max(phi0c)
rho0c = rho0c/np.max(rho0c)

xg,yg = np.meshgrid(phi0c, rho0c)

plt.figure()
plt.scatter(xg, yg)
plt.xlabel('phi0c')
plt.ylabel('rho0c')
plt.title('normalized')
plt.show()

