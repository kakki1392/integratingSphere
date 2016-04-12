import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

plt.figure(1)

s = np.linspace(0.0, 250, 500)
pdfs = np.zeros((500,6))
pdfs[:,0] = s
for i in range(5,10):
        myFile = "data/emptySphereR1.0Rho0.9" + str(i) + ".dat"
        data = np.loadtxt(myFile)
        pathLength = data[:,1]
        kde = stats.gaussian_kde(pathLength)
        kdepdf = kde.evaluate(s)
        pdfs[:,i-4] = kdepdf
        plt.plot(s,kdepdf)
np.savetxt('pdfs.txt', pdfs, delimiter=' ')
plt.show()

#hist1, edges1 = np.histogram(collision1, 80, range=(0,150), density=True)
#hist2, edges2 = np.histogram(collision2, 80, range=(0,150), density=True)
#plt.plot(edges1[0:-1], hist1, 'r', edges2[0:-1], hist2, 'b')
#plt.xlim([0,150])
#plt.show()
