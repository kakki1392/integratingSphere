import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

plt.figure(1)

s = np.linspace(0.0, 400, 300)
file1 = "data/emptySphereR1.0Rho0.95.dat"
file2 = "data/emptySphereR1.25Rho0.95.dat"
file3 = "data/emptySphereR1.5Rho0.95.dat"

data1 = np.loadtxt(file1)
data2 = np.loadtxt(file2)
data3 = np.loadtxt(file3)

path1 = data1[:,1]
path2 = data2[:,1]
path3 = data3[:,1]

mean1 = np.mean(path1)
mean2 = np.mean(path2)
mean3 = np.mean(path3)

std1 = np.std(path1, ddof=1)
std2 = np.std(path2, ddof=1)
std3 = np.std(path3, ddof=1)

print("R=1.0, mean: " + str(mean1) + ", std: " + str(std1) +"\n")
print("R=1.25, mean: " + str(mean2) + ", std: " + str(std2) +"\n")
print("R=1.5, mean: " + str(mean3) + ", std: " + str(std3) +"\n")



kde1 = stats.gaussian_kde(path1);
kde2 = stats.gaussian_kde(path2);
kde3 = stats.gaussian_kde(path3);

kdepdf1 = kde1.evaluate(s) 
kdepdf2 = kde2.evaluate(s) 
kdepdf3 = kde3.evaluate(s) 

plt.plot(s, kdepdf1, 'r', s, kdepdf2, 'b', s, kdepdf3, 'g')
plt.show()
