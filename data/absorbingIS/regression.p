import numpy as np
from scipy import stats
data = np.loadtxt("absorbing.txt")
print(data[0,:])
slope, intercept, r_value, p_value, std_err = stats.linregress(data[:,0], data[:,2])
print slope
