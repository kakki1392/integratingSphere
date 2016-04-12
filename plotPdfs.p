import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

plt.figure(1)

pdfs = np.loadtxt("pdfs.txt")
s = pdfs[:,0]

plt.plot(s, pdfs[:,1], label='test1')
plt.plot(s, pdfs[:,2], label='test2')
plt.plot(s, pdfs[:,3], label='test3')
plt.plot(s, pdfs[:,4], label='test4')
plt.plot(s, pdfs[:,5], label='test5')

plt.legend()
plt.show()

