from pylab import*

data = np.loadtxt("q1_Integration.txt")
data = data.transpose()

plt.semilogx(data[0], data[1], 'ro-', label='Uniform Dist')
plt.semilogx(data[0], data[2], 'bo-', label='Exponential Dist')
plt.semilogx(data[0],np.ones(len(data[0]))*0.74682, 'k--', label='Exact')

plt.legend()
plt.grid()
plt.xlabel("N")
plt.ylabel("Integration Value")
plt.xlim(data[0][0]/2, data[0][-1]*5)
plt.show()
