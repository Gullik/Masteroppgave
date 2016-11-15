import numpy as np
import pylab as plt

x = (1, (2)**3, (3)**3, (4)**3)
y = (1.28, 11.76, 28.9, 54.65)


plt.figure()
plt.plot(x,y, "*")
plt.plot(x,y)
plt.ylabel("Time (s)")
plt.xlabel("#Processors")
plt.text(2**3,10, "$128^3$")

plt.title("Cores vs Time")
plt.grid()
plt.savefig("")
plt.show()
