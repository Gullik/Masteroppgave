import numpy as np
import pylab as plt

x = (1, (2)**3, (3)**3, (4)**3)
y = (1.28, 11.76, 28.9, 54.65)


plt.figure()
plt.plot(x,y, "*")
plt.plot(x,y)
plt.ylabel("Time (s)")
plt.xlabel("#Processors")
plt.text(2**3-3,13, "$64^3$")
plt.text(3**3-3,30, "$96^3$")
plt.text(4**3-3,50, "$128^3$")
plt.title("Cores vs Time")
plt.grid()
plt.savefig("scalingMG")
plt.show()
