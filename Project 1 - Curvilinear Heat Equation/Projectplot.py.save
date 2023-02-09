import matplotlib.pyplot as plt
import math

L2X = [0.35287, 0.26399, 0.193936]
# L2Y = [2.55850487372923e-02, 6.32110849680624e-03, 1.56924708511860e-03]
grid = [1/10, 1/20, 1/40]

slopex = (math.log(L2X[-1]) - math.log(L2X[0])) / (math.log(grid[-1]) - math.log(grid[0]))
# slopey = (math.log(L2Y[-1]) - math.log(L2Y[0])) / (math.log(grid[-1]) - math.log(grid[0]))
print("X Order of Accuracy is equal to " + str(slopex))
# print("Y Order of Accuracy is equal to " + str(slopey))
plt.loglog(grid, L2X, label='L2X Norm Vs. Grid Size')
# plt.loglog(grid, L2Y, label='L2Y Norm Vs. Grid Size')
plt.grid(which='both')
plt.legend(loc='best')
plt.title("L2 Norm Vs. Grid Size")
plt.xlabel('Grid Size')
plt.ylabel("L2 Norm")
plt.show()
