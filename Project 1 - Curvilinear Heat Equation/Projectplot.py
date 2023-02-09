import matplotlib.pyplot as plt
import math

# L2X = [0.35287, 0.26399, 0.193936]
# L2Y = [2.55850487372923e-02, 6.32110849680624e-03, 1.56924708511860e-03]
# L2DXDETA = [9.038810282333e-03, 2.576592491063e-03, 6.831326554923e-04]
# L2DXDXI = [5.682364126323e-02, 2.823597735774e-02, 1.407873016788e-02]
# L2DYDETA = [9.874171399161e-03, 2.683241860139e-03, 6.965737007848e-04]
# L2DYDXI = [1.108703697596e-02, 4.057991669475e-03, 1.705857812837e-03]
L2DXDETA = [4.938263900465e-03, 1.246543322903e-03, 3.122956535715e-04]
L2DXDXI = [5.707798228544e-02, 2.822170368627e-02, 1.407418046097e-02]
L2DYDETA = [4.938310503551e-03, 1.246541830174e-03, 3.122947972526e-04]
L2DYDXI = [7.465451523880e-03, 3.299502580964e-03, 1.585317886082e-03]
grid = [1/10, 1/20, 1/40]

slopeDXDETA = (math.log(L2DXDETA[-1]) - math.log(L2DXDETA[0])) / (math.log(grid[-1]) - math.log(grid[0]))
slopeL2DXDXI = (math.log(L2DXDXI[-1]) - math.log(L2DXDXI[0])) / (math.log(grid[-1]) - math.log(grid[0]))
slopeL2DYDETA = (math.log(L2DYDETA[-1]) - math.log(L2DYDETA[0])) / (math.log(grid[-1]) - math.log(grid[0]))
slopeL2DYDXI = (math.log(L2DYDXI[-1]) - math.log(L2DYDXI[0])) / (math.log(grid[-1]) - math.log(grid[0]))



# slopey = (math.log(L2Y[-1]) - math.log(L2Y[0])) / (math.log(grid[-1]) - math.log(grid[0]))
print("L2DXDETA Order of Accuracy is equal to " + str(slopeDXDETA))
print("L2DXDXI Order of Accuracy is equal to " + str(slopeL2DXDXI))
print("L2DYDETA Order of Accuracy is equal to " + str(slopeL2DYDETA))
print("L2DYDXI Order of Accuracy is equal to " + str(slopeL2DYDXI))
# print("Y Order of Accuracy is equal to " + str(slopey))
plt.loglog(grid, L2DXDETA, label='L2DXDETA Norm Vs. Grid Size')
plt.loglog(grid, L2DXDXI, label='L2DXDXI Norm Vs. Grid Size')
plt.loglog(grid, L2DYDETA, label='L2DYDETA Norm Vs. Grid Size')
plt.loglog(grid, L2DYDXI, label='L2DYDXI Norm Vs. Grid Size')
# plt.loglog(grid, L2Y, label='L2Y Norm Vs. Grid Size')
plt.grid(which='both')
plt.legend(loc='best')
plt.title("L2 Norm Vs. Grid Size")
plt.xlabel('Grid Size')
plt.ylabel("L2 Norm")
plt.show()
