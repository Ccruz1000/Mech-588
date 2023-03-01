import matplotlib.pyplot as plt
import math

dx_dxi = [5.707798228544e-02, 2.822170368627e-02, 1.407418046097e-02]
dx_deta = [4.938263900313e-03, 1.246543322775e-03, 3.122956536154e-04]
dy_dxi = [3.224174737130e-02, 1.574351892013e-02, 7.826952764783e-03]
dy_deta = [4.938310503551e-03, 1.246541830174e-03, 3.122947972526e-04]
grid = [1/16, 1/32, 1/64]

slopedx_dxi = (math.log(dx_dxi[-1]) - math.log(dx_dxi[0])) / (math.log(grid[-1]) - math.log(grid[0]))
slopedx_deta = (math.log(dx_deta[-1]) - math.log(dx_deta[0])) / (math.log(grid[-1]) - math.log(grid[0]))
slopedy_dxi = (math.log(dy_dxi[-1]) - math.log(dy_dxi[0])) / (math.log(grid[-1]) - math.log(grid[0]))
slopedy_deta = (math.log(dy_deta[-1]) - math.log(dy_deta[0])) / (math.log(grid[-1]) - math.log(grid[0]))



# slopey = (math.log(L2Y[-1]) - math.log(L2Y[0])) / (math.log(grid[-1]) - math.log(grid[0]))
print("dx_dxi Order of Accuracy is equal to " + str(slopedx_dxi))
print("dx_deta Order of Accuracy is equal to " + str(slopedx_deta))
print("dy_dxi Order of Accuracy is equal to " + str(slopedy_dxi))
print("dy_deta Order of Accuracy is equal to " + str(slopedy_deta))
# print("Y Order of Accuracy is equal to " + str(slopey))
plt.loglog(grid, dx_dxi, label='dx_dxi Norm Vs. Grid Size')
plt.loglog(grid, dx_deta, label='dx_deta Norm Vs. Grid Size')
plt.loglog(grid, dy_dxi, label='dy_dxi Norm Vs. Grid Size')
plt.loglog(grid, dy_deta, label='dy_deta Norm Vs. Grid Size')
# plt.loglog(grid, L2Y, label='L2Y Norm Vs. Grid Size')
plt.grid(which='both')
plt.legend(loc='best')
plt.title("L2 Norm Vs. Grid Size")
plt.xlabel('Grid Size')
plt.ylabel("L2 Norm")
plt.show()
