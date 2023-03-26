import matplotlib.pyplot as plt
import math

ABS2 = [3.687916700571e-01, 7.849697758724e-01]
Midpoint = [3.687916700554e-01, 7.849675203488e-01]
# RK2 = [3.687916700554e-01, 7.849675203488e-01]
RK2 = [3.688876084236e-01, 7.879862040321e-01]
grid = [36/math.sqrt(32), 36/math.sqrt(2068)]

slopeABS2 = (math.log(ABS2[-1]) - math.log(ABS2[0])) / (math.log(grid[-1]) - math.log(grid[0]))
slopeMidpoint = (math.log(Midpoint[-1]) - math.log(Midpoint[0])) / (math.log(grid[-1]) - math.log(grid[0]))
slopeRK2 = (math.log(RK2[-1]) - math.log(RK2[0])) / (math.log(grid[-1]) - math.log(grid[0]))

print("ABS2 Order of Accuracy is equal to " + str(slopeABS2))
print("Midpoint Order of Accuracy is equal to " + str(slopeMidpoint))
print("RK2 Order of Accuracy is equal to " + str(slopeRK2))

plt.loglog(grid, ABS2, label='ABS2')
plt.loglog(grid, Midpoint, label='Midpoint')
plt.loglog(grid, RK2, label='RK2')

# plt.loglog(grid, L2Y, label='L2Y Norm Vs. Grid Size')
plt.grid(which='both')
plt.legend(loc='best')
plt.title("L2 Norm Vs. Grid Size")
plt.xlabel('Grid Size')
plt.ylabel("L2 Norm")
plt.show()
