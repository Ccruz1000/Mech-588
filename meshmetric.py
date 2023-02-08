#Import Packages
import numpy as np
import matplotlib.pyplot as plt

#********************************************************************
# Curvillinear Mesh Programming Assignment Mesh Metric Generation
# Mech 588 - Advanced CFD UBC Mechanical Engineering Winter 2023
# Christian Rowsell (40131393)
#********************************************************************

# User defined functions

# Define function to calculate xi and eta given x and y
def meshmetric(xi, eta):
    x = xi + np.sin((eta/10) * np.pi)
    y = eta + np.cos((xi/10) * np.pi)
    return x, y

# Define boundaries
miny = 0
minx = 0
maxx = 5
maxy = 5

# Define number of divisions
nx = 11
ny = 11

# Define X and Y coordinates
eta_list = np.linspace(miny, maxy, ny)
xi_list = np.linspace(minx, maxx, nx)
xi, eta = np.meshgrid(xi_list, eta_list)
print(xi_list)

# Calculate xi and eta coordinates
x, y = meshmetric(xi, eta)

plt.plot(x, y)
plt.plot(y,x)
plt.ylim(miny, maxy)
plt.xlim(minx, maxx)
plt.show()
