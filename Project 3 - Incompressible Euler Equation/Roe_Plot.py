# Import packages
import numpy as np
import matplotlib.pyplot as plt
import os

# User defined packages

# Get path to retrieve data and save results
path = os.getcwd()
save_path = path + '\Results'

# Initialize variables
datafile = save_path + '\Final_cond.txt'  # open initial conditions file
data = np.loadtxt(datafile, dtype=float, skiprows=1, delimiter=' ')

# Organize data
x = data[:, 0]
rho = data[:, 1]
u = data[:, 2]
P = data[:, 3]
E = data[:, 4]

# Plot initial condition data
fig = plt.subplots()

numplot = 4
axisfont = 16
titlefont = 18

# Plot rho
ax1 = plt.subplot(numplot, 1, 1)
plt.plot(x, rho, 'forestgreen')
plt.ylabel('ρ\n', fontsize=axisfont)
plt.grid(True)
plt.tick_params(axis='x', bottom=False, labelbottom=False)
plt.xlim(x[0], x[-1])
plt.title("Initial Conditions", fontsize=titlefont)

# Plot u
ax2 = plt.subplot(numplot, 1, 2)
plt.plot(x, u, 'red')
plt.ylabel('u', fontsize=axisfont)
plt.grid(True)
plt.tick_params(axis='x', bottom=False, labelbottom=False)
plt.xlim(x[0], x[-1])

# Plot P
ax3 = plt.subplot(numplot, 1, 3)
plt.plot(x, P, 'blue')
plt.ylabel('P\n', fontsize=axisfont)
plt.grid(True)
plt.tick_params(axis='x', bottom=False, labelbottom=False)
plt.xlim(x[0], x[-1])

# Plot E
ax4 = plt.subplot(numplot, 1, 4)
plt.plot(x, E, 'orange')
plt.ylabel('E\n', fontsize=axisfont)
plt.grid(True)
plt.xlabel('x', fontsize=axisfont)
plt.xlim(x[0], x[-1])

# Adjust plots
plt.subplots_adjust(left=0.15)
plt.subplots_adjust(bottom=0.12)
plt.subplots_adjust(top=0.90)

# Save and show plots
plt.savefig(save_path + "\Initial_conditions.png")
plt.show()
plt.close()




