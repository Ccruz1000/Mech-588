# Import packages
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import time

# User defined packages

# Plot a single timestep
def Plot_Single(datafile, title):
    data = np.loadtxt(datafile, dtype=float, skiprows=1, delimiter=' ')

    # Organize data
    x = data[:, 0]
    rho = data[:, 1]
    u = data[:, 2]
    P = data[:, 3]
    E = data[:, 4]

    # Plot condition data
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
    plt.title(title, fontsize=titlefont)

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

def Plot_Timesteps(datafile, title):
    # Read files
    files = glob.glob(datafile + '\*.txt')
    i = 0  # Select how often to plot
    numplot = 5  # Choose how many timesteps to plot
    plotcheck = int(len(files) / numplot)
    for file in files:
        if (i % plotcheck == 0):
            data = np.loadtxt(file, dtype=float, skiprows=1, delimiter=' ')

            # Organize data
            x = data[:, 0]
            rho = data[:, 1]
            u = data[:, 2]
            P = data[:, 3]
            E = data[:, 4]

            # Plot condition data
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
            plt.title(title + str(i), fontsize=titlefont)

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

            plt.show()
            plt.close()
        i += 1

# Get path to retrieve data and save results

path = os.getcwd()
save_path = path + '\Results'

# Define Steger filenames
Stegertimepath = path + '\StegerTimesteps'
initialStegerdata = save_path + '\Init_Steger_cond.txt'  # open initial conditions file
finalStegerdata = save_path + '\Final_Steger_cond.txt' # Open final conditions file

# Define Roe filenames
Roetimepath = path + '\RoeTimesteps'
initialRoedata = save_path + '\Init_Roe_cond.txt'  # open initial conditions file
finalRoedata = save_path + '\Final_Roe_cond.txt' # Open final conditions file

# # Plot Steger
Plot_Single(initialStegerdata, "Steger Initial Conditions")
Plot_Single(finalStegerdata, "Steger Final Conditions")
# Plot_Timesteps(Stegertimepath, "Steger Timestep ")

# Plot Roe
Plot_Single(initialRoedata, "Roe Initial Conditions")
Plot_Single(finalRoedata, "Roe Final Conditions")
# Plot_Timesteps(Roetimepath, "Roe Timestep ")





