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

# Created interactive plot of many timeseries
# def Plot_TimeSeries(timepath, title):
#     # Create list of files for the time series
#     files = glob.glob(timepath + '\*.txt')
#     data = np.loadtxt(files[0], dtype=float, skiprows=1, delimiter=' ')
#
#     # Create initial plot
#     # Organize data
#     x = data[:, 0]
#     rho = data[:, 1]
#     u = data[:, 2]
#     P = data[:, 3]
#     E = data[:, 4]
#
#     # # run GUI event loop
#     # plt.ion()
#
#     # Plot condition data
#     fig = plt.subplots()
#
#     numplot = 4
#     axisfont = 16
#     titlefont = 18
#
#     # Plot rho
#     ax1 = plt.subplot(numplot, 1, 1)
#     rholine = plt.plot(x, rho, 'forestgreen')
#     plt.ylabel('ρ\n', fontsize=axisfont)
#     plt.grid(True)
#     plt.tick_params(axis='x', bottom=False, labelbottom=False)
#     plt.xlim(x[0], x[-1])
#     plt.title(title, fontsize=titlefont)
#
#     # Plot u
#     ax2 = plt.subplot(numplot, 1, 2)
#     uline = plt.plot(x, u, 'red')
#     plt.ylabel('u', fontsize=axisfont)
#     plt.grid(True)
#     plt.tick_params(axis='x', bottom=False, labelbottom=False)
#     plt.xlim(x[0], x[-1])
#
#     # Plot P
#     ax3 = plt.plot(numplot, 1, 3)
#     Pline = plt.subplot(x, P, 'blue')
#     plt.ylabel('P\n', fontsize=axisfont)
#     plt.grid(True)
#     plt.tick_params(axis='x', bottom=False, labelbottom=False)
#     plt.xlim(x[0], x[-1])
#
#     # Plot E
#     ax4 = plt.subplot(numplot, 1, 4)
#     Eline = plt.plot(x, E, 'orange')
#     plt.ylabel('E\n', fontsize=axisfont)
#     plt.grid(True)
#     plt.xlabel('x', fontsize=axisfont)
#     plt.xlim(x[0], x[-1])
#
#     # Adjust plots
#     plt.subplots_adjust(left=0.15)
#     plt.subplots_adjust(bottom=0.12)
#     plt.subplots_adjust(top=0.90)
#
#     # Create subplots
#     for file in files:
#
#         # Open timestep data
#         newdata = np.loadtxt(file, dtype=float, skiprows=1, delimiter=' ')
#
#         # Organize data
#         new_rho = newdata[:, 1]
#         new_u = newdata[:, 2]
#         new_P = newdata[:, 3]
#         new_E = newdata[:, 4]
#
#         # Update plot values
#         rholine.set_xdata(x)
#         rholine.set_ydata(new_rho)
#
#         uline.set_xdata(x)
#         uline.set_ydata(new_u)
#
#         Pline.set_xdata(x)
#         Pline.set_ydata(new_P)
#
#         Eline.set_xdata(x)
#         Eline.set_ydata(new_rho)
#
#         # Draw updated values
#         fig.canvas.draw()
#         fig.canvas.flush_events()
#
#         time.sleep(0.1)
#



# Get path to retrieve data and save results

path = os.getcwd()
save_path = path + '\Results'
timepath = path + '\StegerTimesteps'
initialdata = save_path + '\Init_Steger_cond.txt'  # open initial conditions file
finaldata = save_path + '\Final_Steger_cond.txt' # Open final conditions file

Plot_Single(initialdata, "Initial Conditions")
Plot_Single(finaldata, "Final Conditions")
# Plot_TimeSeries(timepath, "Timestep")

files = glob.glob(timepath + '\*.txt')
i = 0 # Select how often to plot
numplot = 5  # Choose how many timesteps to plot
# plotcheck = int(len(files) / numplot)
plotcheck = 1

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
        plt.title("Timestep", fontsize=titlefont)

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
