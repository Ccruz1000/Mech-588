# Import packages
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import time

# User defined packages

# Plot a single timestep
def Plot_Single(datafile, title, savetitle):
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
    plt.savefig(save_path + savetitle)
    # plt.show()
    plt.close()

# Plot each individual timestep
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

           #  plt.show()
            plt.close()
        i += 1

# Plot all 3 final conditions on the same graph
def Plot_Results(title, Roe, Steger, Lax):
    Roedata = np.loadtxt(Roe, dtype=float, skiprows=1, delimiter=' ')
    StegerData = np.loadtxt(Steger, dtype=float, skiprows=1, delimiter=' ')
    LaxData = np.loadtxt(Lax, dtype=float, skiprows=1, delimiter=' ')
    # Organize data
    x = Roedata[:, 0]  # Same for all 3
    rho = [Roedata[:, 1], StegerData[:, 1], LaxData[:, 1]]
    u = [Roedata[:, 1], StegerData[:, 1], LaxData[:, 1]]
    P = [Roedata[:, 2], StegerData[:, 2], LaxData[:, 2]]
    E = [Roedata[:, 3], StegerData[:, 3], LaxData[:, 3]]

    # Plot colors
    roecolor = 'g'
    stegercolor = 'b'
    laxcolor = 'r'

    # Plot condition data
    fig = plt.subplots()

    numplot = 4
    axisfont = 16
    titlefont = 18

    # Plot rho
    ax1 = plt.subplot(numplot, 1, 1)
    plt.plot(x, rho[0], c=roecolor, label="Roe")  # Plot Roe
    plt.plot(x, rho[1], c=stegercolor, label="Steger-Warming")  # Plot Steger
    plt.plot(x, rho[2], c=laxcolor, label="Lax-Friedrich")  # Plot LAX
    plt.ylabel('ρ\n', fontsize=axisfont)
    plt.grid(True)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)
    plt.xlim(x[0], x[-1])
    plt.legend(bbox_to_anchor=(0.82, 0.59))

    plt.title(title, fontsize=titlefont)

    # Plot u
    ax2 = plt.subplot(numplot, 1, 2)
    plt.plot(x, u[0], c=roecolor, label="Roe")  # Plot Roe
    plt.plot(x, u[1], c=stegercolor, label="Steger-Warming")  # Plot Steger
    plt.plot(x, u[2], c=laxcolor, label="Lax-Friedrich")  # Plot LAX
    plt.ylabel('u', fontsize=axisfont)
    plt.grid(True)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)
    plt.xlim(x[0], x[-1])

    # Plot P
    ax3 = plt.subplot(numplot, 1, 3)
    plt.plot(x, P[0], c=roecolor, label="Roe")  # Plot Roe
    plt.plot(x, P[1], c=stegercolor, label="Steger-Warming")  # Plot Steger
    plt.plot(x, P[2], c=laxcolor, label="Lax-Friedrich")  # Plot LAX
    plt.ylabel('P\n', fontsize=axisfont)
    plt.grid(True)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)
    plt.xlim(x[0], x[-1])

    # Plot E
    ax4 = plt.subplot(numplot, 1, 4)
    plt.plot(x, E[0], c=roecolor, label="Roe")  # Plot Roe
    plt.plot(x, E[1], c=stegercolor, label="Steger-Warming")  # Plot Steger
    plt.plot(x, E[2], c=laxcolor, label="Lax-Friedrich")  # Plot LAX
    plt.ylabel('E\n', fontsize=axisfont)
    plt.grid(True)
    plt.xlabel('x', fontsize=axisfont)
    plt.xlim(x[0], x[-1])

    # Adjust plots
    plt.subplots_adjust(left=0.15)
    plt.subplots_adjust(bottom=0.12)
    plt.subplots_adjust(top=0.90)

    # Save and show plots
    plt.savefig(save_path + "\Final_Results.png")
    plt.show()
    plt.close()


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

# Define Lax-Friedrich filenames
Laxtimepath = path + '\LaxTimesteps'
initialLaxdata = save_path + '\Init_Lax_cond.txt'  # open initial conditions file
finalLaxdata = save_path + '\Final_Lax_cond.txt' # Open final conditions file

# # Plot Steger
Plot_Single(initialStegerdata, "Steger Initial Conditions", "/Init_Steger_cond.png")
Plot_Single(finalStegerdata, "Steger Final Conditions", "/Init_Steger_cond.png")
# Plot_Timesteps(Stegertimepath, "Steger Timestep ")

# Plot Roe
Plot_Single(initialRoedata, "Roe Initial Conditions", "/Init_Roe_cond.png")
Plot_Single(finalRoedata, "Roe Final Conditions", "/Final_Roe_cond.png")
# Plot_Timesteps(Roetimepath, "Roe Timestep ")

# Plot Lax
Plot_Single(initialLaxdata, "Strong Shock Wave Initial Conditions", "/Init_Lax_cond.png")
Plot_Single(finalLaxdata, "Lax Final Conditions", "/Final_Lax_cond.png")
# Plot_Timesteps(Roetimepath, "Roe Timestep ")

# Plot all one same graph
Plot_Results("Strong Shock Wave Test", finalRoedata, finalStegerdata, finalLaxdata)





