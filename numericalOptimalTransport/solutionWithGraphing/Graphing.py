from matplotlib import pyplot as plt
from scipy.integrate import quad
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import rc
import imageio
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from scipy.integrate import quad
import numpy as np
from scipy.integrate import quad

def compareIntegrals(solution, f, g):
    '''
    Takes in the solution list, and the functions for f and g, and
    plots the integrals that they are comparing.
    '''
    xValues = np.linspace(-1,1,len(solution))
    estUx = [-1]
    dx = xValues[1] - xValues[0]
    for i in range(len(solution) - 1):
        if i != 0:
            u_x = (solution[i+1] - solution[i - 1]) / (2 * dx)
            estUx.append(u_x)
    estUx.append(1)

    gList = []
    fList = []
    for i in range(len(xValues)):
        fVal, err = quad(f,-1, xValues[i])
        gVal, err = quad(g,-1,estUx[i])
        fList.append(fVal)
        gList.append(gVal)
    diffList = []
    for i in range(len(fList)):
        diffList.append(abs(fList[i] - gList[i]))
    plt.plot(xValues, fList, label='Int F Values', linestyle=":", color="r")
    plt.plot(xValues, gList, label='Int G Values', linestyle=":", color="blue")
    plt.title('Comparing the Integrals of f and g')
    print(max(diffList))

def plotFirstDer(xValues, graphing):
    '''
    Plot the first derivative and return a list of the first derivatives
    '''
    dx = xValues[1] - xValues[0]
    fd = []
    for index in range(0, len(graphing)):
        estUx = [-1]
        for i in range(len(graphing[index]) - 1):
            if i != 0:
                u_x = (graphing[index][i+1] - graphing[index][i - 1]) / (2 * dx)
                estUx.append(u_x)
        estUx.append(1)
        fd.append(estUx)

    n_lines = len(fd)
    c = np.arange(0, n_lines)

    norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.rainbow)
    cmap.set_array([])

    fig, ax = plt.subplots(dpi = 100)
    for i in range(len(fd)):
        ax.plot(xValues, fd[i], c=cmap.to_rgba(i + 1))
    cbar = fig.colorbar(cmap, ticks=[0, int(len(fd)/2), len(fd)-1])
    ax.set(xlabel=r'X', ylabel=r'$\displaystyle u_{x}$' )
    cbar.set_label(r'Time Steps (1000)')
    plt.grid()
    plt.show()
    return fd

def plotSecDer(secondD, f, xValues, fd ):
    n_lines = len(secondD)
    c = np.arange(0, n_lines)

    norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.rainbow)
    cmap.set_array([])

    real = []
    ux = fd[len(fd) - 1]
    for i in range(len(xValues)):
        real.append(2 * f(xValues[i]))

    fig, ax = plt.subplots(dpi=100)
    for i in range(len(secondD)):
        ax.plot(xValues, secondD[i], c=cmap.to_rgba(i + 1))
    cbar = fig.colorbar(cmap, ticks=[0, int(len(secondD) / 2), len(secondD) - 1])
    ax.set(xlabel=r'X', ylabel=r'$\displaystyle u_{xx}$')
    cbar.set_label(r'Time Steps (1000)')
    plt.grid()
    plt.show()

def plot_for_offset(dxarray, norm, cmap, howmany, xValues, u_x):
    '''
    Used for the gif animation
    '''
    cmap.set_array([])
    fig, ax = plt.subplots(dpi=100)
    for i in range(0, howmany):
           ax.plot(xValues, dxarray[i], c=cmap.to_rgba(i + 1))
    ax.grid()
    ax.set(xlabel=r'X', ylabel=r'$\displaystyle u_x$',
           title=r'Estimated $\displaystyle u_x$ Over Time, $\displaystyle \frac{9}{20}x + \frac{1}{2} \rightarrow \frac{1}{4}(x+2) $')
    cbar = fig.colorbar(cmap, ticks=[0, int(len(u_x) / 2), len(u_x) - 1])
    cbar.set_label(r'Time Steps (10000)')

    # IMPORTANT ANIMATION CODE HERE
    # Used to keep the limits constant
    ax.set_ylim(-1.1, 1.1)

    # Used to return the plot as an image array
    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image

def gifPlot(u_x):
    n_lines = len(u_x)
    c = np.arange(0, n_lines)

    norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.rainbow)

    frames = []
    for i in range(len(u_x)):
        frame = plot_for_offset(u_x, norm, cmap, i)
        frames.append(frame)
        plt.close('all')

    kwargs_write = {'fps': 10.0, 'quantizer': 'nq'}
    imageio.mimsave('./nine20th.gif', frames, fps=10)