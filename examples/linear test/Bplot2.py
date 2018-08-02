

'''
    Provides function that will plot the reconstructed magnetic field caluctlated from
    the alogrithm in recconstruct.py
    '''

import math
import sys
import os.path
import deflect
import matplotlib as mpl
mpl.use('Agg') # Headless plotting (avoids python-tk GUI requirement)
from matplotlib.ticker import FixedLocator
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import re

# log magnitude
def magnetic_field(Bxpath, Bypath):
    '''
        The goal is to calculate the magnetic field per bin
        
        Parameters
        ----------
        Br (2D array of (x,y)): B Field per (x,y)
        
        Returns
        -------
        BrMag (2D array): Log reconstructed B Field per bin
        '''
    magB = np.log10(np.sqrt(Bxpath**2 + Bypath**2))
    return BrMag
font = {'family': 'serif',
    'color':  'black',
        'weight': 'normal',
        'size': 32,
    }
# Load the plasma coordinates.
X = np.loadtxt('solve/plasma_x.txt', delimiter=',')
Y = np.loadtxt('solve/plasma_y.txt', delimiter=',')

# Checkpoint interval.
interval = 10000

# Plasma parameters.
ri = 1.3 # Distance from proton source to plasma.
li = 0.2 # Distance across plasma.
rs = 20 # Distance from plasma to screen.
v = 5.24e9 # Velocity of protons.

stretch = 13.1 / 10.2
#plt.rc('text', usetex=True)

#################################################

# Plot all of the solve steps for this problem, saving them as
# images/magBpath######.png.
for i in range(400):
    print ("Step " + str(i*interval))
    phix = np.loadtxt('solve/phix'+str(interval*i)+'.txt', delimiter=',')
    phiy = np.loadtxt('solve/phiy'+str(interval*i)+'.txt', delimiter=',')
    phin = np.loadtxt('solve/phin'+str(interval*i)+'.txt', delimiter=',')
    wBx, wBy = deflect.reconstruct(ri, li, rs, v, X, Y, phix, phiy)
    Bxpath, Bypath = deflect.magpath(wBx, wBy)
    

    fig = plt.figure()
    fig.set_figwidth(6.0 * stretch)
    fig.set_figheight(6.0)
    ax = fig.add_subplot(1, 1, 1)
    BMag = magnetic_field(Bxpath, Bypath)
    vmin = BMag.min()
    vmax = BMag.max()
    norm = mpl.colors.Normalize(vmin= vmin,vmax= vmax)
    
    fig, ax = plt.subplots()
    strm = ax.streamplot(X[:, 0], Y[0, :], Bxpath, Bypath, color=BMag,
                         linewidth=2, cmap=cm.RdYlGn, density=2.0, arrowsize=2.0, norm = norm)
    fig.colorbar(strm.lines)
    
    #################################################
                         
    
    xmin = round(min(X[:, 0]), 1)
    xmax = round(max(X[:, 0]), 1)
    ymin = round(min(Y[0, :]), 1)
    ymax = round(max(Y[0, :]), 1)
    ax.set_xlim(int(xmin) - 0.5, int(xmax) + 0.5)
    ax.set_ylim(int(ymin) - 0.5, int(ymax) + 0.5)
    ax.set(title= typeName + ": Log " + r" $B_\perp$ Projection (G cm)",
           ylabel=r"Y (cm)", xlabel=r"X (cm)")
    ax.tick_params(labelsize='large')
                                                            
    fig.savefig("B_" + str(i*interval) + ".png", format='png')

