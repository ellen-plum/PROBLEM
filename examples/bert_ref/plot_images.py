import numpy as np
import matplotlib.pyplot as plt
import problem.deflect
"""
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
"""
"""
# Plot all of the solve steps for this problem, saving them as
# images/magBpath######.png.
for i in range(400):
    print ("Step " + str(i*interval))
    phix = np.loadtxt('solve/phix'+str(interval*i)+'.txt', delimiter=',')
    phiy = np.loadtxt('solve/phiy'+str(interval*i)+'.txt', delimiter=',')
    phin = np.loadtxt('solve/phin'+str(interval*i)+'.txt', delimiter=',')
    wBx, wBy = problem.deflect.reconstruct(ri, li, rs, v, X, Y, phix, phiy)
    Bxpath, Bypath = problem.deflect.magpath(wBx, wBy)
    
    fig, ax = plt.subplots()
    
    magB = np.sqrt(Bxpath**2 + Bypath**2) # make the axis adjust according to magnetic field.......
    maxB = np.amax(magB)
    ticks = np.linspace(0,maxB,12)
    cs = ax.pcolorfast(X[:,0], Y[0,:], magB, vmin=0, vmax=maxB)
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('y (cm)')
    cb = plt.colorbar(cs, ax=ax,ticks=ticks)

    ax.set_title('Calculated Magnitude Path-Int Magnetic Field (G*cm) Step: '
                 '{}'.format(interval*i))

    fig.savefig('images/magBpath'+str(interval*i).zfill(6)+'.png')
    print('Step: {}'.format(i))

"""
import math
import sys
import os.path
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
    magB = np.sqrt(Bxpath**2 + Bypath**2)/10**3 # convert from Gcm to Tmm
    return magB
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
ri = .9 # Distance from proton source to plasma.
li = 0.2 # Distance across plasma.
rs = 15.3 # Distance from plasma to screen.
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
    wBx, wBy = problem.deflect.reconstruct(ri, li, rs, v, X, Y, phix, phiy)
    Bxpath, Bypath = problem.deflect.magpath(wBx, wBy)
    
    
    fig = plt.figure()
    fig.set_figwidth(6.0 * stretch)
    fig.set_figheight(6.0)
    ax = fig.add_subplot(1, 1, 1)
    BMag = magnetic_field(Bxpath, Bypath)
    vmin = BMag.min()
    vmax = BMag.max()
    norm = mpl.colors.Normalize(vmin= vmin,vmax= vmax)
    
    fig, ax = plt.subplots() # streamplot takes in (y, x) :(
    strm = ax.streamplot(Y [0, :], X[:, 0], Bypath, Bxpath,color=BMag,
                         linewidth=2, cmap=cm.RdYlGn, density=2.0, arrowsize=2.0, norm = norm)
    fig.colorbar(strm.lines)
    #################################################
    xmin = min(X[:, 0])
    xmax = max(X[:, 0])
    ymin = min(Y[0, :])
    ymax = max(Y[0, :])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set(title= " Bert (3MeV) " + r" $B_\perp$ Projection (T mm)",
           ylabel=r"Y (cm)", xlabel=r"X (cm)")
    ax.tick_params(labelsize='large')
                         
    fig.savefig("B_" + str(i*interval) + ".png", format='png')


