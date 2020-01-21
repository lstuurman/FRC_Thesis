# create 3d bar plot of all sigma,omega values 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import glob
import pickle
import numpy as np
import re

def load_data():
    omegas = np.array([])
    sigmas = np.array([])
    sig_files = glob.glob('../data/exp1/sigma*')
    om_files = glob.glob('../data/exp1/omega*')
    sig_files.sort()
    om_files.sort()
    for i in range(len(om_files)):
        # extract values from dict and then reshape into flat list : 
        omega = np.array(list(pickle.load(open(om_files[i], 'rb')).values()))
        sigma = np.array(list(pickle.load(open(sig_files[i], 'rb')).values()))
        omegas = np.append(omegas,omega.flatten())
        sigmas = np.append(sigmas,sigma.flatten())
    return omegas,sigmas

def scatter_plot():
    omegas,sigmas = load_data()
    plt.figure(figsize=(15,15))
    plt.scatter(sigmas,omegas, c = 'm',alpha = .5,s = 40)
    plt.scatter(6.7,-.27, c = 'g' , s = 80)
    plt.show()

def plot_per_type(xwindow,ywindow):
    colors = ['c','g','y','0.75']
    sig_files = glob.glob('../data/exp1/sigma*')
    om_files = glob.glob('../data/exp1/omega*')
    sig_files.sort()
    om_files.sort()
    plt.figure(figsize=(15,15))
    for i in range(len(sig_files)):
        # extract values from dict and then reshape into flat list : 
        omega = np.array(list(pickle.load(open(om_files[i], 'rb')).values()))
        sigma = np.array(list(pickle.load(open(sig_files[i], 'rb')).values()))
        name = re.split(("sigma"),sig_files[i])[-1][:-4]
        plt.scatter(sigma.flatten(),omega.flatten(),c = colors[i], s = 40,alpha = .5, label = name)
    plt.scatter(6.7,-.27, c = 'r' , s = 80)
    plt.xlim(xwindow)
    plt.ylim(ywindow)
    plt.ylabel('$\omega$')
    plt.xlabel('$\sigma$')
    plt.title("Omega and Sigma values for differen types of generated random graphs")
    plt.legend()
    plt.show()


def plot3d():
    # load data : 
    # code adapted from https://matplotlib.org/3.1.0/gallery/mplot3d/hist3d.html
    omegas,sigmas = load_data()
    print(len(omegas))
    print(len(sigmas))
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(projection = '3d')
    # create bins for histogram (frequency counts)
    hist,xedges,yedges = np.histogram2d(omegas,sigmas, bins = 10, range = [[-1,1],[1,200]])

    # position for bars : 
    xpos,ypos = np.meshgrid(xedges[:-1] + .25, yedges[:-1] + .25, indexing = 'ij')
    xpos = xpos.ravel()     # flatten arrays
    ypos = ypos.ravel() 
    zpos = 0     

    # Construct arrays with the dimensions for the 16 bars.
    dx = dy = 0.5 * np.ones_like(zpos)
    dz = hist.ravel()
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')
    plt.show()

if __name__ == "__main__":
    #plot3d()
    #scatter_plot()
    plot_per_type((1,20),(-.5,.5))
# ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')

# # This import registers the 3D projection, but is otherwise unused.
# from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

# import matplotlib.pyplot as plt
# import numpy as np

# # Fixing random state for reproducibility
# np.random.seed(19680801)


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# x, y = np.random.rand(2, 100) * 4
# hist, xedges, yedges = np.histogram2d(x, y, bins=4, range=[[0, 4], [0, 4]])

# # Construct arrays for the anchor positions of the 16 bars.
# xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25, indexing="ij")
# xpos = xpos.ravel()
# ypos = ypos.ravel()
# zpos = 0

# # Construct arrays with the dimensions for the 16 bars.
# dx = dy = 0.5 * np.ones_like(zpos)
# dz = hist.ravel()

# ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')

# plt.show()