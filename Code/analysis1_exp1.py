# create 3d bar plot of all sigma,omega values 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import glob
import pickle
import numpy as np
import re
import seaborn as sns
from matplotlib.patches import Ellipse
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
sns.set_style('darkgrid')

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
    colors = ['c','g','y','0.2','b','chocolate']
    sig_files = glob.glob('../data/exp1/sigma*')
    om_files = glob.glob('../data/exp1/omega*')
    sig_files.sort()
    om_files.sort()
    all_sigmas = np.array([])
    all_omegas = np.array([])
    plt.figure(figsize=(15,15))
    for i in range(len(sig_files)):
        # extract values from dict and then reshape into flat list : 
        omega = np.array(list(pickle.load(open(om_files[i], 'rb')).values())).flatten()
        sigma = np.array(list(pickle.load(open(sig_files[i], 'rb')).values())).flatten()
        name = re.split(("sigma"),sig_files[i])[-1][:-4]
        plt.scatter(sigma,omega,c = colors[i], s = 30,alpha = .3, label = name)
        all_sigmas = np.append(all_sigmas,sigma)
        all_omegas = np.append(all_omegas,omega)

    # plot Novkovic with circle around
    plt.scatter(6.7,-.27, c = 'r' , s = 80, alpha = .5, label = 'Novkovic')
    x_nov,y_nov = 6.7,-.27
    #print(all_omegas.shape,all_sigmas.shape)

    xr =  .25 * np.std(all_sigmas) 
    yr =  .25 * np.std(all_omegas)
    print(4*xr,4*yr)
    # PLOT SQUARE :
    #plt.plot([x_nov - xr,x_nov+xr],[y_nov-yr,y_nov-yr],c = 'r')
    #plt.plot([x_nov - xr,x_nov+xr],[y_nov+yr,y_nov+yr],c = 'r')
    #plt.plot([x_nov - xr,x_nov - xr],[y_nov+yr,y_nov-yr],c = 'r')
    #plt.plot([x_nov + xr,x_nov+xr],[y_nov-yr,y_nov+yr],c = 'r')
    # legend en sstuff
    plt.xlim((x_nov-xr,x_nov+xr))
    plt.ylim(y_nov-yr,y_nov+yr)
    plt.ylabel('$\omega$')
    plt.xlabel('$\sigma$')
    plt.legend()
    #plt.show()


def plot_frc_like():
    colors = ['c','g','y','0.2','b','chocolate']
    sig_files = glob.glob('../data/exp1/sigma*')
    om_files = glob.glob('../data/exp1/omega*')
    sig_files.sort()
    om_files.sort()
    plt.figure(figsize=(15,15))

    for i in range(len(sig_files)):
        # extract values from dict and then reshape into flat list : 
        omega = pickle.load(open(om_files[i], 'rb'))
        sigma = pickle.load(open(sig_files[i], 'rb'))
        for key in sigma.keys():
            #print(key)
            # check for keys with 176N and .25 V/E
            if '176N_0.25V/E' in key:
                sig = np.array(sigma[key]).flatten()
                om = np.array(omega[key]).flatten()
                plt.scatter(sig,om,c = colors[i], s = 40,alpha = .3)
        # plot nothing with same collor for legend
        name = re.split(("sigma"),sig_files[i])[-1][:-4]
        plt.scatter(sig[0],om[0],c = colors[i], s = 40,alpha = .3, label = name)

    # plot Novkovic with circle around
    plt.scatter(6.7,-.27, c = 'r' , s = 80, alpha = .5, label = 'Novkovic')

    plt.ylabel('$\omega$')
    plt.xlabel('$\sigma$')
    plt.title("Omega and Sigma values for graphs with N = 176 and E +- 685")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    #plot3d()
    #scatter_plot()
    plot_per_type((0,150),(-1.2,1.2))
    #plot_averages()
    #plot_frc_like()



# def plot3d():
#     # load data : 
#     # code adapted from https://matplotlib.org/3.1.0/gallery/mplot3d/hist3d.html
#     omegas,sigmas = load_data()
#     print(len(omegas))
#     print(len(sigmas))
#     fig = plt.figure(figsize = (15,15))
#     ax = fig.add_subplot(projection = '3d')
#     # create bins for histogram (frequency counts)
#     hist,xedges,yedges = np.histogram2d(omegas,sigmas, bins = 10, range = [[-1,1],[1,200]])

#     # position for bars : 
#     xpos,ypos = np.meshgrid(xedges[:-1] + .25, yedges[:-1] + .25, indexing = 'ij')
#     xpos = xpos.ravel()     # flatten arrays
#     ypos = ypos.ravel() 
#     zpos = 0     

#     # Construct arrays with the dimensions for the 16 bars.
#     dx = dy = 0.5 * np.ones_like(zpos)
#     dz = hist.ravel()
#     ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')
#     plt.show()



# def plot_averages():
#      # plot average and variance per graph type and parameter set :
#     colors = ['c','g','y','0.2','b']
#     sig_files = glob.glob('../data/exp1/sigma*')
#     om_files = glob.glob('../data/exp1/omega*')
#     sig_files.sort()
#     om_files.sort()
#     plt.figure(figsize=(15,15))
#     # save all values to use std for circle : 
#     # all_sigs = []
#     # all_oms = []
#     count = 0
#     for i in range(len(sig_files)):
#         # extract values from dict and then reshape into flat list : 
#         omega = pickle.load(open(om_files[i], 'rb'))
#         sigma = pickle.load(open(sig_files[i], 'rb'))
#         name = re.split(("sigma"),sig_files[i])[-1][:-4]
#         # compute averages per parameter set :
#         for key,svalue in sigma.items():
#             ovalue = omega[key]
#             s_average ,o_average = np.average(svalue),np.average(ovalue)
#             s_var , o_var = np.var(svalue) , np.var(ovalue)
#             # plot meann : 
#             plt.scatter(s_average,o_average,c = colors[i], s = 20,alpha = .8)
#             circle1 = Ellipse((s_average,o_average),s_var,o_var,color= colors[i],alpha = .2)
#             plt.gcf().gca().add_artist(circle1)
#             count += 1
#     print('total number of points ; ',count)
#         # all_oms.append(omega)
#         # all_sigs.append(sigma)

#     # plot Novkovic with circle around
#     plt.scatter(6.7,-.27, c = 'r' , s = 80, alpha = .5, label = 'Novkovic')
#     # circle : 
#     # theta = np.linspace(0,2*np.pi,100)
#     # xr = xwindow[-1] - xwindow[0]
#     # yr = ywindow[-1] - ywindow[0]
#     # #xr = np.std(all_sigs) 
#     # #yr = np.std(all_oms)

#     # x =  .05*xr * np.cos(theta) + 6.7 #1/ywindow[-1] *
#     # y = .05*yr * np.sin(theta) - .27
#     # plt.plot(x,y, c = 'r')
#     plt.xlim((0,50))
#     # plt.ylim(ywindow)
#     plt.ylabel('$\omega$')
#     plt.xlabel('$\sigma$')
#     plt.title("Omega and Sigma averages with variance for all created graphs")
#     plt.legend()
#     plt.show()
