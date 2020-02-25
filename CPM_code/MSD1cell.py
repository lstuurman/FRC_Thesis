import numpy as np 
import cpm
import matplotlib.pyplot as plt
from CPM_helpers1 import *
import math
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import time

def compute_MSD(displacement):
    # range of sliding windows from 1 to steps - 100
    delta_t = np.arange(1,len(displacement) - 100)
    MSD = np.zeros(len(delta_t))
    for i,dt in enumerate(delta_t):
        # convolve to get displacements : 
        moving_av = np.convolve(displacement, np.ones(dt), 'valid')
        MSD[i] = np.average(np.square(moving_av))
        #MSD[i] = np.sqrt(np.average(np.square(moving_av)))

    # plot : 
    plt.plot(delta_t,MSD)
    plt.title('MSD of lonely cell')
    plt.ylabel('MSD')
    plt.xlabel('$\delta t$')
    plt.savefig('testdat/MSD2.png')
    plt.show()
    
    return delta_t, MSD

def compute_autocorrelaton(displacement):
    # compute angles between vectors with distance dt
    delta_t = np.arange(1,len(displacement) - 100)
    vecs = [displacement[i+1] - displacement[i] for i in range(len(displacement) - 1)]
    angles = []
    for dt in delta_t:
        phis = [np.dot(vecs[i],vecs[i+dt])/(norm(vecs[i]) * norm(vecs[i+dt])) for i in range(len(vecs) - dt)]
        angles.append(np.average(phis))

    plt.plot(delta_t,angles)
    plt.title('Correlation of angles in cell track')
    plt.ylabel('Angle phi')
    plt.xlabel('$\delta t$')
    #plt.legend()
    plt.savefig('testdat/autocorr2.png')
    plt.show()
    return delta_t,angles


def compute_AutoCorr_WRONG(angles):
    r_vals = []
    p_vals = []
    cos_cor = []
    #print(angles)
    angles = [x for x in angles if np.isfinite(x)]
    delta_t = np.arange(1,len(angles) - 100)
    #angles = angles[np.isfinite(angles)]
    #print(len(angles))
    for dt in delta_t:
        corrs = pearsonr(angles[:len(angles) - dt], angles[dt:])
        corrs2= [np.dot(angles[i],angles[i + dt]) for i in range(len(angles) - dt)]
        r_vals.append(corrs[0])
        p_vals.append(corrs[1])
        cos_cor.append(np.average(corrs2))
        #print(corrs)
    #print(r_vals)
    #print(p_vals)
    plt.plot(delta_t,r_vals,label = "R value")
    plt.plot(delta_t,p_vals,'--',label = "p-value")
    #plt.plot(delta_t,cos_cor, label = 'Cosine Angles')
    plt.title('Autocorrelation of angles in cell track')
    plt.ylabel('Pearson r')
    plt.xlabel('$\delta t$')
    plt.legend()
    plt.savefig('testdat/autocorr2.png')
    plt.show()

    return delta_t,r_vals, p_vals

def furth(t,D,P):
    return 6*D*(t - P*(1 - math.e**(-t/P)))


def fit_Motilty(delta_t,MSD):
    popt,_ = curve_fit(furth,delta_t,MSD)
    print("Motility coefficient : ",popt[0])
    print("Persistance parameter : ",popt[1])
    return popt

def MSD_lonelycell():
    # setup : 
    begin_pos = [128,128,128]
    begin_pos = np.random.randint(0,256,size=3)
    simulation = single_cell_setup1(256,begin_pos,1000,5000,FRC=False)
    # run : 
    volume_track,cell_track = run_sim_1cell(simulation,2000)
    cell_track = handle_boundaries(cell_track,pr = False)
    plot_celltrack(cell_track)
    print('Percentage scanned volume : ', scanned_volume(volume_track))
    displ, angles = analyse_track(cell_track)
    # MSD
    delta_t,MSD = compute_MSD(displ)
    # autocorr : 
    t,auto_corr,pvalues = compute_AutoCorr_WRONG(angles)
    t,new_angles = compute_autocorrelaton(displ)
    
    popt = fit_Motilty(delta_t,MSD)


if __name__ == "__main__":
    t1 = time.time()
    MSD_lonelycell()
    print('computing time : ', time.time() - t1)

