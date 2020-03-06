import numpy as np 
import glob
from CPM_helpers1 import *
from MSD1cell import *
import pandas as pd
import re 
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt


def get_motility(track,plot = False):
    # cell track is list of 3d coordinates.
    if plot:
        plot_celltrack(track)
    #vol = scanned_volume(volume_track)
    displ, angles = analyse_track(track)
    # MSD
    delta_t,MSD = compute_MSD(displ,plot = plot)
    # autocorr : 
    _,_,_,slope,p = compute_AutoCorr_WRONG(angles,plot = plot)
    #t,new_angles = compute_autocorrelaton(displ)
    
    popt = fit_Motilty(delta_t,MSD)

    return popt[0],popt[1],slope,p

def roseplot(track):
    # get list of angles : 
    _,angles = analyse_track(track)

    n_numbers = 100
    bins_number = 16  # the [0, 360) interval will be subdivided into this
    # number of equal bins
    bins = np.linspace(0.0, 2 * np.pi, bins_number + 1)
    #angles = 2 * np.pi * np.random.rand(n_numbers)
    n, _, _ = plt.hist(angles, bins)

    #plt.clf()
    width = 2 * np.pi / bins_number
    ax = plt.subplot(1, 1, 1, projection='polar')
    bars = ax.bar(bins[:bins_number], n, width=width, bottom=0.0)
    for bar in bars:
        bar.set_alpha(0.5)
    plt.show()

def new_auto(cell_track): 
    averages = []
    for dt in range(0,300):
        # angles per dt: 
        cosines = []
        for i in range(len(cell_track) - 2 - dt):
            point1 = cell_track[i]
            point2 = cell_track[i + 1]
            point3 = cell_track[i + dt]
            point4 = cell_track[i + dt + 1]
            if list(point2) == list(point4):
                pass
            v1 = point2 - point1
            v2 = point4 - point3
            #cos = np.clip(np.dot(v1,v2)/(norm(v1) * norm(v2)),-1,1)
            cos = np.dot(v1,v2)/(norm(v1) * norm(v2))
            cosines.append(cos)
        if len(cosines) > 100: #minimum 100 mcs with actual displacement
            averages.append(np.average(cosines))
    return averages


def build_df(files):
    data_dict = {'Motility':[],'Persistance':[],'AutoSlope':[],'P-val':[]
    ,'Lambda':[],'Max_act':[]}
    numbers = '[-+]?\d*\.\d+|\d+'
    for f in files:
        track = np.loadtxt(f)
        m,p,slope,p_val = get_motility(track)
        data_dict['Motility'].append(m)
        data_dict['Persistance'].append(p)
        data_dict['AutoSlope'].append(slope)
        data_dict['P-val'].append(p_val)
        l,Max = f.split('MAX')
        data_dict['Lambda'].append(int(re.findall(numbers,l)[1]))
        data_dict['Max_act'].append(int(re.findall(numbers,Max)[0][:-1]))
    
    df = pd.DataFrame(data_dict)
    df.to_csv('../results/CPM_nofrc2.csv')

if __name__ == "__main__":
    files = glob.glob('../results/nofrc2/*.txt')
    files.sort()
    # build_df(files)
    # test
    f = files[18]
    track = np.loadtxt(f)
    roseplot(track)
    auto_cor = new_auto(track)
    plt.clf()
    plt.plot(auto_cor)
    plt.show()
    #print(get_motility(track,plot = True))