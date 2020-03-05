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
        data_dict['Lambda'].append(int(re.findall(numbers,l)[0]))
        data_dict['Max_act'].append(int(re.findall(numbers,Max)[0][:-1]))
    
    df = pd.DataFrame(data_dict)
    df.to_csv('../results/CPM_nofrc1.csv')

if __name__ == "__main__":
    files = glob.glob('../results/first_cpm/*.txt')
    files.sort()
    #build_df(files)
    # test
    f = files[40]
    track = np.loadtxt(f)
    print(get_motility(track,plot = True))