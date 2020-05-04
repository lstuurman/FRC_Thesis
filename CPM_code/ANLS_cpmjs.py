import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import re

import sys
sys.path.insert(0,'../CPM_code')
from MSD1cell import *
from OrderNpersist import * 

def ANLSIS():
    # create data frame with fields : 
    #|| act/prfdir | density | speed | perisitance | order ||
    # data :
    files = glob.glob('/home/lau/Desktop/Thesis Stuff/cpmjsDensity/*')
    num_ptrn = '[-+]?\d*\.\d+|\d+'
    files.sort()
    # loop over files and compute average speed | peristance | order
    rows = []
    for f in files:
        print(f)
        data = pd.read_csv(f, sep="\t", names=['time','id','ctype', 'x', 'y', 'z'])   
        # create list of arrays of cell tracks 
        ids = data.id.unique()
        tracks = [data[data.id == ID].to_numpy()[:,-3:] for ID in ids]
        print(len(ids),len(tracks))
        print(tracks[0].shape)
        
        # speed and order: 
        vec_tracks = np.array([to_vecs(t) for t in tracks])
        speed = np.mean([norm(v) for v in vec_tracks])
        order = Order_tracks(vec_tracks)
        # persistance : 
        autocors = [new_auto(t) for t in tracks]
        persist = Persist_tracks(autocors)
        print(speed,persist,order)
        # params : 
        sim_type = f.split('_')[0].split('/')[-1]
        density = density = float(re.findall(num_ptrn,f.split('_')[1])[0])
        rows.append([sim_type,density,speed,persist,order])
        
    df = pd.DataFrame(data= rows, columns = ['Model','Density','speed','peristance','order'])
    df.to_csv('results_dens.csv')

ANLSIS()