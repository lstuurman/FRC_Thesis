import pandas as pd
import numpy as np
import glob
import re
from numpy.linalg import norm
import sys
sys.path.insert(0,'../../CPM_code')
sys.path.insert(0,'../../Code')
#from MSD1cell import *
from OrderNpersist import to_vecs,Order_tracks2 ,order_radius,Persist_tracks 
from analyse_tracks import new_auto

def handle_boundaries(cell_track,pr = False):
    # look for boundary crossings in any
    # of the coordinates
    cell_track2 = cell_track.copy()
    for i in range(len(cell_track) - 1):
        dif = np.subtract(cell_track[i],cell_track[i+1])
        for j,coordinate in enumerate(dif):
            if coordinate > 50:
                # went over boundary from 256 -> 
                cell_track2[:i + 1,j] -= 50

            elif coordinate < -50:
                cell_track2[:i + 1,j] += 50
    return cell_track2

def ANLSIS():
    # create data frame with fields : 
    #|| act/prfdir | density | speed | perisitance | order ||
    # data :
    files = glob.glob('../../data/cpmjs/density/*')
    num_ptrn = '[-+]?\d*\.\d+|\d+'
    files.sort()
    # loop over files and compute average speed | peristance | order
    rows = []
    print(files)
    for f in files:
        print(f)
        data = pd.read_csv(f, sep="\t", names=['time','id','ctype', 'x', 'y', 'z'])   
        # create list of arrays of cell tracks 
        ids = data.id.unique()
        tracks = [data[data.id == ID].to_numpy()[:,-3:] for ID in ids]
        # fix boundaries : 
        tracks = [handle_boundaries(t) for t in tracks]
        print(len(ids),len(tracks))
        print(tracks[0].shape)
        
        # speed and order: vec_tracks = np.array([[to_vecs(t) for t in track] for track in tracks])
        vec_tracks = np.array([to_vecs(t) for t in tracks])
        # for i,v in enumerate(vec_tracks[0]):
        #     print(tracks[0][i],tracks[0][i+1])
        #     print(v,norm(v))
        #     exit()
        speed = np.mean([[norm(v) for v in vec_track] for vec_track in vec_tracks])
        order = Order_tracks2(vec_tracks)
        lcl_order = order_radius(tracks,20)
        # persistance : 
        autocors = [new_auto(t) for t in tracks]
        persist = Persist_tracks(autocors)
        print(speed,np.average(persist),np.average(order),np.average(lcl_order))
        # params : 
        sim_type = f.split('_')[0].split('/')[-1]
        density = density = float(re.findall(num_ptrn,f.split('_')[1])[0])
        rows.append([sim_type,density,speed,np.average(persist),np.average(order),np.average(lcl_order)])
    df = pd.DataFrame(data= rows, columns = ['Model','Density','speed','peristance','order','lcl_order'])
    df.to_csv('Dens3.csv')

ANLSIS()
