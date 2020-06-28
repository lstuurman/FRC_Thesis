import pandas as pd
import numpy as np
import glob
from numpy.linalg import norm
import re
from scipy.spatial.distance import euclidean
import sys
import matplotlib.pyplot as plt
sys.path.insert(0,'../')
from OrderNpersist import Persist_tracks,to_vecs
from analyse_tracks import new_auto,auto_cor_pooled
from ANLS_DENS import OrderAt_T,local_order
from random import sample
import time
from itertools import product

def handle_boundaries(cell_track,pr = False):
    # look for boundary crossings in any
    # of the coordinates
    cell_track2 = cell_track.copy()
    for i in range(len(cell_track) - 1):
        dif = np.subtract(cell_track[i],cell_track[i+1])
        for j,coordinate in enumerate(dif):
            if coordinate > 32:
                # went over boundary from 64 -> 0
                cell_track2[:i + 1,j] -= 64

            elif coordinate < -32:
                # form 0 -> 64

                cell_track2[:i + 1,j] += 64
    return cell_track2

def wrap(deltas, dimension):
    deltas[deltas < -(dimension//2)] = deltas[deltas < -(dimension//2)] + dimension
    deltas[deltas > (dimension//2)] = deltas[deltas > (dimension//2)] - dimension
    return deltas

def normalize(deltas):
    return deltas / np.sqrt(np.sum(deltas ** 2, 1))[:, None]

def autocor(arr):
    dimension = 64
    #arr = np.array(df[["x","y","z"]])
    deltas = arr[1:] - arr[:-1]
    deltas = wrap(deltas, dimension)
    deltas = normalize(deltas)
    dot_for_dt = []
    for delta_x in range(100):
        dots = np.sum(deltas[delta_x:] * deltas[:len(deltas)-delta_x], 1)
        dot_for_dt.append(np.mean(dots))
    return dot_for_dt


def build_csv(path):
    """ Path is folder containing all celltracks of all max-lambda combinations.
    """
    # find unique param pairs in folder :
    num_ptrn = '[-+]?\d*\.\d+|\d+'
    all_files = glob.glob(path)
    param_strings = [f.split('_')[4] for f in all_files]
    params1 = set([tuple(re.findall(num_ptrn,s)) for s in param_strings])
    iters = [str(i) for i in range(2)]
    params2 = list(product(list(params1),iters))
    #print(params2)
    params2 = [(p[0][0],p[0][1],p[1]) for p in params2 if type(p[0]) == tuple] + list(params1)
    #print(params2)
    rows = []
    deviation_rows = []
    for i,prms in enumerate(params2):
        print(prms)
        t1 = time.time()
        files = [f for f in all_files if '_' + prms[0] + 'M' in f and 'X' + prms[1] + '_' in f]
        if len(prms) == 3:
            files = [f for f in files if 'iter'+prms[2] in f]
            itr = int(prms[2])
        else:
            files = [f for f in files if 'iter' not in f]
            itr = 2
        #files = [f for f in files if int(re.findall(num_ptrn,f)[-1]) > 1]
        files.sort(reverse = True,key = lambda x: int(re.findall(num_ptrn,x)[-1]))
        # only take first 100 cells for faster analysis :
        #files = files[:100]
        # extract tracks from files :
        tracks = [np.loadtxt(f) for f in files]
        tracks = [handle_boundaries(t) for t in tracks]
        print('number of cells for paramset : ',len(tracks))

        vec_tracks = np.array([to_vecs(t) for t in tracks])
        speeds = [[norm(v) for v in vec_track] for vec_track in vec_tracks]
        speeds = [item for sublist in speeds for item in sublist]
        speed = np.mean(speeds)
        std_speed = np.std(speeds)
        print('Speed calculated')
        #peristance : 
        # pooled persistance : 
        t2 = time.time()
        dts = 100
        #deltas = arr[1:] - arr[:-1]
        #deltas = wrap(deltas, dimension)
        #deltas = normalize(deltas)
        #dot_for_dt = []
        #for delta_x in range(200):
        #    dots = np.sum(deltas[delta_x:] * deltas[:len(deltas)-delta_x], 1)
        #    dot_for_dt.append(np.mean(dots))
        #dots_for_dts = [[] for i in range(dts)]
        #pooled_dts = [auto_cor_pooled(t,dots_for_dts,dts) for t in tracks]
        #averages2 = [np.mean(i) for i in dots_for_dts]
        try:
            averages = [autocor(t) for t in tracks]
        except:
            continue
        t3 = time.time()
        #pers2 = Persist_tracks([averages2])
        pooled_pers = Persist_tracks(averages)
        t4 = time.time()
        #print(np.mean(pooled_pers),pers2)
        #if len(pooled_pers) == 0:
        #    pooled_pers = np.nan
        #else:
            #pooled_pers = pooled_pers[0]

        rows.append([prms[0],prms[1],itr,speed,std_speed,np.mean(pooled_pers),np.std(pooled_pers)])
        print(rows[-1])
        print('computing time :',time.time() - t1)
        print('dots per dts time : ',t3 - t2)
        print('half time :',t4 - t3)
    df1 = pd.DataFrame(data = rows,
        columns = ['Lambda', 'Max_act','iter','speed','speed_std','persistance','std_persit'])
    df1.to_csv('FITFULL/ACT2_all.csv')


if __name__ == "__main__":
    build_csv('../../data/FITFULL_ACT_FRC/thin64_1/*')
