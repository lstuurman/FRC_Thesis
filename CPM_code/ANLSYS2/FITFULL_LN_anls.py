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

def build_csv(path):
    """ Path is folder containing all celltracks of all max-lambda combinations.
    """
    # find unique param pairs in folder :
    num_ptrn = '[-+]?\d*\.\d+|\d+'
    all_files = glob.glob(path)
    param_strings = [f.split('_')[4] for f in all_files]
    params = set([tuple(re.findall(num_ptrn,s)) for s in param_strings])

    rows = []
    deviation_rows = []
    for i,prms in enumerate(params):
        print(prms)
        files = [f for f in all_files if '_' + prms[0] + 'M' in f and 'X' + prms[1] + '_' in f]
        #files = [f for f in files if int(re.findall(num_ptrn,f)[-1]) > 1]
        files.sort(reverse = True,key = lambda x: int(re.findall(num_ptrn,x)[-1]))
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
        dts = 100
        dots_for_dts = [[] for i in range(dts)]
        pooled_dts = [auto_cor_pooled(t,dots_for_dts,dts) for t in tracks]
        averages = [np.mean(i) for i in dots_for_dts]
        pooled_pers = Persist_tracks([averages])
        if len(pooled_pers) == 0:
            pooled_pers = np.nan
        else:
            pooled_pers = pooled_pers[0]

        rows.append([prms[0],prms[1],speed,std_speed,pooled_pers])
        print(rows[-1])

    df1 = pd.DataFrame(data = rows,
        columns = ['Lambda', 'Max_act','speed','speed_std','persistance'])
    df1.to_csv('FITFULL/ACT1.csv')


if __name__ == "__main__":
    build_csv('../../data/FITFULL_ACT_FRC/thin44_1/*')