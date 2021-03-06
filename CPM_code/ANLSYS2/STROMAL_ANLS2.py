import pandas as pd
import numpy as np
import glob
from numpy.linalg import norm
import re
from scipy.spatial.distance import euclidean
import sys
import matplotlib.pyplot as plt
sys.path.insert(0,'../')
from OrderNpersist import Persist_tracks,Order_tracks2,order_radius,to_vecs
from analyse_tracks import new_auto,auto_cor_pooled
from ANLS_DENS import OrderAt_T,local_order
import time

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
    gtypes = ['ER','BA', 'WS', 'PW','GM']
    folder = glob.glob(path)
    files = glob.glob(path+ '/*')
    params = set([s.split(re.findall(num_ptrn,s)[-2])[-1] for s in files])
    params = {p for p in params if p != '.txt'}
    print(params)
    ind_rows = []
    global_rows = []
    for gt in params:
        t1 = time.time()
        files = glob.glob(path+ '/*' + gt)
        #print(files)
        print(gt)
        Type = gt[:2]
        itr = gt[2]
        print(Type,itr)
        files.sort(reverse = True,key = lambda x: int(re.findall(num_ptrn,x)[-1]))
        # extract tracks from files :
        tracks = [np.loadtxt(f) for f in files]
        print('unfiltered cells ',len(tracks))
        tracks = [t for t in tracks if len(t) == 200]

        # handle boundaries ; 
        tracks2 = [handle_boundaries(t) for t in tracks]
        vec_tracks = np.array([to_vecs(t) for t in tracks])

        ordr2,std_ordr2 = OrderAt_T(vec_tracks)

        dts = 100
        dots_for_dts = [[] for i in range(dts)]
        _ = [auto_cor_pooled(t,dots_for_dts,dts) for t in tracks]
        averages = [np.mean(i) for i in dots_for_dts]
        pooled_pers = Persist_tracks([averages])
        if len(pooled_pers) == 0:
            pooled_pers = np.nan
        else:
            pooled_pers = pooled_pers[0]


        global_rows.append([Type,int(itr),pooled_pers,ordr2,std_ordr2,lcl_ordr,std_lcl])
        print('number of cells for paramset : ',len(tracks))
        # speed : 
        vec_tracks = np.array([to_vecs(t) for t in tracks])
        speeds = [np.average([norm(v) for v in vec_track]) for vec_track in vec_tracks]
        print('Speed calculated')
        #peristance : 
        autocors = [] #[new_auto(t) for t in tracks]
        for ti,t in enumerate(tracks):
            autocors.append(new_auto(t))
            #print(ti)
        half_times = Persist_tracks(autocors)
        print('Length half Times : ',len(half_times))
        # save individual
        for i,ht in enumerate(half_times):
            ind_rows.append([Type,itr,i,speeds[i],ht])
        print('Persistance calculated')

#        global_rows.append([dens,pooled_pers,ordr1,ordr2,std_ordr2,lcl_ordr,std_lcl])
        print(np.average(speeds),global_rows[-1])
        print(ind_rows[-1])
        print('cumputed '+ gt + ' in',time.time() - t1)

    df1 = pd.DataFrame(data = ind_rows,
        columns = ['type','iter','cell_id','speed','persist'])
    df1.to_csv('STROMAL/PRFDR_ind.csv')
    df2 = pd.DataFrame(data = global_rows,
        columns = ['type','iter','pooled_persist','global_order','std_glbl_order','lcl_order','std_lcl'])
    df2.to_csv('STROMAL/PRFDR_glob.csv')

if __name__ == "__main__":
    build_csv('../../data/STROMAL_PRFDR')
