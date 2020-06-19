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
    folders = glob.glob(path)
    param_strings = [f.split('_')[3] for f in folders]
    params = set([re.findall(num_ptrn,s)[0] for s in param_strings])

    ind_rows = []
    global_rows = []
    for foldr in folders:
        files = glob.glob(foldr+ '/*')
        #print(files)
        print(foldr)
        dens = re.findall(num_ptrn,foldr.split('_')[-1])[0]
        print(dens)
        files.sort(reverse = True,key = lambda x: int(re.findall(num_ptrn,x)[-1]))
        # extract tracks from files :
        tracks = [np.loadtxt(f) for f in files]
        vec_tracks = np.array([to_vecs(t) for t in tracks])
        # lcl order :
        lcl_ordrs = local_order(tracks,vec_tracks,4)
        lcl_ordr = np.average(lcl_ordrs)
        std_lcl = np.std(lcl_ordrs)
        
        # handle boundaries ; 
        tracks = [handle_boundaries(t) for t in tracks]
        vec_tracks = np.array([to_vecs(t) for t in tracks])

        ordr1 = Global_order(vec_tracks)
        ordr2,std_ordr2 = OrderAt_T(vec_tracks)

        dts = 100
        dots_for_dts = [[] for i in range(dts)]
        pooled_dts = [auto_cor_pooled(t,dots_for_dts,dts) for t in tracks]
        averages = [np.mean(i) for i in dots_for_dts]
        pooled_pers = Persist_tracks([averages])
        if len(pooled_pers) == 0:
            pooled_pers = np.nan
        else:
            pooled_pers = pooled_pers[0]


        global_rows.append([dens,pooled_pers,ordr1,ordr2,std_ordr2,lcl_ordr,std_lcl])
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
            ind_rows.append([dens,i,speeds[i],ht])
        print('Persistance calculated')

#        global_rows.append([dens,pooled_pers,ordr1,ordr2,std_ordr2,lcl_ordr,std_lcl])
        print(np.average(speeds),global_rows[-1])
        print(ind_rows[-1])

    df1 = pd.DataFrame(data = ind_rows,
        columns = ['dens','cell_id','speed','persist'])
    df1.to_csv('density/DENS_PRFDR_ind.csv')
    df2 = pd.DataFrame(data = global_rows,
        columns = ['density','pooled_persist','global_order','sum_order','std_sum_order','lcl_order','std_lcl'])
    df2.to_csv('density/DENS_PRFDR_global.csv')

if __name__ == "__main__":
    build_csv('../../data/increase_DENS_PRFDR/*')