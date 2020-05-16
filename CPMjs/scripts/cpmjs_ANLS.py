import pandas as pd
import numpy as np
import glob
from numpy.linalg import norm
from itertools import product
import re

import sys
sys.path.insert(0,'../../CPM_code/')
from OrderNpersist import Persist_tracks,Order_tracks2,order_radius,to_vecs
from analyse_tracks import new_auto

def handle_boundaries(cell_track,pr = False):
    # look for boundary crossings in any
    # of the coordinates
    cell_track2 = cell_track.copy()
    for i in range(len(cell_track) - 1):
        dif = np.subtract(cell_track[i],cell_track[i+1])
        for j,coordinate in enumerate(dif):
            if coordinate > 64:
                # went over boundary from 256 -> 
                cell_track2[:i + 1,j] -= 64

            elif coordinate < -64:
                cell_track2[:i + 1,j] += 64
    return cell_track2


def Global_order(vec_tracks):
    # sum vec instead of calculating angles : 
    orders = 0
    for vtrack in vec_tracks:
        for v in vtrack:
            orders += v
    return norm(orders)

def OrderAt_T(vec_tracks):
    # find smallest track : 
    min_length = min([len(t) for t in vec_tracks])
    
    # reshape vectors to square matrix : 
    vt2 = np.zeros((len(vec_tracks),min_length,3))
    for i,v in enumerate(vec_tracks):
        vt2[i] = v[:min_length]
                   
    print(vt2.shape)
    orders = []
    for i in range(min_length):
        vecs_at_t = vt2[:,i]
        ordr = 0
        for v in vecs_at_t:
            ordr += v
        orders.append(ordr)
    print('lenght of order shoulde be min_lenght. Min_lenght',min_length,'==',len(orders))
    av_ordr = np.average([norm(v) for v in orders])
    std_ordr = np.std([norm(v) for v in orders])
        
    return av_ordr,std_ordr

def local_order(tracks,vec_tracks,bins = 5):
    # find smallest track : 
    min_length = min([len(t) for t in vec_tracks])
    flat_tracks = [item for sublist in tracks for item in sublist.flatten()]
    max_d = max(flat_tracks) # end of domain 
    # reshape vectors to square matrix : 
    vt2 = np.zeros((len(vec_tracks),min_length,3))
    for i,v in enumerate(vec_tracks):
        vt2[i] = v[:min_length]
    
    # go through vectors per timestep :
    orders = []
    for i in range(min_length):
        print('Time in track:',i)
        # bin according to position :
        vecs_at_t = vt2[:,i]
        pos = [tracks[j][i] for j in range(len(vecs_at_t))]
        binned_pos = np.digitize(pos,np.linspace(10,max_d,bins))
        # loop over bins : 
        for bn in np.unique(binned_pos,axis = 0):
            # find bins vectors in same bins : 
            indcs = [all(x) for x in binned_pos == bn]
            vecs_n_bin = vecs_at_t[indcs]
            # consider at least 2 cells in one bin :
            if len(vecs_n_bin) < 2:
                continue
            # sum 
            norm_bin = norm(sum(vecs_n_bin))
            orders.append(norm_bin)
    return orders

def build_csv(path):
    """ Path is folder containing all celltracks of all Persist -lambda combinations.
    """
    # find unique param pairs in folder :
    num_ptrn = '[-+]?\d*\.\d+|\d+'
    all_files = glob.glob(path)
    #params = [tuple(re.findall(num_ptrn,s)) for f in all_files]

    rows = []
    deviation_rows = []
    for i,f in enumerate(all_files):
        data = pd.read_csv(f, sep="\t", names=['time','id','ctype', 'x', 'y', 'z'])
        ids = data.id.unique()
        tracks = [data[data.id == ID].to_numpy()[:,-3:] for ID in ids]
        prms = re.findall(num_ptrn,f.split('/')[-1])
        print(prms) 
        # fix boundaries : 
        tracks = [handle_boundaries(t) for t in tracks]
        print('number of cells for paramset : ',len(tracks))
        if len(tracks) == 0:
            print('no cell track??',f)
            continue
        # speed : 
        vec_tracks = np.array([to_vecs(t) for t in tracks])
        speeds = [[norm(v) for v in vec_track] for vec_track in vec_tracks]
        speeds = [item for sublist in speeds for item in sublist]
        speed = np.mean(speeds)
        std_speed = np.std(speeds)
        print('Speed calculated')
        #peristance : 
        autocors = [] #[new_auto(t) for t in tracks]
        for ti,t in enumerate(tracks):
            autocors.append(new_auto(t))
            print(ti)
        half_times = Persist_tracks(autocors)
        ht = np.average(half_times)
        std_ht = np.std(half_times)
        print('Persistance calculated')
        #order : 
        ordr1 = Global_order(vec_tracks)
        ordr2,std_ordr2 = OrderAt_T(vec_tracks)
        print('Global order calculated')
        lcl_ordrs = local_order(tracks,vec_tracks)
        lcl_ordr = np.average(lcl_ordrs)
        std_lcl = np.std(lcl_ordrs)

        rows.append([prms[0],prms[1],speed,ht,ordr1,ordr2,lcl_ordr])
        deviation_rows.append([prms[0],prms[1],std_speed,std_ht,std_ordr2,std_lcl])
        print(rows[-1])

    df1 = pd.DataFrame(data = rows,
        columns = ['Lambda', 'Persist','speed','persistance','sum_order','global_order','lcl_order'])
    df1.to_csv('V500_PRFDR_single.csv')
    df2 = pd.DataFrame(data = deviation_rows,
        columns = ['Lambda', 'Perist','speed','persistance','global_order','lcl_order'])
    df2.to_csv('V500_PRFDR_single_std.csv')

if __name__ == "__main__":
    build_csv('../../data/cpmjs/single_PRFDR_V500/*')
