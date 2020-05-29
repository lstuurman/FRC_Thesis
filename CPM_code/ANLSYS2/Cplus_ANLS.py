import pandas as pd
import numpy as np
import glob
from numpy.linalg import norm
from itertools import product
import re
from scipy.spatial.distance import euclidean
import sys
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 
sys.path.insert(0,'../')
from OrderNpersist import Persist_tracks,Order_tracks2,order_radius,to_vecs
from analyse_tracks import new_auto
from CPM_helpers1 import plot_celltrack

def save_track(cell_track,prms1,prms2,i):
    # plot path of center of mass  :
    x = [c[0] for c in cell_track]
    y = [c[1] for c in cell_track]
    z = [c[2] for c in cell_track]

    #fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x,y,z)
    ax.scatter3D(x[0],y[0],z[0],label = 'start')
    ax.scatter3D(x[-1],y[-1],z[-1],label = 'end')
    ax.legend()
    plt.savefig('../../results/ACT_CHECKS/tracks4/singleLowE' + prms1 + prms2 + '_' + str(i) + '.png')

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
    """ Path is folder containing all celltracks of all max-lambda combinations.
    """
    # find unique param pairs in folder :
    num_ptrn = '[-+]?\d*\.\d+|\d+'
    all_files = glob.glob(path)
    param_strings = [f.split('_')[3] for f in all_files]
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
        print('number of cells for paramset : ',len(tracks))
        #for iter,t in enumerate(tracks):
        #    save_track(t,prms[0],prms[1],iter)
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
            #plt.plot(autocors[-1],range(len(autocors[-1])))
            #plt.savefig('../../results/ACT_CHECKS/AC4/' + prms[0] + prms[1] + 'Single500LowE_' + str(ti) +'.png')
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
        columns = ['Lambda', 'Max_act','speed','persistance','sum_order','global_order','lcl_order'])
    df1.to_csv('single_ACTZOOM5.csv')
    df2 = pd.DataFrame(data = deviation_rows,
        columns = ['Lambda', 'Max_act','speed','persistance','global_order','lcl_order'])
    df2.to_csv('single_ACTZOOM5_std.csv')

if __name__ == "__main__":
    build_csv('../../data/FIT_speedy2/single_ZOOM5/*')
