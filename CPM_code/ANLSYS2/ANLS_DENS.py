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
from analyse_tracks import new_auto,auto_cor_pooled
from CPM_helpers1 import plot_celltrack

def handle_boundaries(cell_track,pr = False):
    # look for boundary crossings in any
    # of the coordinates
    cell_track2 = cell_track.copy()
    for i in range(len(cell_track) - 1):
        dif = np.subtract(cell_track[i],cell_track[i+1])
        for j,coordinate in enumerate(dif):
            if coordinate > 16:
                # went over boundary from 256 -> 0
                if pr:
                    print('Jumped from :',cell_track[i],'to :',cell_track[i+1])
                    print('Adding ',256, ' to rest of cell track') #cell_track[i,j]
                    print('changed axis : ',j)
                    print('Old coordinat : ',cell_track[i])

                cell_track2[:i + 1,j] -= 32

                if pr:
                    print('New coordinate : ',cell_track[i])
                    print(i,j)

            elif coordinate < -16:
                # form 0 -> 256
                if pr:
                    print('Jumped from :',cell_track[i],'to :',cell_track[i+1])
                    print('Adding ', 256, ' to previous of cell track') 
                    print('Old coordinat : ',cell_track[i])

                cell_track2[:i + 1,j] += 32

                if pr:
                    print('New coordinate : ',cell_track[i])
                    print(i,j)
    return cell_track2


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
                   
    #print(vt2.shape)
    orders = []
    for i in range(min_length):
        vecs_at_t = vt2[:,i]
        ordr = 0
        for v in vecs_at_t:
            ordr += v
        orders.append(ordr)
    #print('lenght of order shoulde be min_lenght. Min_lenght',min_length,'==',len(orders))
    av_ordr = np.average([norm(v) for v in orders])
    std_ordr = np.std([norm(v) for v in orders])
        
    return av_ordr,std_ordr

def local_order1(tracks,vec_tracks,bins = 5):
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
        #print('Time in track:',i)
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

def local_order(tracks,vec_tracks,bins = 4):
    # find smallest track : 
    min_length = min([len(t) for t in vec_tracks])
    print('smallest track : ',min_length)
    flat_tracks = [item for sublist in tracks for item in sublist.flatten()]
    max_d = max(flat_tracks) # end of domain : 
    # reshape vectors to square matrix : 
    vt2 = np.zeros((len(vec_tracks),min_length,3))
    for i,v in enumerate(vec_tracks):
        vt2[i] = v[:min_length]
    
    # go through vectors per timestep :
    orders = []
    for i in range(min_length):
        # bin according to position :
        vecs_at_t = vt2[:,i]
        pos = [tracks[j][i] for j in range(len(vecs_at_t))]
        binned_pos = np.digitize(pos,np.linspace(10,max_d,bins))
        # loop over bins : 
        for bn in np.unique(binned_pos,axis = 0):
            # find bins vectors in same bins : 
            indcs = [all(x) for x in binned_pos == bn] #binned_pos.all(bn)
            vecs_n_bin = vecs_at_t[indcs]
            # filter crossing of periodic boundary : 
            vecs_n_bin = [v for v in vecs_n_bin if norm(v) < 16]
            # consider at least 2 cells in one bin :
            if len(vecs_n_bin) < 2:
                continue
            # sum 
            norm_bin = norm(sum(vecs_n_bin))
            orders.append(norm_bin/len(vecs_n_bin))
    return orders

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
        #order : 
#        ordr1 = Global_order(vec_tracks)
#        ordr2,std_ordr2 = OrderAt_T(vec_tracks)
        print('Global order calculated')
#        lcl_ordrs = local_order(tracks,vec_tracks)
#        lcl_ordr = np.average(lcl_ordrs)
#        std_lcl = np.std(lcl_ordrs)
        
#        dts = 100
#        dots_for_dts = [[] for i in range(dts)]
#        pooled_dts = [auto_cor_pooled(t,dots_for_dts,dts) for t in tracks]
#        averages = [np.mean(i) for i in dots_for_dts]
#        pooled_pers = Persist_tracks([averages])
#        if len(pooled_pers) == 0:
#            pooled_pers = np.nan
#        else:
#            pooled_pers = pooled_pers[0]


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
