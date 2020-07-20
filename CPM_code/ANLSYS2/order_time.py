import pandas as pd
import numpy as np
from numpy.linalg import norm
import glob
import re
import sys
sys.path.insert(0,'../')
from OrderNpersist import to_vecs
import time
from STROMAL_ANLS import handle_boundaries

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
    #av_ordr = [norm(v) for v in orders]
        
    return orders

def order_time(path):
    # regex for finding correct files:
    num_ptrn = '[-+]?\d*\.\d+|\d+'
    gtypes = ['ER','BA', 'WS', 'PW','GM']
    folder = glob.glob(path)
    files = glob.glob(path+ '/*')
    params = set([s.split(re.findall(num_ptrn,s)[-2])[-1] for s in files])
    params = {p for p in params if p != '.txt'}
    n_params = len(list(params))
    print(n_params)
    #data_dict = {'time':[],'iter':[],'type':[],'x':[],'y':[],'z':[]}
    data_rows = []
    for gt in params:
        t1 = time.time()
        files = glob.glob(path+ '/*' + gt)
        print(len(files))
        print(gt)
        Type = gt[:2]
        itr = gt[2]
        tracks = [np.loadtxt(f) for f in files]
        tracks = [handle_boundaries(t) for t in tracks if len(t) == 200] # filter out old short tryout files
        vec_tracks = np.array([to_vecs(t) for t in tracks])
        n_cells = len(vec_tracks)
        orders = OrderAt_T(vec_tracks)
        for i,ordr in enumerate(orders):
            data_rows.append([i,itr,Type,ordr[0],ordr[1],ordr[2],tracks[i][0],tracks[i][1],tracks[i][2],n_cells])

    df = pd.DataFrame(data = data_rows,columns = ['time','iter','type','v_x','v_y','v_z','x','y','z','n_cells'])
    df.to_csv('STROMAL/2PRFDR_ordr_T.csv')

if __name__ == "__main__":
    order_time('../../data/STROMAL_PRFDR3')


