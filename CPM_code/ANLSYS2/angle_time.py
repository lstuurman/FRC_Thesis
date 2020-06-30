import pandas as pd
import numpy as np
from numpy.linalg import norm
import glob
import re
import sys
sys.path.insert(0,'../')
from OrderNpersist import to_vecs
import time

def save_angles(path):
    # regex for finding correct files:
    num_ptrn = '[-+]?\d*\.\d+|\d+'
    gtypes = ['ER','BA', 'WS', 'PW','GM']
    folder = glob.glob(path)
    files = glob.glob(path+ '/*')
    params = set([s.split(re.findall(num_ptrn,s)[-2])[-1] for s in files])
    params = {p for p in params if p != '.txt'}
    n_params = len(list(params))
    print(n_params)
    #data_dict = {'time':[],'iter':[],'cell_id':[],'type':[],'angle':[]}
    data_rows = []
    relative_rows = []
    for gt in params:
        t1 = time.time()
        files = glob.glob(path+ '/*' + gt)
        print(len(files))
        print(gt)
        Type = gt[:2]
        itr = gt[2]
        tracks = [np.loadtxt(f) for f in files]
        tracks = [t for t in tracks if len(t) == 200] # filter out old short tryout files
        vec_tracks = np.array([to_vecs(t) for t in tracks])
        #print('vec_tracks.shape : ',vec_tracks.shape)
        # loop over vector tracks per cell per timestep:
        ids = range(len(vec_tracks))
        times = range(len(vec_tracks[0])-1)
        for id in ids:
            for t in times:
                v1 = vec_tracks[id,t]
                v2 = vec_tracks[id,t+1]
                cos = np.dot(v1,v2)/(norm(v1) * norm(v2))
                angle = np.arccos(np.clip(cos,-1,1))
                data_rows.append([t,int(itr),id,Type,np.degrees(angle)])
                # angel relative to vector (1,1,1)
                cos2 = np.dot(v1,(1,1,1))/(norm(v1) * norm((1,1,1)))
                angle2 = np.arccos(np.clip(cos2,-1,1))
                relative_rows.append([t,int(itr),id,Type,np.degrees(angle2)])

        print(data_rows[-1])
        print(relative_rows[-1])
        #break
        print('paramset in :',time.time() - t1)
    df = pd.DataFrame(data = data_rows,columns = ['time','iter','cell_id','type','angle'])
    df.to_csv('STROMAL/PRFDR_angles.csv')
    df2 = pd.DataFrame(data = relative_rows,columns = ['time','iter','cell_id','type','angle'])
    df2.to_csv('STROMAL/PRFDR_relangles.csv')

if __name__ == "__main__":
    save_angles('../../data/STROMAL_PRFDR')
