import numpy as np 
import pandas as pd
from numpy.linalg import norm
from ANLS_DENS import Global_order
from scipy.ndimage import distance_transform_edt as edt
from STROMAL_ANLS import handle_boundaries
import sys 
sys.path.insert(0,'../')
from OrderNpersist import to_vecs
import pickle

#### Quantify stream formation at spefic sites

def in_radius(point,cntr,radius):
    if norm(point - cntr) < radius:
        return True
    else:
        return False

def find_cntrs(FRC):
    # negative : 
    FRC = FRC == 0
    dist_matrix = edt(FRC)
    # find centers :
    big_indeces = np.argpartition(dist_matrix,kth = 10,axis = None)[-10:]
    big_vals = np.take_along_axis(dist_matrix, big_indeces, axis=None) 
    cntrs = []
    for val in set(big_vals):
        cords = np.where(dist_matrix == val)
        cntrs += list(zip(*cords))

    return np.array(cntrs)

def edge_centers(FRC):
    pass


def stream_order(tracks,vec_tracks,centers,radius):
    
    # go through vectors per timestep :
    order_at_cntrs = np.zeros((len(vec_tracks),len(centers)))
    for i in range(len(vec_tracks)):
        # bin according to position :
        vecs_at_t = vec_tracks[:,i]
        vec_pos = [(tracks[j][i],vecs_at_t[j]) for j in range(len(vecs_at_t))]

        # loop over bins : 
        for cntr in centers:
            VP_at_cntr = [pair[1] for pair in vec_pos if in_radius(pair[0],cntr,radius)]
            ordr = Global_order(VP_at_cntr)
            order_at_cntrs[i,cntr] = ordr


    return order_at_cntrs


def streams(path):
    # find cell tracks  of some relevant files : 
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
    for gt in params[:1]:
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
        print(len(tracks))
        tracks2 = [handle_boundaries(t) for t in tracks]
        vec_tracks = np.array([to_vecs(t) for t in tracks2])

        frc_file = '../data/FRCs/' + str(itr) + Type + '64_diam3.pkl'
        #frc_gfile =  '../data/FRCs/GRAPH'+ str(itr) + Type + '64_diam3.pkl'
        frc = pickle.load(open(frc_file,'rb'))
        #graph = pickle.load(open(frc_gfile,'rb'))
        cntrs = find_cntrs(frc)
        cntr_ordrs = stream_order(tracks,vec_tracks,cntrs,4)


    df = pd.DataFrame(data = cntr_ordrs)
    df.to_csv('STROMAL/cntr_ordrs.csv')



    #     # find smallest track : 
    # min_length = min([len(t) for t in vec_tracks])
    # print('smallest track : ',min_length)
    # flat_tracks = [item for sublist in tracks for item in sublist.flatten()]
    # max_d = max(flat_tracks) # end of domain : 
    # # reshape vectors to square matrix : 
    # vt2 = np.zeros((len(vec_tracks),min_length,3))
    # for i,v in enumerate(vec_tracks):
    #     vt2[i] = v[:min_length]