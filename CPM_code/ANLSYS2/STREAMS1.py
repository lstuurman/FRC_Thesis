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
import glob 
import re
#### Quantify stream formation at spefic sites

def in_radius(point,cntr,radius):
    if norm(point - cntr) < radius:
        return True
    else:
        return False



def fill_circles(M):
    # pad with 4s to prevent circles over boundary:
    M = np.pad(M,1,'constant', constant_values = 1)
    
    width,height,depth = M.shape
    # blueprint for circular masks :
    x_axis = np.arange(width)
    y_axis = np.arange(height) 
    z_axis = np.arange(depth)

    # negative of image to calculate maximim dists
    negative = M == 0.
    dist_matrix = edt(negative)
    max_dist = int(max(dist_matrix.flatten()))
    
    # cut of boundaries : 
    cut_off = int(.5 * max_dist)
    dists = dist_matrix[cut_off:-cut_off,cut_off:-cut_off,cut_off:-cut_off]
    cx,cy,cz = np.where(dists == dists.max())
    cx += cut_off
    cy += cut_off
    cz += cut_off
    #print(cx,cy,cz)
    max_radius = int(dist_matrix[cx[0],cy[0],cz[0]])#[0]
    
    # list for color values of circles : 
    colors = np.linspace(1,3,max_radius + 1)[::-1]
    # save radii used : 
    cntr_list = [[cx[0],cy[0],cz[0]]]
    # frames for animation :
    #dist_frames = []
    #frames = []
    
    while len(cntr_list) < 10:
        #print([cx[0],cy[0],cz[0]])
        #print(max_radius)
        # place largest possible circle at possition with max distance 
        mask = (x_axis[np.newaxis,:,np.newaxis]-cy[0])**2 + (y_axis[:,np.newaxis,np.newaxis]-cx[0])**2 + (z_axis[np.newaxis,np.newaxis,:] -cz[0])**2 < max_radius**2
        indeces = np.where(mask)
        # set mask to nice collor value
        mask = mask  * colors[max_radius]
        #print(indeces)
        if sum(M[indeces]) != 0.: # Check if circle actually fits
            max_radius -= 1 
        else:
            # update :
            M = M + mask
            negative = M == 0.
            dist_matrix = edt(negative)
            #radii_list.append(max_radius)

            # find new circle : 
            max_dist = int(max(dist_matrix.flatten()))
            # cut of boundaries : 
            cut_off = int( .5 *  max_dist)#.5 *
            dists = dist_matrix[cut_off:-cut_off,cut_off:-cut_off,cut_off:-cut_off]
            cx,cy,cz = np.where(dists == dists.max())
            cx += cut_off
            cy += cut_off
            cz += cut_off
            max_radius = int(dist_matrix[cx[0],cy[0],cz[0]])
            cntr_list.append([cx[0],cy[0],cz[0]])
            
    return cntr_list

def edge_centers(FRC):
    pass


def stream_order(tracks,vec_tracks,centers,radius):
    print(vec_tracks.shape)
    # go through vectors per timestep :
    order_at_cntrs = np.zeros((len(vec_tracks[0]),len(centers)))
    for i in range(len(vec_tracks[0])):
        # bin according to position :
        vecs_at_t = vec_tracks[:,i]
        vec_pos = [(tracks[j][i],vecs_at_t[j]) for j in range(len(vecs_at_t))]

        # loop over bins : 
        for j,cntr in enumerate(centers):
            #print(cntr)
            VP_at_cntr = [pair[1] for pair in vec_pos if in_radius(pair[0],cntr,radius)]
            ordr = Global_order(VP_at_cntr)
            order_at_cntrs[i,j] = ordr


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
    for gt in list(params)[:1]:
        #t1 = time.time()
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

        frc_file = '../../data/FRCs/' + str(itr) + Type + '64_diam3.pkl'
        #frc_gfile =  '../../data/FRCs/GRAPH'+ str(itr) + Type + '64_diam3.pkl'
        frc = pickle.load(open(frc_file,'rb'))
        #graph = pickle.load(open(frc_gfile,'rb'))
        #cntrs = find_cntrs(frc)
        cntrs = fill_circles(frc)
        #print(cntrs)
        #print(len(cntrs))
        cntr_ordrs = stream_order(tracks,vec_tracks,cntrs,4)
        print(cntr_ordrs.shape)
        print(np.mean(cntr_ordrs))

    df = pd.DataFrame(data = cntr_ordrs)
    df.to_csv('STROMAL/cntr_ordrs.csv')


if __name__ == "__main__":
    streams('../../data/STROMAL_PRFDR')

    #     # find smallest track : 
    # min_length = min([len(t) for t in vec_tracks])
    # print('smallest track : ',min_length)
    # flat_tracks = [item for sublist in tracks for item in sublist.flatten()]
    # max_d = max(flat_tracks) # end of domain : 
    # # reshape vectors to square matrix : 
    # vt2 = np.zeros((len(vec_tracks),min_length,3))
    # for i,v in enumerate(vec_tracks):
    #     vt2[i] = v[:min_length]

#def find_cntrs(FRC):
#    # negative : 
#    FRC = FRC == 0
#    dist_matrix = edt(FRC)
#    # find centers :
#    #big_indeces = np.argpartition(dist_matrix,kth = 10,axis = None)[-10:]
#    #big_vals = np.take_along_axis(dist_matrix, big_indeces, axis=None)[:-1]
#    flat_edt = dist_matrix.flatten()
#    #print(flat_edt.shape)
#    big_vals = flat_edt[np.argsort(flat_edt)]
#    vls_at_cntrs = np.array(sorted(list(set(big_vals)))[::-1])
#    cntrs = []
#    count = 0
#    #print(vls_at_cntrs[:10],vls_at_cntrs[-10:])
#    #print('lenght biggest values ; ',vls_at_cntrs)
#    for val in vls_at_cntrs:
#        cords = np.where(dist_matrix == val)
#        cntr_list = list(zip(*cords))
#        if len(cntrs) > 0:
#            #for c in cntrs:
#            for c_ in cntr_list:
#                if not all([in_radius(np.array(c_),np.array(c),8) for c in cntrs]):
#                    cntrs.append(c_)
#                    break
#                #break
#            #break
#                #cntrs += [c_ for c_ in cntr_list if not in_radius(np.array(c_),np.array(c),8)]
#        elif len(cntrs) > 10:
#            break
#        else:
#            cntrs.append(tuple(cntr_list[0]))
#        print(val)
#        print(cntrs)
#    for val in vls_at_cntrs[:10]:
#        print('distance ',val)
#        cords = np.where(dist_matrix == val)
#        cntrs += list(zip(*cords))
#        print(list(zip(*cords)))
#    return np.array(cntrs)

