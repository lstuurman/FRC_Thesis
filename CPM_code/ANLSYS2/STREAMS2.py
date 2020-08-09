import numpy as np 
import pandas as pd
from numpy.linalg import norm
from ANLS_DENS import Global_order
from scipy.ndimage import distance_transform_edt as edt
from STROMAL_ANLS import handle_boundaries
from scipy.spatial import distance
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
    radii_list = [max_radius]

    while max_radius > 5:
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
            radii_list.append(max_radius)

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
            
    return radii_list,cntr_list

def edge_centers(g):
    # create tuples of dist,cntr
    dsts_cntrs = []
    for n1,n2 in g.edges():
        x1,y1,z1 = g.nodes[n1]['pos']
        x2,y2,z2 = g.nodes[n2]['pos']
        d = distance.euclidean((x1,y1,z1),(x2,y2,z2))
        c = [(x1 + x2)/2,(y1 + y2)/2,(z1 + z2)/2]
        dsts_cntrs.append((d,c))
    
    # order on max dist : 
    dsts_cntrs.sort(key = lambda x:x[0])
    return [c[0] for c in dsts_cntrs],[c[1] for c in dsts_cntrs]

def edge_centers_adjstd(g):
    
    dsts_cntrs = []
    for n1,n2 in g.edges():
        x1,y1,z1 = g.nodes[n1]['pos']
        x2,y2,z2 = g.nodes[n2]['pos']
        coords = [x1,y1,z1,x2,y2,z2]
        # truncate coordinates to be withing bounds :
        #coords = coords - 96 
        #print(coords)
        for i in range(len(coords)):
            coords[i] = coords[i] - 96
            if coords[i] < 0:
                coords[i] = 0
            if coords[i] > 64:
                coords[i] = 64

        x1,y1,z1,x2,y2,z2 = coords
        #print(coords)
        d = distance.euclidean((x1,y1,z1),(x2,y2,z2))
        c = [(x1 + x2)/2,(y1 + y2)/2,(z1 + z2)/2]
        #print(d,c)
        #exit()
        dsts_cntrs.append((d,c))
    #print(dsts_cntrs)
    # order on max dist : 
    dsts_cntrs.sort(key = lambda x:x[0])
    #print(dsts_cntrs[-10:])
    return [c[0] for c in dsts_cntrs],[c[1] for c in dsts_cntrs] 
 


def stream_order(tracks,vec_tracks,centers,radius):
    #print(vec_tracks.shape)
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
    edge_rows = []
    gap_rows = []
    for gt in list(params):
        #t1 = time.time()
        files = glob.glob(path+ '/*' + gt)
        #print(files)
        print(gt)
        Type = gt[:-5]
        if Type == 'OFRC':
            continue
        itr = gt[-5]
        print(Type,itr)
        files.sort(reverse = True,key = lambda x: int(re.findall(num_ptrn,x)[-1]))
        # extract tracks from files :
        tracks = [np.loadtxt(f) for f in files]
        print('unfiltered cells ',len(tracks))
        tracks = [t for t in tracks if len(t) == 200]
        print(len(tracks))
        tracks2 = [handle_boundaries(t) for t in tracks]
        vec_tracks = np.array([to_vecs(t) for t in tracks2])

        frc_file = '../../data/FRCs/' + Type + '64_diam3.pkl'
        frc_gfile =  '../../data/FRCs/GRAPH' + Type + '64_diam3.pkl'
        frc = pickle.load(open(frc_file,'rb'))
        g = pickle.load(open(frc_gfile,'rb'))
        if Type == 'GM' or Type == 'WS':
            edge_lengths,edge_cntrs = edge_centers(g)
        else:
            edge_lengths,edge_cntrs = edge_centers_adjstd(g)
        #print(edge_cntrs)
        radii,cntrs = fill_circles(frc)
        #print(cntrs)
        #print(len(cntrs))
        # make edge lists same size as gap list : 
        edge_lengths = edge_lengths[:len(cntrs)]
        edge_cntrs = edge_cntrs[:len(cntrs)]
        print(len(cntrs),len(radii),len(edge_cntrs),len(edge_lengths))
        gap_ordrs = stream_order(tracks,vec_tracks,cntrs,5)
        edge_ordrs = stream_order(tracks,vec_tracks,edge_cntrs,5)
        for t,gaps in enumerate(gap_ordrs):
            for c,gp in enumerate(gaps):
                gap_rows.append([t,c,cntrs[c],radii[c],gp,Type,itr])
                edge_rows.append([t,c,edge_cntrs[c],edge_lengths[c],edge_ordrs[t,c],Type,itr])
        #print(edge_ordrs.shape)
        print(len(cntrs))
        print(np.mean(gap_ordrs))
        print(np.mean(edge_ordrs))

    df = pd.DataFrame(data = gap_rows,columns = ['time','center','coords','radius','order','type','iter'])
    df.to_csv('STROMAL/gap_ordrs3.csv')
    df = pd.DataFrame(data = edge_rows,columns = ['time','center','coords','length','order','type','iter'])
    df.to_csv('STROMAL/edge_ordrs3.csv')



if __name__ == "__main__":
    streams('../../data/STROMAL_ACT3')

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


