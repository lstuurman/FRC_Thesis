import numpy as np 
import glob
#from MSD1cell import *
import pandas as pd
import re 
import pickle
from numpy.linalg import norm
from itertools import product,permutations,combinations
from analyse_tracks import new_auto

def to_vecs(cell_track):
    vecs = []
    for i in range(80):
        point1 = cell_track[i]
        point2 = cell_track[i + 1]
        v1 = np.array(point2 - point1)
        vecs.append(v1)
    return np.array(vecs)

def Order(files):
    #files = glob.glob(path)
    tracks = [np.loadtxt(f) for f in files]
    vec_tracks = np.array([to_vecs(t) for t in tracks])
    speeds = [norm(v) for v in vec_tracks]

    # list of order 
    orders = []
    for i in range(len(vec_tracks[0])):
        vecs_at_t = vec_tracks[:,i]
        # ignore 0 vectors : (nonmoving cells)
        vecs_at_t = [t for t in vecs_at_t if norm(t) != 0]

        pairs = combinations(vecs_at_t,2)
        cosines = [np.dot(v1,v2)/(norm(v1) * norm(v2)) for (v1,v2) in pairs] #if list(v1) != list(v2)
        orders += cosines
    # check for empties : 
    if len(orders) == 0 or len(speeds) == 0:
        print('Speeds : ',speeds)
    return np.mean(orders),np.mean(speeds)
    
def Order_tracks(vec_tracks):
    # list of order 
    orders = []
    for i in range(len(vec_tracks[0])):
        vecs_at_t = vec_tracks[:,i]
        # ignore 0 vectors : (nonmoving cells)
        vecs_at_t = [t for t in vecs_at_t if norm(t) != 0]

        pairs = combinations(vecs_at_t,2)
        cosines = [np.dot(v1,v2)/(norm(v1) * norm(v2)) for (v1,v2) in pairs] #if list(v1) != list(v2)
        orders.append(np.average(cosines))

    return np.mean(orders)

def Persist_tracks(autocors):

    half_lives = []
    # find half time auto corr : 
    for ac in autocors:
        if len(ac) == 0:
            # weird fucking non moving cell or something..
            continue
        ac = np.array(ac)
        #print(ac[0])
        max_indeces = np.where(ac > .5)[0]
        min_indeces = np.where(ac < .5)[0]
        if len(max_indeces) == 0 or len(min_indeces) == 0:
            print(ac)
            continue
        dt_0 = max(np.where(ac > .5)[0])
        dt_1 = min(np.where(ac < .5)[0])
        cos_0 = ac[dt_0]
        cos_1 = ac[dt_1]
        t_half = dt_0 + (.5 - cos_0)/(cos_1 - cos_0) * (dt_1 - dt_0)
        #print(cos_0,cos_1,t_half)
        half_lives.append(t_half)
    if len(half_lives) == 0:
        print(ac)
    return np.mean(half_lives)

def Persist(files):
    #files = glob.glob(path)
    tracks = [np.loadtxt(f) for f in files]
    autocors = [new_auto(t) for t in tracks]

    half_lives = []
    # find half time auto corr : 
    for ac in autocors:
        if len(ac) == 0:
            # weird fucking non moving cell or something..
            continue
        ac = np.array(ac)
        # print(np.where(ac > .5))
        # print(np.where(ac > .5)[0])
        # print(ac[0])
        max_indeces = np.where(ac > .5)[0]
        min_indeces = np.where(ac < .5)[0]
        if len(max_indeces) == 0 or len(min_indeces) == 0:
            print(ac)
            continue
        dt_0 = max(np.where(ac > .5)[0])
        dt_1 = min(np.where(ac < .5)[0])
        cos_0 = ac[dt_0]
        cos_1 = ac[dt_1]
        t_half = dt_0 + (.5 - cos_0)/(cos_1 - cos_0) * (dt_1 - dt_0)
        # print(dt_0,dt_1)
        # print(cos_0,cos_1)
        # print(t_half)
        half_lives.append(t_half)
    if len(half_lives) == 0:
        print(ac)
    return np.mean(half_lives)

def build_csv(path):
    """ Path is folder containing all celltracks of all max-lambda combinations.
    """
    # find unique param pairs in folder :
    num_ptrn = '[-+]?\d*\.\d+|\d+'
    all_files = glob.glob(path)
    param_strings = [f.split('_')[2] for f in all_files]
    params = set([tuple(re.findall(num_ptrn,s)) for s in param_strings])
    
    rows = []
    for i,prms in enumerate(params):
        files = [f for f in all_files if prms[0] + 'M' in f and prms[1] + '_' in f][20:30]
        print(len(files))
        ordr,speed = Order(files)
        prst = Persist(files)
        rows.append([i,prms[0],prms[1],speed,prst,ordr])
        if i > 5:
            break

    df = pd.DataFrame(data = rows,
        columns = ['index', 'Lambda', 'Max_act','speed','persistance','order'])
    df.to_csv('result_per_param.csv')
    



if __name__ == "__main__":
    build_csv('../data/full_ln1/*')