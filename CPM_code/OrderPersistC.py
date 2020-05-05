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


