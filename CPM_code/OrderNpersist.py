import numpy as np 
import glob
from MSD1cell import *
import pandas as pd
import re 
import pickle
from numpy.linalg import norm
from itertools import product,permutations
from analyse_tracks import new_auto

def to_vecs(cell_track):
    vecs = []
    for i in range(80):
        point1 = cell_track[i]
        point2 = cell_track[i + 1]
        v1 = np.array(point2 - point1)
        vecs.append(v1)
    return np.array(vecs)

def Order(path):
    files = glob.glob(path)
    tracks = [np.loadtxt(f) for f in files]
    vec_tracks = np.array([to_vecs(t) for t in tracks])

    # list of order 
    orders = []
    for i in range(len(vec_tracks[0])):
        vecs_at_t = vec_tracks[:,i]
        pairs = permutations(vecs_at_t)
        cosines = [np.dot(v1,v2)/(norm(v1) * norm(v2)) for (v1,v2) in pairs] #if list(v1) != list(v2)
        angels_per_t += cosines

    return np.mean(orders)
    

def Persist(path):
    files = glob.glob(path)
    tracks = [np.loadtxt(f) for f in files]
    autocors = [new_auto(t) for t in tracks]

    half_lives = []
    # find half time auto corr : 
    for ac in autocors:
        dt_0 = max(np.where(ac > .5)[0])
        dt_1 = min(np.where(ac < .5)[0])
        cos_0 = ac[dt_0]
        cos_1 = ac[dt_1]
        t_half = dt_0 + (.5 - cos_0)/(cos_1 - cos_0) * (dt_1 - dt_0)
        print(dt_0,dt_1)
        print(cos_0,cos_1)
        print(t_half)
        break
        half_lives.append(t_half)

    return np.mean(half_lives)

