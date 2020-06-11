import networkx as nx
import numpy as np
import networkx as nx 
import pickle
import matplotlib.pyplot as plt 
import glob
from scipy import ndimage
import pandas as pd
from scipy.spatial import distance
import sys
sys.path.insert(0,'../Code')
from Generate_graphs import ER_p_value,WS_K_value,M_power_cluster,E_M_relation
from graph_drawing_algs import fruchterman_reingold,normalize_postions
from Bresenheim import nodesInCube,fill_cube,adjust_thickness


def nodesInCube(g,positions,dim):
    # input is 3d list of [x,y,z] postitions
    dim = dim - 1
    xyz = np.array(positions)
    positions = xyz * dim
    positions = positions.astype(int)
    # replace node data:
    for n in g.nodes(data = True):
        n[1]['pos'] = positions[n[0]]
    return g,positions

def distances(dim,g):
    ### Fill cube of dimxdimxdim with points calculated by bresehheim based on stufff###
    #cube_draw = np.zeros((dim,dim,dim)) #import cv2
    distances = []
    for n1,n2 in g.edges():
        x1,y1,z1 = g.nodes[n1]['pos']
        x2,y2,z2 = g.nodes[n2]['pos']
        d = distance.euclidean((x1,y1,z1),(x2,y2,z2))
        distances.append(d)
        #points_on_line = Bresenham3D(x1,y1,z1,x2,y2,z2)
        #print(sys.getsizeof(points_on_line))
        #for p in points_on_line:
         #   cube_draw[p] = 1
    #np.savetxt('../results/Graphs_stuff/edge_lengths.txt',distances)
    return distances

def net_to_cube_adapted(g_type = 'ER',dim = 256):
    # generate graph of given type:
    print('Generating graph')
    if g_type == 'ER':
        p = ER_p_value(700,.25)
        g = nx.erdos_renyi_graph(700,p)#7000
        g = fruchterman_reingold(g,20)
    elif g_type == 'BA':
        g = nx.barabasi_albert_graph(700,4)#700
        g = fruchterman_reingold(g,20)
    elif g_type == 'WS':
        k = WS_K_value(1000,.25)
        p = 0.027825594022071243
        g = nx.watts_strogatz_graph(1000,k,p)#1000
        g = fruchterman_reingold(g,20)
    elif g_type == 'PW':
        m = M_power_cluster(700,.25)
        p = 0.666666666666666
        g = nx.powerlaw_cluster_graph(700,m,p)# 700
        g = fruchterman_reingold(g,20)
    elif g_type == 'GM':
        r = 20/dim # 20microns
        g = nx.random_geometric_graph(4500,r,dim = 3)#4500
    # extract positions from nx graph object : 
    positions = []
    for n,data in g.nodes(data = True):
        positions.append(data['pos'])
    positions = normalize_postions(positions)

    #bin possitions and place in grid :
    print('Filling cube with graph') 
    if g_type == 'GM' or g_type == 'WS':
        g,_ = nodesInCube(g,positions,dim)
        return g
    else:
        g,_ = nodesInCube(g,positions,4 * dim)
        return g

def save_cubes():
    g_types = ['ER','BA', 'WS', 'PW', 'GM']
    for graph_type in g_types: 
        g = net_to_cube_adapted(g_type = graph_type,dim = 256)
        dsts = distances(256,g)
        fname = '../results/Graphs_stuff/edge_lengths_' + graph_type +'.txt'
        np.savetxt(fname,dsts)
save_cubes()

#g = nx.random_geometric_graph(4500,20/256,dim = 3)
#positions = []
#for n,data in g.nodes(data = True):
#    positions.append(data['pos'])

#g,_ = nodesInCube(g,positions,256)
#distances(256,g)
