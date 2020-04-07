import networkx as nx
import numpy as np
import networkx as nx 
import pickle
import matplotlib.pyplot as plt 
import glob
from scipy import ndimage
import pandas as pd
from scipy.spatial import distance

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
    np.savetxt('../results/Graphs_stuff/edge_lengths.txt',distances)

g = nx.random_geometric_graph(4500,20/256,dim = 3)
positions = []
for n,data in g.nodes(data = True):
    positions.append(data['pos'])

g,_ = nodesInCube(g,positions,256)
distances(256,g)
