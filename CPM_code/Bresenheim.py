import networkx as nx 
import numpy as np 
from mayavi import mlab
#from Code.ThreeDdraw import draw_plotly
import sys
sys.path.insert(0,'../Code')
# Add the Code folder path to the sys.path list
# sys.path.append('/../Code/')
from ThreeDdraw import draw_plotly


### FUNCTIONS FOR PLACING GRAPH IN CPM ###

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


# Python3 code for generating points on a 3-D line  
# using Bresenham's Algorithm 
### ADOPTTED FROM ###
#https://www.geeksforgeeks.org/bresenhams-algorithm-for-3-d-line-drawing/
  
def Bresenham3D(x1, y1, z1, x2, y2, z2): 
    ListOfPoints = [] 
    ListOfPoints.append((x1, y1, z1)) 
    dx = abs(x2 - x1) 
    dy = abs(y2 - y1) 
    dz = abs(z2 - z1) 
    if (x2 > x1): 
        xs = 1
    else: 
        xs = -1
    if (y2 > y1): 
        ys = 1
    else: 
        ys = -1
    if (z2 > z1): 
        zs = 1
    else: 
        zs = -1
  
    # Driving axis is X-axis" 
    if (dx >= dy and dx >= dz):         
        p1 = 2 * dy - dx 
        p2 = 2 * dz - dx 
        while (x1 != x2): 
            x1 += xs 
            if (p1 >= 0): 
                y1 += ys 
                p1 -= 2 * dx 
            if (p2 >= 0): 
                z1 += zs 
                p2 -= 2 * dx 
            p1 += 2 * dy 
            p2 += 2 * dz 
            ListOfPoints.append((x1, y1, z1))
            #print(sys.getsizeof(ListOfPoints))
  
    # Driving axis is Y-axis" 
    elif (dy >= dx and dy >= dz):        
        p1 = 2 * dx - dy 
        p2 = 2 * dz - dy 
        while (y1 != y2): 
            y1 += ys 
            if (p1 >= 0): 
                x1 += xs 
                p1 -= 2 * dy 
            if (p2 >= 0): 
                z1 += zs 
                p2 -= 2 * dy 
            p1 += 2 * dx 
            p2 += 2 * dz 
            ListOfPoints.append((x1, y1, z1)) 
            #print(sys.getsizeof(ListOfPoints))
  
    # Driving axis is Z-axis" 
    else:         
        p1 = 2 * dy - dz 
        p2 = 2 * dx - dz 
        while (z1 != z2): 
            z1 += zs 
            if (p1 >= 0): 
                y1 += ys 
                p1 -= 2 * dz 
            if (p2 >= 0): 
                x1 += xs 
                p2 -= 2 * dz 
            p1 += 2 * dy 
            p2 += 2 * dx 
            ListOfPoints.append((x1, y1, z1)) 
            #print(sys.getsizeof(ListOfPoints))
    return ListOfPoints 

### END ADOPTED CODE ###

def fill_cube(dim,g):
    ### Fill cube of dimxdimxdim with points calculated by bresehheim based on stufff###
    cube_draw = np.zeros((dim,dim,dim))#import cv2
    for n1,n2 in g.edges():
        x1,y1,z1 = g.nodes[n1]['pos']
        x2,y2,z2 = g.nodes[n2]['pos']
        points_on_line = Bresenham3D(x1,y1,z1,x2,y2,z2)
        #print(sys.getsizeof(points_on_line))
        for p in points_on_line:
            cube_draw[p] = 1
    return cube_draw

def adjust_thickness(cube,t):
    # create padding for boundaries : 
    pad_shape = tuple(d+2*t for d in cube.shape)
    pad = np.zeros(pad_shape)
    pad[:cube.shape[0],:cube.shape[1],:cube.shape[2]] = cube
    
    # go over gridpoints part of graph: 
    indeces = np.where(pad==1)
    for triplet in zip(*indeces):
        x,y,z = triplet
        # fill moores neigbourhood
        pad[x-t:x+t,y-t:y+t,z-t:z+t] = 1
    
    # cuttof padding : 
    cube = pad[:cube.shape[0],:cube.shape[1],:cube.shape[2]]
    return cube

def test(D,plotly = False,mayavi = False):
    # show how code works : 
    # set dimension for grid : 
    D = 256
    r = 20/256# 20microns
    # create graphs with positions : 
    g = nx.random_geometric_graph(4500,r,dim = 3)
    # nice 3d visualazition in plotly : 
    matrix = nx.to_numpy_matrix(g)
    positions = []
    for n,data in g.nodes(data = True):
        positions.append(data['pos'])
    if plotly:
        draw_plotly(positions,matrix)

    # bin possitions and place in grid : 
    g,_ = nodesInCube(g,positions,D)
    c = fill_cube(D,g)
    #adjust thickness so that it fills +- 17% of cube : 
    thicckk = adjust_thickness(c,2)
    if mayavi:
        #mlab.clf()
        mlab.options.offscreen = True
        mlab.contour3d(thicckk)
        #mlab.show()
        mlab.savefig('test.png')
    print('percentage of volume occupied by frc',(np.sum(thicckk)/D**3)*100)
    return thicckk

if __name__ == "__main__":
    test(256,mayavi=True)
