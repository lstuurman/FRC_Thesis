import numpy as np
import networkx as nx 
import pickle
import matplotlib.pyplot as plt 
import glob
from scipy import ndimage

import sys
sys.path.insert(0,'../Code')
from Generate_graphs import ER_p_value,WS_K_value,M_power_cluster,E_M_relation
from graph_drawing_algs import fruchterman_reingold
from Bresenheim import nodesInCube,fill_cube,adjust_thickness


##script for gap analysis of multiple networks

# generate networks in cube : 
def net_to_cube(g_type = 'ER',dim = 256):
    # generate graph of given type:
    print('Generating graph')
    if g_type == 'ER':
        p = ER_p_value(4500,.25)
        g = nx.erdos_renyi_graph(4500,p)
        g = fruchterman_reingold(g,25)
    elif g_type == 'BA':
        g = nx.barabasi_albert_graph(4500,4)
        g = fruchterman_reingold(g,25)
    elif g_type == 'WS':
        k = WS_K_value(4500,.25)
        p = 0.027825594022071243
        g = nx.watts_strogatz_graph(4500,k,p)
        g = fruchterman_reingold(g,25)
    elif g_type == 'PW':
        m = M_power_cluster(4500,.25)
        p = 0.666666666666666
        g = nx.powerlaw_cluster_graph(4500,m,p)
        g = fruchterman_reingold(g,4500)
    elif g_type == 'GM':
        r = 20/256 # 20microns
        g = nx.random_geometric_graph(4500,r,dim = 3)

    # extract positions from nx graph object : 
    positions = []
    for n,data in g.nodes(data = True):
        positions.append(data['pos'])

    # bin possitions and place in grid :
    print('Filling cube with graph') 
    g,_ = nodesInCube(g,positions,dim)
    c = fill_cube(dim,g)
    # adjust thickness so that it fills +- 17% of cube : 
    thicckk = adjust_thickness(c,2)
    return thicckk

def save_cubes():
    g_types = ['ER', 'BA', 'WS', 'PW', 'GM']
    for graph_type in g_types: 
        cube = net_to_cube(g_type = graph_type,dim = 256)
        fname = '../data/cubes/' + graph_type + '.pkl'
        pickle.dump(cube,open(fname,'wb'))
        print('created cube with ',graph_type, ' graph')

def take_random_slices(cube,scaling_f = 4):
    # random position for slice 
    index = np.random.randint(0,252)
    axes = np.random.randint(3)
    #  rotate and sum 4 layers of the cut : 
    cut = np.rot90(cube,axes)[index:index + 4]
    compressed_cut = np.sum(cut,axis = 0)
    # enlarge for better resolution of circles : 
    M = np.kron(compressed_cut, np.ones((scaling_f,scaling_f)))
    #plt.imshow(M)
    #plt.show()
    return M


def fill_circles(M):
    # pad with 4s to prevent circles over boundary:
    M = np.pad(M,1,'constant', constant_values = 4)
    
    width,height = M.shape
    # blueprint for circular masks :
    x_axis = np.arange(width)
    y_axis = np.arange(height)
    
    # negative of image to calculate maximim dists
    negative = M == 0.
    dist_matrix = ndimage.distance_transform_edt(negative)
    max_dist = int(max(dist_matrix.flatten()))
    
    # cut of boundaries : 
    cut_off = int(.5 * max_dist)
    dists = dist_matrix[cut_off:-cut_off,cut_off:-cut_off]
    cx,cy = np.where(dists == dists.max())
    cx += cut_off
    cy += cut_off
    max_radius = int(dist_matrix[cx,cy][0])
    
    # list for color values of circles : 
    colors = np.linspace(1,3,max_radius + 1)[::-1]
    # save radii used : 
    radii_list = []
    # frames for animation :
    dist_frames = []
    frames = []
    
    while max_radius > 7:
        # place largest possible circle at possition with max distance 
        mask = (x_axis[np.newaxis,:]-cy[0])**2 + (y_axis[:,np.newaxis]-cx[0])**2 < max_radius**2
        indeces = np.where(mask)
        # set mask to nice collor value
        mask = mask  * colors[max_radius]

        if sum(M[indeces]) != 0.: # Check if circle actually fits
            max_radius -= 1 
        else:
            # update :
            M = M + mask
            negative = M == 0.
            dist_matrix = ndimage.distance_transform_edt(negative)
            radii_list.append(max_radius)

            # find new circle : 
            max_dist = int(max(dist_matrix.flatten()))
            # cut of boundaries : 
            cut_off = int( .5 *  max_dist)#.5 *
            dists = dist_matrix[cut_off:-cut_off,cut_off:-cut_off]
            cx,cy = np.where(dists == dists.max())
            cx += cut_off
            cy += cut_off
            max_radius = int(dist_matrix[cx,cy][0])
            
    return M,radii_list

def main():
    #save_cubes()
    files = glob.glob('../data/cubes/*.pkl')
    for f in files:
        cube = pickle.load(open(f,'rb'))
        radii_data = []
        gtype = f[-6:-4]
        for i in range(5):
            M = take_random_slices(cube)
            M,radii = fill_circles(M)
            radii_data.append(radii)
            plt.imshow(M)
            #plt.show()
            plt.savefig('../data/cubes' + gtype + str(i) + '.png')
            print('../data/cubes' + gtype + str(i) + '.png')
            print('finished gap analysis : ',gtype,' ' + str(i))
        dfile = '../data/cubes/radii' + gtype + '.txt'
        np.savetxt(dfile,np.array(radii_data))



if __name__ == "__main__":
    main()

