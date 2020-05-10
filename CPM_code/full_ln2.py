import cpm
import numpy as np
from CPM_helpers1 import *
from MSD1cell import *
import pickle
import pandas as pd
import time
import random

def handle_boundaries2(cell_track,pr = False):
    # look for boundary crossings in any
    # of the coordinates
    cell_track2 = cell_track.copy()
    for i in range(len(cell_track) - 1):
        dif = np.subtract(cell_track[i],cell_track[i+1])
        for j,coordinate in enumerate(dif):
            if coordinate > 64:
                # went over boundary from 256 -> 0
                if pr:
                    print('Jumped from :',cell_track[i],'to :',cell_track[i+1])
                    print('Adding ',256, ' to rest of cell track') #cell_track[i,j]
                    print('changed axis : ',j)
                    print('Old coordinat : ',cell_track[i])

                cell_track2[:i + 1,j] -= 64

                if pr:
                    print('New coordinate : ',cell_track[i])
                    print(i,j)
                
            elif coordinate < -64:
                # form 0 -> 256
                if pr:
                    print('Jumped from :',cell_track[i],'to :',cell_track[i+1])
                    print('Adding ', 256, ' to previous of cell track') 
                    print('Old coordinat : ',cell_track[i])

                cell_track2[:i + 1,j] += 64

                if pr:
                    print('New coordinate : ',cell_track[i])
                    print(i,j)
    return cell_track2



def setup_frc(D):
    # show how code works : 
    # set dimension for grid : 
    D = 64
    r = 20/64# 20microns
    # create graphs with positions : 
    g = nx.random_geometric_graph(80,r,dim = 3)
    # nice 3d visualazition in plotly : 
    #matrix = nx.to_numpy_matrix(g)
    positions = []
    for n,data in g.nodes(data = True):
        positions.append(data['pos'])

    # bin possitions and place in grid : 
    g,_ = nodesInCube(g,positions,D)
    c = fill_cube(D,g)
    cube = adjust_thickness(c,2)
    #adjust thickness so that it fills +- 17% of cube : 
    # show crossection ; 


    # mlab.contour3d(cube)
    # mlab.show()

    
    # slicee = cube[32:32+4]
    # plt.imshow(np.sum(slicee,axis = 0))
    # plt.show()

    # slicee = cube[32:32+4,:32,:32]
    # plt.imshow(np.sum(slicee,axis = 0))
    # plt.show()

    print('percentage of volume occupied by frc',(np.sum(cube)/D**3)*100)
    return cube

def setup(l_act,m_act):
        ### SET UP CPM ###
    # params from Nino: multiplicated the adhesion engergies by 10
    # and because lambda of .1 not possible here. 
    # params suitable for single cell in empty space
    dimension = 64
    number_of_types = 3
    temperature = 70

    # initialize : 

    simulation = cpm.Cpm(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 2,target_area = 1000, lambda_area=250)
    simulation.set_constraints(cell_type = 2, lambda_perimeter = 20, target_perimeter = 5400)#8600
    simulation.set_constraints(cell_type = 2, lambda_act = l_act, max_act = m_act) # 2500, max_act = 42
    # adhesion ; 
    simulation.set_constraints(cell_type = 1,other_cell_type = 2,adhesion = -50)
    simulation.set_constraints(cell_type = 2,other_cell_type = 2,adhesion = 100)
    simulation.set_constraints(cell_type = 2,other_cell_type = 0,adhesion = 0)


    print('Creating FRC')
    frc_in_cube = setup_frc(64)
    print('Loading into simulation : ')
    simulation.initialize_from_array(frc_in_cube,1)
    print('Done')

    ### fill cube with cells
    # number of cells :
    #frc_in_cube = simulation.get_state() // 2**24 == 1 
    free_voxels = 64**3  - np.count_nonzero(frc_in_cube)
    n_cells = np.floor(1 * (free_voxels/1000))
    # sample random positions : 
    free_indeces = np.where(frc_in_cube == 0.)
    possible_seeds = np.array(list(zip(free_indeces[0],free_indeces[1],free_indeces[2])))
    indx = np.random.randint(len(possible_seeds), size = int(n_cells))
    seeds = possible_seeds[indx,:]

    for c in seeds:
        simulation.add_cell(c[0],c[1],c[2],2)
    
    print('number of cells : ', len(seeds))

    s = simulation.get_state()

    celltypes = s % 2**24
    t = []
    n_cells = np.unique(celltypes)

    # little warmup run 
    simulation.run(100)

    for n in n_cells:
        print(n)
        celltypes = s % 2**24 == n
        size = np.sum(celltypes)
        t.append(size)
        print(n,size)

    # free voxels :

    return simulation

def runsim(simulation,steps):
    state = simulation.get_state() 

    ids = state % 2**24
    n_cells = np.unique(ids)

    iters = int(steps) # /10
    cofmass_track = np.zeros((len(n_cells)+1,iters,3))
    print(cofmass_track.shape)
    #cell_sizes = []
    t0 = time.time()
    for i in range(iters):
        simulation.run(10)
        cell_sizes = []
        for n in n_cells:
            cell = state % 2**24 == n
            # check cell_size : 
            size = np.sum(cell)
            cell_sizes.append(size)
            if size < 500:
                #print('to small')
                continue
            cofmass_track[n-1,i] = np.array(real_cofmass(cell,64,pr = False))
            #print(cofmass_track[n-1,i])
            #print(cofmass_track[n-2,i])
        if i == 0:
            t1 = time.time() - t0
            print('expected computing time : ',t1 * steps)
        print('iteration : ',i,cofmass_track[1,i])
        #print('number of small cells ',np.sum([1 for i in cell_sizes if i < 100]))
        #print(len(cofmass_track))
    print(cofmass_track.shape)
    return cofmass_track #,cell_sizes

if __name__ == "__main__":
    sim = setup(2000,1000)
    tracks = runsim(sim,1000)
    # for i,track in enumerate(tracks):
    #     fname = '../data/full_ln/frc_track_80dens' + str(i) + '.txt'
    #     np.savetxt(fname,track)
    #tracks = tracks.reshape((len(tracks) * 20,3))
    cell_track = tracks[1]
    plot_celltrack(cell_track)
    cell_track = handle_boundaries2(cell_track)
    plot_celltrack(cell_track)
    np.savetxt('testdat/testtrack_frc.txt',tracks[1])
