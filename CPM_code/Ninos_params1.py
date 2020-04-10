import cpm
import numpy as np
from CPM_helpers1 import *
from MSD1cell import *
import pickle
import pandas as pd
import time
import random

def setup_frc(D):
    # show how code works : 
    # set dimension for grid : 
    D = 64
    r = 20/64# 20microns
    # create graphs with positions : 
    g = nx.random_geometric_graph(80,r,dim = 3)
    # nice 3d visualazition in plotly : 
    matrix = nx.to_numpy_matrix(g)
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

def setup():
        ### SET UP CPM ###
    # params from Nino: multiplicated the adhesion engergies by 10
    # and because lambda of .1 not possible here. 
    # params suitable for single cell in empty space
    dimension = 64
    number_of_types = 3
    temperature = 200

    # initialize : 

    simulation = cpm.Cpm(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 2,target_area = 1000, lambda_area=250)
    simulation.set_constraints(cell_type = 2, lambda_perimeter = 2, target_perimeter = 500)#8600
    simulation.set_constraints(cell_type = 2, lambda_act = 20000, max_act = 80) # 2500, max_act = 42
    # adhesion ; 
    simulation.set_constraints(cell_type = 1,other_cell_type = 2,adhesion = -50)
    simulation.set_constraints(cell_type = 2,other_cell_type = 2,adhesion = 100)
    simulation.set_constraints(cell_type = 2,other_cell_type = 0,adhesion = 0)
    #simulation.initialize_from_array(cube_with_type,1)


    #print('Creating FRC')
    #frc_in_cube = setup_frc(64)
    #print('Loading into simulation : ')
    #simulation.initialize_from_array(frc_in_cube,1)
    #print('Done')

    ### fill cube with cells
    # number of cells :
    frc_in_cube = simulation.get_state() // 2**24 == 1 
    free_voxels = 64**3  - np.count_nonzero(frc_in_cube)
    n_cells = np.floor(.8 * (free_voxels/1000))
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
    n_cells = len(np.unique(celltypes))
    # print(n_cells)
    # exit()
    print('all cells full : ',n_cells * 700)
    # little warmup run 
    while np.sum(t) < n_cells * 700:
        cellids = s % 2**24
        t = []
        for n in range(2,n_cells + 2):
            t.append(np.sum(cellids == n))
        print(np.sum(t))
        simulation.run(10)
        #print(np.sum(t))
    print(t)
    return simulation

def runsim(simulation,steps):
    state = simulation.get_state() 

    ids = state % 2**24
    n_cells = len(np.unique(ids))

    iters = int(steps) # /10
    cofmass_track = np.zeros((n_cells,iters,3))
    #cell_sizes = []
    t0 = time.time()
    for i in range(iters):
        simulation.run(10)
        cell_sizes = []
        for n in range(2,n_cells + 1):
            cell = state % 2**24 == n
            # check cell_size : 
            size = np.sum(cell)
            cell_sizes.append(size)
            if size < 500:
                #print('to small')
                continue
            cofmass_track[n-2,i] = np.array(real_cofmass(cell, pr = False))
        if i == 0:
            t1 = time.time() - t0
            print('expected computing time : ',t1 * steps)
        print('iteration : ',i)
        print('number of small cells ',np.sum([1 for i in cell_sizes if i < 100]))
    return cofmass_track #,cell_sizes

if __name__ == "__main__":
    sim = setup()
    tracks = runsim(sim,500)
    for i,track in enumerate(tracks):
        fname = '../data/full_ln/frc_track_80dens' + str(i) + '.txt'
        np.savetxt(fname,track)
    #np.savetxt('../data/full_ln/frc_sizes.txt',sizes)
