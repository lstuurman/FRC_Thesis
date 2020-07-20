import cpm
import numpy as np
from CPM_helpers1 import real_cofmass

import pickle
import time
import random

from itertools import product
from multiprocessing import Pool
import multiprocessing
from Bresenheim import nodesInCube,fill_cube,adjust_thickness

def return_bounds(dim):
    cube = np.zeros((dim,dim,dim))
    #set face of cube to 1:
    cube[0,:,:] = 1
    cube[-1,:,:] = 1
    cube[:,0,:] = 1
    cube[:,-1,:] = 1
    cube[:,:,0] = 1
    cube[:,:,-1] = 1
    return cube.astype(np.uint32)

def setup(g_type):
        ### SET UP CPM ###
    dimension = 256
    number_of_types = 4
    temperature = 200

    # initialize : 

    simulation = cpm.Cpm3d(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 2,target_area = 150, lambda_area=25)
    simulation.set_constraints(cell_type = 2, lambda_perimeter = .2, target_perimeter = 1200) #8600
    simulation.set_constraints(cell_type = 2, lambda_act = 1200, max_act = 40) # 2500, max_act = 42
    # adhesion ; 
    simulation.set_constraints(cell_type = 2,other_cell_type = 1,adhesion = -5)
    simulation.set_constraints(cell_type = 2,other_cell_type = 2,adhesion = 10)
    simulation.set_constraints(cell_type = 2,other_cell_type = 0,adhesion = 0)
    simulation.set_constraints(cell_type = 2,other_cell_type = 3,adhesion = 10)

    simulation.set_constraints(cell_type = 1, fixed=True)
    simulation.set_constraints(cell_type = 3,fixed = True)
    try:
        #print('Creating FRC')
        dfile = '../data/FRCs/256_' + g_type + '.pkl'
        frc_in_cube = pickle.load(open(dfile,'rb')).astype(np.uint32)
        print('Loading into simulation : ')
        bounds = return_bounds(dimension) 
        frc_in_cube[bounds == 1] = 0
        simulation.initialize_from_array(frc_in_cube,1)
        simulation.initialize_from_array(bounds,3)
        print('Done')
    except:
        #pass
        print('No FRC')
        # fill boundaries with boundary cell :
        bounds = return_bounds(dimension) 
        simulation.initialize_from_array(bounds,3)
    print('boundaries loaded')
    ### fill cube with cells
    # number of cells :
    frc_in_cube = simulation.get_state() // 2**24 == 1 
    free_voxels = 256**3  - np.count_nonzero(frc_in_cube)
    n_cells = np.floor(1 * (free_voxels/150))
    print(n_cells)
    # sample random positions : 
    free_indeces = np.where(frc_in_cube == 0)
    possible_seeds = np.array(list(zip(free_indeces[0],free_indeces[1],free_indeces[2])))
    indx = np.random.randint(len(possible_seeds), size = int(n_cells))
    seeds = possible_seeds[indx,:]
    for c in seeds:
        simulation.add_cell(c[0],c[1],c[2],2)
    
    print('number of cells seeded: ', len(seeds))
    s = simulation.get_state()
    celltypes = s % 2**24
    t = []
    n_cells = np.unique(celltypes)
    print(n_cells)
    # little warmup run
    print('warming up')
    t1 = time.time() 
    simulation.run(50)
    print('warm up in : ',time.time() - t1)

    for n in n_cells:
        #print(n)
        celltypes = s % 2**24 == n
        size = np.sum(celltypes)
        t.append(size)
        #print(n,size)

    print('free voxels : ',256**3 - sum(t))
    return simulation


def runsim(simulation,steps):
    state = simulation.get_state() 

    ids = state % 2**24
    n_cells = np.unique(ids)

    iters = int(steps) # /10
    cofmass_track = np.zeros((len(n_cells)-3,iters,3))
    print(cofmass_track.shape)
    print(len(n_cells))
    print(n_cells[2:])
    print([i for i in range(n_cells[-1]) if i not in n_cells])
    t0 = time.time()
    for i in range(iters):
        simulation.run(40)
        #cell_sizes = []
        for ind,n in enumerate(n_cells[3:]):
            cell = state % 2**24 == n
            # check cell_size : 
            #size = np.sum(cell)
            #print(size)
            #cell_sizes.append(size)
            #if size <100:
                #print('to small')
            #    continue
            cofmass_track[ind,i] = np.array(real_cofmass(cell,64,pr = False))

        if i == 0:
            t1 = time.time() - t0
            print('expected computing time : ',t1 * steps)
        print(ind,i)
    return cofmass_track


def run_grid_point(params):
    gtype,iter = params
    t1 = time.time()
    print(params)
    #lambda_act,max_act = params
    sim = setup(gtype)
    # run : 
    cell_track = runsim(sim,200)
    for i,track in enumerate(cell_track):
        fname = 'CELL' + str(i) +  gtype + str(iter)
        np.savetxt('../data/STROMAL_ACT4/'+fname+'.txt',track)
    print('computed : ',params, 'in ',time.time() - t1)

if __name__ == "__main__":

    jobs = []
    inpts = [('GM',0),('GM',1),('NOFRC',0),('NOFRC',1)]
    
    print(len(inpts))
    for inpt in inpts:
        run_grid_point(inpt)
        #p = multiprocessing.Process(target = run_grid_point,args = (inpt,))
        #jobs.append(p)
        #p.start()
    #for p in jobs:
    #    p.join()
