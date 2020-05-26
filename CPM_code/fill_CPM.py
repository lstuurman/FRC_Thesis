import cpm
import numpy as np
import pandas as pd
import time
import random
from itertools import product
from multiprocessing import Pool
from CPM_helpers1 import real_cofmass
import os

def handle_boundaries(cell_track,pr = False):
    # look for boundary crossings in any
    # of the coordinates
    cell_track2 = cell_track.copy()
    for i in range(len(cell_track) - 1):
        dif = np.subtract(cell_track[i],cell_track[i+1])
        for j,coordinate in enumerate(dif):
            if coordinate > 16:
                # went over boundary from 256 -> 0
                if pr:
                    print('Jumped from :',cell_track[i],'to :',cell_track[i+1])
                    print('Adding ',256, ' to rest of cell track') #cell_track[i,j]
                    print('changed axis : ',j)
                    print('Old coordinat : ',cell_track[i])

                cell_track2[:i + 1,j] -= 32

                if pr:
                    print('New coordinate : ',cell_track[i])
                    print(i,j)

            elif coordinate < -16:
                # form 0 -> 256
                if pr:
                    print('Jumped from :',cell_track[i],'to :',cell_track[i+1])
                    print('Adding ', 256, ' to previous of cell track') 
                    print('Old coordinat : ',cell_track[i])

                cell_track2[:i + 1,j] += 32

                if pr:
                    print('New coordinate : ',cell_track[i])
                    print(i,j)
    return cell_track2


def setup(dens):
        ### SET UP CPM ###
    # params from Nino: multiplicated the adhesion engergies by 10
    # and because lambda of .1 not possible here. 
    # params suitable for single cell in empty space
    dimension = 32
    number_of_types = 2
    temperature = 20

    # initialize : 

    simulation = cpm.Cpm3d(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 1,target_area = 150, lambda_area=25)
    simulation.set_constraints(cell_type = 1, lambda_perimeter = .2, target_perimeter = 1500) #8600
    simulation.set_constraints(cell_type = 1, lambda_persistence = 3000, persistence_diffusion = .76,persistence_time = 15) # 2500, max_act = 42
    # adhesion ; 
    #simulation.set_constraints(cell_type = 1,other_cell_type = 2,adhesion = -5)
    simulation.set_constraints(cell_type = 1,other_cell_type = 1,adhesion = 10)
    simulation.set_constraints(cell_type = 0,other_cell_type = 1,adhesion = 0)

    ### fill cube with cells
    # number of cells :
    frc_in_cube = simulation.get_state() // 2**24 == 1 
    free_voxels = 32**3  - np.count_nonzero(frc_in_cube)
    n_cells = np.floor(dens * (free_voxels/150))
    # sample random positions : 
    free_indeces = np.where(frc_in_cube == 0.)
    possible_seeds = np.array(list(zip(free_indeces[0],free_indeces[1],free_indeces[2])))
    indx = np.random.randint(len(possible_seeds), size = int(n_cells))
    seeds = possible_seeds[indx,:]

    for c in seeds:
        simulation.add_cell(c[0],c[1],c[2],1)
    
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
    cofmass_track = np.zeros((len(n_cells),iters,3))
    print(cofmass_track.shape)
    #cell_sizes = []
    t0 = time.time()
    for i in range(iters):
        simulation.run(20)
        cell_sizes = []
        #centers = simulation.get_centroids()
        #print(centers)
        for n in n_cells:
            cell = state % 2**24 == n
            # check cell_size : 
            size = np.sum(cell)
            #print(size)
            cell_sizes.append(size)
            if size <100:
                #print('to small')
                continue
            cofmass_track[n,i] = np.array(real_cofmass(cell,32,pr = False))
            #print(cofmass_track[n,i])
            #print(n,size)
            #print(cofmass_track[n-2,i])
        if i == 0:
            t1 = time.time() - t0
            print('expected computing time : ',t1 * steps)
    #print(cofmass_track[-1])
    #print(cofmass_track.shape)
    return cofmass_track #,cell_sizes

def run_grid_point(density):
    t1 = time.time()
    os.mkdir("150V_DENS"  +str(density))
    #lambda_act,max_act = params
    # iterate 5 times : 
    cell_tracks = []
    for _ in range(1):
        sim = setup(density)
        # run : 
        cell_track = runsim(sim,500)
        #for t in cell_track:
            #cell_tracks.append(t)
        cell_tracks.append(cell_track[0:])
    for i,track in enumerate(cell_tracks):
        newtrack = handle_boundaries(track)
        #cell_tracks[i] = newtrack
        fname = "150V_DENS"  +str(density) + "/CELL_" + str(i)
        np.savetxt('../data/increase_DENS_PRFDR/'+fname+'.txt',newtrack)
    #print('computed : ',params, 'in ',time.time() - t1)


def gridsearch():
    ### runn multi T-cell simulations for with different combinations of params
    ### to find some good parameters for cell track autocorrelation
    # input : 
    #l_act = np.linspace(1000,5000,num=10,dtype=int)
    #l_act = np.array([500,750,1000,2000,3000,4000,5000])
    l_act = np.array([50,100,200,300,400,500,600,700,800,900,1000,2500,5000,10000,20000])
    #max_act = np.array([10,50,75,100,150,200,500])
    #max_act = np.linspace(1000,5000,num = 5,dtype=int)
    max_act = np.array([0.1,0.2,0.3,0.4,.5,.6,.7,.8,.9,1.0])
    inputs = [(x[0],x[1]) for x in product(l_act,max_act)]
    # run in parallel : 
    cpus = 10 #.cpu_count() - 15
    print('Using ',cpus,'cores')
    p = Pool(cpus)
    output = np.array(p.map(run_grid_point,inputs))
    p.close()
    p.join()

if __name__ == "__main__":
    #sim = setup(2000,20)
    # run : 
    #cell_track = runsim(sim,500)

    gridsearch()
