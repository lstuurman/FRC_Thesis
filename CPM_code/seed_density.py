import numpy as np 
import pandas as pd 
import glob 
import cpm  
import random
### SCRIPT FOR EXPERIMENT TO SEE RELATION BETWEEN CELLS SEEDED AND ACTUAL DENSITY

def seed_cpm():
    # standard cpm setup : 
    dimension = 64
    number_of_types = 2
    temperature = 70

    # initialize : 

    simulation = cpm.Cpm(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 2,target_area = 150, lambda_area=250)
    simulation.set_constraints(cell_type = 2, lambda_perimeter = 20, target_perimeter = 1400)#8600
    # adhesion ; 
    simulation.set_constraints(cell_type = 2,other_cell_type = 1,adhesion = 100)
    simulation.set_constraints(cell_type = 2,other_cell_type = 0,adhesion = 0)

    free_voxels = 64**3
    max_cells = free_voxels/150
    # data to save : 
    rows = []
    # iteratively seed cells : 
    for n in range(max_cells):
        # get sim state and occupied voxels: 
        full_voxels = simulation.get_state() // 2**24 == 2
        free_indeces = np.where(full_voxels == 0.)
        # list possible empty spaces
        possible_seeds = np.array(list(zip(free_indeces[0],free_indeces[1],free_indeces[2])))
        # indx = np.random.randint(len(possible_seeds), size = int(n_cells))
        # seeds = possible_seeds[indx,:]
        seed = random.sample(possible_seeds,1)
        # add T cell : 
        simulation.add_cell(seed[0],seed[1],seed[2],2)
        # burnin sim : 
        sim.run(200)

        # save data : 
        rows.append([n,len(possible_seeds)])

    df = pd.DataFrame(data = rows, columns= ['Seeded_cells','Free_voxels'])
    df.to_csv('increasing_pressure.csv')