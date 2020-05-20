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
    temperature = 20

    # initialize : 

    simulation = cpm.Cpm(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 1,target_area = 150, lambda_area=25)
    simulation.set_constraints(cell_type = 1, lambda_perimeter = 20, target_perimeter = 1400)#8600
    # adhesion ; 
    simulation.set_constraints(cell_type = 1,other_cell_type = 1,adhesion = 10)
    simulation.set_constraints(cell_type = 1,other_cell_type = 0,adhesion = 0)

    free_voxels = 64**3
    max_cells = free_voxels/150
    # data to save : 
    rows = []
    size_matrix = []
    # iteratively seed cells : 
    for n in range(int(max_cells)):
        # get sim state and occupied voxels: 
        full_voxels = simulation.get_state() // 2**24 == 1
        free_indeces = np.where(full_voxels == 0.)
        # list possible empty spaces
        possible_seeds = np.array(list(zip(free_indeces[0],free_indeces[1],free_indeces[2])))
        if len(possible_seeds) < 100:
            rows.append([n*10,len(possible_seeds)])
            break
        indx = np.random.randint(len(possible_seeds), size = 10)
        #else:
        #    indx = np.random.randint(len(possible_seeds),size = 1)

        seed = possible_seeds[indx,:]
        #print(seed)
        #seed = random.sample(possible_seeds,1)
        # add T cells :
        for c in seed:
            simulation.add_cell(c[0],c[1],c[2],1)
        # burnin sim : 
        simulation.run(50)
        # cell sizes : 
        sizes = []
        cell_ids = simulation.get_state() % 2**24
        n_cells = np.unique(cell_ids)

        for i_d in n_cells:
            celltypes = simulation.get_state % 2**24 == i_d
            size = np.sum(celltypes)
            sizes.append(size)
        # save data : 
        rows.append([n*10,len(possible_seeds),sum(sizes)])
        size_matrix.append(sizes)

        print(rows[-1])
        print(sizes)
    df = pd.DataFrame(data = rows, columns= ['Seeded_cells','Free_voxels','Cell_volume'])
    df.to_csv('increasing_pressureV150.csv')
    nparr = np.array(size_matrix)
    np.savetxt('size_matrix.txt',nparr)

seed_cpm()

