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
    simulation.set_constraints(cell_type = 2,target_area = 500, lambda_area=25)
    simulation.set_constraints(cell_type = 2, lambda_perimeter = 20, target_perimeter = 3500)#8600
    #simulation.set_constraints(cell_type = 2, lambda_act = 2000, max_act = 100)
    # adhesion ; 
    simulation.set_constraints(cell_type = 2,other_cell_type = 1,adhesion = 100)
    simulation.set_constraints(cell_type = 2,other_cell_type = 0,adhesion = 0)

    free_voxels = 64**3
    max_cells = free_voxels/500
    # data to save : 
    rows = []
    # iteratively seed cells : 
    for n in range(int(max_cells)):
        # get sim state and occupied voxels: 
        full_voxels = simulation.get_state() // 2**24 == 2
        free_indeces = np.where(full_voxels == 0.)
        # list possible empty spaces
        possible_seeds = np.array(list(zip(free_indeces[0],free_indeces[1],free_indeces[2])))
        if len(possible_seeds) < 500:
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
            simulation.add_cell(c[0],c[1],c[2],2)
        # burnin sim : 
        simulation.run(50)

        # save data : 
        rows.append([n*10,len(possible_seeds)])
        print(rows[-1])
    df = pd.DataFrame(data = rows, columns= ['Seeded_cells','Free_voxels'])
    df.to_csv('increasing_pressureV500.csv')

seed_cpm()

