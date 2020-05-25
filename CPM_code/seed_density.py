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

    simulation = cpm.Cpm3d(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 1,target_area = 150, lambda_area=25)
    simulation.set_constraints(cell_type = 1, lambda_perimeter = .2, target_perimeter = 1500)#8600
    simulation.set_constraints(cell_type = 1, lambda_act = 2000,max_act = 10)
    # adhesion ; 
    simulation.set_constraints(cell_type = 1,other_cell_type = 1,adhesion = 10)
    simulation.set_constraints(cell_type = 1,other_cell_type = 0,adhesion = 0)

    free_voxels = 64**3
    max_cells = free_voxels/150
    # data to save : 
    rows = []
    size_matrix = []
    # Also write to file : 
    mat_file = open('matrix_moving.txt','w')
    row_file = open('rows_moving.txt','w')
    s = simulation.get_state()
    # iteratively seed cells : 
    for n in range(int(max_cells)):
        # get sim state and occupied voxels: 
        full_voxels = s // 2**24 == 2
        free_indeces = np.where(full_voxels == 0.)
        # list possible empty spaces
        possible_seeds = np.array(list(zip(free_indeces[0],free_indeces[1],free_indeces[2])))
        if len(possible_seeds) < 100:
            rows.append([n*10,len(possible_seeds)])
            break
        indx = np.random.randint(len(possible_seeds), size = 50)
        #else:
        #    indx = np.random.randint(len(possible_seeds),size = 1)

        seed = possible_seeds[indx,:]
        #print(seed)
        #seed = random.sample(possible_seeds,1)
        # add T cells :
        for c in seed:
            simulation.add_cell(c[0],c[1],c[2],1)
        # burnin sim : 
        simulation.run(100)
        # cell sizes : 
        sizes = []
        cell_ids = s % 2**24
        n_cells = np.unique(cell_ids)

        for i_d in n_cells[1:]:
            celltypes = s % 2**24 == i_d
            size = np.sum(celltypes)
            sizes.append(size)
        # save data : 
        rows.append([n*10,len(possible_seeds),sum(sizes)])
        size_matrix.append(sizes)
        row_to_write = str(n*10) + '\t' + str(len(possible_seeds)) + '\t' + str(sum(sizes)) + '\n'
        row_file.write(row_to_write)
        mat_file.write("\t".join(str(x) for x in sizes) + '\n')
        print(rows[-1])
        print(sizes)
    mat_file.close()
    row_file.close()
    df = pd.DataFrame(data = rows, columns= ['Seeded_cells','Free_voxels','Cell_volume'])
    df.to_csv('increasing_pressureV150_moving.csv')
    # make matrix square :
    length = max(map(len,size_matrix))
    print(size_matrix)
    size_matrix = np.array([xi + [0]*(length - len(xi)) for xi in size_matrix])
    #nparr = np.array(size_matrix)
    np.savetxt('size_matrix_moving.txt',size_matrix)

seed_cpm()

