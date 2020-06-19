import cpm
import numpy as np
from CPM_helpers1 import real_cofmass
#from CPM_helpers1 import *
#from MSD1cell import *
import pickle
import time
import random
import networkx as nx
from itertools import product
from multiprocessing import Pool
from Bresenheim import nodesInCube,fill_cube,adjust_thickness
# Fit act model fully pupulated with FRC structure : 

def setup_frc(D):
    # set dimension for grid : 
    D = 64
    r = 20/64# 20microns

    # create graphs with positions : 
    g = nx.random_geometric_graph(170,r,dim = 3)
    positions = []
    for n,data in g.nodes(data = True):
        positions.append(data['pos'])

    # bin possitions and place in grid : 
    g,_ = nodesInCube(g,positions,D)
    c = fill_cube(D,g)
    cube = adjust_thickness(c,1)


    print('percentage of volume occupied by frc',(np.sum(cube)/D**3)*100)

    # save as pickle :
    pickle.dump(cube.astype(np.uint32),open('../data/FRCs/GM64_diam3.pkl','wb'))


def setup(l_act,m_act):
        ### SET UP CPM ###
    dimension = 64
    number_of_types = 3
    temperature = 200

    # initialize : 

    simulation = cpm.Cpm3d(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 2,target_area = 150, lambda_area=25)
    simulation.set_constraints(cell_type = 2, lambda_perimeter = .2, target_perimeter = 1200) #8600
    simulation.set_constraints(cell_type = 1, lambda_persistence = 2000, persistence_diffusion = .87,persistence_time = 15)
    #simulation.set_constraints(cell_type = 2, lambda_act = int(l_act), max_act = int(m_act)) # 2500, max_act = 42
    # adhesion ; 
    simulation.set_constraints(cell_type = 1,other_cell_type = 2,adhesion = -5)
    simulation.set_constraints(cell_type = 1,other_cell_type = 1,adhesion = 10)
    simulation.set_constraints(cell_type = 0,other_cell_type = 1,adhesion = 0)


    #print('Creating FRC')
    frc_in_cube = pickle.load(open('../data/FRCs/GM64_diam3.pkl','rb'))
    #print('Loading into simulation : ')
    simulation.initialize_from_array(frc_in_cube,1)
    #print('Done')

    ### fill cube with cells
    # number of cells :
    frc_in_cube = simulation.get_state() // 2**24 == 1 
    free_voxels = 64**3  - np.count_nonzero(frc_in_cube)
    n_cells = np.floor(1 * (free_voxels/150))

    # sample random positions : 
    free_indeces = np.where(frc_in_cube == 0.)
    possible_seeds = np.array(list(zip(free_indeces[0],free_indeces[1],free_indeces[2])))
    indx = np.random.randint(len(possible_seeds), size = int(n_cells))
    seeds = possible_seeds[indx,:]
    for c in seeds:
        simulation.add_cell(c[0],c[1],c[2],2)
    
    #print('number of cells : ', len(seeds))
    s = simulation.get_state()

    celltypes = s % 2**24
    t = []
    n_cells = np.unique(celltypes)

    # little warmup run 
    simulation.run(50)

    for n in n_cells:
        #print(n)
        celltypes = s % 2**24 == n
        size = np.sum(celltypes)
        t.append(size)
        #print(n,size)

    print('free voxels : ',64**3 - sum(t))
    return simulation


def runsim(simulation,steps):
    state = simulation.get_state() 

    ids = state % 2**24
    n_cells = np.unique(ids)

    iters = int(steps) # /10
    cofmass_track = np.zeros((len(n_cells)-2,iters,3))
    print(cofmass_track.shape)
    print(len(n_cells))
    print(n_cells[2:])
    print([i for i in range(n_cells[-1]) if i not in n_cells])
    t0 = time.time()
    for i in range(iters):
        simulation.run(50)
        cell_sizes = []
        for ind,n in enumerate(n_cells[2:]):
            cell = state % 2**24 == n
            # check cell_size : 
            size = np.sum(cell)
            #print(size)
            cell_sizes.append(size)
            if size <100:
                #print('to small')
                continue
            cofmass_track[ind,i] = np.array(real_cofmass(cell,64,pr = False))

        if i == 0:
            t1 = time.time() - t0
            print('expected computing time : ',t1 * steps)
        print(ind,i)
    return cofmass_track

def run_grid_point(params):
    t1 = time.time()
    print(params)
    lambda_act,max_act = params
    sim = setup(lambda_act,max_act)
    # run : 
    cell_track = runsim(sim,200)
    for i,track in enumerate(cell_track):
        fname = 'LAMBDA_'+str(lambda_act) +'MAX'+str(max_act)+'_' + str(i)
        np.savetxt('../data/FITFULL_ACT_PRFDR/thin64_1/CELL'+fname+'.txt',track)
    print('computed : ',params, 'in ',time.time() - t1)


def gridsearch():
    ### runn multi T-cell simulations for with different combinations of params
    ### to find some good parameters for cell track autocorrelation
    l_act = np.linspace(2000,4000,11)
    max_act = np.linspace(10,100,10)
    inputs = [(x[0],x[1]) for x in product(l_act,max_act)]
    
    #for inp in inputs[-1:]:
    #    run_grid_point(inp)

    # run in parallel : 
    cpus = 12 #.cpu_count() - 15
    print('Using ',cpus,'cores')
    p = Pool(cpus)
    output = np.array(p.map(run_grid_point,inputs))
    p.close()
    #p.join()

if __name__ == "__main__":
    #sim = setup(2000,20)
    # run : 
    #cell_track = runsim(sim,500)

    gridsearch()
    #setup_frc(64)


