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
from Generate_graphs import ER_p_value,WS_K_value,M_power_cluster,E_M_relation
from graph_drawing_algs import fruchterman_reingold,normalize_postions

# Fit act model fully pupulated with FRC structure : 

def setup_frc(D,g_type):
    # set dimension for grid : 
    D = 64
    r = 20/64# 20microns

    print('Generating graph')
    if g_type == 'ER':
        p = ER_p_value(270,.25)
        g = nx.erdos_renyi_graph(270,p)#7000
        g = fruchterman_reingold(g,20)
    elif g_type == 'BA':
        g = nx.barabasi_albert_graph(220,4)#700
        g = fruchterman_reingold(g,20)
    elif g_type == 'WS':
        k = WS_K_value(600,.25)
        p = 0.027825594022071243
        g = nx.watts_strogatz_graph(600,k,p)#1000
        g = fruchterman_reingold(g,20)
    elif g_type == 'PW':
        m = M_power_cluster(310,.25)
        p = 0.666666666666666
        g = nx.powerlaw_cluster_graph(310,m,p)# 700
        g = fruchterman_reingold(g,20)
    elif g_type == 'GM':
        r = 20/dim # 20microns
        g = nx.random_geometric_graph(170,r,dim = 3)#4500
    
    positions = []
    for n,data in g.nodes(data = True):
        positions.append(data['pos'])

    positions = normalize_postions(positions)
    print('Filling cube with graph') 
    if g_type == 'GM' or g_type == 'WS':
        g,_ = nodesInCube(g,positions,D)
        c = fill_cube(D,g)
        #adjust thickness so that it fills +- 17% of cube : 
        cube = adjust_thickness(c,1)
    else:
        g,_ = nodesInCube(g,positions,4 * D)
        c = fill_cube(4 * D,g)
        #adjust thickness so that it fills +- 17% of cube : 
        thicckk = adjust_thickness(c,1)
        cube = thicckk[96:-96,96:-96,96:-96]


    print('percentage of volume occupied by frc',(np.sum(cube)/D**3)*100)

    # save as pickle :
    dfile = '../data/FRCs/' + g_type + '64_diam3.pkl'
    pickle.dump(cube.astype(np.uint32),open(dfile,'wb'))


def setup(g_type):
        ### SET UP CPM ###
    dimension = 64
    number_of_types = 3
    temperature = 200

    # initialize : 

    simulation = cpm.Cpm3d(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 2,target_area = 150, lambda_area=25)
    simulation.set_constraints(cell_type = 2, lambda_perimeter = .2, target_perimeter = 1200) #8600
    #simulation.set_constraints(cell_type = 2, lambda_act = 3000, max_act = 40) # 2500, max_act = 42
    simulation.set_constraints(cell_type = 2, lambda_persistence = 600, persistence_diffusion = .9,persistence_time = 15)
    # adhesion ; 
    simulation.set_constraints(cell_type = 2,other_cell_type = 1,adhesion = -5)
    simulation.set_constraints(cell_type = 2,other_cell_type = 2,adhesion = 10)
    simulation.set_constraints(cell_type = 2,other_cell_type = 0,adhesion = 0)

    simulation.set_constraints(cell_type = 1, fixed=1)
    #print('Creating FRC')
    try:
        dfile = '../data/FRCs/' + g_type + '64_diam3.pkl'
        frc_in_cube = pickle.load(open(dfile,'rb'))
        #print('Loading into simulation : ')
        simulation.initialize_from_array(frc_in_cube,1)
        #print('Done')
    except:
        pass

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
        simulation.run(100)
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

def run_grid_point(gtype):
    t1 = time.time()
    #print(params)
    #lambda_act,max_act = params
    for iter in range(3,5):
        sim = setup(gtype)
    # run : 
        cell_track = runsim(sim,200)
        for i,track in enumerate(cell_track):
            fname = 'CELL' + str(i) +  gtype[1:] + str(iter)
            np.savetxt('../data/STROMAL_PRFDR3/'+fname+'.txt',track)
    print('computed : ',params, 'in ',time.time() - t1)


def gridsearch():
    ### runn multi T-cell simulations for with different combinations of params
    ### to find some good parameters for cell track autocorrelation
    #l_act = np.linspace(2000,4000,11)
    #max_act = np.linspace(10,100,10)
    #inputs = [(x[0],x[1]) for x in product(l_act,max_act)]
    inputs = ['5WS','2PW','0GM'] # '5ER','4BA','NOFRC'
    #for inp in inputs[-1:]:
    #    run_grid_point(inp)

    # run in parallel : 
    cpus = 10 #.cpu_count() - 15
    print('Using ',cpus,'cores')
    p = Pool(cpus)
    output = np.array(p.map(run_grid_point,inputs))
    p.close()
    #p.join()

if __name__ == "__main__":
    #sim = setup(2000,20)
    # run : 
    #cell_track = runsim(sim,500)
    #g_types = ['ER', 'WS'] # 'BA'
    #for graph_type in g_types: 
    #    setup_frc(64,graph_type)
    run_grid_point('0GM')
    #setup_frc(64)

