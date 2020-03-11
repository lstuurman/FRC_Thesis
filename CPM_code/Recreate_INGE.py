import cpm
import numpy as np 
from CPM_helpers1 import real_cofmass
import pickle

def INGE_setup(dim = 1028,position = (512,512,512),l_act = 40,m_act = 50):
    ### SET UP CPM ###
    # params from Inge: multiplicated the adhesion engergies by 10 
    # Sample every 5 steps: 
    # and because lambda of .1 not possible here. 
    # params suitable for single cell in empty space
    dimension = dim
    number_of_types = 3
    temperature = 7

    # initialize : 

    simulation = cpm.Cpm(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 2,target_area = 1800, lambda_area=25)
    simulation.set_constraints(cell_type = 2, lambda_perimeter = 0.01, target_perimeter = 8600)#8600
    simulation.set_constraints(cell_type = 2, lambda_act = l_act, max_act = m_act) # 2500, max_act = 42
    # adhesion ; 
    simulation.set_constraints(cell_type = 2,other_cell_type = 2,adhesion = 15)
    simulation.set_constraints(cell_type = 2,other_cell_type = 0,adhesion = 5)

    # add cell :
    simulation.add_cell(position[0],position[1],position[2],2)

    # burnin time of 500 stpes : 
    simulation.run(500)

    return simulation

def run_inge(simulation,steps = 1000):
    state = simulation.get_state() 
    act_state = simulation.get_act_state()
    dimension = state.shape[0]

    iters = int(steps) # /10
    volume_track = np.zeros(tuple([dimension] * 3))
    cofmass_track = np.empty((iters,3))
    for i in range(iters):
        simulation.run(5)
        cell = state // 2**24 == 2
        ids = state % 2**24 == 2
        n_cells = np.unique(ids)
        if len(n_cells) > 2:
            print('cell split up..')
        cofmass_track[i] = np.array(real_cofmass(cell, pr = False))

        if i%100 == 0:
            volume_track = volume_track + act_state

    # save : 
    np.savetxt('testdat/INGE_track.txt',cofmass_track)
    pickle.dump(volume_track,open('INGE_volume.pkl','wb'))

if __name__ == "__main__":
    sim = INGE_setup()
    run_inge(sim)