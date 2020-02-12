### Test CPM with and without FRC ###

import cpm
from CPM_helpers1 import *
from Bresenheim import *

def cell_types(s_state):
    return len(np.unique(s_state))

def tcell_indeces(s_state):
    return np.where(s_state == 2)

def nt_cells(s_state):
    i = np.where(s_state == 2)
    return len(np.unique(s_state[i]))

def test_without_frc():
    simulation = single_cell_setup1()
    s = simulation.get_state()
    cellids = s % 2**24
    celltypes = s // 2**24
    # print number of cells and sizes of cells before running : 
    n_types = cell_types(celltypes)
    t_indeces = np.where(celltypes == 2)
    n_tcells = nt_cells(celltypes)
    vol_tcells = np.sum(celltypes != False)
    print('\nbefore running inital MC steps : ')
    print('number of cell types : ', n_types)
    print('number of t cells : ', n_tcells)
    print('volume of t cells : ', vol_tcells)

    # run
    simulation.run(10)
    print('\nafter running 10 inital MC steps : ')
    print('number of cell types : ', cell_types(s // 2**24))
    print('number of t cells : ', nt_cells(s // 2**24))
    print('volume of t cells : ', np.sum(s // 2**24 != False))

def test_with_frc():
    simulation = single_cell_setup1()
    frc = test(256,plotly=False,mayavi=False)
    frc = np.copy(frc.astype(np.uint32))
    frc[frc==1] += 2**24
    simulation.initialize_from_array(frc,0)
    s = simulation.get_state()
    cellids = s % 2**24
    celltypes = s // 2**24
    # print number of cells and sizes of cells before running : 
    n_types = cell_types(celltypes)
    t_indeces = np.where(celltypes == 2)
    n_tcells = nt_cells(celltypes)
    vol_tcells = np.sum(celltypes != False)
    print('\nbefore running inital MC steps : ')
    print('number of cell types : ', n_types)
    print('number of t cells : ', n_tcells)
    print('volume of t cells : ', vol_tcells)
    exit()
    # run
    simulation.run(10)
    print('\nafter running 10 inital MC steps : ')
    print('number of cell types : ', cell_types(s // 2**24))
    print('number of t cells : ', nt_cells(s // 2**24))
    print('volume of t cells : ', np.sum(s // 2**24 != False))

#test_without_frc()
test_with_frc()