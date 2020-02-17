import numpy as np 
import cpm
import matplotlib.pyplot as plt
from CPM_helpers1 import *
# from scipy.ndimage.measurements import center_of_mass
# from scipy.ndimage import label
# from scipy.spatial.distance import euclidean
# from numpy.linalg import norm
# from mpl_toolkits import mplot3d

def compute_MSD(displacement):
    # range of sliding windows from 1 to steps - 100
    delta_t = np.arange(1,len(displacement) - 10)
    MSD = np.zeros(len(delta_t))
    for i,dt in enumerate(delta_t):
        # convolve to get displacements : 
        moving_av = np.convolve(displacement, np.ones(dt), 'valid')
        MSD[i] = np.average(np.square(moving_av))
    return delta_t, MSD

def MSD_lonelycell():
    # setup : 
    simulation = single_cell_setup1()
    # run : 
    volume_track,cell_track = run_sim_1cell(simulation,500)
    cell_track = handle_boundaries(cell_track)
    print('Percentage scanned volume : ', scanned_volume(volume_track))
    # analyse
    cell_track = handle_boundaries(cell_track)
    print(cell_track)
    displ, angles = analyse_track(cell_track)
    delta_t,MSD = compute_MSD(displ)
    plt.plot(delta_t,MSD)
    plt.show()
    plot_celltrack(cell_track)

if __name__ == "__main__":
    MSD_lonelycell()