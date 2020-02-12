import numpy as np
from Bresenheim import *
import cpm
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import center_of_mass
from scipy.spatial.distance import euclidean
from numpy.linalg import norm
from mpl_toolkits import mplot3d

def single_cell_setup1():
    ### SET UP CPM ###
    # params from Inge: multiplicated the adhesion engergies by 10
    # and because lambda of .1 not possible here. 
    dimension = 256
    number_of_types = 3
    temperature = 7

    # initialize : 

    simulation = cpm.Cpm(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 2,target_area = 1800, lambda_area=250)
    simulation.set_constraints(cell_type = 2, lambda_perimeter = 100, target_perimeter = 1000)
    simulation.set_constraints(cell_type = 2, lambda_act = 200, max_act = 42) # 160,40
    # adhesion ; 
    simulation.set_constraints(cell_type = 1,other_cell_type = 2,adhesion = 150)
    simulation.set_constraints(cell_type = 2,other_cell_type = 2,adhesion = 150)
    simulation.set_constraints(cell_type = 2,other_cell_type = 0,adhesion = 50)
    #simulation.initialize_from_array(cube_with_type,1)

    simulation.add_cell(128,128,128,2)
    return simulation

def run_sim(simulation,steps):
    cell_state = simulation.get_state()
    act_state = simulation.get_act_state()
    dimension = cell_state.shape[0]

    volume_track = np.zeros(tuple([dimension] * 3))
    cofmass_track = []

    iters = int(steps/10)
    for i in range(iters):
        simulation.run(10)
        cofmass_track.append(center_of_mass(cell_state))
        if i%10 == 0:
            volume_track = volume_track + act_state

    return volume_track,cofmass_track

def scanned_volume(volume_track):
    dim = volume_track.shape[0]
    return len(np.where(volume_track != 0.)[0]) / dim**3

def analyse_track(cell_track):
    # total displacement : 
    tot_displ = 0
    # angles : 
    angles= []
    for i in range(len(cell_track) - 1):
        point1 = np.array(cell_track[i])
        point2 = np.array(cell_track[i + 1])
        tot_displ += euclidean(point1,point2)

        if i != len(cell_track) - 2:
            point3 = np.array(cell_track[i + 2])
            v1 = point2 - point1
            v2 = point3 - point2
            cos = np.dot(v1,v2)/(norm(v1) * norm(v2))
            angles.append(np.arccos(np.clip(cos,-1,1)))

    # # autocorrelation : 
    # corr = np.correlate(cell_track, cell_track, mode='full')
    # auto_correlation =  corr[corr.size/2:]

    return tot_displ,angles#,auto_correlation

def plot_celltrack(cell_track):
    # plot path of center of mass  :
    x = [c[0] for c in cell_track]
    y = [c[1] for c in cell_track]
    z = [c[2] for c in cell_track]

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x,y,z)
    ax.scatter3D(x[0],y[0],z[0],label = 'start')
    ax.scatter3D(x[-1],y[-1],z[-1],label = 'end')
    ax.legend()
    fig.show()

if __name__ == "__main__":
    simulation = single_cell_setup1()
    volume_track,cell_track = run_sim(simulation,1000)
    percentage_scanned = scanned_volume(volume_track)
    displ, angles = analyse_track(cell_track)
    plot_celltrack(cell_track)