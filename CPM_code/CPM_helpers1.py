import numpy as np
from Bresenheim import *
import cpm
import matplotlib.pyplot as plt
from scipy.ndimage.measurements import center_of_mass
from scipy.ndimage import label , generate_binary_structure
from scipy.spatial.distance import euclidean
from numpy.linalg import norm
from mpl_toolkits import mplot3d
from mayavi import mlab

#source env/bin/activate

def single_cell_setup1(dim,position,l_act,m_act,FRC = False):
    ### SET UP CPM ###
    # params from Inge: multiplicated the adhesion engergies by 10
    # and because lambda of .1 not possible here. 
    # params suitable for single cell in empty space
    dimension = dim
    number_of_types = 3
    temperature = 7

    # initialize : 

    simulation = cpm.Cpm(dimension, number_of_types, temperature)
    # LAmbdas ; 
    simulation.set_constraints(cell_type = 2,target_area = 150, lambda_area=250)
    simulation.set_constraints(cell_type = 2, lambda_perimeter = 20, target_perimeter = 1500)#8600
    simulation.set_constraints(cell_type = 2, lambda_act = l_act, max_act = m_act) # 2500, max_act = 42
    # adhesion ; 
    simulation.set_constraints(cell_type = 1,other_cell_type = 2,adhesion = 1)
    simulation.set_constraints(cell_type = 2,other_cell_type = 2,adhesion = 150)
    simulation.set_constraints(cell_type = 2,other_cell_type = 0,adhesion = 5)
    #simulation.initialize_from_array(cube_with_type,1)
    if FRC:
        print('Creating FRC')
        frc_in_cube = test(dimension,plotly=False,mayavi=False)
        print('Loading into simulation : ')
        simulation.initialize_from_array(frc_in_cube,1)
        print('Done')

    #simulation.add_cell(128,128,128,2)
    simulation.add_cell(position[0],position[1],position[2],2)
    #simulation.add_cell(256,256,256,2)
    
    # run untill cell has approx target area : 
    s = simulation.get_state()
    celltypes = s // 2**24
    t = celltypes == 2
    
    while np.sum(t) < 140:
        celltypes = s // 2**24
        t = celltypes == 2
        simulation.run(1)
        print(np.sum(t))
    return simulation

def dist_to_axis(x,ref):
    # return distance between point x and reference
    distances = []
    for i,axis in enumerate(x):
        d = axis - ref[i]
        if d < 128 and d > -128:
            distances.append(d)
        elif d < 0:
            # other patch on the other side of box
            distances.append(d + 256)
        else : 
            distances.append(d - 256)
    return distances



def real_cofmass(cell,pr = False):
    """ cell = simulation.get_state() % 2**24 == id """
    neighborhood =  generate_binary_structure(3,4)
    labels, num_features = label(cell, structure=neighborhood)
     # compute sizes of pieces of cell: (bigger then 1)
    sizes = [(np.count_nonzero(labels == i),i) for i in range(1,num_features + 1) if np.count_nonzero(labels == i) > 1]
    # Just one patch : 
    if len(sizes) == 1:
        return center_of_mass(cell)

    # multiple patches : 
    sorted_sizes = sorted(sizes,key = lambda x: x[0])[::-1]
    big_labels = [x[1] for x in sorted_sizes]
    masses = center_of_mass(cell,labels = labels, index = big_labels)
    biggest_patch = np.array(masses[0]) #* (sorted_sizes[0][0] / total_size)
    biggest_size = sorted_sizes[0][0]

    if pr:
        print(sorted_sizes)
        print(masses)
    # take distances between boundaries of every axis : 
    dists = list(map(dist_to_axis,masses[1:],[biggest_patch for i in range(len(masses) - 1)]))
    # loop over distances and add or destract from main blob :

    for i,d in enumerate(dists):
        print('boundary')
        biggest_size += sorted_sizes[i + 1][0]
        if pr:
            print('Biggest patch at : ', biggest_patch)
            print('Added coordate :', masses[i + 1])
            print('With difference : ',  np.array(d))
            print('with weight : ', (sorted_sizes[i + 1][0] / biggest_size)) #biggest_size
            print(np.array(d) * (sorted_sizes[i + 1][0] / biggest_size))# biggest_size
        
        biggest_patch += np.array(d) * (sorted_sizes[i + 1][0] / biggest_size)# + biggest_size +
        
        # correct negative values : 
        for i in range(3):
            if biggest_patch[i] < 0:
                biggest_patch[i] += 256
            elif biggest_patch[i] > 256.:
                biggest_patch[i] -= 256
        if pr:
            print('NEW biggest patch : ', biggest_patch)
    #biggest_patch = [x if x>0 and x < 256. elif x < 0 + 256 else x - 256 for x in biggest_patch]
    return biggest_patch

def run_sim_1cell(simulation,steps):
    state = simulation.get_state() 
    act_state = simulation.get_act_state()
    dimension = state.shape[0]

    iters = int(steps/10) # /10
    volume_track = np.zeros(tuple([dimension] * 3))
    cofmass_track = np.empty((iters,3))
    for i in range(iters):
        simulation.run(10)
        cell = state // 2**24 == 2
        ids = state % 2**24 == 2
        n_cells = np.unique(ids)
        if len(n_cells) > 2:
            print('cell split up..')
        cofmass_track[i] = np.array(real_cofmass(cell, pr = False))

        print(cofmass_track[i])
        if i%100 == 0:
            volume_track = volume_track + act_state
    print(cofmass_track.shape)
    return volume_track,cofmass_track

def scanned_volume(volume_track):
    dim = volume_track.shape[0]
    return len(np.where(volume_track != 0.)[0]) / dim**3

def handle_boundaries(cell_track,pr = False):
    # look for boundary crossings in any
    # of the coordinates
    cell_track2 = cell_track.copy()
    for i in range(len(cell_track) - 1):
        dif = np.subtract(cell_track[i],cell_track[i+1])
        for j,coordinate in enumerate(dif):
            if coordinate > 128:
                # went over boundary from 256 -> 0
                if pr:
                    print('Jumped from :',cell_track[i],'to :',cell_track[i+1])
                    print('Adding ',256, ' to rest of cell track') #cell_track[i,j]
                    print('changed axis : ',j)
                    print('Old coordinat : ',cell_track[i])

                cell_track2[:i + 1,j] -= 256

                if pr:
                    print('New coordinate : ',cell_track[i])
                    print(i,j)
                
            elif coordinate < -128:
                # form 0 -> 256
                if pr:
                    print('Jumped from :',cell_track[i],'to :',cell_track[i+1])
                    print('Adding ', 256, ' to previous of cell track') 
                    print('Old coordinat : ',cell_track[i])

                cell_track2[:i + 1,j] += 256

                if pr:
                    print('New coordinate : ',cell_track[i])
                    print(i,j)
    return cell_track2
            

def analyse_track(cell_track):
    # total displacement : 
    disp_per_t = []
    # angles : 
    angles= []
    for i in range(len(cell_track) - 1):
        point1 = cell_track[i]
        point2 = cell_track[i + 1]
        disp_per_t.append(euclidean(point1,point2))

        if i != len(cell_track) - 2:
            point3 = cell_track[i + 2]
            v1 = point2 - point1
            v2 = point3 - point2
            cos = np.clip(np.dot(v1,v2)/(norm(v1) * norm(v2)),-1,1)
            #angles.append(np.arccos(np.clip(cos,-1,1)))
            print(cos)
            #if np.isfinite(cos):
            angles.append(cos)

    # # autocorrelation : 
    # corr = np.correlate(cell_track, cell_track, mode='full')
    # auto_correlation =  corr[corr.size/2:]

    return disp_per_t,angles#,auto_correlation

def plot_celltrack(cell_track):
    # plot path of center of mass  :
    x = [c[0] for c in cell_track]
    y = [c[1] for c in cell_track]
    z = [c[2] for c in cell_track]

    #fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x,y,z)
    ax.scatter3D(x[0],y[0],z[0],label = 'start')
    ax.scatter3D(x[-1],y[-1],z[-1],label = 'end')
    ax.legend()
    plt.show()

if __name__ == "__main__":
    simulation = single_cell_setup1(256,2500,42,FRC=False)
    volume_track,cell_track = run_sim_1cell(simulation,500)
    # percentage_scanned = scanned_volume(volume_track)
    # displ, angles = analyse_track(cell_track)
    plot_celltrack(cell_track)
    cell_track = handle_boundaries(cell_track)
    plot_celltrack(cell_track)


    # if any(c < 250 and c > 10 for c in center_of_mass(labels)): #el
    #     # cell not near boundary -> can be regarded as one cell
    #     #sizes = [np.count_nonzero(labels == i) for i in range(1,num_features + 1)]
    #     max_indeces = np.where(sizes == np.amax(sizes))[0]
    #     index = np.random.choice(max_indeces,1,replace = False)
    #     # set smaller patches to zero :
    #     labels[labels != index+1] = 0
    #     return center_of_mass(labels)


# def coffmass2(cell):
#     total_size = np.sum(cell)
#     labels,num_features = label(cell)
#     sizes = [(np.count_nonzero(labels == i),i) for i in range(1,num_features + 1)]
#     masses = center_of_mass(cell,labels = labels)
#     # take distances between boundaries of every axis : 
#     dists = list(map(dist_to_axis,masses))

#     # go over patches : 
#     cofmass = []
#     for i in sizes:
#         label = i[0]
#         size = i[0]
#         centroid = center_of_mass(cell,labels = labels, index = big_labels)


        # if i > 0:
            #d = euclidean(cofmass_track[i],cofmass_track[i - 1])
            # if d > 10.:
            #     real_cofmass(simulation.get_state() % 2**24 == 1, pr = True)
            #     mlab.clf()
            #     mlab.contour3d(simulation.get_state())
            #     mlab.show()