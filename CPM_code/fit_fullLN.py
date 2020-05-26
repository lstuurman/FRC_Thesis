import numpy as np
import cpm
from CPM_helpers1 import *
from MSD1cell import *
from itertools import product
from multiprocessing import Pool
import os
import pickle
import pandas as pd
import time
from full_ln2 import * 

def handle_boundaries2(cell_track,pr = False):
    # look for boundary crossings in any
    # of the coordinates
    cell_track2 = cell_track.copy()
    for i in range(len(cell_track) - 1):
        dif = np.subtract(cell_track[i],cell_track[i+1])
        for j,coordinate in enumerate(dif):
            if coordinate > 32:
                # went over boundary from 256 -> 0
                if pr:
                    print('Jumped from :',cell_track[i],'to :',cell_track[i+1])
                    print('Adding ',256, ' to rest of cell track') #cell_track[i,j]
                    print('changed axis : ',j)
                    print('Old coordinat : ',cell_track[i])

                cell_track2[:i + 1,j] -= 64

                if pr:
                    print('New coordinate : ',cell_track[i])
                    print(i,j)
                
            elif coordinate < -32:
                # form 0 -> 256
                if pr:
                    print('Jumped from :',cell_track[i],'to :',cell_track[i+1])
                    print('Adding ', 256, ' to previous of cell track') 
                    print('Old coordinat : ',cell_track[i])

                cell_track2[:i + 1,j] += 64

                if pr:
                    print('New coordinate : ',cell_track[i])
                    print(i,j)
    return cell_track2

def run_sim(params):
    t1 = time.time()
    lambda_act,max_act = params
    # iterate 5 times : 
    cell_tracks = []
    for _ in range(10):
        sim = setup(lambda_act,int(max_act))
        # run : 
        cell_track = runsim(sim,500)
        #for t in cell_track:
            #cell_tracks.append(t)
        cell_tracks.append(cell_track[-1])
    for i,track in enumerate(cell_tracks):
        newtrack = handle_boundaries2(track)
        #cell_tracks[i] = newtrack
        fname = 'LAMBDA_'+str(lambda_act) +'MAX'+str(max_act)+'_' + str(i)
        np.savetxt('../data/FIT_speedy2/150V_singleZOOM3/CELL'+fname+'.txt',newtrack)
    print('computed : ',params, 'in ',time.time() - t1)

    #return data


def gridsearch():
    ### runn multi T-cell simulations for with different combinations of params
    ### to find some good parameters for cell track autocorrelation
    # input : 
    #l_act = np.linspace(1000,5000,num=10,dtype=int)
    #l_act = np.array([500,750,1000,2000,3000,4000,5000])
    #l_act = np.array([500,1000,2500,5000,10000])
    #l_act = np.array([50,100,200,300,400,500,600,700,800,900,1000,2500,5000,10000,20000])
    l_act = np.linspace(50,90,21)
    #max_act = np.array([10,50,75,100,150,200,500])
    #max_act = np.linspace(1000,5000,num = 5,dtype=int)
    max_act = np.linspace(100,1000,10)
    print(l_act,max_act)
    inputs = [(x[0],x[1]) for x in product(l_act,max_act)]
    # run in parallel : 
    cpus = 10 #.cpu_count() - 15
    print('Using ',cpus,'cores')
    p = Pool(cpus)
    output = np.array(p.map(run_sim,inputs))
    p.close()
    p.join()
    # save data :
   #pickle.dump(output,open('testdat/raw_nofrc.pkl','wb'))

# data = #{'Motility':[],'Persistance':[],'Scanned_volume':[],'Lambda_act':[],'Max_act':[]}
    # for i,tup in enumerate(inputs):
    #     data = output[i]
    #     data['Scanned_volume'].append(np.average(data[:,2]))
    #     data['Motility'].append(np.average(data[:,0]))
    #     data['Persistance'].append(np.average(data[:,1]))
    #     data['Lambda_act'].append(tup[0])
    #     data['Max_act'].append(tup[1])
		
	
	# # save data as csv :
    # df = pd.DataFrame(data)
    # df.to_csv('testdat/no_frc1.csv')

if __name__ == "__main__":
    gridsearch()
    #run_sim((2000,1500))
    #l_act = np.linspace(1000,10000,num=10,dtype=int)
    #max_act = np.linspace(50,1000,num = 10,dtype=int)
    #inputs = [(x[0],x[1]) for x in product(l_act,max_act)]
    #run_sim(inputs[0])
