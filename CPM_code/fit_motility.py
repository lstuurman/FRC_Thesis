import numpy as np
import cpm
from CPM_helpers1 import *
from MSD1cell import *
from itertools import product
from multiprocessing import Pool
import os
import pickle
import pandas as pd

def run_sim(params):
    lambda_act,max_act = params
    # setup : 
    #begin_pos = np.random.randint(0,256,size=3)
    begin_pos = [128,128,128]
    # data to return : 
    data = []
    print(params)
    for i in range(3):
        simulation = single_cell_setup1(256,begin_pos,lambda_act,max_act,FRC=False)
        # run : 
        volume_track,cell_track = run_sim_1cell(simulation,2000)
        cell_track = handle_boundaries(cell_track,pr = False)
        # scanned volume : 
        vol = scanned_volume(volume_track)
        print('Percentage scanned volume : ', vol)
        displ, angles = analyse_track(cell_track)
        # MSD
        delta_t,MSD = compute_MSD(displ)
        # autocorr : 
        t,auto_corr,pvalues = compute_AutoCorr_WRONG(angles)
        t,new_angles = compute_autocorrelaton(displ)
        
        popt = fit_Motilty(delta_t,MSD)
        data.append([popt[0],popt[1],vol,volume_track,cell_track,auto_corr,pvalues,new_angles])

    return data

def gridsearch():
    ### runn single T-cell simulations for with different combinations of params
    ### to find some good parameters for cell motility and persistance
    # input : 
    l_act = np.geomspace(500,10000,num=10,dtype=int)
    max_act = np.geomspace(10,5000,num = 10,dtype=int)
    inputs = [(x[0],x[1]) for x in product(l_act,max_act)]
    # run in parallel : 
    cpus = os.cpu_count() - 1
    print('Using ',cpus,'cores')
    p = Pool(cpus)
    output = np.array(p.map(run_sim,inputs))
    p.close()
    p.join()
    # save data :
    pickle.dump(output,open('testdat/raw_nofrc.pkl'))

    data = {'Motility':[],'Persistance':[],'Scanned_volume':[],'Lambda_act':[],'Max_act':[]}
    for i,tup in enumerate(inputs):
        data = output[i]
        data['Scanned_volume'].append(np.average(data[:,2]))
        data['Motility'].append(np.average(data[:,0]))
        data['Persistance'].append(np.average(data[:,1]))
        data['Lambda_act'].append(tup[0])
        data['Max_act'].append(tup[1])
		
	
	# save data as csv :
    df = pd.DataFrame(data)
    df.to_csv('testdat/no_frc1.csv')

if __name__ == "__main__":
    gridsearch()