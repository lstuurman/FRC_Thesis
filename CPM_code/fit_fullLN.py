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

def run_sim(params):
    t1 = time.time()
    lambda_act,max_act = params
    # data to return : 
    data = []
    for i in range(1):
        sim = setup(lambda_act,max_act)
        # run : 
        cell_track = runsim(sim,500)
        cell_track = handle_boundaries(cell_track,pr = False)

        # popt = fit_Motilty(delta_t,MSD)
        fname = 'LAMBDA_'+str(lambda_act) +'MAX'+str(max_act)+str(i)
        np.savetxt('../data/full_LN2/CELL'+fname+'.txt',cell_track)
    print('computed : ',params, 'in ',time.time() - t1)

    #return data


def gridsearch():
    ### runn multi T-cell simulations for with different combinations of params
    ### to find some good parameters for cell track autocorrelation
    # input : 
    l_act = np.linspace(1000,10000,num=10,dtype=int)
    max_act = np.linspace(50,1000,num = 10,dtype=int)
    inputs = [(x[0],x[1]) for x in product(l_act,max_act)]
    # run in parallel : 
    cpus = os.cpu_count() - 20
    print('Using ',cpus,'cores')
    p = Pool(cpus)
    output = np.array(p.map(run_sim,inputs))
    p.close()
    p.join()
    # save data :
    #pickle.dump(output,open('testdat/raw_nofrc.pkl','wb'))

    # data = {'Motility':[],'Persistance':[],'Scanned_volume':[],'Lambda_act':[],'Max_act':[]}
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
