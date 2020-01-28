### 
import numpy as np
import networkx as nx
import pickle as pkl
import scipy.optimize
import time
from find_EV_relation_geometric_graph import return_radius
from find_EV_relation_geometric_graph import polynomial
from numpy import linalg
from Generate_graphs import *
import multiprocessing

# create functions for every type of graphs that saves the graph :

def save_ER(VE_ratios,N_nodes,repetitions):
    ER_data = gen_ER_graphs(VE_ratios,N_nodes,repetitions)
    pkl.dump(ER_data,open('../data/exp1/exp1_ER.pkl','wb'))

def save_WS(VE_ratios,N_nodes,repetitions):
    WS_data = gen_WS_graphs(VE_ratios,N_nodes,repetitions)
    pkl.dump(WS_data,open('../data/exp1/exp1_WS.pkl','wb'))

def save_power(VE_ratios,N_nodes,repetitions):
    power_data = gen_power_cluster_graphs(VE_ratios,N_nodes,repetitions)
    pkl.dump(power_data,open('../data/exp1/exp1_power.pkl','wb'))

def save_geom(VE_ratios,N_nodes,repetitions):
    geom_data = generate_geometric_graphs(VE_ratios,N_nodes,repetitions)
    pkl.dump(geom_data,open('../data/exp1/exp1_geom.pkl','wb'))

def save_BA(VE_ratios,N_nodes,repetitions):
    BA_data = gen_BA_graphs(VE_ratios,N_nodes,repetitions)
    pkl.dump(BA_data,open('../data/exp1/exp1_BA.pkl','wb'))

if __name__ == "__main__":

    global repetitions
    repetitions = 10
    VE_ratios = np.linspace(.125,.375,5)
    N_nodes = [100,176,500,1000]

    proces_1 = multiprocessing.Process(name = 'ER', target = save_ER,args = (VE_ratios,N_nodes,repetitions))
    proces_2 = multiprocessing.Process(name = 'WS', target = save_WS,args = (VE_ratios,N_nodes,repetitions))
    proces_3 = multiprocessing.Process(name = 'power', target = save_power,args = (VE_ratios,N_nodes,repetitions))
    proces_4 = multiprocessing.Process(name = 'geom', target = save_geom,args = (VE_ratios,N_nodes,repetitions))
    proces_5 = multiprocessing.Process(name = 'BA', target = save_BA,args = (VE_ratios,N_nodes,repetitions))

    proces_1.start()
    proces_2.start()
    proces_3.start()
    proces_4.start()
    proces_5.start()
