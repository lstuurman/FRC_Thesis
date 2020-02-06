# script to compute omega and sigma small world measurs of graphs 
from Generate_graphs import ER_p_value
from Generate_graphs import WS_K_value
from multiprocessing import Process,Manager,Pool
import pickle
import networkx as nx
import numpy as np
from itertools import product

def gen_equiv_random(N,ve_ratio):
    # generate 1000 random graphs, equivalent to the ones used in the experiments
    p = ER_p_value(N,ve_ratio)
    stats  = [] # list of tuples with (average_clustering,average_path) for every generated network. 
    count  = 0
    while count < 1000:
        g = nx.erdos_renyi_graph(N,p)
        if nx.is_connected(g):
            # average clustering and path length : 
            Cl = nx.average_clustering(g)
            L = nx.average_shortest_path_length(g)
            stats.append((Cl,L))
            count += 1
            #print(count)
    print(N,g.number_of_edges())
    
    return stats



def save_random_graphs(N_list, VE_list):
    data_dict = {}
    for N in N_list:
        for ve in VE_ratios:
            stats = gen_equiv_random(N,ve)
            key = str(N) + 'N_' + str(ve) + 'V/E'
            data_dict[key] = stats
    
    pickle.dump(data_dict,open('../data/helper_data/equivalent_random_stats.pkl','wb'))

def generate_lattice(N_list,VE_list):
    # generate one lattice per N/veratio
    data_dict = {}
    for N in N_list:
        for ve in VE_list:
            K = WS_K_value(N,ve)
            g = nx.watts_strogatz_graph(N,K,0)
            Cl = nx.average_clustering(g)
            key = str(N) + 'N_' + str(ve) + 'V/E'
            data_dict[key] = Cl
    
    pickle.dump(data_dict,open('../data/helper_data/equivalent_lattice_stats.pkl','wb'))
# /home/lau/GIT/FRC_Thesis/data/helper_data

if __name__ == "__main__":
    VE_ratios = np.linspace(.125,.375,5)
    N_nodes = [100,176,500,1000]
    product = product(N_nodes,VE_ratios)
    keys = [str(x[0]) + 'N_' + str(x[1]) + 'V/E' for x in product]
    inputs = [(x[0],x[1]) for x in product]

    p = Pool(5)
    outputs  = p.map(gen_equiv_random,inputs)
    p.close()
    p.join()
    data_dict = {}
    for i in range(len(keys)):
        print(keys[i])
        print(outputs[i])
        data_dict[keys[i]] = outputs[i]
    
    pickle.dump(data_dict,open('../data/helper_data/equivalent_random_stats.pkl','wb'))

    # for key in keys:
    #     print(key)
    #     print('\n')
    # print(len(keys))
    # save_random_graphs(N_nodes,VE_ratios)
    generate_lattice(N_nodes,VE_ratios)
