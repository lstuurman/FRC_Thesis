# script to compute omega and sigma small world measurs of graphs 
from Generate_graphs import ER_p_value
from Generate_graphs import WS_K_value
from multiprocessing import Process,Manager,Pool
import pickle
import networkx as nx
import numpy as np
import time
from itertools import product

def ER_p_value2(N,E):
    # helper function to find p parameter resulting in certian V/E ratio
    n_edges = N * (N/E)**-1
    return 2*n_edges / (N * (N-1))

def gen_equiv_random(tup):
    E,N = tup
    # generate 1000 random graphs, equivalent to the ones used in the experiments
    p = ER_p_value2(N,E)
    p2 = ER_p_value(N,N/E)
    stats  = [] # list of tuples with (average_clustering,average_path) for every generated network. 
    count  = 0
    edges = []
    print(tup)
    print(p2,p)
    while count < 1:
        #print(tup)
        #print(p)
        #print(p2)
        g = nx.fast_gnp_random_graph(N,p)
        #print(nx.is_connected(g))
        #print(g.number_of_edges())
        #exit()
        if nx.is_connected(g):
            # average clustering and path length : 
            Cl = nx.average_clustering(g)
            L = nx.average_shortest_path_length(g)
            stats.append((Cl,L))
            count += 1
            #print(count)
            edges.append(g.number_of_edges())
            print('HIT')
    #print(N,g.number_of_edges())
    f = '../data/helper_data/'+str(N)+'N_'+str(E)+'E.txt'
    np.savetxt(f,edges)
    #print(tup)
    #print(N,g.number_of_edges())
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
    #VE_ratios = np.linspace(.125,.375,5)
    #N_nodes = [100,176,500,1000]
    #product = list(product(N_nodes,VE_ratios))
    product = np.loadtxt('../results/Graphs_stuff/unique_NEtupes')
    keys = [str(int(x[0])) + 'N_' + str(int(x[1])) + 'E' for x in product]
    inputs = [(int(x[0]),int(x[1])) for x in product]
    t0 = time.time()
    p = Pool(10)
    outputs  = p.map(gen_equiv_random,inputs)
    p.close()
    p.join()
    data_dict = {}
    print('Generated random graphs in : ', time.time() - t0)
    for i in range(len(keys)):
        data_dict[keys[i]] = outputs[i]

    pickle.dump(data_dict,open('../data/helper_data/equivalent_random_stats2.pkl','wb'))

    # for key in keys:
    #     print(key)
    #     print('\n')
    # print(len(keys))
    # save_random_graphs(N_nodes,VE_ratios)
    generate_lattice(N_nodes,VE_ratios)
