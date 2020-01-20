# script to compute omega and sigma small world measurs of graphs 
from Generate_graphs import ER_p_value
from Generate_graphs import WS_K_value
import pickle
import networkx as nx
import numpy as np

def gen_equiv_random(N,ve_ratio):
    # generate 1000 random graphs, equivalent to the ones used in the experiments
    p = ER_p_value(N,ve_ratio)
    stats  = [] # list of tuples with (average_clustering,average_path) for every generated network. 
    count  = 0
    while count < 10:
        g = nx.erdos_renyi_graph(N,p)
        if nx.is_connected(g):
            # average clustering and path length : 
            Cl = nx.average_clustering(g)
            L = nx.average_shortest_path_length(g)
            stats.append((Cl,L))
            count += 1
            print(count)
    print(N,g.number_of_edges())
    
    return stats

def save_random_graphs(N_list, VE_list):
    data_dict = {}
    for N in N_list:
        for ve in VE_ratios:
            stats = gen_equiv_random(N,ve)
            key = str(N) + 'N_' + str(ve) + 'V/E'
            data_dict[key] = stats
    
    pickle.dump(data_dict,open('comparison_graphs/equivalent_random_stats.pkl','wb'))

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
    
    pickle.dump(data_dict,open('comparison_graphs/equivalent_lattice_stats.pkl','wb'))



if __name__ == "__main__":
    VE_ratios = np.linspace(.125,.375,5)
    N_nodes = [100,176,500]
    save_random_graphs(N_nodes,VE_ratios)
    generate_lattice(N_nodes,VE_ratios)
