### AIM : Generate different graphs
###     : calc basic stats : degrees, path_lengths, clustering
###     : to safe   -> Number of edges
###                    Number of nodes  
###                    list of degrees
###                    list of path_lengths
###                    list of clustering coefficients of all graphs

import numpy as np
import networkx as nx
import pickle as pkl
import scipy.optimize
import time
from find_EV_relation_geometric_graph import return_radius
#from find_EV_relation_geometric_graph import polynomial
from numpy import linalg


def get_metrics(g):
    # start with some quick and metrics
    # degree per node : 
    node_degree_pairs = list(g.degree())
    degrees = list(list(zip(*node_degree_pairs))[1]) # take only a list of degrees
    clustering_coefs  = list(nx.clustering(g).values())
    av_clustering = nx.average_clustering(g)
    paths = nx.shortest_path(g)
    #paths_list = [list(i.values()) for i in paths.values()] # convert weird dict to list countaining paths
    path_lengths = []
    for key1,val1 in paths.items():
        for key2,val2 in val1.items():
            if key1 != key2:
                path_lengths.append(len(val2))
    # for i in range(len(paths_list)):
    #     for j in range(len(paths_list[0])):
    #         if i != j:
    #             path_lengths.append(len(paths_list[i][j]))

    # check if we have correct number of paths ;
    n = g.number_of_nodes()
    expected = .5*(n*(n-1))
    av_pl = nx.average_shortest_path_length(g)
    #print(len(paths) == expected)
    #print(len(path_lengths),expected)
    # betweennes centrality : 
    centrality = list(nx.betweenness_centrality(g).values())
    # eigenvalues : 
    vals,vecs = linalg.eigh(nx.to_numpy_matrix(g))
    max_eig = max(vals)
    print(g.number_of_edges())
    # also return adjacency matrix to regenerate graph
    return [degrees,clustering_coefs,path_lengths,av_clustering,av_pl,centrality,max_eig,nx.to_numpy_matrix(g),g.number_of_edges()]

def ER_p_value(N,ve_ratio):
    # helper function to find p parameter resulting in certian V/E ratio
    n_edges = N * ve_ratio**-1
    return 2*n_edges / (N * (N-1))


def gen_ER_graphs(ratios,N_nodes,repetitions):
    """
    Function to generate wide range of ER networks and gather statistics
    data : dictionary where keys are 'xN_yV/E' with x en y specific N and V/E ratio values
           every entry in dict is matrix of [[degrees,clustering_coefs,path_lengths],[degrees,clustering_coefs,path_lengths],... etc}
           So 3 dimensional matrix 
    """
    data = {}
    for N in N_nodes:
        for r in ratios:
            # initialize field in dict 
            data_key = str(N) + 'N_' + str(r) + 'V/E'
            data[data_key] = []

            count = 0
            while count < repetitions * 10:
                p = ER_p_value(N,r)
                ER_g = nx.erdos_renyi_graph(N,p)
                if nx.is_connected(ER_g):
                    stats = get_metrics(ER_g)
                    data[data_key].append(stats)
                    count += 1
    return data

# nextttt watts-Str0gatzz
#g = nx.watts_strogatz_graph(176, 8, 0)
# params K and P
def WS_K_value(N,ve_ratio):
    # Helper function to find proper K value for N-V/E pair
    n_edges = N * ve_ratio**-1
    return int(2*n_edges/N)

def gen_WS_graphs(ratios,N_nodes,repetitions):
    # Take p_values logaritmically spaced:
    p_vals = np.logspace(-2,0,10)

    data = {}
    for N in N_nodes:
        for r in ratios:
            for p in p_vals:
                # initialize field in dict 
                data_key = str(N) + 'N_' + str(r) + 'V/E ' + str(p) + 'P'
                data[data_key] = []

                for iter in range(repetitions):
                    K = WS_K_value(N,r)
                    WS_g = nx.watts_strogatz_graph(N,K,p)
                    stats = get_metrics(WS_g)
                    data[data_key].append(stats)
    return data

def radius_geoGraph(N,ve_ratio):
    # https://mathoverflow.net/questions/124579/mean-minimum-distance-for-n-random-points-on-a-unit-square-plane
    # https://www.youtube.com/watch?v=i4VqXRRXi68

    n_edges = N * ve_ratio**-1
    return return_radius(N,n_edges)
    
def generate_geometric_graphs(ratios,N_nodes,repetitions):
    data = {}
    for N in N_nodes:
        for r in ratios:
            data_key = str(N) + 'N_' + str(r) + 'V/E '
            data[data_key] = []
            count = 0 
            #for iter in range(repetitions*10):
            while count < 1:#repetitions*10
                radius = radius_geoGraph(N,r)
                geo_g = nx.random_geometric_graph(N,radius, dim=2, p=2)
                if nx.is_connected(geo_g):
                    print(radius)
                    print(N)
                    stats = get_metrics(geo_g)
                    data[data_key].append(stats)
                    count += 1
    return data


def E_M_relation(x,N,E):
    # K² - KN + E = 0
    return x**2 - N*x + E


def M_power_cluster(N,ve_ratio):
    # to find correct value for K  resulting in approximatly correct density
    # we need to find the roots of K² - KN + E = 0
    n_edges = N * ve_ratio**-1
    sol = scipy.optimize.root(E_M_relation,10,args = (N,n_edges))
    return int(sol.x)

def gen_power_cluster_graphs(ratios,N_nodes,repetitions):
    # Take p_values between 0 and 1:
    p_vals = np.linspace(0,1,10)

    data = {}
    for N in N_nodes:
        for r in ratios:
            for p in p_vals:
                # initialize field in dict 
                data_key = str(N) + 'N_' + str(r) + 'V/E ' + str(p) + 'P'
                data[data_key] = []

                for iter in range(repetitions):
                    M = M_power_cluster(N,r)
                    power_g = nx.powerlaw_cluster_graph(N,M,p)
                    stats = get_metrics(power_g)
                    data[data_key].append(stats)
    return data

def gen_BA_graphs(ratios,N_nodes,repetitions):
    data = {}
    for N in N_nodes:
        for r in ratios:
            # initialize field in dict 
            data_key = str(N) + 'N_' + str(r) + 'V/E '
            data[data_key] = []

            for iter in range(repetitions * 10):
                #M = M_power_cluster(N,r) 
                # r**-1 gives you the approximate value for M number of edges to add for every node
                BA_g = nx.barabasi_albert_graph(N,int(r**-1))
                #power_g = nx.powerlaw_cluster_graph(N,M,p)
                stats = get_metrics(BA_g)
                data[data_key].append(stats)
    return data


if __name__ == "__main__":
    global repetitions
    repetitions = 10
    VE_ratios = np.linspace(.125,.375,5)
    N_nodes = [100,176,500,1000]

    #t = time.time()
    #ER_data = gen_ER_graphs(VE_ratios,N_nodes,repetitions)
    #pkl.dump(ER_data,open('../data/exp1/exp1_ER.pkl','wb'))
    #print('ER graphs generated in : ',time.time() - t)
    #WS_data = gen_WS_graphs(VE_ratios,N_nodes,repetitions)
    #pkl.dump(WS_data,open('../data/exp1/exp1_WS.pkl','wb'))
    #print('WS graphs generated after : ',time.time() - t)
    #power_data = gen_power_cluster_graphs(VE_ratios,N_nodes,repetitions)
    #pkl.dump(power_data,open('../data/exp1/exp1_power.pkl','wb'))
    #print('Clustered powerlaw graphs generated after : ',time.time() - t)
    geom_data = generate_geometric_graphs(VE_ratios,N_nodes,repetitions)
    pkl.dump(geom_data,open('../data/exp1/exp1_geom.pkl','wb'))
    #print('random geometric graphs generated after : ',time.time() - t)
    #BA_data = gen_BA_graphs(VE_ratios,N_nodes,repetitions)
    #pkl.dump(BA_data,open('../data/exp1/exp1_BA.pkl','wb'))
    #print('BA graphs generated after : ',time.time() - t)


    
    
    
    



#date = re.search("([0-9]{4}\-[0-9]{2}\-[0-9]{2})"



    # # test helper functions : 
    # print('expected P_val erdos_renyii : ' , 0.04448051948)
    # print('helper result : ', ER_p_value(176,.25))
    # print('expected K : ',8)
    # print('helper result : ',WS_K_value(176,.25))
    # print('expected k value' , 4)
    # print('helper_result : ', M_power_cluster(400,.25))













# class Graph_Generator:
#     def __init__(self,g_type,n):
#         self.gtype = g_type
#         self.graph_types = [nx.erdos_renyi_graph,nx.watts_strogatz_graph,]
#         self.n_nodes = n_nodes

#     def generate(self):
#         # iterate over number of parameters : 
#         for p 
#         # determine graph type : 

# clustering = nx.average_clustering(g)
# shortest_path = nx.average_shortest_path_length(g)
# print('clustering : ',clustering)
# print('shortest path : ',shortest_path)
# # slow metrics for small world ness : 
# sigma = nx.sigma(g,niter = 5)
# omega = nx.omega(g,niter = 5)
# print('small-world metrics: s = ',sigma,'w',omega)
# return [clustering,shortest_path,sigma,omega]
