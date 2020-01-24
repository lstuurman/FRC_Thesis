# script to check results of my own sigma and omega calculations : 
import pickle
import networkx as nx

def check():
    fname = '/home/lau/GIT/FRC_Thesis/data/exp1/filtered_graphs/graphs_data.pkl'
    data = pickle.load(open(fname,'rb'))

    for graph in data:
        name = graph[0]
        stats = graph[1]
        matrix = stats[-1]
        g = nx.from_numpy_matrix(matrix)
        print(name)
        sigma,omega = nx.sigma(g,niter=10) , nx.omega(g,niter=10)
        print(sigma,omega)

check()