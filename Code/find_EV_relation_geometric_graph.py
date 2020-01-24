### Script to empirically find relation between number of edges and radius R

import networkx as nx 
import numpy as np 
import numpy.polynomial
import matplotlib.pyplot as plt 
import pickle

def find_ER(N,r):
    # generate x random geometric graphs and return E
    E_s = []
    for i in range(5):
        g = nx.random_geometric_graph(N,r, dim=2, pos=None, p=2, seed=None)
        E_s.append(g.number_of_edges())
    return E_s

def scatter_plot(N,r_list):
    data = []
    r_data = []
    for r in r_list:
        E = find_ER(N,r)
        data += E
        r_data += len(E)*[r]
    
    # plt.figure(figsize=(10,10))
    # plt.scatter(r_data,data)
    # plt.title('Relation between radius and number of edges in /n random geometric graph')
    # plt.ylabel('Number of Edges')
    # plt.xlabel('Radius')
    # plt.show()

    # save data : 
    pickle.dump([data,r_data],open('../data/helper_data/' + str(N) + 'edge_data.pkl', 'wb'))

    return data,r_data

def polynomial(coefs,x):
    d  = len(coefs)
    return sum([x**(d - i) * coefs[i-1] for i in range(1,d)])

def return_radius(N,E):
    f = pickle.load(open('../data/helper_data/' + str(N) + 'edge_data.pkl', 'rb'))
    x = f[1]
    y = f[0]
    degr = 3
    coefs = np.polyfit(x,y,degr)
    return polynomial(coefs,E)

def fit_polynomial(N,degr):
    f = pickle.load(open('../data/helper_data/' + str(N) + 'edge_data.pkl', 'rb'))
    x = f[1]
    y = f[0]
    coefs = np.polyfit(x,y,degr)
    print(coefs)
    # plot found polynomial over data : 
    x_fit = np.linspace(x[0],x[-1],100)
    y_fit = [polynomial(coefs,i) for i in x_fit]
    
    plt.figure(figsize=(10,10))
    plt.scatter(x,y, label = 'DATA')
    plt.plot(x_fit,y_fit,c = 'r',label = 'FIT')
    plt.title('Relation between radius and number of edges in \n random geometric graph')
    plt.ylabel('Number of Edges')
    plt.xlabel('Radius')
    plt.legend()
    plt.savefig('../data/helper_data/EV_relation_geometric' + str(N) + str(x[-1]) + '.png')
    plt.show()


if __name__ == '__main__':
    degree = 4
    N_list = [100,176,500]
    for N in N_list:
        r_list = np.linspace(0,1,100)
        Edges,radius = scatter_plot(N,r_list)
        fit_polynomial(N,degree)