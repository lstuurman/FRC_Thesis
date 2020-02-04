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
        g = nx.random_geometric_graph(N,r, dim=3, pos=None, p=2, seed=None)
        E_s.append(g.number_of_edges())
    return E_s

def scatter_plot(N,r_list):
    data = []
    r_data = []
    for r in r_list:
        E = find_ER(N,r)
        data += E
        r_data += len(E)*[r]
    print(r_data)
    
    # plt.figure(figsize=(10,10))
    # plt.scatter(r_data,data)
    # plt.title('Relation between radius and number of edges in /n random geometric graph')
    # plt.ylabel('Number of Edges')
    # plt.xlabel('Radius')
    # plt.show()

    # save data : 
    pickle.dump([data,r_data],open('../data/helper_data/' + str(N) + 'edge_data3d.pkl', 'wb'))

    return data,r_data

# def polynomial(coefs,x):
#     d  = len(coefs)
#     p = np.poly1d(coefs)
#     return p(x)
    #return sum([x**(d - i) * coefs[i-1] for i in range(1,d)])
    #return sum([x**(i) * coefs[i] for i in range(len(coefs))])
    #return sum(x**(i) * coefs[i] for i in range(len(coefs)))

def return_radius3d(N,n_edges):
    f = pickle.load(open('../data/helper_data/' + str(N) + 'edge_data3d.pkl', 'rb'))
    r = f[1]
    E = f[0]
    degr = 3
    coefs = np.polyfit(E,r,degr)
    p = np.poly1d(coefs)
    return p(n_edges)

def fit_polynomial(N,degr):
    f = pickle.load(open('../data/helper_data/' + str(N) + 'edge_data3d.pkl', 'rb'))
    r = f[1]
    E = f[0]
    coefs = np.polyfit(E,r,degr)
    print(coefs)
    # plot found polynomial over data : 
    x_fit = np.linspace(min(E),max(E),100)
    print(E[0],E[-1])
    print(r[0],r[-1])
    polynomial = np.poly1d(coefs)
    y_fit = [polynomial(i) for i in x_fit]
    print(y_fit[0],y_fit[-1])
    plt.figure(figsize=(10,10))
    plt.scatter(E,r, label = 'DATA')
    plt.plot(x_fit,y_fit,c = 'r',label = 'FIT')
    plt.title('Relation between radius and number of edges in \n random geometric graph')
    plt.xlabel('Number of Edges')
    plt.ylabel('Radius')
    plt.legend()
    plt.savefig('../data/helper_data/EV_relation_geometric3d' + str(N) + str(r[-1]) + '.png')
    plt.show()


if __name__ == '__main__':
    degree = 3
    N_list = [100,176,500,1000]
    for N in N_list:
        r_list = np.linspace(0.05,.4,100)
        Edges,radius = scatter_plot(N,r_list)
        fit_polynomial(N,degree)