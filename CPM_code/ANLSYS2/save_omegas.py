import networkx as nx #import smallworld
import pickle
import numpy as np

def calc_omegas(path):
    graph = pickle.load(open(path,'rb'))
    return nx.sigma(graph,niter = 100,nrand = 50)

if __name__ == "__main__":
    gnames = ['5WS','2PW','0GM']
    lines = '' 
    for g in gnames:
        filename = '../../data/FRCs/GRAPH'+ g + '64_diam3.pkl'
        om = calc_omegas(filename)
        lines += g + '\t' + str(om) + '\n'
        print(lines)

    with open('STROMAL_sigma.txt3','w') as f:
        f.write(lines)
