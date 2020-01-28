### Script to find graphs with certain eucleadean destance to novkovic 
import networkx as nx 
import matplotlib.pyplot as plt 
import numpy as np
import glob
import pickle
import re

def find_sigmas(dist):
    # Find graphs with sigma vallues as novokovich and save in list of ['graph_type','parameter_set',index]
    NOV_sigma = 6.7
    files = glob.glob('../data/exp1/sigma*.pkl')
    points_list = []
    sigmas = []
    for f in files:
        data = pickle.load(open(f,'rb'))
        # go over param_keys : 
        for key,value in data.items():
            for i in range(len(value)):
                if abs(NOV_sigma - value[i]) < dist:
                    name = re.split(("sigma"),f)[-1][:-4]
                    points_list.append((name,key,i))
                    sigmas.append(value[i])
    return sigmas,points_list


def find_omegas(dist):
    # same function as above but for omegas
    NOV_omega = -.27
    files = glob.glob('../data/exp1/omega*.pkl')
    points_list = []
    omegas = []
    for f in files:
        data = pickle.load(open(f,'rb'))
        # go over param_keys : 
        for key,value in data.items():
            for i in range(len(value)):
                if abs(NOV_omega - value[i]) < dist:
                    name = re.split(("omega"),f)[-1][:-4]
                    points_list.append((name,key,i))
                    omegas.append(value[i])
    return omegas,points_list



def similar_graphs(dist_sig,dist_om):
    # overlapp between similar sigmas and omegas : 
    o,omegas = find_omegas(dist_om)
    s,sigmas = find_sigmas(dist_sig)
    print(len(omegas))
    print(len(sigmas))
    intersection = list(set(sigmas) & set(omegas))
    print(len(intersection))

    #filtered_graphs = find_sigmas(dist)
    graph_data = []
    # look voor original data : 
    for hit in intersection:
        file_name = '../data/exp1/exp1_' + hit[0] + '.pkl'
        key = hit[1]
        index = hit[2]
        sigma,omega = s[sigmas.index(hit)], o[omegas.index(hit)]

        # open file : 
        data = pickle.load(open(file_name, 'rb'))
        graph_data.append((hit[0] + key + str(index), data[key][index],(sigma,omega)))
        print(hit)

    return graph_data

def plot_filtered_graphs(fname):
    # set fontsize for plots : 
    plt.rcParams.update({'font.size': 22})

    graphs = pickle.load(open(fname,'rb'))

    for graph in graphs:
        name = graph[0]
        stats = graph[1]

        # use saved adj matrix to regenerate graph : 
        matrix = stats[-1]
        g = nx.from_numpy_matrix(matrix)
        # plot : 
        pos = nx.kamada_kawai_layout(g)
        # check if sigma and omega measures are correct : 
        plt.figure(figsize=(15,15))
        plt.title(str(graph[-1]))# + str(corr_stats)
        print(str(graph[-1]))#,str(corr_stats)
        nx.draw(g,pos)
        name = '/home/lau/GIT/FRC_Thesis/results/NOV_similar_graphs/' + name.replace("/","") + '.png'
        print(name)
        plt.tight_layout()
        plt.savefig(name)


if __name__ == "__main__":
    graphs = similar_graphs(1,.05)
    # save graphs 
    fname = '../data/exp1/filtered_graphs/graphs_data.pkl'
    fname = '/home/lau/GIT/FRC_Thesis/data/exp1/filtered_graphs/graphs_data.pkl'
    pickle.dump(graphs,open(fname,'wb'))
    plot_filtered_graphs(fname)











# def find_omegas(dist):
#     # same function as above but for omegas
    
    
#     points_list = []
#     for f in files:
#         data = pickle.load(open(f,'rb'))
#         # go over param_keys : 
#         for key,value in data.items():
#             for i in range(len(value)):
#                 if abs(NOV_omega - value[i]) < dist:
#                     name = re.split(("omega"),f)[-1][:-4]
#                     points_list.append((name,key,i))
#     return points_list


# def find_sigmas(dist):
#     # find points around novkovic return key of parameter for model, and index in format ('graph_type','param_key',i)
#     # NOV_sigma = 6.7
#     # NOV_omega = -.27
#     NOV = np.array((6.7,-.27))
#     sig_files = glob.glob('../data/exp1/sigma*.pkl')
#     om_files = glob.glob('../data/exp1/omega*.pkl')
#     sig_files.sort()
#     om_files.sort()
#     print(sig_files)
#     print(om_files)
#     points_list = []
#     for i in range(len(sig_files)):
#         sig_data = pickle.load(open(sig_files[i],'rb'))
#         om_data = pickle.load(open(om_files[i],'rb'))
#         # go over param_keys : 
#         for key,value in sig_data.items():
#             for j in range(len(value)):
#                 point = np.array((value[j],om_data[key][j]))
#                 distance = np.linalg.norm(point - NOV)
#                 if distance < dist:
#                     name = re.split(("sigma"),sig_files[i])[-1][:-4]
#                     points_list.append((name,key,i,point))
#     print(len(points_list))
#     return points_list