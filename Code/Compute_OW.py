# script to compute omega and sigma stats of generated graphs
import pickle
import numpy as np
import re
import glob
import pandas as pd

def extract_averages(path):
    data = pickle.load(open(path,'rb'))
    av_dict = {}
    for key,value in data.items():
        av_dict[key] = [[x[3] for x in value],[x[4] for x in value]]
        #av_dict[key] = [[np.average(x[1]) for x in value],[np.average(x[2]) for x in value]]
    return av_dict

def extract_specific_stats(path,stat_index1,stat_index2):
    data = pickle.load(open(path,'rb'))
    av_dict = {}
    for key,value in data.items():
        av_dict[key] = [[np.average(x[stat_index1]) for x in value],[np.average(x[stat_index2]) for x in value]]
    return av_dict

def comparison_averages(path):
    data = pickle.load(open(path,'rb'))
    av_dict = {}
    for key,value in data.items():
        av_dict[key] = [[np.average(x[0]) for x in value],[np.average(x[1]) for x in value]]
    return av_dict

def compute_sigma(data,random,lattice):#
    # data is dictionary with average cluster coeffs and path lengths for generated graphs
    # random is a dict containing equivalant random graphs clustering and path lengths
    sigma = {}
    omega = {}
    for key,value in data.items():
        # make shure we have a key that exists in random data (no extra params)
        end = key.find('V/E')
        short_key = key[:end + 3]
        sigma[key] = []     #initialize dict field
        omega[key] = []
        # search for equivalent stats: 
        CL_rand = random[short_key][0]
        Pl_rand = random[short_key][1]
        CL_latt = lattice[short_key]
        CLs_data = value[0]
        Pls_data = value[1]

        # iterate over all generated graphs with specific parameters : N/E
        for i in range(len(CLs_data)):
            #for i in range(len(Cl)):
            CL = CLs_data[i]
            Pl = Pls_data[i]
            sigma_estimates = []
            omega_estimates = []
            # compute sigma for every equivalent graph : 
            for j in range(len(CL_rand)):
                sigma_estimates.append((CL/CL_rand[j])/(Pl/Pl_rand[j]))
                omega_estimates.append(Pl_rand[j]/Pl - CL/CL_latt)

            # # take average over all networks of this parameter : 
            sigma[key].append(np.average(sigma_estimates))
            omega[key].append(np.average(omega_estimates))
            # or don't take average : 
            # sigma[key].append(sigma_estimates)
            # omega[key].append(omega_estimates)
        
    return sigma,omega

def OW_allgraphs():
    # Comparison graphs
    random = comparison_averages('../data/helper_data/equivalent_random_stats.pkl')
    lattice = pickle.load(open('../data/helper_data/equivalent_lattice_stats.pkl','rb'))
    graphs = glob.glob('../data/exp1/exp1_*.pkl')
    # to dicts because some models have more parameters and thus more values
    averages1 = {}
    averages2 = {}
    names = []
    for graph_file in graphs:
        print(graph_file)
        # don't plot averages. but raw ow_stats of graphs : 
        graph = extract_averages(graph_file)
        #print(pickle.load(open(graph_file,'rb')))
        #print(graph)
        #pass
        sigmas,omegas = compute_sigma(graph,random,lattice)

        test_sigma = np.array(list(sigmas.values())).flatten()
        test_omega = np.array(list(omegas.values())).flatten()
        #print(len(list(sigmas.values())), len(list(omegas.values())))
        print(len(test_omega))
        print(len(test_sigma))
        # save raw sigmas and omegas : 
        graph_name = re.split("_",graph_file)[-1]
        names.append(graph_name[:-4])
        pickle.dump(sigmas,open('../data/exp1/sigma' + graph_name,'wb'))
        pickle.dump(omegas,open('../data/exp1/omega' + graph_name,'wb'))
        # save averages to write as csv
        sig_values = [np.average(x) for x in sigmas.values()]
        om_values = [np.average(x) for x in omegas.values()]
        data = list(zip(sig_values,om_values))
        if len(data) > 15:
            averages1[graph_name] = data
        else:
            averages2[graph_name] = data
        #averages[graph_name[:-4]] = [list(zip(sig_values,om_values))]
    
    # convert to pandas and save as csv
    df1 = pd.DataFrame(averages1)#, columns = names
    #df1 = df1.transpose()
    df1.to_csv('../data/exp1/Omega_Sigma1.csv')
    # dict 2: 
    df2 = pd.DataFrame(averages2)#, columns = names
    #df2 = df2.transpose()
    df2.to_csv('../data/exp1/Omega_Sigma2.csv')

if __name__ == "__main__":
    OW_allgraphs()

    data = extract_averages('../data/exp1/exp1_ER.pkl')
    random = comparison_averages('../data/helper_data/equivalent_random_stats.pkl')
    lattice = pickle.load(open('../data/helper_data/equivalent_lattice_stats.pkl','rb'))

    #data = comp_averages('data/exp1_ER.pkl')
    sigma,omega = compute_sigma(data,random,lattice)
    # for key,val in sigma.items():
    #     print(key)
    #     print('sigma : ',val,' | omega : ',omega[key])
            
