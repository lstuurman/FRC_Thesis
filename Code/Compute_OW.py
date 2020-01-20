# script to compute omega and sigma stats of generated graphs
import pickle
import numpy as np
import re

def extract_averages(path):
    data = pickle.load(open(path,'rb'))
    av_dict = {}
    for key,value in data.items():
        av_dict[key] = [[np.average(x[1]) for x in value],[np.average(x[2]) for x in value]]
    return av_dict

def comparison_averages(path):
    data = pickle.load(open(path,'rb'))
    av_dict = {}
    for key,value in data.items():
        av_dict[key] = [[np.average(x[0]) for x in value],[np.average(x[1]) for x in value]]
    return av_dict

def compute_sigma(data,random,lattice):
    # data is dictionary with average cluster coeffs and path lengths for generated graphs
    # random is a dict containing equivalant random graphs clustering and path lengths
    sigma = {}
    for key,value in data.items():
        # make shure we have a key that exists in random data (no extra params)
        end = key.find('V/E')
        short_key = key[:end + 3]
        sigma[key] = []     #initialize dict field
        # search for equivalent stats: 
        CL_rand = random[short_key][0]
        Pl_rand = random[short_key][1]
        CLs_data = value[0]
        Pls_data = value[1]

        # iterate over all generated graphs with specific parameters : N/
        for i in range(len(value)):
            #for i in range(len(Cl)):
            CL = CLs_data[i]
            Pl = Pls_data[i]
            sigma_estimates = []
            # compute sigma for every equivalent graph : 
            for j in range(len(CL_rand)):
                sigma_estimates.append((CL/CL_rand[j])/(Pl/Pl_rand[j]))
        
            sigma[key].append(np.average(sigma_estimates))


if __name__ == "__main__":
    data = extract_averages('data/exp1_ER.pkl')
    random = comparison_averages('comparison_graphs/equivalent_random_stats.pkl')
    #data = comp_averages('data/exp1_ER.pkl')
    sigmas = compute_sigma(data,random)
    for key,val in sigmas.items():
        print(key)
        print(val)
            