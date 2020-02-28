# Plot other statistics of saved networks : 

import networkx as nx 
import numpy as np 
import pickle 
import matplotlib.pyplot as plt 
from Compute_OW import extract_averages
import glob
import seaborn as sns
import pandas as pd
import re
sns.set(style="darkgrid")


def edges_to_df():
    # plot of clustering coeffient and path length
    # 2 figures 
    #   5 subplots with points that are paremeter sets
    #   size of points would be variance
    # per type :
    files = glob.glob('../data/exp1/exp1_*.pkl')
    # lists as blue print for dataframe : 
    n_edges = []
    std_edges = []
    param_set = []
    graph_type = []
    for f in files:
        data = extract_averages(f)
        for key,value in data.items():
            end = key.find('V/E')
            short_key = key[:end + 3]
            graph_name = re.split("_",f)[-1][:-4]
            param_set.append(short_key)
            graph_type.append(graph_name)

            n_edges.append(np.average(value))
            std_edges.append(np.std(value))

    
    df = pd.DataFrame.from_records([n_edges,std_edges,param_set,graph_type])
    df = df.T
    df.columns = ['n_edges','std_edges','params','type']

    df.to_csv('../results/edges.csv')


edges_to_df()
