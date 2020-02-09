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


def to_df():
    # plot of clustering coeffient and path length
    # 2 figures 
    #   5 subplots with points that are paremeter sets
    #   size of points would be variance
    # per type :
    files = glob.glob('/home/lau/GIT/FRC_Thesis/data/exp1/exp1_*.pkl')
    # lists as blue print for dataframe : 
    clustering = []
    path = []
    std_clustering = []
    std_path = []
    n_edges = []
    param_set = []
    graph_type = []
    for f in files:
        data = extract_averages(f)
        for key,value in data.items():
            end = key.find('V/E')
            short_key = key[:end + 3]
            graph_name = re.split("_",f)[-1][:-4]
            clustering.append(np.average(value[0]))
            path.append(np.average(value[1]))
            std_clustering.append(np.std(value[0]))
            std_path.append(np.std(value[1]))
            param_set.append(short_key)
            graph_type.append(graph_name)
    
    df = pd.DataFrame.from_records([clustering,path,std_clustering,std_path,param_set,graph_type])
    df = df.T
    df.columns = ['clustering','path_length','std_clustering','std_path','params','type']
    # add cumulative variance :
    norm_clust =  ((df['std_clustering'] - df['std_clustering'].min())/(df['std_clustering'].max() - df['std_clustering'].min()))
    norm_pl = ((df['std_path'] - df['std_path'].min())/(df['std_path'].max() - df['std_path'].min()))
    df['CUMUL_STD'] = norm_clust + norm_pl

    df.to_csv('/home/lau/GIT/FRC_Thesis/data/exp1/CL_P.csv')


def scatter():
    data = pd.read_csv('/home/lau/GIT/FRC_Thesis/data/exp1/CL_P.csv')
    #g = sns.FacetGrid(data, col="type", hue="type",sharey = False, sharex = False)
    #g.map(sns.relplot,x = "path_length", y = "clustering",size = "std_clustring", data = data,alpha=.5)
    #g.add_legend()
    # g = sns.relplot(x="total_bill", y="tip",hue="day", col="time", data=tips)
    g = sns.relplot(x = "path_length", y = "clustering",size = "CUMUL_STD", data = data,col = "type",sizes = (15,200),hue = "type",alpha=.5,facet_kws={'sharey': False, 'sharex': False})
    leg = g._legend
    # truncate legend texts:
    for t in leg.texts:
        if len(t.get_text()) > 9:
            t.set_text(t.get_text()[:4])
    leg.set_bbox_to_anchor([1., 0.7])  # coordinates of lower left of bounding box
    plt.show()

def histogram():
    data = pd.read_csv('/home/lau/GIT/FRC_Thesis/data/exp1/CL_P.csv')
    g = sns.FacetGrid(data, col="type", hue="type",sharey = False, sharex = False)
    g.map(plt.hist,"clustering",alpha = .5,bins = 20)
    g.add_legend()
    plt.show()
    plt.close()
    g2 = sns.FacetGrid(data, col="type", hue="type",sharey = False, sharex = False)
    g2.map(plt.hist,"path_length",alpha = .5)
    g2.add_legend()
    plt.show()

#to_df()
scatter()
histogram()