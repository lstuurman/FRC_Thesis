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
    files = glob.glob('../data/exp1/exp1_*.pkl')
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
            n_edges.append(value[-1])
    
    df = pd.DataFrame.from_records([clustering,path,std_clustering,std_path,n_edges,param_set,graph_type])
    df = df.T
    df.columns = ['clustering','path_length','std_clustering','std_path','n_edges','params','type']
    # add cumulative variance :
    norm_clust =  ((df['std_clustering'] - df['std_clustering'].min())/(df['std_clustering'].max() - df['std_clustering'].min()))
    norm_pl = ((df['std_path'] - df['std_path'].min())/(df['std_path'].max() - df['std_path'].min()))
    df['CUMUL_STD'] = norm_clust + norm_pl

    df.to_csv('../data/exp1/CL_P.csv')

def to_df1():
    # plot of clustering coeffient and path length
    # 2 figures 
    #   5 subplots with points that are paremeter sets
    #   size of points would be variance
    # per type :
    sig_files = glob.glob('../data/exp1/sigma*.pkl')
    og_files = glob.glob('../data/exp1/omega*.pkl')
    sig_files.sort()
    og_files.sort()
    # lists as blue print for dataframe : 
    sigmas = []
    omegas = []
    std_sigmas = []
    std_omegas = []
    param_set = []
    graph_type = []
    for i in range(len(sig_files)):
        #print(sig_files[i])
        sig_data = pickle.load(open(sig_files[i],'rb'))
        #print(data)
        #sig_data = extract_averages(sig_files[i])
        om_data = pickle.load(open(og_files[i],'rb'))
        for key,value in sig_data.items():
            end = key.find('V/E')
            short_key = key[:end + 3]
            graph_name = re.split("sigma",sig_files[i])[-1][:-4]
            #print(np.average[value])
            sigmas.append(np.average(value))
            omegas.append(np.average(om_data[key]))
            std_sigmas.append(np.std(value))
            std_omegas.append(np.std(om_data[key]))
            param_set.append(short_key)
            graph_type.append(graph_name)

    df = pd.DataFrame.from_records([sigmas,omegas,std_sigmas,std_omegas,param_set,graph_type])
    df = df.T
    print(df.head())
    df.columns = ['sigma','omega','std_sigma','std_omega','params','type']
    # add cumulative variance :
    norm_clust =  ((df['std_sigma'] - df['std_sigma'].min())/(df['std_sigma'].max() - df['std_sigma'].min()))
    norm_pl = ((df['std_omega'] - df['std_omega'].min())/(df['std_omega'].max() - df['std_omega'].min()))
    df['CUMUL_STD'] = norm_clust + norm_pl

    df.to_csv('../data/exp1/omega_sigma.csv')

def scatter():
    data = pd.read_csv('../results/omega_sigma.csv')
    g = sns.FacetGrid(data, col="type",col_wrap=3, hue="type",sharey = False, sharex = False)#,size = "CUMUL_STD"
    g.map(sns.regplot,"sigma", "omega",fit_reg = False, scatter_kws={'s': data["CUMUL_STD"] * 100,'alpha':.5})#,size = "CUMUL_STD"
    g.add_legend()
    # g = sns.relplot(x="total_bill", y="tip",hue="day", col="time", data=tips)
    #g = sns.relplot(x = "sigma", y = "omega",size = "CUMUL_STD", data = data,col = "type",sizes = (15,200),hue = "type",alpha=.5,facet_kws={'sharey': False, 'sharex': False})
    #leg = g._legend
    # truncate legend texts:
    # for t in leg.texts:
    #     if len(t.get_text()) > 9:
    #         t.set_text(t.get_text()[:4])
    # leg.set_bbox_to_anchor([1., 0.7])  # coordinates of lower left of bounding box
    plt.show()

def histogram():
    data = pd.read_csv('../data/exp1/CL_P.csv')
    g = sns.FacetGrid(data, col="type", hue="type",sharey = False, sharex = False)
    g.map(plt.hist,"clustering",alpha = .5,bins = 20)
    g.add_legend()
    plt.show()
    plt.close()
    g2 = sns.FacetGrid(data, col="type", hue="type",sharey = False, sharex = False)
    g2.map(plt.hist,"path_length",alpha = .5)
    g2.add_legend()
    plt.show()

if __name__ == "__main__":
    #to_df()
    #to_df1()
    # #to_df()
    scatter()
    # histogram()
