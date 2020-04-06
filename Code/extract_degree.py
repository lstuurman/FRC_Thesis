
import numpy as np
import re
import glob
import pandas as pd
import pickle

def extract_degrees(path):
    files = glob.glob(path)
    degrees = []
    for f in files:
        print('loading data ',f)
        pattern = '176N_0.3125V/E '
        dict = pickle.load(open(f,'rb'))
        name = f.split('_')[-1][:-4]
        for key in dict.keys():
            if pattern in key:
                data = dict[key]
                print('extracting data')
                for i,network in enumerate(data):
                    for d in network[0]:
                        degrees.append([d,i,name])
    #np.savetxt('../data/exp1/degrees.txt',degrees)
    df = pd.DataFrame(degrees)
    df.columns = ['Degree','iter','type']
    print(df.head())
    df.to_csv('../results/degrees.csv')


extract_degrees('../data/exp1/exp1*.pkl')
