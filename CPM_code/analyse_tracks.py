import numpy as np 
import glob
#from CPM_helpers1 import *
from MSD1cell import *
import pandas as pd
import re 
#from mpl_toolkits import mplot3d
#import matplotlib.pyplot as plt
#import seaborn as sns
import pickle
from numpy.linalg import norm
#from vis_frc import visualize_frc
#sns.set_style('darkgrid')


def get_motility(track,plot = False):
    # cell track is list of 3d coordinates.
    if plot:
        plot_celltrack(track)
    #vol = scanned_volume(volume_track)
    displ, angles = analyse_track(track)
    # MSD
    delta_t,MSD = compute_MSD(displ,plot = plot)
    # autocorr : 
    _,_,_,slope,p = compute_AutoCorr_WRONG(angles,plot = plot)
    #t,new_angles = compute_autocorrelaton(displ)
    
    popt = fit_Motilty(delta_t,MSD)

    return popt[0],popt[1],slope,p

def roseplot(track):
    # get list of angles : 
    _,angles = analyse_track(track)

    n_numbers = 100
    bins_number = 16  # the [0, 360) interval will be subdivided into this
    # number of equal bins
    bins = np.linspace(0.0, 2 * np.pi, bins_number + 1)
    #angles = 2 * np.pi * np.random.rand(n_numbers)
    n, _, _ = plt.hist(angles, bins)

    #plt.clf()
    width = 2 * np.pi / bins_number
    ax = plt.subplot(1, 1, 1, projection='polar')
    bars = ax.bar(bins[:bins_number], n, width=width, bottom=0.0)
    for bar in bars:
        bar.set_alpha(0.5)
    plt.show()

def new_auto(cell_track): 
    averages = []
    for dt in range(0,300):
        # angles per dt: 
        cosines = []
        for i in range(len(cell_track) - 1 - dt):
            point1 = cell_track[i]
            point2 = cell_track[i + 1]
            point3 = cell_track[i + dt]
            point4 = cell_track[i + dt + 1]

            v1 = point2 - point1
            v2 = point4 - point3
            # check for zero fectors to prevent division by zero: 
            if (list(point2) == list(point4) or norm(v1) == 0) or norm(v2) == 0:
                #print('non moving cell')
                continue
            #cos = np.clip(np.dot(v1,v2)/(norm(v1) * norm(v2)),-1,1)
            cos = np.dot(v1,v2)/(norm(v1) * norm(v2))
            cosines.append(cos)
        if len(cosines) > 100: #minimum 100 mcs with actual displacement
            averages.append(np.average(cosines))
    return averages

def all_autos(path):
    files = glob.glob(path)
    print(len(files))
    files.sort()
    num_pattern = '[-+]?\d*\.\d+|\d+'
    #uniq_names = set([name[:-5] for name in files])
    rows = []
    for id_,name in enumerate(files):
        #f_names = glob.glob(name + '*.txt')
        # params : 
        l,Max= name.split('MAX')
        print(l,Max)
        _lambda = int(re.findall(num_pattern,l)[1])
        _max = int(re.findall(num_pattern,Max)[0])
        cell_id = int(re.findall(num_pattern,Max)[1])
        track = np.loadtxt(name)
        correlations = new_auto(track)
        #tracks = [np.loadtxt(f) for f in f_names]
        #correlations = [new_auto(t) for t in tracks]
        #min_l = min([len(x) for x in correlations])
        # truncate on smallest track :
        #correlations = [x[:min_l] for x in correlations]
        # average all 3 correlations: 
        #av_cors = np.average(correlations,axis = 0)

        # append rows to create nice dataframe
        for i,x in enumerate(correlations):
            rows.append([id_,cell_id,i,x,_lambda,_max])
        print(name)
        print(_lambda, ':',_max,':',cell_id)
        #break

    df = pd.DataFrame(data = rows,columns = ['id','cell_id','dt','auto correlation','lambda','max act'])
    df.to_csv('../results/autocorr_fullLN_frc2.csv')
    # # plot :::
    # sns.lineplot(x = 'dt', y = 'auto correlation', hue = 'lambda', data = df)#, style = 'max act' ,
    # plt.show()

def all_autos_average(path):
    files = glob.glob(path+'*.txt')
    print(len(files))
    files.sort()
    num_pattern = '[-+]?\d*\.\d+|\d+'
    uniq_names = set([name[:-5] for name in files])
    rows = []
    for id,name in enumerate(uniq_names):
        f_names = glob.glob(name + '*.txt')
        # params : 
        l,Max = name.split('MAX')
        _lambda = int(re.findall(num_pattern,l)[1])
        _max = int(re.findall(num_pattern,Max)[0][:-1])
        # other data ; 
        tracks = [np.loadtxt(f) for f in f_names]
        correlations = [new_auto(t) for t in tracks]
        min_l = min([len(x) for x in correlations])
        # truncate on smallest track :
        correlations = [x[:min_l] for x in correlations]
        # average all 3 correlations: 
        av_cors = np.average(correlations,axis = 0)

        # append rows to create nice dataframe
        for i,x in enumerate(av_cors):
            rows.append([id,i,x,_lambda,_max])
        print(name)

    df = pd.DataFrame(data = rows,columns = ['id','dt','auto correlation','lambda','max act'])
    df.to_csv('../results/autocorrelation_average_nofrc1.csv')

def get_volume(path):
    # test 
    files = glob.glob(path+'*.pkl')
    files.sort()
    print(files)
    volumes_dict = {'FILE':[],'VOLUME':[]}
    for f in files:
        try: 
            track = pickle.load(open(f,'rb'))
            #print(np.where(track != 0.))
            vol = scanned_volume(track)
            print(vol)
            volumes_dict['FILE'].append(f[14:])
            volumes_dict['VOLUME'].append(vol)
            #break #visualize_frc(track)
        except:
            pass
    
    volume_df = pd.DataFrame.from_dict(volumes_dict)
    print(volume_df.head())
    # to text file also :
    volume_df.to_csv('../results/volume_frc1.csv')
    np.savetext('../results/vulume_frc1.txt',volume_df.values,fm = '%d')
    #volume_df.to_csv('../results/volume_frc1.csv')
    
def build_df(files):
    data_dict = {'Cellid': [],'Motility':[],'Persistance':[],'AutoSlope':[],'P-val':[]
    ,'Lambda':[],'Max_act':[],'Speed':[]}
    numbers = '[-+]?\d*\.\d+|\d+'
    for f in files:
        track = np.loadtxt(f)
        displ,_ = analyse_track(track)
        data_dict['Speed'].append(sum(displ)/len(displ))
        m,p,slope,p_val = get_motility(track)
        data_dict['Motility'].append(m)
        data_dict['Persistance'].append(p)
        data_dict['AutoSlope'].append(slope)
        data_dict['P-val'].append(p_val)
        l,Max = f.split('MAX')
        data_dict['Lambda'].append(int(re.findall(numbers,l)[1]))
        data_dict['Max_act'].append(int(re.findall(numbers,Max)[0]))
        data_dict['Cellid'].append(int(re.findall(numbers,Max)[1]))
        print(f)

    df = pd.DataFrame(data_dict)
    df.to_csv('../results/fullCPM_frc1.csv')

if __name__ == "__main__":
    path = '../data/full_LN3/*.txt'
    all_autos(path)
    #get_volume(path)
    #files = glob.glob(path)
    #build_df(files)

    # files = glob.glob(path)
    # files.sort()
    # build_df(files)

    #scanned = get_volume(path)
    # test
    # f = files[18]
    # track = np.loadtxt(f)
    # roseplot(track)
    # auto_cor = new_auto(track)
    # plt.clf()
    # plt.plot(auto_cor)
    # plt.show()
    #print(get_motility(track,plot = True))
