import cpm
import numpy as np
from full_ln2 import setup_frc,setup
import pickle
import pandas as pd
import time
import random
import matplotlib.pyplot as plt
import glob
import matplotlib.image as mpimg
import matplotlib.animation as animation
import imageio
import re 

def simulate():
    sim = setup(2000,800)
    state = sim.get_act_state()
    for i in range(500):
        print('starting simulation')
        sim.run(10)
        # write state slice to disc: 
        cube = state % 2 ** 24
        cut = np.sum(cube[32:36],axis = 0)
        # enlarge : 
        frame = np.kron(cut,np.ones((4,4)))
        frame[frame == 0] = 1000
        plt.imshow(frame)
        plt.savefig('../data/animation/2000_800_actfr' + str(i) + '.png')


    
def animate(i):
    files = glob.glob('/home/lau/Desktop/Thesis Stuff/anim_data/2000_800_fr*')
    files.sort()
    im.set_array(mpimg.imread(files[i]))
    return [im]

def create_anim():
    fps = 5
    nSeconds = 100
    files = glob.glob('/home/lau/Desktop/Thesis Stuff/anim_data/2000_800_fr*')
    files.sort()
    #print(files)
    fig = plt.figure(figsize=(20,20))
    f1 = mpimg.imread(files[0])
    global im
    im = plt.imshow(f1)

    anim = animation.FuncAnimation(fig,animate) #,frames = fps * nSeconds,interval= 200

    anim.save('full_ln_anim.mp4')

def create_anim2():
    filenames = glob.glob('/home/lau/Desktop/Thesis Stuff/anim_data/2000_800_actfr*')
    filenames.sort(key=lambda var:[int(x) if x.isdigit() else x for x in re.findall(r'[^0-9]|[0-9]+', var)])
    with imageio.get_writer('full_lnact_2000_800_anim.mp4', mode='I') as writer:
        for filename in filenames:
            print(filename)
            image = imageio.imread(filename)
            writer.append_data(image)

if __name__ == '__main__':
#    simulate()
    #im = None
    create_anim2()
