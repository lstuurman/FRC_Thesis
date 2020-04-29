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

def simulate():
    sim = setup(2000,800)
    state = sim.get_act_state()
    for i in range(10):
        print('starting simulation')
        sim.run(10)
        # write state slice to disc: 
        cube = state % 2 ** 24
        cut = np.sum(cube[32:36],axis = 0)
        # enlarge : 
        frame = np.kron(cut,np.ones((4,4)))
        frame[frame == 0] = 1000
        plt.imshow(frame)
        plt.savefig('../data/animation/fr' + str(i) + '.png')


    
def animate(i):
    files = glob.glob('../data/animation/fr*')
    im.set_array(mpimg.imread(files[i]))
    return [im]

def create_anim():
    fps = 30
    nSeconds = 10
    files = glob.glob('../data/animation/fr*')

    fig = plt.figure(figsize=(20,20))
    f1 = mpimg.imread(files[0])
    global im
    im = plt.imshow(f1)

    anim = animation.FuncAnimation(fig,animate,
        frames = 10,interval= 1000/fps)

    anim.save('full_ln_anim.mp4')

if __name__ == '__main__':
    simulate()
    # im = None
    # create_anim()