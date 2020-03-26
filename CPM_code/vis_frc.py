import plotly.graph_objects as go
import numpy as np
import networkx as nx
import sys
import pickle
sys.path.insert(0,'../Code')
from mayavi import mlab
from ThreeDdraw import draw_plotly
from Bresenheim import *
from CPM_helpers1 import *
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

# def visualize_frc2(cube):
#     X,Y,Z = np.mgrid[0:256:256j,0:256:256j,0:256:256j]
#     fig = go.Figure(data = go.Isosurface(
#         x = X.flatten(),
#         y = Y.flatten(),
#         z = Z.flatten(),
#         value = cube.flatten(),
#         isomin = .2,
#         isomax = .7,
#         opacity = .1,
#         surface_count = 10))
#     plot(fig)#,auto_open = False,filename = 'testdat/nofrc1.html'

def visualize_frc(frc):
    fig = go.Figure(data=[go.Mesh3d(x=frc[0], y=frc[1], z=frc[2],
                   alphahull=5,
                   opacity=0.4,
                   color='cyan')])
    fig.show()

if __name__ == '__main__':
    #cube = test(256)
    # sim = single_cell_setup1(256,(128,128,128),5000,500,FRC = True)
    # state = sim.get_state()
    # frc = np.where(state // 2**24 == 1)
    frc = pickle.load(open('../results/VOLLAMBDA_1893MAX6290.pkl','rb'))
    print(len(np.where(frc !=0)[0]))
    mlab.contour3d(frc)
    mlab.show()
    #visualize_frc(frc)
