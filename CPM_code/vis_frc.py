import plotly.graph_objects as go
import numpy as np
import networkx as nx
import sys
sys.path.insert(0,'../Code')
from ThreeDdraw import draw_plotly
from Bresenheim import *
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

def visualize_frc(cube):
    X,Y,Z = np.mgrid[0:256:256j,0:256:256j,0:256:256j]
    fig = go.Figure(data = go.Isosurface(
        x = X.flatten(),
        y = Y.flatten(),
        z = Z.flatten(),
        value = cube.flatten(),
        isomin = .2,
        isomax = .7,
        opacity = .1,
        surface_count = 25))
    plot(fig,auto_open = False,filename = 'testdat/frc1.html')

if __name__ == '__main__':
    cube = test(256)
    visualize_frc(cube)
