### Come up with some drawing functions : 
#import igraph
import plotly.express as px
import plotly as py
import plotly.graph_objs as go
import pickle
import networkx as nx 
import numpy as np 


def draw_plotly(pos,matrix):
    """ pos should be a list of [[x,y,z],...,[x_n,y_n,z_n]] positions of the nodes
    matrix should be a adjacency matrix of the graph to draw : 
     """
    matrix = np.array(matrix)
    g = nx.from_numpy_matrix(matrix)
    xyz = np.array(pos)
    x,y,z = xyz.T
    degrees = [d[-1]/100 for d in g.degree()]

    x_edges = []
    y_edges = []
    z_edges = []

    for edge in g.edges():
        x0, y0, z0 = pos[edge[0]] 
        x1, y1, z1 = pos[edge[1]]
        x_edges.append(x0)
        x_edges.append(x1)
        x_edges.append(None)
        y_edges.append(y0)
        y_edges.append(y1)
        y_edges.append(None)
        z_edges.append(z0)
        z_edges.append(z1)
        z_edges.append(None)
    
    print(len(list(g.edges())))
    print(len(y_edges))
    
    trace1=go.Scatter3d(x=x_edges,
               y=y_edges,
               z=z_edges,
               mode='lines',
               line=dict(color='rgb(125,125,125)', width=1.5)
               )
    
    trace2=go.Scatter3d(x=x,
               y=y,
               z=z,
               mode='markers',
               marker=dict(symbol='circle',
                             size=7,
                             color=degrees,
                             colorscale='Viridis',
                             line=dict(color='rgb(50,50,50)', width=0.5)
                             )
               )

    axis=dict(showbackground=False,
            showline=False,
            zeroline=False,
            showgrid=False,
            showticklabels=False,
            title=''
            )

    layout = go.Layout(
            title="3D Representation of randomg geometric graph",
            width=1000,
            height=1000,
            showlegend=False,
            scene=dict(
                xaxis=dict(axis),
                yaxis=dict(axis),
                zaxis=dict(axis),
            ))

    data=[trace1,trace2]
    fig=go.Figure(data=data, layout=layout)
    fig.show()

def kamada_draw():
    fname = '/home/lau/GIT/FRC_Thesis/data/exp1/filtered_graphs/graphs_data.pkl'
    graphs = pickle.load(open(fname,'rb'))

    for graph in graphs[:1]:
        name = graph[0]
        stats = graph[1]
        print(name)
        # use saved adj matrix to regenerate graph : 
        matrix = stats[-1].tolist()
        #graph = igraph.Graph.Adjacency(matrix)
        g = igraph.Graph.Adjacency(matrix)
        pos = g.layout_kamada_kawai()    
        igraph.plot(g, layout = pos)
        draw_plotly(pos,matrix)


# kamada_draw()
