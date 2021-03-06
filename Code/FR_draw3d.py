### FRUCHTERMANN REINGOlD (slow) ###

# HELPER FUNCTIONS : 
def optimal_dist(C,N):
    area = 1
    C = .2
    return C * (area/N)**(1/3)

def dist(a,b):
    """ a and b are both points in 3 d space np.array([x,y,z])
    """
    return numpy.linalg.norm(a-b)

def f_attract(distance,optimal_dist):
    return distance**2/optimal_dist

def f_repulse(distance,optimal_dist):
    C = .2
    return  C * (optimal_dist**2) / distance # maybe a minus ..

def fruchterman_reingold(g,num_iters):
    """ g = networkx graph 
        assume network to be drawn in 3d unit cube
    """
    C = .2
    opt_dist = optimal_dist(C,g.number_of_nodes())

    #opt_dist = np.array([optim_dist] * 3)
    temp = 1/10
    # assign inintial random possitions :
    for _,n in g.nodes(data = True):
        n['pos'] = np.random.rand(3)

    # main loop:
    for i in range(num_iters):
        # calculate repulsive forces : 
        for i1,n1 in g.nodes(data = True):
            n1['disp'] = 0
            for i2,n2 in g.nodes(data = True):
                if i1 != i2:
                    d = n1['pos'] - n2['pos']
                    #n1['disp'] = f_repulse(d,opt_dist)
                    n1['disp'] = n1['disp'] + (d/norm(d)) * f_repulse(norm(d),opt_dist)

        # calculate attractive forces between edges : 
        for E in g.edges(data = True):
            n1,n2 = g.nodes[E[0]],g.nodes[E[1]]
            d = n1['pos'] - n2['pos']
            #n1['disp'] = f_attract(d,opt_dist)
            #n2['disp'] = f_attract(norm(d),opt_dist)
            n1['disp'] = n1['disp'] - (d/norm(d)) * f_attract(norm(d),opt_dist)
            n2['disp'] = n2['disp'] + (d/norm(d)) * f_attract(norm(d),opt_dist)
        
        # limit displacement with temperature
        for _,n in g.nodes(data = True):
            final_displacement = (n['disp']/norm(n['disp'])) * min([norm(n['disp']),temp])
            # check if optimal distance is already achieved :
            if np.isnan(final_displacement).any():
                pass
            else:
                n['pos'] = n['pos'] + final_displacement
            # check if position is not out of bounds : 
            for j in range(3):
                if n['pos'][j] > 1.:
                    n['pos'][j] = 1.
                elif n['pos'][j] < 0.:
                    n['pos'][j] = 0.
                elif np.nan in n['pos']:
                     n['pos'] = np.random.rand(3)

        temp += - 1/(10*num_iters)
        #print(temp)
        print(i)
    return g


def test_3d():
    g = nx.random_geometric_graph(200,.24, dim = 3)
    print(g.number_of_edges())
    print(nx.is_connected(g))
    while not nx.is_connected(g):
      g = nx.random_geometric_graph(200,.25, dim = 3)
    #g = nx.erdos_renyi_graph(170,.15)
    print(nx.is_connected(g))
    matrix = nx.to_numpy_matrix(g)
    g = fruchterman_reingold(g,100)
    positions = []
    for n,data in g.nodes(data = True):
        positions.append(data['pos'])
    
    draw_plotly(positions,matrix)
    return positions