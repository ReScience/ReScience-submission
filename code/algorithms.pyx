# -*- coding: utf-8 -*-
#
#    Copyright 2017: Frank Stollmeier
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
# cython: boundscheck=False

import numpy as np
cimport numpy as np
import networkx as nx
import random


# Wrapper classes for random number generation 

cdef extern from "<random>" namespace "std":
    
    cdef cppclass mt19937:
        mt19937()
        mt19937(unsigned int seed)
        void seed(unsigned int seed)
        
    cdef cppclass uniform_real_distribution[T]:
        uniform_real_distribution()
        uniform_real_distribution(T a, T b)
        T operator()(mt19937 gen)

    cdef cppclass uniform_int_distribution[T]:
        uniform_int_distribution()
        uniform_int_distribution(T a, T b)
        T operator()(mt19937 gen)

cdef mt19937 rng = mt19937(1)


# Functions to generate networkx-graphs


def create_square_lattice_graph4(N1,N2, periodic = True):
    '''Create a regular square lattice graph with N1*N2 nodes and k=4 edges per node.'''
    g = nx.grid_graph(dim=[N1,N2],periodic=periodic)
    return g


def create_square_lattice_graph6(N1,N2, periodic = True):
    '''Create a regular square lattice graph with N1*N2 nodes and k=6 edges per node.'''
    g = nx.grid_graph(dim=[N1,N2],periodic=periodic)
    for n in g.nodes():
        g.add_edge(n,((n[0]+1)%N1,(n[1]+1)%N2))
    return g


def create_triangular_lattice_graph(N1,N2):
    '''Create a regular triangular lettice graph with N1*N2 nodes and k=8 edges per node.'''
    g = nx.grid_graph(dim=[N1,N2],periodic=True)
    for i in range(N1):
        for j in range(N2):
            g.add_edge((i,j),((i+1)%N1,(j+1)%N2))
            g.add_edge(((i+1)%N1,j),(i,(j+1)%N2))
    return g


def create_hexagonal_lattice_graph(N1,N2):
    '''Create a regular hexagonal lattice graph with int(N1*N2*2/3.0) nodes and k=3 edges per node.
    Example parameters:
    N1=10 and N2=15 will produce a grid with N=100 nodes.
    N1=25 and N2=30 will produce a grid with N=500 nodes. 
    Note:
    For certain parameters it happens that some nodes on the boundary of the grid have two or four neighbors. The function will print a warning in these cases.
    '''
    if N1<3 or N2<3 or N1%3==2 or N2%3==2:
        #raise ValueError("These values are not allowed.")
        print("Warning: With these values N1 and N2, there may be some nodes on the boundary of the lattice with 2 or 4 edges.\
        If you prefer all nodes having exactly three edges, choose values which do not satisfy N1<3 or N2<3 or N1%3==2 or N2%3==2.")
    g = create_square_lattice_graph6(N1,N2, periodic=True)
    for i in range(N1):
        for j in range(N2):
            if (i+j)%3==0:
                g.remove_node((i%N1,j%N2))
            if i==0 and j%3==2:
                g.remove_edge((i,j),(N1-1,j-1))
            if i%3==2 and j==0:
                u,v = (i,j), (i-1,N2-1)
                if g.has_edge(u,v):
                    g.remove_edge(u,v)
    return g


def create_connected_random_graph(N,k):
    '''Create a connected random graph with N nodes and k edges per node on average.'''
    #Algorithm: 
    #First, add nodes and connect them to already connected nodes to ensure that the graph is connected.
    #Second, add edges between randomly chosen pairs of nodes until the desired number of edges is present.
    g = nx.Graph()
    g.add_node(0)
    nodes = [0]
    for i in range(1,N):
        g.add_node(i)
        g.add_edge(i,random.choice(nodes))
        nodes.append(i)
    n_edges = np.array([len(g.edges(n)) for n in g.nodes()])
    while np.mean(n_edges) <= k:
        n1,n2 = random.randint(0,N-1), random.randint(0,N-1)
        if n1 != n2 and not g.has_edge(n1,n2):
            g.add_edge(n1,n2)
            n_edges[n1] = n_edges[n1] + 1
            n_edges[n2] = n_edges[n2] + 1
    return g


def convert_networkx_graph_to_adjacency_list(G):
    '''Convert a networkx graph G to an adjacency list stored in a numpy array of shape (N,N+1) 
    in which row i stores the indices of the neighbors of node i followed by -1 values until the end of the row.
    '''
    adj_list = - np.ones((G.number_of_nodes(),G.number_of_nodes()+1),dtype=np.int)
    index_dict = dict([(n,i) for i,n in enumerate(G.nodes_iter())])
    for node in G.nodes_iter():
        edges = G.edges(node)
        adj_list[index_dict[node],:len(edges)] = [index_dict[e[1]] for e in edges]  
    return adj_list





# Function to measure fixation probabilities

def measure_fixation_probabilities(graph_generator, graph_generator_arguments, w, c, n_graphs, n_runs, benefit_to_cost_ratios, dview = None, iteration_type = 1):
    ''' For each value in benefit_to_cost_ratios, do n_graphs*n_runs simulations to measure the fixation probabilities.
    
    Parameters:
    graph_generator: Function which returns a graph as an adjacency list. This adjacency list should be a numpy array of shape (N,N+1) in which row i stores the indices of the neighbors of node i followed by -1 values until the end of the row.
    graph_generator_arguments: List of arguments passed to the graph_generator function.
    w:  Selection strength.
    c:  Costs.
    n_graphs: The number of graphs generated with graph_generator.
    n_runs: The number of simulations on each graph.
    benefit_to_cost_ratios: List of b/c ratios.
    dview (optional): ipyparallel.client.view.DirectView object for parallelization. The default value is None, meaning no parallelization.
    iteration_type (optional): Either 1 (death-birth process) or 0 (imitation process). Default is 1. 
    
    Returns:
    Fixation probabilities as a numpy-array of the same length as benefit_to_cost_ratios.
    '''
    fixation_probabilities = []
    seeds = np.linspace(1, 2**32-1, n_graphs, dtype=np.uint32)
    for bc_ratio in benefit_to_cost_ratios:
        b = c * bc_ratio
        if dview is None:
            results = [measure_fixation_probability(graph_generator(*graph_generator_arguments), n_runs, b, c, w, s, iteration_type) for s in seeds]
        else:
            dview.push(dict(graph_generator=graph_generator, graph_generator_arguments=graph_generator_arguments, n_runs=n_runs, b=b, c=c, w=w, iteration_type=iteration_type))
            dview.scatter('seeds',seeds)
            dview.execute('results = [measure_fixation_probability(graph_generator(*graph_generator_arguments), n_runs, b, c, w, s, iteration_type) for s in seeds]')
            results = dview.gather('results')
        fixation_probabilities.append(np.mean(results))
    return np.array(fixation_probabilities) 


def measure_fixation_probability(np.ndarray[np.int64_t, ndim=2] graph_nd, int n_runs, float b, float c, float w, unsigned int seed, int iteration_type):
    '''Measure fixation probability by running many simulations starting with a random node as the initial cooperator. 
    
    Parameters:
    graph_nd: The graph given as a numpy array of shape (N,N+1) in which row i stores the indices of the neighbors of node i followed by -1 values until the end of the row.
    n_runs: Number of simulations.
    b: Benefit value.
    c: Cost value.
    w: Selection strength.
    seed: Seed for the random number generator.
    iteration_type: Either 1 (death-birth process) or 0 (imitation process).
    
    Returns:
    The ratio of the number of simulations that ended with a fixation of cooperators.
    '''
    cdef long int [:,:] graph = graph_nd
    cdef int N = graph.shape[0]
    cdef int first_cooperator_index
    rng.seed(seed)
    cdef uniform_int_distribution[int] dist = uniform_int_distribution[int](0,N-1)
    cdef int C_fixations = 0
    cdef int i = 0
    for i in range(n_runs):
        first_cooperator_index = dist(rng)
        if run_until_fixation(graph, N, first_cooperator_index, b, c, w, iteration_type):
            C_fixations += 1
    return C_fixations / float(n_runs)




# Algorithms to simulate games on graphs

cdef bint run_until_fixation(long int [:,:] graph, int N, int cooperator_index, float b, float c, float w, int iteration_type):
    '''Run a simulation until all nodes have the same strategy. In the initial state all nodes defect except the one that is specified with cooperator_index.
    
    Parameters:
    graph: The graph given as a numpy array of shape (N,N+1) in which row i stores the indices of the neighbors of node i followed by -1 values until the end of the row.
    N: Size of the graph.
    cooperator_index: Index of the initial cooperator.
    b: Benefit value.
    c: Cost value.
    w: Selection strength.
    iteration_type: Either 1 (death-birth process) or 0 (imitation process).
    
    Returns:
    True or false, depending on whether the fixation of cooperators was successful.  
    '''
    cdef int i = 0
    cdef int j = 0
    cdef uniform_int_distribution[int] dist1 = uniform_int_distribution[int](0,N-1)
    cdef uniform_real_distribution[double] dist2 = uniform_real_distribution[double](0.0,1.0)

    #Initialization
    cdef long int [:] strategies = np.zeros(N, dtype=np.int)
    strategies[cooperator_index] = 1
    cdef int neighbors_index
    cdef long int [:] adjacent_cooperators = np.zeros(N, dtype=np.int)
    cdef long int [:] adjacent_defectors = np.zeros(N, dtype=np.int)
    cdef int adj_cooperators
    cdef int adj_defectors
    for i in range(N):
        adj_cooperators = 0
        adj_defectors = 0
        j = 0
        neighbors_index = graph[i,j]
        while neighbors_index >= 0:
            if strategies[neighbors_index]:
                adj_cooperators += 1
            else:
                adj_defectors += 1
            j += 1
            neighbors_index = graph[i,j]
        adjacent_cooperators[i] = adj_cooperators
        adjacent_defectors[i] = adj_defectors
    
    #Simulation
    cdef int n_cooperators = 1
    cdef int counter = 0
    cdef int update_node_index
    cdef double random_number
    if iteration_type: #death-birth
        iteration = iteration_DB
    else:              #imitation
        iteration = iteration_IM
    while n_cooperators > 0 and n_cooperators < N:
        update_node_index = dist1(rng)
        random_number = dist2(rng)
        n_cooperators += iteration(graph, strategies, adjacent_cooperators, adjacent_defectors, update_node_index, random_number, b, c, w)
        counter += 1
    return n_cooperators == N




cdef int iteration_DB(long int [:,:] graph, long int [:] strategies, long int [:] adjacent_cooperators, long int [:] adjacent_defectors, int update_node_index, float random_number, float b, float c, float w) nogil:
    '''Do one death-birth update.
    
    Parameters:
    graph: The graph given as a numpy array of shape (N,N+1) in which row i stores the indices of the neighbors of node i followed by -1 values until the end of the row.
    strategies: Integer array of shape (N) which stores the strategy of each node. If node i is a cooperator(defector), then strategies[i] is 1(0).
    adjacent_cooperators: Integer array of shape (N) which stores the number of adjacent cooperators of each node.
    adjacent_defectors: Integer array of shape (N) which stores the number of adjacent defectors of each node.
    update_node_index: Index of the node where the death and the birth event happen.
    random_number: A float value between 0 and 1.
    b: Benefit value.
    c: Cost value.
    w: Selection strength.
        
    Returns:
    -1, 0 or +1, depending on whether the number of cooperators decreased by one, stayed constant, or increased by one.
    '''
    cdef int previous_strategy = strategies[update_node_index]
    
    #Sum up fitness of adjacent neighbors
    cdef float F_C = 0
    cdef float F_D = 0
    cdef float f
    cdef int neighbors_index = graph[update_node_index,0]
    cdef int i = 0
    while neighbors_index >= 0:
        f = fitness(strategies[neighbors_index], adjacent_cooperators[neighbors_index], adjacent_defectors[neighbors_index], b, c, w)
        if strategies[neighbors_index]:
            F_C += f
        else:
            F_D += f
        i += 1
        neighbors_index = graph[update_node_index,i]
    
    #Birth step
    cdef int diff_n_cooperators = 0
    if random_number < F_C / (F_C+F_D):
        strategies[update_node_index] = 1
        if not previous_strategy:
            diff_n_cooperators = +1
    else:
        strategies[update_node_index] = 0
        if previous_strategy:
            diff_n_cooperators = -1
    
    #Update the stored information about the neighbors of each node
    if diff_n_cooperators != 0:
        i = 0
        while True:
            neighbors_index = graph[update_node_index,i]
            if neighbors_index < 0:
                break
            adjacent_cooperators[neighbors_index] += diff_n_cooperators
            adjacent_defectors[neighbors_index] -= diff_n_cooperators
            i += 1
    
    return diff_n_cooperators


cdef int iteration_IM(long int [:,:] graph, long int [:] strategies, long int [:] adjacent_cooperators, long int [:] adjacent_defectors, int update_node_index, float random_number, float b, float c, float w) nogil:
    '''Do one imitation step.
    
    Parameters:
    graph: The graph given as a numpy array of shape (N,N+1) in which row i stores the indices of the neighbors of node i followed by -1 values until the end of the row.
    strategies: Integer array of shape (N) which stores the strategy of each node. If node i is a cooperator(defector), then strategies[i] is 1(0).
    adjacent_cooperators: Integer array of shape (N) which stores the number of adjacent cooperators of each node.
    adjacent_defectors: Integer array of shape (N) which stores the number of adjacent defectors of each node.
    update_node_index: Index of the node where the death and the birth event happen.
    random_number: A float value between 0 and 1.
    b: Benefit value.
    c: Cost value.
    w: Selection strength.
        
    Returns:
    -1, 0 or +1, depending on whether the number of cooperators decreased by one, stayed constant, or increased by one.
    '''
    cdef int previous_strategy = strategies[update_node_index]
    
    #Sum up fitness of adjacent neighbors
    cdef float F_C = 0
    cdef float F_D = 0
    cdef float f
    cdef int neighbors_index = graph[update_node_index,0]
    cdef int i = 0
    while neighbors_index >= 0:
        f = fitness(strategies[neighbors_index], adjacent_cooperators[neighbors_index], adjacent_defectors[neighbors_index], b, c, w)
        if strategies[neighbors_index]:
            F_C += f
        else:
            F_D += f
        i += 1
        neighbors_index = graph[update_node_index,i]
    
    #Imitation
    cdef int diff_n_cooperators = 0
    cdef float f0 = fitness(strategies[update_node_index], adjacent_cooperators[update_node_index], adjacent_defectors[update_node_index], b, c, w)
    if not previous_strategy and random_number < F_C / (F_C+F_D+f0):
        strategies[update_node_index] = 1
        diff_n_cooperators = +1
    elif previous_strategy and random_number < F_D / (F_C+F_D+f0):
        strategies[update_node_index] = 0
        diff_n_cooperators = -1
    
    #Update the stored information about the neighbors of each node
    if diff_n_cooperators != 0:
        i = 0
        while True:
            neighbors_index = graph[update_node_index,i]
            if neighbors_index < 0:
                break
            adjacent_cooperators[neighbors_index] += diff_n_cooperators
            adjacent_defectors[neighbors_index] -= diff_n_cooperators
            i += 1
    return diff_n_cooperators


cdef inline float fitness(int strategy, int n_cooperators, int n_defectors, float b, float c, float w) nogil:
    '''Calculate the fitness of a node. 
    
    Parameters:
    strategy: Integer value, 1 (0) if the node is a cooperator (defector).
    n_cooperators: The number of adjacent nodes that cooperate.
    n_defectors: The number of adjacent nodes that defect.
    b: Benefit value.
    c: Cost value.
    w: Selection strength.
    '''
    cdef float payoff
    if strategy:
        payoff = b*n_cooperators - c*(n_cooperators + n_defectors)
    else:
        payoff = b*n_cooperators
    return 1 - w + w*payoff

