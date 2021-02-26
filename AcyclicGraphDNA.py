import numpy as np
import itertools

import networkx as nx
import pymetis

from graphviz import Digraph

import re


def sign(node, edge):
    if node == edge[0]:
        return -1 #output edge
    elif node == edge[1]:
        return 1 #input edge
    else:
        return 0 # node not in edge
    
def get_edges_for_node(node, edges):
    node_edges = []
    for edge in edges:
        if node in edge:
            node_edges.append(edge)
    return node_edges

#!This transformation could be used when Longest path=Hamiltonian path and only for acyclic graphs!
def to_qubo(adjacency_matrix, A=1):
    all_edges = np.argwhere(adjacency_matrix==1).tolist()
    dim = len(all_edges)
    nodes_quantity = len(adjacency_matrix)
    M = np.zeros((dim, dim))
    
    for node in range(nodes_quantity):
        node_edges = get_edges_for_node(node, all_edges)
        node_edges_plus = [edge for edge in node_edges if sign(node, edge)==1]
        node_edges_minus = [edge for edge in node_edges if sign(node, edge)==-1]

        # interaction of edges, sign depends on node direction: -1 output, +1 input
        for edges in itertools.combinations(node_edges_plus, 2):      
            i, j = all_edges.index(edges[0]), all_edges.index(edges[1])
            M[i][j] += sign(node, edges[0])*sign(node, edges[1])*A
            M[j][i] += sign(node, edges[0])*sign(node, edges[1])*A
            
        for edges in itertools.combinations(node_edges_minus, 2):      
            i, j = all_edges.index(edges[0]), all_edges.index(edges[1])
            M[i][j] += sign(node, edges[0])*sign(node, edges[1])*A
            M[j][i] += sign(node, edges[0])*sign(node, edges[1])*A

        # diagonal part
        for edge in node_edges_plus:
            i = all_edges.index(edge)
            M[i][i] -= A
        # diagonal part
        for edge in node_edges_minus:
            i = all_edges.index(edge)
            M[i][i] -= A
    
    return M


def load_file(filename):
    fd = open(filename , 'r')

    lines = []
    segments = []
    links = []
    containments = []
    for line in fd:
        if line[0] == 'S':
            segments.append(line[1:])
        if line[0] == 'L':
            links.append(line[1:])
        if line[0] == 'C':
            containments.append(line[1:])
        lines.append(line)

    print("segments: %s, links: %s, containments: %s, lines: %s" %(len(segments), len(links), len(containments), len(lines)))
    
    return segments, links, containments


# get all nodes
def get_nodes(segments):
    nodes = []
    for segment in segments:
        name = re.search('[!-)+-<>-~][!-~]*', segment)
        nodes.append(name[0]+'+')
        nodes.append(name[0]+'-')
    print('\nnodes quantity:', len(nodes))
    nodes_labels = dict(enumerate(nodes))
    return nodes, nodes_labels


# get all links
def get_edges(links, nodes):
    dim = len(nodes)
    adjacency_matrix = np.zeros(shape=(dim, dim))
    swap = lambda sign: '+' if sign=='-' else '-'
    
    links_edges = {}
    for link in links:
        l = re.findall('[!-)+-<>-~][!-~]*', link)
        node1 = l[0]+l[1]
        node2 = l[2]+l[3]
        if node1 in nodes:
            i, j = nodes.index(node1), nodes.index(node2)
            adjacency_matrix[i, j] = 1
            links_edges[(i, j)] = [node1, node2]

            #adding strand
            node1 = l[2]+swap(l[3])
            node2 = l[0]+swap(l[1])
            i, j = nodes.index(node1), nodes.index(node2)
            adjacency_matrix[i, j] = 1
            links_edges[(i, j)] = [node1, node2]
            
    return adjacency_matrix, links_edges


def get_initial_graph(segments, links):
    nodes, nodes_labels = get_nodes(segments)
    adjacency_matrix, links_edges = get_edges(links, nodes)
    
    edges = list(links_edges.values())
    graph = nx.DiGraph()
    graph.add_edges_from(edges)
    
    return graph, adjacency_matrix, nodes, nodes_labels, links_edges


def get_part_graph(nodes, links):
    G = nx.DiGraph()

    for link in links.values():
        node1 = link[0]
        node2 = link[1]
        if (node1 in nodes) and (node2 in nodes):
            G.add_edge(node1, node2)

    return G


def get_parted_graphs(part_map, node_labels, edges):
    part_graphs = []
    
    for part, part_nodes in part_map.items():
        G = nx.DiGraph()

        for edge in edges:
            node1 = edge[0]
            node2 = edge[1]
            
            labels = [node_labels[i] for i in part_nodes]
            if (node1 in labels) and (node2 in labels):
                G.add_edge(node1, node2)

        part_graphs.append(G)     
    
    return part_graphs


def get_parted_graphs(part_map, nodes_labels, edges):
    part_graphs = []
    
    for part, part_nodes in part_map.items():
        G = nx.DiGraph()

        for edge in edges:
            node1 = edge[0]
            node2 = edge[1]
            
            labels = [nodes_labels[i] for i in part_nodes]
            if (node1 in labels) and (node2 in labels):
                G.add_edge(node1, node2)

        part_graphs.append(G)     
    
    return part_graphs


def part_graph(cluster_number, adjacency_matrix, log=False):
    
    adjacency_list = adjacency_matrix_to_adjacency_list(adjacency_matrix)
    cut_count, part_vert = pymetis.part_graph(cluster_number, adjacency=adjacency_list)

    # Find nodes in each partition:
    part_map = {}
    for i, p in enumerate(set(part_vert)):
        part_map[p] = np.argwhere(np.array(part_vert) == p).ravel()

    if log:
        for part, nodes in part_map.items():
            print('part %s, nodes quantity: %s' %(part, len(nodes)))
            print(nodes)

    return part_map, cut_count, part_vert


def draw_graph(node_labels, solve_edges, edges):
    size = len(node_labels)
    node_cnt = len(node_labels)
    g = Digraph() 
    
    for i in range(size):
        g.node(node_labels[i])

    for edge in edges:
        if edge not in solve_edges:
            g.edge(node_labels[edge[0]], node_labels[edge[1]], color = 'grey')

    for edge in solve_edges:
        g.node(node_labels[edge[0]], style='filled', fillcolor='white')
        g.node(node_labels[edge[1]], style='filled', fillcolor='white')

        g.edge(node_labels[edge[0]], node_labels[edge[1]], color = 'lightskyblue', label = '', fontsize='20', fontcolor='red')

    return g

def get_parted_graphs_graphiz(nodes, node_labels, edges):
    part_graphs = []
    colors = ['yellow', 'red', 'green', 'orange', 'brown']
    for part, part_nodes in nodes.items():
        G = Digraph() 

        for i in part_nodes:
            G.node(node_labels[i], style='filled', fillcolor=colors[part])
        for edge in edges:
            node1 = edge[0]
            node2 = edge[1]
            labels = [node_labels[i] for i in part_nodes]
            if (node1 in labels) and (node2 in labels):
                G.edge(node1, node2, color='black')

        part_graphs.append(G)     
    
    return part_graphs


def adjacency_matrix_to_adjacency_list(adjacency_matrix):
    adjacency_list = []
    indices = np.argwhere(adjacency_matrix==1)
    
    for node in range(len(adjacency_matrix)):
        adjacency_list.append([])
        for index in indices:
            if index[0] == node:
                adjacency_list[node].append(index[1])
    adjacency_list = [np.asarray(el) for el in adjacency_list]
    return adjacency_list    


def complement_to_undirected_graph(graph=None, edges=None):
    if edges == None:
        edges = graph.edges()
    undirected_edges = []
    for edge in edges:
        undirected_edges.append(edge)
        edge_reverse = [edge[1], edge[0]]
        undirected_edges.append(edge_reverse)

    graph = nx.DiGraph()
    graph.add_edges_from(undirected_edges)
        
    return graph


def get_strand_graph(graph):
    undirected_graph = complement_to_undirected_graph(graph)
    adjacency_matrix = nx.to_numpy_matrix(undirected_graph)
    nodes_labels = list(undirected_graph.nodes())

    #part graph into two
    part_map, _ , _ = part_graph(cluster_number=2, adjacency_matrix=adjacency_matrix)

    #take one strand
    strand_graph = get_parted_graphs(part_map, nodes_labels, edges=graph.edges)[0]
    return strand_graph


def get_source_and_sink(nodes, edges, log=False):
    node_degrees = {}
    for node in nodes:
        input_degree = 0
        output_degree = 0
        for edge in edges:
            if node==edge[1]:
                input_degree += 1
            elif node==edge[0]:
                output_degree += 1

            node_degrees[node] = [input_degree, output_degree]

    result = []
    for node, degree in node_degrees.items():
        if degree[0] == 0:
            if log:
                print('node %s has degree 1 in input, index %s' %(node, nodes.index(node)))
            result.append(node)
        elif degree[1] == 0:
            if log:
                print('node %s has degree 1 in output, index %s' %(node, nodes.index(node)))
            result.append(node)
    return result


def get_unconnected_nodes(adjacency_matrix):
    idx_row = np.where(~adjacency_matrix.any(axis=0))[0]
    idx_columns = np.where(~adjacency_matrix.any(axis=1))[0]
    print('rows to delete:', idx_row)
    print('columns to delete:', idx_columns)
    idx = np.intersect1d(idx_row, idx_columns)
    print('indices to delete', idx)
    
    return idx


def drop_isolated_nodes_and_update(adjacency_matrix, nodes, links):
 
    idx = get_unconnected_nodes(adjacency_matrix)
    print("idices of unconnected nodes", idx)

    for i in np.flip(idx.flatten()):
        print("node %s with index %s deleted" %(nodes[i], i))
        del nodes[i]

    nodes_labels = dict(enumerate(nodes))
    adjacency_matrix = np.delete(adjacency_matrix, idx, axis=1) # drop zero rows
    adjacency_matrix = np.delete(adjacency_matrix, idx, axis=0) # drop zero columns

    print('adjacency_matrix.shape:', adjacency_matrix.shape)
    #print(nodes_labels)

    adjacency_matrix, links_edges = get_edges(links, nodes)
    print('Links quantity: ', len(links_edges))
    
    return nodes, nodes_labels, links_edges


def draw_solution(graphs, strand_graph, strand_graph_undirected, part_map, solve_edges, show_all_edges=True):
    bridge_nodes = []
    for pg in graphs:
        bridge_nodes += get_source_and_sink(list(pg.nodes), list(pg.edges))
    bridge_nodes

    node_labels = list(strand_graph_undirected.nodes)
    edges = strand_graph.edges
    bridge_eges = []
    
    for node1 in bridge_nodes:
        for node2 in bridge_nodes:
            if (node1, node2) in edges:
                bridge_eges.append((node1, node2))

    graph = Digraph() 
    colors = ['yellow', 'red', 'green', 'orange', 'brown']
    for part, nodes in part_map.items():
        for i in nodes:
            graph.node(node_labels[i], style='filled', fillcolor=colors[part])

    if show_all_edges:
        for edge in edges:
            if edge not in solve_edges and edge not in bridge_eges:
                graph.edge(edge[0], edge[1], color='black')

    for edge in solve_edges:
        graph.edge(edge[0], edge[1], color = 'lightskyblue', label = '', fontsize='20', fontcolor='red')

    for node1 in bridge_nodes:
        for node2 in bridge_nodes:
            if (node1, node2) in edges:
                graph.edge(node1, node2, color = 'red', label = '', fontsize='20', fontcolor='red')    
    
    return graph
