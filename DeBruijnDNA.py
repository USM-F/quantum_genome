import numpy as np
from graphviz import Digraph

def draw_graph(adj, node_labels, path_spins = [0] * 10, kmer_len = 3, suffix_len = 2):
    size = len(adj)
    node_cnt = len(adj)
    g = Digraph() 
    
    def get_pos(node_idx):
        sub_list = path_spins[node_idx * node_cnt: (node_idx + 1) * node_cnt]
        if 1 in sub_list:
            return sub_list.index(1)
        else:
            return 0
    nodes = {}
    for i in range(size):
        g.node(node_labels[i])
        for j in range(size):
            idx1 = get_pos(i)
            idx2 = get_pos(j)
            c = "lightskyblue"
            l = ""
            if (idx2 == (idx1 + 1)):
                c = "black"
                l = str(idx2)
            if adj[i][j] == 1:
                g.edge(node_labels[i], node_labels[j], color = c, label = l, fontsize='20', fontcolor='red')
                if l!="":
                    nodes[int(l)] = [node_labels[i], node_labels[j]]

    return g, nodes
        
'''
Random De Brujin graph
Input:
 . Nucleotide sequence (optional) : "ACTC"
Output:
 . adjacency matrix: [[0,1,0],[0,0,1],[0,0,0]]
 . node_labels: {0 : "ACT", 1 : "CTC"}
'''
def make_debr(seq = "", seq_len = 10, kmer_len = 3, suffix_len = 2):
    if seq == "":
        seq = random_seq(seq_len, kmer_len)
    else:
        seq_len = len(seq)

    nodes = []
    edges = set()
    node_labels = {}
    
    prefix_map = {}
    suffix_map = {}
    
    def add_edge(kmer, prefix):
        if prefix in prefix_map:
            for next_kmer in prefix_map[prefix]:
                edges.add((kmer, next_kmer))

    for i in range(0, seq_len - kmer_len + 1, kmer_len-suffix_len):
        kmer = seq[i:i+kmer_len]
        nodes.append(kmer)
        node_labels[len(node_labels)] = kmer
        prefix = kmer[0:suffix_len]
        suffix = kmer[-suffix_len:]
        if prefix not in prefix_map:
            prefix_map[prefix] = set()
        if suffix not in suffix_map:
            suffix_map[suffix] = set()
        prefix_map[prefix].add(kmer)
        suffix_map[suffix].add(kmer)
        
    for i in range(0, seq_len - kmer_len + 1, kmer_len-suffix_len):
        kmer = seq[i:i+kmer_len]
        next_prefix = seq[i+kmer_len-suffix_len:i+kmer_len]
        add_edge(kmer, next_prefix)

    node_cnt = len(nodes)
    adj = np.zeros((node_cnt, node_cnt))
    for edge in edges:
        idx_1 = nodes.index(edge[0])
        idx_2 = nodes.index(edge[1])
        adj[idx_1][idx_2] = 1
        
    return adj, node_labels

'''
Convert adjacency graph into QUBO matrix
Input:
 . adjacency matrix : 
 . ttype : transformation type ("hamiltonian_path" or "longest_path")
Output:
 . QUBO matrix : 
'''
def to_qubo(adj):
    node_cnt = len(adj)
    pos_cnt = len(adj)

    q_size = node_cnt * pos_cnt

    def spin_idx(node_idx, pos):
        return node_idx * node_cnt + pos

    Q = np.zeros((q_size, q_size))
    p = 1
    for i in range(node_cnt):
        for j in range(pos_cnt):
            idx1 = spin_idx(i,j)
            for k in range(pos_cnt):
                idx2 = spin_idx(i,k)
                Q[idx1][idx2] = p
                Q[idx2][idx1] = p

            for k in range(node_cnt):
                idx2 = spin_idx(k,j)
                Q[idx1][idx2] = p
                Q[idx2][idx1] = p
            Q[idx1][idx1] = -p

    for i in range(pos_cnt - 1):
        for j in range(node_cnt):
            idx1 = spin_idx(j, i)
            for k in range(node_cnt):
                idx2 = spin_idx(k, i + 1)
                if adj[j,k] != 1:
                    Q[idx1][idx2] = p
                    Q[idx2][idx1] = p

    return Q
