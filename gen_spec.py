import re
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

fac_dir = '/home/tim/research/cuebit_col/fac_output/'

def midpoints(x):
    return (x[:-1]+x[1:])/2

def gaussian(x, A, mu, sigma):
    return A * 1/(sigma*(2*math.pi)**.5) * np.exp(-(x-mu)**2/(2*sigma**2))

def read_tr(trans_file):
    with open(trans_file,'r') as f:
        lines = f.readlines()
        for i,line in enumerate(lines):
            nt_re = re.match(r'NTRANS\t= ([0-9]+)',line)
            if nt_re:
                num_ens = int(nt_re[1])
                start = i
                break
        tr_arr = np.array([line.strip().split() for line in lines[start+4:]])

    tr_arr = tr_arr[:,[0,2,4,6]].astype(float)
    # upper_index, lower_ind, energy, A

    tr_arr = np.column_stack((tr_arr,tr_arr[:,3]))
    for i,line in enumerate(tr_arr):
        tr_arr[i,4] = line[4]/np.sum(tr_arr[(tr_arr[:,0]==line[0]),3])
    tree = {int(path[0]):tr_arr[(tr_arr[:,0]==path[0]),:] for path in tr_arr}
    non_leaf_nodes = np.unique(tr_arr[:,:2],axis=0).astype(int)
    non_leaf_nodes = set([node[0] for node in non_leaf_nodes])

    return tree, non_leaf_nodes

def cx_xss(file, en):
    with open(file, 'r') as f:
        lines = f.readlines()
        for i,line in enumerate(lines):
            ne_re = re.match(r'NE0\t= ([0-9]+)',line)
            if ne_re:
                num_ens = int(ne_re[1])
                start = i
                break

        ens = np.array([float(line.strip()) for line in lines[start+1:start+1+num_ens]])
        states = np.array([line.strip().split() for line in lines[start+1+num_ens::num_ens+1]])
        en_ind = np.argmin(abs(ens-en))
        xss = np.array([line.strip().split()[1] for line in lines[start+2+num_ens+en_ind::num_ens+1]])
    return ens[en_ind], states, xss

def branch(start_node, start_node_prob, tree, non_leaf_nodes):
        curr_nodes = [[start_node, start_node, 0, 0, start_node_prob]]
        output = []
        while True:
            next_nodes = []
            for _, node, _, _, prob in curr_nodes:
                if node not in non_leaf_nodes: continue
                subtree = tree[node]
                to_append = subtree.copy()
                to_append[:, 4] *= prob
                to_append = to_append.tolist()
                output += to_append
                next_nodes += to_append
            curr_nodes = next_nodes
            if len(curr_nodes) == 0:
                break
        return np.array(output)

def gen_spec(start_state, pop, tree, non_leaf_nodes, sigma, x):
    spec = branch(start_state,pop,tree,non_leaf_nodes)
    y = np.zeros(len(x))
    for line in spec: # upper_i, lower_i, energy, A, b_ratio
        y += gaussian(x, line[4], line[2], sigma)
    return y

file = 'o8_co2'
fac_en, cx_states, xss = cx_xss(f'{fac_dir}{file}/{file}.cx',35800)
tree, non_leaf_nodes = read_tr(f'{fac_dir}{file}/{file}.tr')
x = np.linspace(500,900,5000)
y = np.zeros(len(x))
for cx_state, xs in zip(cx_states[:,0], xss):
    y = y + gen_spec(int(cx_state), float(xs), tree, non_leaf_nodes,2.5,x)

fig, ax = plt.subplots(figsize=(8.09*2,5*2))
plt.plot(x,y/np.max(y))
plt.xlabel('Energy [eV]')
plt.ylabel('Arb. Intensity')
plt.title('O$^{8+}$ on CO$_2$'+f'   {fac_en:.2e} eV')
plt.tight_layout()
plt.minorticks_on()
plt.show()