#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib.pyplot as plt

intersected_file = sys.argv[1]
non_intersected_file = sys.argv[2]

state_to_color = {-1: (0/255, 0/255, 0/255),
                  0: ( 255/255, 255/255, 255/255),
                  1: ( 64/255, 142/255, 39/255),
                  2: ( 144/255, 144/255, 144/255),
                  3: ( 45/255, 57/255, 222/255),
                  4: ( 252/255, 250/255, 108/255),
                  5: ( 247/255, 229/255, 77/255),
                  6: ( 248/255, 223/255, 181/255),
                  7: ( 196/255, 79/255, 243/255),
                  8: ( 65/255, 145/255, 40/255),
                  9: ( 192/255, 87/255, 152/255),
                  10: ( 225/255, 48/255, 34/255),
                  11: ( 235/255, 78/255, 51/255),
                  12: ( 219/255, 132/255, 53/255),
                  13: ( 147/255, 43/255, 236/255),
                  14: ( 127/255, 185/255, 56/255),
                  15: ( 180/255, 37/255, 47/255),
                  16: ( 86/255, 90/255, 204/255),
                  17: ( 89/255, 146/255, 77/255),
                  18: ( 232/255, 50/255, 35/255),
                  19: ( 235/255, 129/255, 49/255),
                  20: ( 0/255, 26/255, 210/255),
                  21: ( 226/255, 52/255, 34/255),
                  22: ( 0/255, 21/255, 191/255),
                  23: ( 225/255, 147/255, 52/255),
                  24: ( 198/255, 42/255, 80/255),
                  25: ( 42/255, 99/255, 24/255),
                  26: ( 191/255, 64/255, 178/255)}

state_to_label = {0: 'quiescent',
                  1: 'transcribed',
                  2: 'heterochromatin',
                  3: 'polycomb',
                  4: 'enhancer-like',
                  5: 'active enhancer-like',
                  6: 'active',
                  7: 'CTCF',
                  8: 'transcribed enhancer-like',
                  9: 'nuclease accessible enhancer-like',
                  10: 'active & nuclease accessible promoter-like',
                  11: 'promoter- & enhancer-like',
                  12: 'active & nuclease accessible enhancer-like',
                  13: 'nuclease accessible CTCF',
                  14: 'active transcribed enhancer-like',
                  15: 'nuclease accessible promoter-like',
                  16: 'polycomb heterochromatin',
                  17: 'transcribed heterochromatin',
                  18: 'active promoter-like',
                  19: 'active promoter- & enhancer-like',
                  20: 'bivalent enhancer-like',
                  21: 'active nuclease accessible promoter- & enhancer-like',
                  22: 'bivalent promoter- & enhancer-like',
                  23: 'active transcribed nuclease accessible promoter- & enhancer-like',
                  24: 'active nuclease accessible CTCF & promoter- & enhancer-like',
                  25: 'transcribed CTCF',
                  26: 'nuclease accessible transcirbed CTCF & enhancer-like'}

no_state_color = ( 0/255, 0/255, 0/255)

cts = []
PL_to_state = {}
for line in open(intersected_file):
    fields = line.strip('\r\n').split('\t')
    PL = (fields[0], fields[1], fields[2])
    ct = fields[3].replace('/project/vision/Data/IDEAS/ideasVisionV20p8Seg', '').replace('.bed', '')
    if ct not in cts:
        cts.append(ct)
    state_color = int(fields[7])
    if PL not in PL_to_state:
        PL_to_state[PL] = {}
    if ct not in PL_to_state[PL]:
        PL_to_state[PL][ct] = state_color
for line in open(non_intersected_file):
    fields = line.strip('\r\n').split('\t')
    PL = (fields[0], fields[1], fields[2])
    PL_to_state[PL] = {}
    for ct in cts:
        PL_to_state[PL][ct] = -1

PL_states_arr = np.zeros((len(PL_to_state.keys()), len(cts)))
PL_states_c_arr = np.zeros((len(PL_to_state.keys()), len(cts)), dtype=np.object)
for i, PL in enumerate(PL_to_state.keys()):
    for j, ct in enumerate(cts):
        if ct in PL_to_state[PL]:
            PL_states_arr[i,j] = PL_to_state[PL][ct]
            PL_states_c_arr[i,j] = state_to_color[PL_to_state[PL][ct]]
        else:
            PL_states_arr[i,j] = -1
            PL_states_c_arr[i,j] = state_to_color[-1]

fig, ax = plt.subplots(figsize=(14,14))
for i in range(PL_states_c_arr.shape[0]):
    ax.scatter(np.arange(PL_states_c_arr.shape[1]), np.full(PL_states_c_arr.shape[1], i), color=PL_states_c_arr[i], s=0.45)
ax.set_xticks(np.arange(PL_states_c_arr.shape[1]))
ax.set_xticklabels(cts, fontsize=12)
ax.set_yticks([])
ax.set_yticklabels([])
ax.set_xlabel('cell type', fontsize=18)
ax.set_ylabel('Partition Location', fontsize=12)
fig.savefig('pl_state_dot_plot.png')
plt.close(fig)

fig, ax = plt.subplots()

x = (range(27))
new_x = [3*i for i in x]

for i in x:
    ax.bar(new_x[i], np.sum(PL_states_arr==i)/PL_states_arr.flatten().shape[0], color=state_to_color[i], edgecolor='black', width=1)
ax.set_xticks(new_x)
ax.set_xticklabels(np.arange(27), fontsize=8)
ax.set_xlabel('Epigenetic state', fontsize=18)
ax.set_ylim(0,1)
ax.set_yticks(np.arange(0,1.2,0.2))
#ax.set_yticklabels(np.arange(0,1.2,0.2), fontsize=14)
ax.set_ylabel('Fraction of partition locations\nacross 20 cell types', fontsize=18)
plt.tight_layout()
fig.savefig('pl_state_bar.png')
plt.close(fig)
