#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys
import numpy as np
import seaborn as sns
import itertools

""" Usage: python look_at_further_split.py breakdown_of_further_split.txt num_regions_not_further_split.txt further_split_int_train.bed further_split_int_test.bed further_split_int_ref.bed """


num_not_further_split = int(open(sys.argv[2]).readlines()[0])

further_split_file_overview = sys.argv[1]

num_further_split = 0
num_regions_from_further_split = 0
num_regions_per_further = []
together_to_num = {}
for line in open(further_split_file_overview):
	num_further_split += 1
	fields = line.strip('\r\n').split()
	num_regions_from_further_split += int(fields[0])
	num_regions_per_further.append(int(fields[0]))
	together_to_num[fields[1]] = {'total': int(fields[0]), 'num_train': 0, 'num_test': 0, 'num_ref': 0}

def autolabel(rects):
	""" Attach a text label above each bar in *rects*, displaying its height."""
	for rect in rects:
		height = rect.get_height()
		ax.annotate('{}'.format(height), xy=(rect.get_x() + rect.get_width()/2, height), xytext=(0,3), textcoords="offset points", ha='center', va='bottom')


'''bar plot of num regions per further split'''

unique_vals, unique_counts = np.unique(np.sort(num_regions_per_further), return_counts=True)

fig, ax = plt.subplots()
rects = ax.bar(unique_vals, unique_counts)
autolabel(rects)
ax.set_xlabel('Number of offspring regions from each further split parent region', fontsize=12)
ax.set_ylabel('Number of occurrences', fontsize=12)
ax.set_xticks(np.unique(np.sort(num_regions_per_further)))
ax.set_xticklabels(np.unique(np.sort(num_regions_per_further)))
plt.tight_layout()
fig.savefig('bar_new_regions_per_further_split.png')
plt.close(fig)

'''pie chart showing num_not_further_split with num_further_split'''
labels = ['Regions >7 Mbp', 'Regions 4-7 Mbp']
sizes = [num_further_split, num_not_further_split]
explode = [0.1, 0]
colors = ['steelblue', 'limegreen']

fig, ax = plt.subplots()
ax.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
ax.axis('equal')
fig.savefig('pie_num_further_split.png')
plt.close(fig)

'''Are they in the same set?'''

train_intersect_file = sys.argv[3]
test_intersect_file = sys.argv[4]
ref_intersect_file = sys.argv[5]

for file, subkey in zip([train_intersect_file, test_intersect_file, ref_intersect_file], ['num_train', 'num_test', 'num_ref']):
	for line in open(file):
		fields = line.strip('\r\n').split('\t')
		together_to_num[fields[4]][subkey] += 1

where_they_are = np.zeros((len(together_to_num.keys()), 3))
ordered_num_total = np.full((len(together_to_num.keys())), '0')
for i, key in enumerate(together_to_num):
	where_they_are[i, 0] = together_to_num[key]['num_train']/together_to_num[key]['total']
	where_they_are[i, 1] = together_to_num[key]['num_test']/together_to_num[key]['total']
	where_they_are[i, 2] = together_to_num[key]['num_ref']/together_to_num[key]['total']
	ordered_num_total[i] = str(int(together_to_num[key]['total']))

fig, ax = plt.subplots(figsize=(6,16))
ax.set_title('Proportion of offspring regions\nin each set for each parent region', fontsize=18)
im = ax.pcolor(where_they_are, cmap='Greys', vmin=0, vmax=1, edgecolors='black', linewidths=1)
cbar = fig.colorbar(im, ax=ax)
cbar.ax.tick_params(labelsize=18)
tick_marks = np.arange(int(where_they_are.shape[0]))
ax.set_yticks(tick_marks+0.5)
ax.set_yticklabels(tick_marks+1, fontsize=18)
ax.set_xticks(np.arange(3)+0.5)
ax.set_xticklabels(['Train\nset', 'Test\nset', 'Ref\nset'], fontsize=18)
ax.set_ylabel('Parent region', fontsize=18)

fmt = '.2f'
for i,j in itertools.product(range(int(where_they_are.shape[0])), range(int(where_they_are.shape[1]))):
	color='white' if where_they_are[i,j] > 0.5 else 'black'
	if where_they_are[i,j] != 0 and where_they_are[i,j] != 1:
		plt.text(j+0.5, i+0.5, format(where_they_are[i,j],fmt),fontsize=18, horizontalalignment='center',verticalalignment='center', color=color)

for i in range(int(where_they_are.shape[0])):
	if np.amax(where_they_are[i]) != 1:
		plt.text(3.1, i+.5, ordered_num_total[i], fontsize=12, color='C0', horizontalalignment='center', verticalalignment='center')

#plt.text(3.1, i+2.5, 'Number\nof offspring\nregions', fontsize=12, color='C0', horizontalalignment='center', verticalalignment='center')
plt.text(3.1, -1.1, 'Number\nof offspring\nregions',fontsize=12, color='C0', horizontalalignment='center', verticalalignment='center')
fig.savefig('further_split_where_they_are.png')
plt.close(fig)

