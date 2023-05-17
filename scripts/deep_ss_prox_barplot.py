#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 16:35:31 2023

@author: emre
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Create a new DataFrame with the results
results_df = pd.read_csv('mutational_groups_extended.csv',delimiter='\t')

label_counts = pd.Series(results_df['group']).value_counts().sort_index()


ss_count = sum(results_df['group'].str.count('ss_'))
prox_count = sum(results_df['group'].str.count('prox'))
deep_count = sum(results_df['group'].str.count('deep'))
# Count the occurrences of each label
all_count = len(results_df)
# Calculate the total count for normalization

# Split the counts into separate Series for each group and calculate percentages
ss_counts = (label_counts[label_counts.index.str.startswith('ss')] / all_count) * 100
prox_counts = (label_counts[label_counts.index.str.startswith('prox')] / all_count) * 100
deep_counts = (label_counts[label_counts.index.str.startswith('deep')] / all_count) * 100

# Set the fixed bar width
bar_width = 0.0001

# Generate x-coordinates for the bars
x_ss = np.arange(len(ss_counts)) * (bar_width * 2)
x_prox = np.arange(len(prox_counts)) * (bar_width * 1.25)
x_deep = np.arange(len(deep_counts)) * (bar_width * 2)

# Create the bar plots side by side
fig, ax1 = plt.subplots(1, figsize=(3, 5))

ax1.bar(x_ss, ss_counts, color='mediumpurple', width=0.0001)
ax1.set_title('SS')
ax1.set_ylabel('')
#ax1.set_xticks(x_ss)
#ax1.yticks(np.arange(0, 1, step=0.1))  # Set label locations.
ax1.set_ylim([0, 0.1])
plt.rcParams.update({'font.size': 20})
ax1.set_xticklabels('')
plt.show()

# Create the bar plots side by side
fig, ax2 = plt.subplots(1, figsize=(5, 8))
ax2.bar(x_prox, prox_counts, color='mediumaquamarine', width=0.0001)
ax2.set_title('PROX')
#ax2.set_xticks(x_prox)
ax2.set_ylim([0, 0.3])
ax2.set_xticklabels('')
plt.rcParams.update({'font.size': 15})
plt.show()


# Create the bar plots side by side
fig, ax3 = plt.subplots(1, figsize=(16, 20))
ax3.bar(x_deep, deep_counts, color='steelblue', width=bar_width)
ax3.set_title('DEEP')
#ax3.set_xticks(x_deep)
plt.rcParams.update({'font.size': 30})
ax3.set_xticklabels('')

plt.show()



