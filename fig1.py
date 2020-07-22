#!/usr/bin/env python
import pandas as pd, matplotlib.pyplot as plt, seaborn as sns

font = {'weight': 'normal', 'size': 16}
plt.rc('text', usetex=True)
plt.rc('font', **font)


fig = plt.figure(figsize=(7.5,5))
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

palette='hls'
ax.set_prop_cycle(color=sns.color_palette(palette,10)[::-1])
ax2.set_prop_cycle(color=sns.color_palette(palette,10)[::-1])

ax.axvline(6,linewidth=1,color='k',linestyle='--')
ax2.axvline(10,linewidth=1,color='k',linestyle='--')

df = pd.read_csv('clonal_family_sizes.csv')
for d in df.columns[:-1]:
    df[d]
    ax.plot(df[d],label=d)

df2 = pd.read_csv('branch_lengths.csv')
for d in df2.columns[:-1]:
    df[d]
    ax2.plot(df[d],label=d)

ax.set_xlim([1,10**3])
ax.set_ylim([10**(-5),1])
ax.set_xticks([1,10,100,1000])
ax.set_xscale('log')
ax.set_yscale('log')

ax2.set_xlim([1,30])
ax2.set_ylim([0,0.3])
ax2.set_xticks([1,10,20,30])

ax.set_xlabel('Clonal family size')
ax.set_ylabel('Distribution')

ax2.set_xlabel('Branch length')
ax2.set_ylabel('Distribution')

handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels, loc='upper right', frameon=False, title='Donors')

plt.tight_layout()
plt.savefig('fig1.pdf')