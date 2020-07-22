#!/usr/bin/env python
import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib.patches as patches, matplotlib.gridspec as gridspec
import seaborn as sns
from itertools import product
from collections import Counter, defaultdict

font = {'weight': 'normal', 'size': 14}
plt.rc('text', usetex=True)
plt.rc('font', **font)
color = {"A": "#FF1100", "C":"#1122FF", "G": "#00BB00", "T": "#FFFF00"}

def word2number(word): 
    switcher = {'A': 0,'C': 1,'G': 2,'T': 3,}
    number = 0 
    for i in range(5):
        number += switcher.get(word[4-i],0)*4**i
    return number

def read_model(file_context,file_position):
    xs = np.loadtxt(file_position,delimiter='\t',usecols=[0],skiprows=1,dtype=int)
    betas = pd.read_csv(file_position,delimiter='\t',usecols=[1],dtype=float).values.T[0]

    ws = np.loadtxt(file_context,delimiter='\t',usecols=[0],skiprows=1,dtype=str)
    gammas = pd.read_csv(file_context,delimiter='\t',usecols=[1],dtype=float).values.T[0]
    omegas = pd.read_csv(file_context,delimiter='\t',usecols=[2,3,4,5],dtype=float).values.T[0]
    
    return betas[np.argsort(xs)],gammas[np.argsort(ws)],omegas[np.argsort(ws)]

betas,gammas,_ = read_model('models/model_synthetic_s5f_one_context.tsv','models/model_synthetic_s5f_one_position.tsv')


length = 5
fontsize=15
height_rectangles = 15
xs = range(256)
width=1/256


fig = plt.figure(figsize=(15,10))

gs = gridspec.GridSpec(4, 1)
gs.update(wspace=0.025, hspace=0.05)

ax = fig.add_subplot(gs[0])

order_height = [1,2,0,3,-1]  

words = []
i2='A'
for i0,i1,i3,i4 in product(["A", "C", "G", "T"],repeat=4):
    words.append("".join([i4,i1,i2,i0,i3]))

a_mu = [gammas[word2number(w)]-1 for w in words]
a_hot = [int(w[1]=='A' or w[1]=='T') for w in words]

ax.bar([x/256+0.5/256 for x in xs if a_mu[x] >= 0], 
        [h for h in a_mu if h >= 0],
        width=width,
        color="black")
ax.bar([x/256+0.5/256 for x in xs if a_hot[x] > 0], 
        [h for i,h in enumerate(a_mu) if a_hot[i] > 0],
        width=width,
        color="orange")
ax.bar([x/256+0.5/256 for x in xs if a_mu[x] < 0], 
        [5*h for h in a_mu if h < 0],
        width=width,
        color="grey")
ax.bar([x/256+0.5/256 for x in xs if a_hot[x]>0 and a_mu[x] < 0], 
        [5*h for i,h in enumerate(a_mu) if a_hot[i] > 0 and h < 0],
        width=width,
        color="orange",label='WA hotspots', linewidth=0)

for ii, yy in enumerate(order_height):
    for k, x in enumerate(["A", "C", "G", "T"]*(4**ii)):
        rect = patches.Rectangle((k/4**(ii), 
                                 yy/(length+1)*height_rectangles + min(a_mu)
                                 - height_rectangles),
                                 1/4**(ii), height_rectangles/(length+1), 
                                 fill=True,alpha=0.75, 
                                 facecolor=color[x])

        ax.add_artist(rect)
        if(ii < 3):
            rx, ry = rect.get_xy()
            cx = rx + rect.get_width()/2.0
            cy = ry + rect.get_height()/2.0

            ax.annotate(x, (cx, cy), color='k', weight='bold', 
                        fontsize=fontsize, ha='center', va='center',rotation=90)

ax.set_xlim(0,1)
ax.set_ylim(min(a_mu)-7*height_rectangles/6, max(a_mu))

ax.set_ylabel('Mutabilities $\gamma$')

ax.set_xticks([])
ax.set_yticks([-5,0,1,2,3,4])
ax.set_yticklabels(['0','1','','','','5'],rotation = 0)


ax2 = fig.add_subplot(gs[1])

words = []
i2='T'
for i0,i1,i3,i4 in product(["T", "G", "C", "A"],repeat=4):
    words.append("".join([i3,i0,i2,i1,i4]))

t_mu = [gammas[word2number(w)]-1 for w in words]


ax2.bar([x/256+0.5/256 for x in xs if t_mu[x] >= 0], 
        [h for h in t_mu if h >= 0],
        width=width,
        color="black")
ax2.bar([x/256+0.5/256 for x in xs if t_mu[x] < 0], 
        [5*h for h in t_mu if h < 0],
        width=width,
        color="grey")

for ii, yy in enumerate(order_height):
    for k, x in enumerate(["T", "G", "C", "A"]*(4**ii)):
        rect = patches.Rectangle((k/4**(ii), 
                                 yy/(length+1)*height_rectangles + min(a_mu)
                                 - height_rectangles),
                                 1/4**(ii), height_rectangles/(length+1), 
                                 fill=True,alpha=0.75, 
                                 facecolor=color[x])

        ax2.add_artist(rect)
        if(ii < 3):
            rx, ry = rect.get_xy()
            cx = rx + rect.get_width()/2.0
            cy = ry + rect.get_height()/2.0

            ax2.annotate(x, (cx, cy), color='k', weight='bold', 
                        fontsize=fontsize, ha='center', va='center',rotation=270)

ax2.set_xlim(0,1)
ax2.set_ylim(min(a_mu)-7*height_rectangles/6, max(a_mu))

ax2.set_ylabel('Mutabilities $\gamma$')

ax2.set_xticks([])
ax2.set_yticks([-5,0,1,2,3,4])
ax2.set_yticklabels(['0','1','','','','5'],rotation = 0)


ax3 = fig.add_subplot(gs[2])

words = []
i2='C'
for i0,i1,i3,i4 in product(["C", "G", "T", "A"],repeat=4):
    words.append("".join([i4,i1,i2,i0,i3]))

c_mu = [gammas[word2number(w)]-1 for w in words]

c_hot = [int((w[0]=='A' or w[0]=='T') and (w[1]=='A' or w[1]=='G') and (w[3]=='C' or w[3]=='T')) for w in words]
c_cold = [int((w[0]=='C' or w[0]=='G') and (w[1]=='C' or w[1]=='T')) for w in words]

ax3.bar([x/256+0.5/256 for x in xs if c_mu[x] >= 0], 
        [h for h in c_mu if h >= 0],
        width=1./256,
        color="black")
ax3.bar([x/256+0.5/256 for x in xs if c_hot[x] > 0], 
        [h for i,h in enumerate(c_mu) if c_hot[i] > 0],
        width=1./256,
        color="purple")
ax3.bar([x/256+0.5/256 for x in xs if c_cold[x] > 0], 
        [h for i,h in enumerate(c_mu) if c_cold[i] > 0],
        width=1./256,
        color="cyan")
ax3.bar([x/256+0.5/256 for x in xs if c_mu[x] < 0], 
        [5*h for h in c_mu if h < 0],
        width=1./256,
        color="grey")
ax3.bar([x/256+0.5/256 for x in xs if c_hot[x]>0 and c_mu[x] < 0], 
        [5*h for i,h in enumerate(c_mu) if c_hot[i] > 0 and h < 0],
        width=1./256,
        color="purple",label='WRCY hotspots')
ax3.bar([x/256+0.5/256 for x in xs if c_cold[x]>0 and c_mu[x] < 0], 
        [5*h for i,h in enumerate(c_mu) if c_cold[i] > 0 and h < 0],
        width=1./256,
        color="cyan",label='SYC coldspots')


for ii, yy in enumerate(order_height):
    for k, x in enumerate(["C", "G", "T", "A"]*(4**ii)):
        rect = patches.Rectangle((k/4**(ii), 
                                 yy/(length+1)*height_rectangles + min(a_mu)
                                 - height_rectangles),
                                 1/4**(ii), height_rectangles/(length+1), 
                                 fill=True,alpha=0.75, 
                                 facecolor=color[x])

        ax3.add_artist(rect)
        if(ii < 3):
            rx, ry = rect.get_xy()
            cx = rx + rect.get_width()/2.0
            cy = ry + rect.get_height()/2.0

            ax3.annotate(x, (cx, cy), color='k', weight='bold', 
                        fontsize=fontsize, ha='center', va='center',rotation=90)

ax3.set_xlim(0,1)
ax3.set_ylim(min(a_mu)-7*height_rectangles/6, max(a_mu))

ax3.set_ylabel('Mutabilities $\gamma$')

ax3.set_xticks([])
ax3.set_yticks([-5,0,1,2,3,4])
ax3.set_yticklabels(['0','1','','','','5'],rotation = 0)


ax4 = fig.add_subplot(gs[3])

words = []
i2='G'
for i0,i1,i3,i4 in product(["G", "C", "A", "T"],repeat=4):
    words.append("".join([i3,i0,i2,i1,i4]))

g_mu = [gammas[word2number(w)-1]-1 for w in words]

g_hot = [int((w[1]=='A' or w[1]=='G') and (w[3]=='C' or w[3]=='T') and (w[4]=='A' or w[4]=='T')) for w in words]
g_cold = [int((w[3]=='A' or w[3]=='G') and (w[4]=='C' or w[4]=='G')) for w in words]

ax4.bar([x/256+0.5/256 for x in xs if g_mu[x] >= 0], 
        [h for h in g_mu if h >= 0],
        width=1./256,
        color="black")
ax4.bar([x/256+0.5/256 for x in xs if g_hot[x] > 0], 
        [h for i,h in enumerate(g_mu) if g_hot[i] > 0],
        width=1./256,
        color="purple")
ax4.bar([x/256+0.5/256 for x in xs if g_cold[x] > 0], 
        [h for i,h in enumerate(g_mu) if g_cold[i] > 0],
        width=1./256,
        color="cyan")
ax4.bar([x/256+0.5/256 for x in xs if g_mu[x] < 0], 
        [5*h for h in g_mu if h < 0],
        width=1./256,
        color="grey")
ax4.bar([x/256+0.5/256 for x in xs if g_hot[x]>0 and g_mu[x] < 0], 
        [5*h for i,h in enumerate(g_mu) if g_hot[i] > 0 and h < 0],
        width=1./256,
        color="purple",label='WRCY hotspots')
ax4.bar([x/256+0.5/256 for x in xs if g_cold[x]>0 and g_mu[x] < 0], 
        [5*h for i,h in enumerate(g_mu) if g_cold[i] > 0 and h < 0],
        width=1./256,
        color="cyan",label='SYC coldspots')

for ii, yy in enumerate(order_height):
    for k, x in enumerate(["G", "C", "A", "T"]*(4**ii)):
        rect = patches.Rectangle((k/4**(ii), 
                                 yy/(length+1)*height_rectangles + min(a_mu)
                                 - height_rectangles),
                                 1/4**(ii), height_rectangles/(length+1), 
                                 fill=True,alpha=0.75, 
                                 facecolor=color[x])

        ax4.add_artist(rect)
        if(ii < 3):
            rx, ry = rect.get_xy()
            cx = rx + rect.get_width()/2.0
            cy = ry + rect.get_height()/2.0

            ax4.annotate(x, (cx, cy), color='k', weight='bold', 
                        fontsize=fontsize, ha='center', va='center',rotation=270)

ax4.set_xlim(0,1)
ax4.set_ylim(min(a_mu)-7*height_rectangles/6, max(a_mu))

ax4.set_ylabel('Mutabilities $\gamma$')

ax4.set_xticks([])
ax4.set_yticks([-5,0,1,2,3,4])
ax4.set_yticklabels(['0','1','','','','5'],rotation = 0)

#plt.tight_layout()
fig.savefig('/home/spisak/old_shmoof/fig3.pdf')#,bbox_inches='tight')