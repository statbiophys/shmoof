#!/usr/bin/env python
import numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns

font = {'weight': 'normal', 'size': 16}
plt.rc('text', usetex=True)
plt.rc('font', **font)

def read_model(file_context,file_position):
    xs = np.loadtxt(file_position,delimiter='\t',usecols=[0],skiprows=1,dtype=int)
    betas = pd.read_csv(file_position,delimiter='\t',usecols=[1],dtype=float).values.T[0]

    ws = np.loadtxt(file_context,delimiter='\t',usecols=[0],skiprows=1,dtype=str)
    gammas = pd.read_csv(file_context,delimiter='\t',usecols=[1],dtype=float).values.T[0]
    omegas = pd.read_csv(file_context,delimiter='\t',usecols=[2,3,4,5],dtype=float).values.T[0]
    
    return betas[np.argsort(xs)],gammas[np.argsort(ws)],omegas[np.argsort(ws)]

def read_s5f(link_mutability,link_subtitution):
    mutability = pd.read_csv(link_mutability,delimiter=' ',usecols=[0,1],dtype=str)
    ws = mutability['Fivemer'].str.replace("'","")
    gammas = mutability['Mutability'].str.replace("'","").astype(float)
    
    substitution = pd.read_csv(link_subtitution,delimiter=' ',usecols=[0,1,2,3,4],dtype=str)
    ws2 = mutability['Fivemer'].str.replace("'","")
    omegas_a = substitution['A'].str.replace("'","").astype(float)
    omegas_c = substitution['C'].str.replace("'","").astype(float)
    omegas_g = substitution['G'].str.replace("'","").astype(float)
    omegas_t = substitution['T'].str.replace("'","").astype(float)
    omegas = np.array([omegas_a,omegas_c,omegas_g,omegas_t]).T
    
    return np.ones(500),np.array(gammas)[np.argsort(ws)],omegas[np.argsort(ws2)]

def read_std(filename):
    std = pd.read_csv(filename,delimiter='\t',usecols=[1],dtype=float).values.T[0]
    return std

_,s5f,_ = read_s5f('http://clip.med.yale.edu/shm/distribution/Mutability.csv','http://clip.med.yale.edu/shm/distribution/Substitution.csv')

xs = np.arange(0,500,1)
x0 = 60
dx = 100/3.14
sin = np.exp(np.sin((xs-x0)/dx))
cos =  np.exp(np.cos((xs-x0)/dx))
sin[:60] = 0
cos[:60] = 0

betas,gammas,_ = read_model('models/model_synthetic_s5f_one_context.tsv','models/model_synthetic_s5f_one_position.tsv')
betas_std = read_std('models/std_synthetic_s5f_one_position.tsv')
betas_sin,gammas_sin,_ = read_model('models/model_synthetic_s5f_sin_context.tsv','models/model_synthetic_s5f_sin_position.tsv')
betas_sin_std = read_std('models/std_synthetic_s5f_sin_position.tsv')
betas_cos,_,_ = read_model('models/model_synthetic_s5f_cos_context.tsv','models/model_synthetic_s5f_cos_position.tsv')
betas_cos_std = read_std('models/std_synthetic_s5f_cos_position.tsv')

filename = 'histograms/cdr_positions_v5.tsv'
cdr = pd.read_csv(filename,delimiter='\t',usecols=[1],dtype=float).values.T[0]
cdr[:200] = 0

fig = plt.figure(figsize=(15,5))

ax1 = fig.add_subplot(131)

ax1.plot([-0.01,100],[-0.01,100],'k-',linewidth=1)

lns1 = ax1.plot(s5f,gammas/np.mean(gammas),'ro',alpha=0.33,label=r'S5F with $\beta=1$')

ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlim([0.1,10])
ax1.set_ylim([0.1,10])

ax1.set_xlabel('True context mutabilities')
ax1.set_ylabel('Inferred mutabilities $\gamma$')

left, bottom, width, height = [0.1, 0.6, 0.1, 0.3]
ax2 = fig.add_axes([left, bottom, width, height])

ax2.plot([-0.01,100],[-0.01,100],'k-',linewidth=1)

lns2 = ax2.plot(s5f,gammas_sin/np.mean(gammas_sin),'go',alpha=0.33,label=r'S5F with $\beta = \beta^1_x$')

ax2.set_yscale('log')
ax2.set_xscale('log')

ax2.set_xlim([0.1,10])
ax2.set_ylim([0.1,10])

lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc='lower right',frameon=False)


ax3 = fig.add_subplot(132)

ax3.set_xlim([60,400])
ax3.set_ylim([0.1,10])

ax3.set_yscale('log')

ax3.plot(np.arange(0,1000,1),np.ones(1000), 'k-',linewidth=1,label=r'True profile')
ax3.plot(betas, 'r-',label=r'Recovered $\beta = 1$')
ax3.fill_between(range(500),betas+1.96*betas_std, betas-1.96*betas_std,color='red',alpha=0.3,linewidth=0.)

ax3.set_xlabel('Position $x$')
ax3.set_ylabel(r'Inferred mutabilities $\beta$')

handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles[::-1], labels[::-1], loc='upper right',frameon=False)

ax4 = ax3.twinx()

ax4.fill_between(range(500),cdr, np.zeros(500),color='k',alpha=0.1,label='CDR3',linewidth=0.)

ax4.set_ylim([0,1])
ax4.set_yticks([])

ax5 = fig.add_subplot(133)

ax5.plot(range(500), sin/np.mean(sin), 'k-',linewidth=1,label=r'True profiles')
ax5.plot(range(500), cos/np.mean(cos), 'k-',linewidth=1)

ax5.plot(betas_sin/np.mean(betas_sin[70:270]), 'b-',label=r'Recovered $\beta = \beta^2_x$')
ax5.fill_between(range(500),(betas_sin+1.96*betas_sin_std)/np.mean(betas_sin[70:270]), (betas_sin-1.96*betas_sin_std)/np.mean(betas_sin[70:270]),color='b',alpha=0.3,linewidth=0.)

ax5.plot(betas_cos, 'g-',label=r'Recovered $\beta = \beta^1_x$')
ax5.fill_between(range(500),betas_cos+1.96*betas_cos_std, betas_cos-1.96*betas_cos_std,color='g',alpha=0.3,linewidth=0.)

ax5.set_xlabel('Position $x$')

ax5.set_xlim([60,400])
ax5.set_ylim([0.1,10])

ax5.set_yscale('log')

handles, labels = ax5.get_legend_handles_labels()
ax5.legend(handles[::-1], labels[::-1], loc='upper right',frameon=False)

ax6 = ax5.twinx()

ax6.fill_between(range(500),cdr, np.zeros(500),color='k',alpha=0.1,label='CDR3',linewidth=0.)

ax6.set_ylim([0,1])
ax6.set_yticks([])

plt.tight_layout()
fig.savefig('/home/spisak/old_shmoof/fig2.pdf', bbox_inches='tight')
