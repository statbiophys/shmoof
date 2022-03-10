#!/usr/bin/env python
import pandas as pd
import numpy as np
from ete3 import PhyloTree
from tqdm import tqdm
from multiprocessing import Pool,cpu_count
from scipy.optimize import brentq

from utils import *

class ContextPosition():
    def __init__(self,individuals,data_location='oof_data/',threshold=20,seq_length=500):
        
        self.data_location = data_location
        self.seq_length = seq_length
        self.threshold = threshold
        self.exclude_x = []
        self.exclude_w = [0]
        
        align,trees = self.readTrees(individuals)
        self.traverseTrees(align,trees)

    def sequence2glossary(self,sequence,x_position=2): 
        '''
        x_position is the position of the mutating bp in 5-mer context, 
        equals 0,1,2,3,4: from x____ to ____x
        default is the centered position: __x__
        '''
        seq_length = len(sequence) 
        glossary = np.zeros(self.seq_length,dtype=int)
        for i in range(seq_length-4):
            word=sequence[i:(i+5)]        
            glossary[i+x_position] = word2number(word)
        return glossary

    def whereMutations(self,parent,child):
        '''
        Return a pointer to mutated positions
        '''
        indicate = np.zeros(self.seq_length,dtype=bool)
        for x,(nt,old_nt) in enumerate(zip(child,parent)):
            if nt != old_nt and nt != 'N' and old_nt != 'N':
                indicate[x] = 1
        return indicate
            
    def readTrees(self,individuals):
        """ Read trees with corresponding alignment """
        align = pd.DataFrame()
        trees = pd.DataFrame()
        for i in individuals:
            df = pd.read_csv(self.data_location+'{}_alignment.tab'.format(i),sep='\t')
            df['INDIVIDUAL'] = i
            df['FAMILY'] =  df['INDIVIDUAL'].astype(str) + df['FAMILY'].astype(str)
            align = pd.concat([align,df])
            df = pd.read_csv(self.data_location+'{}_trees.tab'.format(i),sep='\t')
            df['INDIVIDUAL'] = i
            df['FAMILY'] =  df['INDIVIDUAL'].astype(str) + df['FAMILY'].astype(str)
            trees = pd.concat([trees,df])
        return align, trees

    def traverseTrees(self,align,trees):
        """ Traverse the tree and extract information from branches """
        glossaries_list = []
        indicators_list  = []
        times_list = []
        for family,local_align in tqdm(align[['FAMILY','NODE','ALIGNMENT']].groupby(['FAMILY'])):
            structure = trees[trees['FAMILY']==family].values[0][1]
            alignment = df2fasta(local_align[['NODE','ALIGNMENT']])
            tree = PhyloTree(structure, alignment, format=1, alg_format='fasta')
            for node in tree.traverse('preorder'):
                if not node.is_leaf() and not node.name=='':
                    glossary = self.sequence2glossary(node.sequence)
                    for child in node.children:
                        # ignore the trunk: germline-to-MRCA branch
                        if not child.name == '0': 
                            indicator = self.whereMutations(child.sequence,node.sequence)
                            nb_of_mutations = sum(indicator)
                            # ignore long branches (optional, specify with self.threshold)     
                            if nb_of_mutations>0 and nb_of_mutations<self.threshold: 
                                glossaries_list.append(glossary)
                                indicators_list.append(indicator)
                                times_list.append(nb_of_mutations/(self.seq_length-sum(glossary == 0)))
        self.indicators = np.array(indicators_list)
        self.glossaries = np.array(glossaries_list)
        self.times = np.array(times_list)
        self.exclude_x.extend(np.where(self.indicators.sum(axis=0)==0))
                              
    def solve_for_x(self,gamma,x):
        '''
        For a given position x update beta given gamma
        '''
        rates_mut = []
        sum_rates_const = 0.0
        for i,delta in enumerate(self.indicators[:,x]):
            w = self.glossaries[i,x]
            if delta==1:
                rates_mut.append(self.times[i]*gamma[w])
            else:
                sum_rates_const += self.times[i]*gamma[w]
        def onevarfunc(beta_x):
            return function_for_x(rates_mut, sum_rates_const, beta_x)

        if onevarfunc(0.01)*onevarfunc(100.0)<0:
            return (x,brentq(onevarfunc, 0.01, 100.0, rtol=0.01))
        else:
            return (x,1.0)
                              
    def update_beta(self,gamma):
        '''
        Update beta given gamma 
        Excluding positions in self.exclude_x
        '''
        beta = np.zeros(self.seq_length,dtype=float)  
        pool = Pool(cpu_count())
        results = [pool.apply_async(self.solve_for_x, args=(gamma,x)) for x in range(self.seq_length)]
        beta_measured = []  
        for res in tqdm(results):
            beta_measured.append(res.get(timeout=120))  
        for x,beta_x in beta_measured:
            beta[x] = beta_x     
        pool.close()
        return beta
                              
    def solve_for_w(self,beta,w):
        '''
        For a given word w update gamma given beta
        '''
        rates_mut = []
        sum_rates_const = 0
        occurences = np.transpose(np.where(self.glossaries==w))
        for i,x in occurences:
            if self.indicators[i,x]==1:
                rates_mut.append(self.times[i]*beta[x])
            else:
                sum_rates_const += self.times[i]*beta[x]
        def onevarfunc(gamma_w):
            return function_for_w(rates_mut,sum_rates_const,gamma_w)

        if onevarfunc(0.01)*onevarfunc(100.0)<0:
            return (w,brentq(onevarfunc, 0.01, 100.0, rtol=0.01))
        else:
            return (w,1.0)
                              
    def update_gamma(self,beta=None):
        '''
        Update gamma given beta 
        Excluding motifs in self.exclude_w
        '''
        if beta is None:
            beta = np.ones(self.seq_length,dtype=float)   
        gamma = np.zeros(1025,dtype=float)   
        pool = Pool(cpu_count())
        results = [pool.apply_async(self.solve_for_w, args=(beta,w)) for w in range(1025) if not w in self.exclude_w] 
        gamma_measured = []
        for res in tqdm(results):
            gamma_measured.append(res.get(timeout=120))      
        for w,gamma_w in gamma_measured:
            gamma[w] = gamma_w
        pool.close()
        return gamma
                              
    def infer(self,nb_of_iterations=10,beta0=None):
        '''
        Maximizing likelihood in self.nb_of_iterations steps 
        '''
        if beta0 is None: beta0=np.ones(self.seq_length,dtype=float)
        gammas = np.zeros((1025,nb_of_iterations),dtype=float)
        betas = np.zeros((self.seq_length,nb_of_iterations),dtype=float)
        for i in range(nb_of_iterations):
            gamma1 = self.update_gamma(beta0)
            beta1 = self.update_beta(gamma1)
            gammas[:,i] = gamma1
            betas[:,i] = beta1    
            beta0 = beta1
        return gammas,betas