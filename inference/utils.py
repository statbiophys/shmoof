#!/usr/bin/env python
import pandas as pd
import numpy as np
import re

switcher = {'A': 0,'C': 1,'G': 2,'T': 3,}
reswitcher = {0: 'A',1 : 'C',2: 'G',3: 'T',}

def df2fasta(df):
    fasta = ""
    for node,align in df.values:
        fasta = fasta + ">{}".format(node) + "\n"
        fasta = fasta + align + "\n"
    return fasta

def only_nts(string):
    match = re.match("^[ACGT]*$", string)
    return match is not None

def word2number(word): 
    '''
    map AAAAA - TTTTT to 1 - 1024
    any 5-mer with N mapped to 0
    ''' 
    if not only_nts(word):
        return 0
    else:
        number = 0 
        for i in range(5):
            number += switcher.get(word[4-i],)*4**i
        return number+1

def number2word(number):
    '''
    word2number(number2word) = word
    '''
    if number == 0:
        return "NNNNN"
    else:
        number -= 1
        rests = []
        word = []
        for i in range(5):
            rests.append(number%4)
            number = number//4
        for r in rests:
            word.append(reswitcher.get(r,'N'))
        return "".join(word[::-1])

def function_for_w(rates,sum_rates_const,gamma):
    '''
    Gamma-dependent function with root gamma*
    Utilizes approximations for small and large gamma*beta
    '''
    summa = 0.0
    for rate in rates:
        if gamma*rate < 0.1:
            summa += 1.0/gamma
        elif gamma*rate > 10.0:
            summa += rate*np.exp(-rate*gamma)
        else:
            summa += rate*(1.0/(np.exp(rate*gamma)-1.0))   
    return summa - sum_rates_const

def function_for_x(rates,sum_rates_const,beta):
    '''
    Beta-dependent function with root beta*
    Utilizes approximations for small and large gamma*beta
    '''
    summa = 0.0
    for rate in rates:
        if beta*rate < 0.1:
            summa += 1.0/beta
        elif beta*rate > 10.0:
            summa += rate*np.exp(-rate*beta)
        else:
            summa += rate*(1.0/(np.exp(rate*beta)-1.0))   
    return summa - sum_rates_const

def howManyMutations(parent,child):
    '''
    return the number of mutations
    '''
    number = 0
    for x,(nt,old_nt) in enumerate(zip(child,parent)):
        if nt != old_nt and nt != 'N' and old_nt != 'N':
            number += 1
    return number