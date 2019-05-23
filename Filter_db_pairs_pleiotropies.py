#!/usr/bin/python3

import pandas as pd
import numpy as np
import requests, sys
from itertools import combinations
from scipy import stats
import matplotlib.pyplot as plt
import pickle
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages
import copy
from scipy.stats import sem, t
from scipy import mean
import re
import os

"""
FILTER PAIRS OF PLEIOTROPIES IN DATABASE
"""

#VARIABLES

pd.options.mode.chained_assignment = None  # default='warn'
PairsDB = pd.read_csv("../results/PairsDB.csv", sep='\t', low_memory=False)#panda creation
catalog = pd.read_csv("../../GWAS_Age_merged.csv", sep='\t', low_memory=False)#panda creation
LD = pd.read_csv('../../../CEU_mod.csv', sep='\t')


#FUNCTIONS

def filter_repeated_pairs(table):
    list_of_pairs = {}
    for i, row in table.iterrows():
        pair_of_SNPs = str(sorted([row['SNPA'], row['SNPB']]))
        pair_of_Diseases = sorted([row['DiseaseA'], row['DiseaseB']])
        if row['DiseaseA'] == row['DiseaseB']:
            table = table.drop(i)
            continue
        if pair_of_SNPs in list_of_pairs.keys():
            if pair_of_Diseases in list_of_pairs[pair_of_SNPs]:
                table = table.drop(i)
            else:
                list_of_pairs[pair_of_SNPs].append(pair_of_Diseases)
        else:
            list_of_pairs[pair_of_SNPs] = []
            list_of_pairs[pair_of_SNPs].append(pair_of_Diseases)
    return table

newPairsDB = filter_repeated_pairs(PairsDB)
newPairsDB.to_csv("../results/filteredPairsDB.csv", sep='\t')
