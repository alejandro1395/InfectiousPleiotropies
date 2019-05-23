#!/usr/bin/python3

import pandas as pd
import numpy as np
import requests, sys
from itertools import combinations
#import seaborn as sns
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
CREATING PAIRS OF PLEIOTROPIES IN DATABASE
"""

#VARIABLES

pd.options.mode.chained_assignment = None  # default='warn'
out_file = "../results/PairsDB.csv"
catalog = pd.read_csv("../../GWAS_Age_merged.csv", sep='\t', low_memory=False)#panda creation
LD = pd.read_csv('../../../CEU_mod.csv', sep='\t')

#First of all let's make a filter to remove rows where the same SNP appears twice

def filter_same_SNP_twice_in_dataset(table):
    table = table[table['SNP1'] != table['SNP2']]
    return table

#Then let's merge the filtered catalog with onsets with the LD info
def merge_GWAS_info_with_linkage_info(gwas_table, linkage):
    PairsDF = pd.DataFrame(columns=['SNPA', 'DiseaseA', 'RiskAllA', 'GroupA',
    'OnsetA', 'POS1', 'SNPB', 'DiseaseB', 'RiskAllB', 'GroupB', 'OnsetB',
    'POS2', 'R2'])
    for i, row in linkage.iterrows():
        df_subset1 = gwas_table.loc[gwas_table['SNPS'] == row['SNP1']].reset_index()
        numDiseases_SNP1 = len(df_subset1)
        df_subset2 = gwas_table.loc[gwas_table['SNPS'] == row['SNP2']].reset_index()
        numDiseases_SNP2 = len(df_subset2)
        if numDiseases_SNP1 >= 1 and numDiseases_SNP2 >= 1:
            for num1 in range(0, numDiseases_SNP1):
                for num2 in range(0, numDiseases_SNP2):
                    PairsDF = PairsDF.append({'CHR': df_subset1.loc[num1, 'CHR_ID'],
                    'SNPA': row['SNP1'],
                    'SNPA': row['SNP1'],
                    'DiseaseA': df_subset1.loc[num1, 'Disease'],
                    'RiskAllA': df_subset1.loc[num1, 'Risk_allele'],
                    'GroupA': df_subset1.loc[num1, 'Group'],
                    'OnsetA': df_subset1.loc[num1, 'Onset'],
                    'POS1': row['POS1'],
                    'SNPB': row['SNP2'],
                    'DiseaseB': df_subset2.loc[num2, 'Disease'],
                    'RiskAllB': df_subset2.loc[num2, 'Risk_allele'],
                    'GroupB': df_subset2.loc[num2, 'Group'],
                    'OnsetB': df_subset2.loc[num2, 'Onset'],
                    'POS2': row['POS2'],
                    'R2': row['R2']},
                    ignore_index=True)
        #print(PairsDF)
    return PairsDF


        #for i in range:
            #df.loc[i] = ['<some value for first>','<some value for second>','<some value for third>']`






LD = filter_same_SNP_twice_in_dataset(LD)
#LD.to_csv(out_file, sep='\t')
PairsDB = merge_GWAS_info_with_linkage_info(catalog, LD)
PairsDB.to_csv(out_file, sep='\t')
