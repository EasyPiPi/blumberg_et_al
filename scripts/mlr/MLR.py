#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: yizhao

multiple linear regression for half-lives and gene features
"""

import os
from os.path import expanduser

import numpy as np
import pandas as pd
from pandas import DataFrame, Series

import matplotlib.pyplot as plt
from matplotlib import cm

import seaborn as sns
import pylab as pl

import statsmodels.formula.api as smf
from sklearn import preprocessing

import itertools


root_dir = expanduser('~/Desktop/github_repo/blumberg_et_al')
output = os.path.join(root_dir, 'output/mlr')
os.makedirs(output, exist_ok=True)

# for gene features, gc percentage is decimal, and length is bp.
gene_feature = pd.read_csv(os.path.join(root_dir, "data/features/human_gene_features.csv"), index_col=0)

# keep the longest transcript
gene_feature = gene_feature.sort_values("len_gene", ascending=False).drop_duplicates(subset="gene_id")
gene_feature.index = gene_feature["gene_id"]

# read published half-life
half_life = pd.read_csv(os.path.join(root_dir, "data/published_half_life/published_half_lives.csv"), index_col = 0)
half_life = half_life.loc[:, ['Schofield_et_al_K562', 'Mele_et_al_K562', 'Wachutka_et_al_K562']]
half_life.columns = ['Schofield_et_al', 'Mele_et_al', 'Wachutka_et_al']

# read blumgberg_et_al
# amit_hl = pd.read_csv(os.path.join(root_dir, "output/half_life/table/half_life_K562_Amit_total_RNA.csv"), index_col="ensembl_gene_id")

# Use ENCODE half-life
amit_hl = pd.read_csv(os.path.join(root_dir, "output/half_life/table/half_life_K562_2014NG_polyA_RNA.csv"), index_col="ensembl_gene_id")
amit_hl = amit_hl[["half_life", "biotype"]]
amit_hl.columns = ["Blumberg_et_al", "biotype"]

half_life = half_life.join(amit_hl)
# log2 transfer for half life 
half_life.loc[:, half_life.columns.str.contains("et_al")] = np.log2(half_life.loc[:, half_life.columns.str.contains("et_al")])

hl_MLR = half_life.join(gene_feature)
hl_MLR.replace([np.inf, -np.inf], np.nan, inplace = True)

# Only protein coding genes
hl_MLR = hl_MLR.loc[hl_MLR.biotype == "protein_coding", :]

# Limit to the same set of genes
hl_MLR = hl_MLR.dropna()

hl_MLR.to_csv(os.path.join(output, "MLR_hl_and_features.csv"))

########################################################
# generate table for R heatmap plot
# get coef and p values from the MLR model

#measurement_name = 'Schofield_et_al'
#hl_df = hl_MLR
#subset_col = sel_col

def MLR_coef_p(measurement_name, hl_df, subset_col):
    
    df = hl_df[[measurement_name] + subset_col]
    df.dropna(inplace = True)
    
    # normalize features
    df.loc[:, sel_col] = preprocessing.scale(df[sel_col])
     
    print("The shape of the df is " + str(df.shape))
    # print(df.head())
    
    factors_string = " + ".join(subset_col)
    whole_formula = measurement_name + " ~ " + factors_string
    
    #print(whole_formula)
    
    MLR_fit = smf.ols(formula= str(whole_formula), data=df).fit()
    
    print(MLR_fit.summary())
    
    params = MLR_fit.params
    params.name = "coefficient"
    
    pvalues = MLR_fit.pvalues
    pvalues.name = "pvalue"
    
    print('R2: ', MLR_fit.rsquared)
    # https://stackoverflow.com/questions/44302099/python-statsmodels-ols-confidence-interval
    conf_interval = MLR_fit.conf_int()
    
    coef_p = pd.concat( [params, pvalues, conf_interval], axis=1 )
    coef_p['Study'] = measurement_name
    
    # drop intercepts
    coef_p.drop(["Intercept"], axis=0, inplace=True)
    
    return coef_p

# loop all the half life measurements and get all the coef and p values
# loop cdna and intron features for spliced genes
sel_col = ['gc_3utr', 'gc_5utr', 'gc_cds', 'gc_intron', 'len_3utr', 'len_5utr', 'len_cds',  'len_intron', 'exonJunDen']

spliced_MLR = DataFrame()

for name in half_life.loc[:, half_life.columns.str.contains("et_al")].columns:
    print(name)
    df = MLR_coef_p(name, hl_df=hl_MLR, subset_col=sel_col)
    spliced_MLR = pd.concat([spliced_MLR, df])

spliced_MLR.index.name = 'covariate'

spliced_MLR.to_csv(os.path.join(output, "MLR_summary_spliced_protein_coding_genes.csv"))

# include SEM results for blumberg et al.
sem_coef = pd.read_csv(os.path.join(root_dir, "output/sem/out/protein_coding_spliced_full_feature_for_same_set_gene_in_MLR_out.csv"))
sem_coef = sem_coef.loc[sem_coef.response == "Half_life" , ["covariate", "est", "ci.lower", "ci.upper"]]
sem_coef.index = sem_coef.covariate

sem_coef.columns = ['variables', 'coefficient', 0, 1]
sem_coef["Study"] = "Blumberg_et_al SEM"

spliced_MLR = pd.concat([spliced_MLR, sem_coef])

# figures
var_dict = {'gc_5utr' : 'G+C 5\'UTR', 'gc_cds': 'G+C cds', 'gc_3utr' : 'G+C 3\'UTR', 'gc_intron' : 'G+C intron', 'len_5utr' : 'len 5\'UTR', 'len_cds' : 'len cds', 'len_3utr' : 'len 3\'UTR', 'len_intron' : 'len intron', 'exonJunDen' : 'spl. junc. den.'}

spliced_MLR["variables"] = spliced_MLR.index
spliced_MLR = spliced_MLR.replace(var_dict)

# https://stackoverflow.com/questions/43159528/error-bars-with-seaborn-and-stripplot
spliced_MLR.sort_values(by = ["Study", "variables"], inplace=True)
spliced_MLR["yerr"] = spliced_MLR['coefficient'] - spliced_MLR[0]

g = sns.pointplot(x="variables", y="coefficient", hue="Study", data=spliced_MLR,                   dodge=0.5, join=False, ci=None, scale = 0.7, palette="Set1")

# g = sns.swarmplot(x="variables", y="coefficient", hue="Study", data=spliced_MLR)

x_coords = []
y_coords = []
for point_pair in g.collections:
    for x, y in point_pair.get_offsets():
        x_coords.append(x)
        y_coords.append(y)

g.errorbar(x_coords, y_coords, yerr=spliced_MLR.yerr, ecolor = [item for item in sns.color_palette(palette="Set1", n_colors=5) for i in range(9)], fmt=' ', zorder=-1)

g.set_xticklabels(labels = pd.unique(spliced_MLR["variables"]), rotation=45, ha="right")
plt.plot(range(0,10), np.repeat(0,10), color='grey', linestyle='dashed', linewidth=2)
plt.tight_layout()
plt.savefig(os.path.join(output, "MLR_spliced_protein_coding_genes.png"))

