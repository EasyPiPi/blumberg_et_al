#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 10:49:30 2019

@author: yixin
"""
import numpy as np
import pandas as pd
from pandas import DataFrame, Series
import os

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import style

import seaborn as sns
import pylab

import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.distributions.empirical_distribution import ECDF

from math import isnan

root_dir = "Desktop/github_repo/blumberg_et_al"
mir_fig_dir = os.path.join(root_dir, "output/miRNA/figure")
mir_tab_dir = os.path.join(root_dir, "output/miRNA/table")
genecat_fig_dir = os.path.join(root_dir, "output/genecat/figure")
os.makedirs(mir_fig_dir, exist_ok=True)
os.makedirs(mir_tab_dir, exist_ok=True)
os.makedirs(genecat_fig_dir, exist_ok=True)

# read in half-life table
all_gene_hl = pd.read_csv(os.path.join(root_dir, "output/half_life/table/half_life_concat_libraries.csv"), index_col = 0)

genetype = pd.read_csv(os.path.join(root_dir, "data/id_mapper/genetype.csv"), index_col = 0)
protein_coding = genetype[(genetype == "protein_coding").values].index

# all_gene_hl = all_gene_hl[['K562_Amit_PR']].dropna()
all_gene_hl = all_gene_hl[['K562_2014NG_PR']].dropna()
all_gene_hl = all_gene_hl[all_gene_hl.index.isin(protein_coding)]

#################################################################
# process miRNA targets
all_targets = pd.read_csv(os.path.join(root_dir, "data/targetscan/Human_release_7_2/Predicted_Targets_Info.default_predictions.txt"), sep = '\t')
all_targets = all_targets[(all_targets['Species ID'] == 9606)]
all_targets.loc[:, 'Gene ID'] = all_targets['Gene ID'].str.split(".").str.get(0).values

miRNA_target_split = all_targets[['miR Family', 'Gene ID']]
miRNA_target_split_gb = miRNA_target_split.groupby("miR Family")
miRNA_target_pieces = dict(list(miRNA_target_split_gb))

#for miRNA, df in miRNA_target_pieces.items():
#
#    print(miRNA)
#    print(df.head)

# compare RMA half life for miRNA targets and non-targets

def target_plot_df(miRNA_name):

    # drop zero values
    df_targets_boxplot = all_gene_hl[all_gene_hl.index.isin(miRNA_target_pieces[miRNA_name]['Gene ID'])]
    df_non_targets_boxplot = all_gene_hl[~all_gene_hl.index.isin(miRNA_target_pieces[miRNA_name]['Gene ID'])]

    print("The number of zero values in targets' half lives is " + str((df_targets_boxplot == 0).sum()))
    df_targets_boxplot = df_targets_boxplot.replace(0, np.nan)

    print("The number of zero values in targets' half lives is " + str((df_non_targets_boxplot == 0).sum()))
    df_non_targets_boxplot = df_non_targets_boxplot.replace(0, np.nan)


    df_targets_boxplot.columns = [miRNA_name + ' targets ']
    df_non_targets_boxplot.columns = [miRNA_name + ' non-targets ']

    print('--------------------------------')

    if not df_targets_boxplot.iloc[:,0].dropna().empty:
        group_1_p = stats.ks_2samp(np.log2(df_targets_boxplot.iloc[:,0].dropna()), np.log2(df_non_targets_boxplot.iloc[:,0].dropna())).pvalue
        target_number = (~df_targets_boxplot.isna()).sum().values
        non_target_number = (~df_non_targets_boxplot.isna()).sum().values

        print('KS test for group1: p = %.3e' % group_1_p)
        print('Target Number: ', target_number)
        print('Non-Target Number:', non_target_number)

        df_targets_plot_long = pd.concat([pd.melt(df_targets_boxplot), pd.melt(df_non_targets_boxplot)], ignore_index=True)

        # Drop all half life which is equal to zero

        # df_targets_plot_long.replace(0, np.nan)
        df_targets_plot_long.value = np.log2(df_targets_plot_long.value)
        df_targets_plot_long.dropna(inplace = True)

        return df_targets_plot_long, group_1_p, target_number, non_target_number

    else:
        return pd.DataFrame(), np.nan, np.nan, np.nan

# CDF

def CDF_plot(df, xlim=(-7,5), color = [cm.tab20(6), cm.tab20(7), cm.tab20(0), cm.tab20(1)]):
    g = sns.FacetGrid(df, hue="variable", height = 5, palette=color, legend_out=False, xlim=xlim, aspect = 1.2)
    # g = g.map(sns.distplot, "value", kde_kws=dict(cumulative=True), hist = False)
    g = g.map(sns.ecdfplot, "value")

    g.add_legend(title = '')
    plt.legend(loc=4)
    plt.xlabel(r'log$_2$$($T$_1$$_/$$_2$$^P$$^R$$)$')
    plt.ylabel('ECDF')
    plt.tight_layout()

#KS_test_target_number_summary = DataFrame(columns = ['miRNA_name', 'group_1_p', 'group_2_p', 'target_number', 'non_target_number'])
KS_test_target_number_summary = DataFrame(columns = ['miRNA_name', 'group_1_p', 'target_number', 'non_target_number'])

#miRNA_name = 'let-7-5p/98-5p'
#miRNA_target_pieces[miRNA_name]['gene_name']

miRNA_target_pieces = {k: miRNA_target_pieces[k] for k in miRNA_target_pieces if not miRNA_target_pieces[k].empty}

i = 0

sns.set(font_scale=1.3)
sns.set_style("ticks")

for miRNA_name in miRNA_target_pieces.keys():

    all_gene_hl[all_gene_hl.index.isin(miRNA_target_pieces[miRNA_name]['Gene ID'])]

    df_targets_plot_long, group_1_p, target_number, non_target_number = target_plot_df(miRNA_name)

    KS_test_target_number_summary.loc[i, 'miRNA_name'] = miRNA_name
    KS_test_target_number_summary.loc[i, 'group_1_p'] = group_1_p
    KS_test_target_number_summary.loc[i, 'target_number'] = target_number
    KS_test_target_number_summary.loc[i, 'non_target_number'] = non_target_number

    if not df_targets_plot_long.empty:
        CDF_plot(df_targets_plot_long, color=None)
        pylab.savefig(os.path.join(mir_fig_dir, ("_").join(miRNA_name.split("/")) + "_targets_half_life_CDF_20190514.pdf"))
        plt.close()

    i = i+1

KS_test_target_number_summary.columns = ['miRNA_name', 'KS_test_p_value', 'target_number', 'non_target_number']

KS_test_target_number_summary.to_csv(os.path.join(mir_tab_dir, "KS_test_summary.csv"))

#####################################################
# for ribosomal proteins and zinc-fingers
ribosomal_proteins = pd.read_csv(os.path.join(root_dir, "data/gene_lists/ribosomal_proteins/HGNC_ribosomal_proteins_20190404.txt"), sep = '\t')
zinc_fingers = pd.read_csv(os.path.join(root_dir, "data/gene_lists/zinc_fingers/HGNC_zinc_fingers_20190404.txt"), sep = '\t')

# define figure function
def gencat_plot_df(gene_name_1, gene_name_2, name_1, name_2, name_3):

    # drop zero values
    df_targets_boxplot_1 = all_gene_hl[all_gene_hl.index.isin(gene_name_1)]
    df_targets_boxplot_2 = all_gene_hl[all_gene_hl.index.isin(gene_name_2)]

    df_non_targets_boxplot = all_gene_hl[~all_gene_hl.index.isin(pd.concat([gene_name_1, gene_name_2]))]

    df_targets_boxplot_1.columns = [name_1]
    df_targets_boxplot_2.columns = [name_2]
    df_non_targets_boxplot.columns = [name_3]

    print('--------------------------------')

    group_1_p = stats.ks_2samp(np.log2(df_targets_boxplot_1.iloc[:,0].dropna()), np.log2(df_non_targets_boxplot.iloc[:,0].dropna())).pvalue
    target_number = (~df_targets_boxplot_1.isna()).sum().values
    non_target_number = (~df_non_targets_boxplot.isna()).sum().values

    print('KS test for ribosomal proteins vs. other genes: p = %.3e' % group_1_p)
    print('Ribosomal protein Number: ', target_number)
    print('Other gene Number:', non_target_number)

    print('--------------------------------')

    group_2_p = stats.ks_2samp(np.log2(df_targets_boxplot_2.iloc[:,0].dropna()), np.log2(df_non_targets_boxplot.iloc[:,0].dropna())).pvalue
    target_number = (~df_targets_boxplot_2.isna()).sum().values
    non_target_number = (~df_non_targets_boxplot.isna()).sum().values

    print('KS test for zinc-finger proteins vs. other genes: p = %.3e' % group_1_p)
    print('Zinc-finger Number: ', target_number)
    print('Other gene Number:', non_target_number)

    df_targets_plot_long = pd.concat([pd.melt(df_targets_boxplot_1), pd.melt(df_targets_boxplot_2), pd.melt(df_non_targets_boxplot)], ignore_index=True)

    # Drop all half life which is equal to zero

    df_targets_plot_long = df_targets_plot_long.replace(0, np.nan)
    df_targets_plot_long.dropna(inplace = True)
    df_targets_plot_long.value = np.log2(df_targets_plot_long.value)
    df_targets_plot_long.dropna(inplace = True)

    return df_targets_plot_long

df_targets_plot_long = gencat_plot_df(ribosomal_proteins['Ensembl gene ID'], zinc_fingers['Ensembl gene ID'], 'Ribosomal proteins', 'Zinc-fingers', 'Others')

#style.available
#sns.set(font_scale=1.3)
#sns.set_style("ticks")
#CDF_plot(df_targets_plot_long, xlim=(-7, 7), color=None)
#pylab.savefig(os.path.join(genecat_fig_dir, "ribosomal_proteins_and_zinc_fingers_proteins_half_life_CDF.pdf"))

#### Use Bootstrap to plot 95% confidence interval ####
#### for ribosomal proteins and zinc-fingers ####
# Step 1 Bootstrap
# https://stackoverflow.com/questions/49328319/random-sampling-from-a-column-several-times-in-python-pandas
# https://stackoverflow.com/questions/36390406/sample-each-group-after-pandas-groupby
# a generator function
def sample_rep(dframe, replicates=None):
    for i in range(replicates):
        # resample half-life within each gene group
        yield dframe.groupby("variable").apply(lambda x: x.sample(frac=1, replace=True)).reset_index(drop=True)

d = pd.concat(sample_rep(df_targets_plot_long, replicates=1000), keys=range(1, 1001), names=["replicate"])

# Step 2 Calculate ECDF
#def ecdf(data):
#    """ Compute ECDF """
#    x = np.sort(data)
#    n = x.size
#    y = np.arange(1, n+1) / n
#    return(x,y)
# https://www.statsmodels.org/stable/generated/statsmodels.distributions.empirical_distribution.ECDF.html?highlight=ecdf#statsmodels.distributions.empirical_distribution.ECDF

def compute_ecdf(df):
    df["ecdf"] = ECDF(df.value)(df.value)
    return df

d = d.groupby(["replicate","variable"]).apply(compute_ecdf)
d = d.reset_index(level="replicate").reset_index(drop=True)

# Step 3 compute Confidential Interval for ECD
def compute_ci(df):
    x, y = np.percentile(df.ecdf, q = [2.5, 97.5])
    z = np.mean(df.ecdf)
    return pd.Series(data = [x, y, z], index = ["lower_bound", "upper_bound", "mean"])

d = d.groupby(["variable", "value"]).apply(compute_ci)
d = d.reset_index()

# Step 4 Customized function to plot the CIs
def plot_ci(category, alpha=0.2):
    plt.fill_between(d.loc[d.variable == category, "value"], d.loc[d.variable == category, "upper_bound"], d.loc[d.variable == category, "lower_bound"], alpha=alpha)

sns.set(font_scale=1.2)
sns.set_style("ticks")

g = sns.relplot(data=d, x="value", y="mean", hue="variable", kind="line", ci=None, height=5, aspect=0.8)
g.set(xlim=(-7, 7))
g._legend.set_title("")
g._legend.set_bbox_to_anchor([1, 0.3])
plot_ci("Others")
plot_ci("Ribosomal proteins")
plot_ci("Zinc-fingers")
plt.xlabel(r'log$_2$$($T$_1$$_/$$_2$$^P$$^R$$)$')
plt.ylabel('ECDF')
plt.tight_layout()
pylab.savefig(os.path.join(genecat_fig_dir, "ribosomal_proteins_and_zinc_fingers_proteins_half_life_ECDF.pdf"))
plt.close()

#### for several selected miRNAs ####
KS_test_target_number_summary = \
    DataFrame(columns = ['miRNA_name', 'group_1_p', \
                         'target_number', 'non_target_number'])

selected_miRs = ["miR-182-5p", "miR-125-5p", "miR-19-3p"]

i = 0

for miRNA_name in selected_miRs:
    print(miRNA_name)
    df_targets_plot_long, group_1_p, target_number, non_target_number = target_plot_df(miRNA_name)
    print("--------------------------------")
    KS_test_target_number_summary.loc[i, 'miRNA_name'] = miRNA_name
    KS_test_target_number_summary.loc[i, 'group_1_p'] = group_1_p
    KS_test_target_number_summary.loc[i, 'target_number'] = target_number
    KS_test_target_number_summary.loc[i, 'non_target_number'] = non_target_number

    # Bootstrap and visualization
    # Step 1 Bootstrap
    d = pd.concat(sample_rep(df_targets_plot_long, replicates=1000), keys=range(1, 1001), names=["replicate"])
    # Step 2 Calculate ECDF
    d = d.groupby(["replicate","variable"]).apply(compute_ecdf)
    d = d.reset_index(level="replicate").reset_index(drop=True)
    # Step 3 compute Confidential Interval for ECD
    d = d.groupby(["variable", "value"]).apply(compute_ci)
    d = d.reset_index()
    # Step 4 visualization
    g = sns.relplot(data=d, x="value", y="mean", hue="variable", kind="line", ci=None, height=5, aspect=0.7)
    g.set(xlim=(-7, 7))
    g._legend.set_title("")
    g._legend.set_bbox_to_anchor([1, 0.3])
    plot_ci(d.variable.unique()[0])
    plot_ci(d.variable.unique()[1])
    plt.xlabel(r'log$_2$$($T$_1$$_/$$_2$$^P$$^R$$)$')
    plt.ylabel('ECDF')
    plt.tight_layout()
    pylab.savefig(os.path.join(mir_fig_dir, ("_").join(miRNA_name.split("/")) + "_targets_half_life_ECDF.pdf"))
    plt.close()

    i = i+1
