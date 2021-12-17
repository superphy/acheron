import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

import os, sys

abx_3_letter_code = ['AMC','AMP','AZM','FOX','TIO','CRO','CHL','CIP','GEN','NAL','STR','FIS',
    'TET','SXT','KAN', 'all']
abx_full_name = ['amoxicillin/clavulanic acid', 'ampicillin', 'azithromycin',
    'cefoxitin', 'ceftiofur', 'ceftriaxone', 'chloramphenicol', 'ciprofloxacin',
    'gentamicin', 'nalidixic acid', 'streptomycin', 'sulfisoxazole', 'tetracycline',
    'trimethoprim/sulfamethoxazole','kanamycin', 'All 15 Antimicrobials']

def label_points(x, y, labels, ax):
    for i, label in enumerate(labels):
        if label in ['CHL', 'NAL', 'AMC']:
            ax.text(x[i]+0.005, y[i]-0.009, label, fontsize = 7)
        else:
            ax.text(x[i]+0.005, y[i]+0.005, label, fontsize = 7)

def acc_vs_features(results, out): #0
    """
    Multifacet graph showing accuracy vs features for all drug, all models
    """
    df = results[(results['train']=='salm_amr') & (results['test']=='none') & (results['validate']=='none') & (results['hyp']=='False')]
    df = df.sort_values(by=['feats'])

    # rename 3 letter codes in df to their full names
    df["attribute"] = df["attribute"].map(dict(zip(abx_3_letter_code,abx_full_name)))

    # build multifacet
    sns.set(style="ticks")
    grid = sns.FacetGrid(df, col = 'attribute',col_order = abx_full_name,hue ="model",hue_order=["XGB","SVM","ANN"], col_wrap = 4, margin_titles=False, legend_out=True)
    grid = (grid.map(plt.plot, "feats", "Within 1 Dilution", alpha=1).add_legend().set_titles("{col_name}"))
    grid.fig.get_children()[-1].set_bbox_to_anchor((0.9,0.1,0,0))
    plt.setp(grid._legend.get_texts(), fontsize='18')
    plt.setp(grid._legend.get_title(), fontsize='20')
    grid.set_ylabels('Accuracy')
    grid.set_xlabels('Number of Features')
    plt.xlim(10,3000)
    plt.ylim(0.8,1)

    if out == 'stdout':
        plt.show()
    else:
        plt.savefig('figures/model_finder_multiface.png', dpi=600)
        plt.clf()

def drug_acc_15x7_bars(results, out): #1
    """
    Clustered bar graph, 15 clusters, each cluster is a dataset matchup
    """
    sns.set(style="darkgrid")
    results['train&test'] = [results['train'][i]+'-->'+results['test'][i] for i in range(len(results.index))]

    df = results[(results['model']=='XGB') & (results['feats']==1000)]

    group = sns.catplot(x='train&test', y = 'Within 1 Dilution',hue = "attribute",hue_order=abx_3_letter_code, data = df, kind = 'bar', legend_out=True)
    plt.xlabel('Dataset Comparison')
    plt.ylabel("Accuracy")
    group._legend.set_title('Antimicrobial')
    """bar_titles = [
    'Trained on GRDI,\nTested on NCBI', '\n\n\nTrained on\nNCBI Ken & Hei,\nTested on GRDI',
    'Trained on GRDI,\nTested on NCBI Ken & Hei','\n\n\nTrained on\nNCBI Ken & Hei\nwith 5-Fold\nCross Validation',
    'Trained on NCBI,\nTested on GRDI','\n\n\nTrained on NCBI\nwith 5-Fold\nCross Validation',
    'Trained on GRDI with\n5-Fold Cross Validation']"""
    bar_titles = ['1','2','3','4']
    for text, label in zip(group._legend.texts, abx_full_name):
        text.set_text(label)
    group.set_xticklabels(bar_titles)
    plt.xticks(rotation=0, fontsize = 7, horizontalalignment='center')
    group.fig.get_children()[-1].set_bbox_to_anchor((1.325,0.6,0,0))
    plt.setp(group._legend.get_title(), fontsize='18')
    #plt.setp(group._legend.get_texts(), fontsize='12')
    #plt.tight_layout(pad = 3)
    if out == 'stdout':
        plt.show()
    else:
        plt.savefig('figures/dataset_clusters.png', dpi=600)
        plt.clf()

    cust_pal = ['fire engine red','water blue', 'bright lime','vibrant purple','cyan','strong pink','dark grass green']
    sns.set_palette(sns.xkcd_palette(cust_pal))
    group = sns.catplot(x='attribute', y = 'Within 1 Dilution',hue = "train&test", order=abx_3_letter_code, data = df, kind = 'bar', legend_out =True)
    plt.xlabel('Antimicrobial')
    plt.ylabel("Accuracy")
    group._legend.set_title('Dataset Comparison')
    """bar_titles = [
    'Trained on GRDI,\nTested on NCBI', '\n\n\nTrained on\nNCBI Ken & Hei,\nTested on GRDI',
    'Trained on GRDI,\nTested on NCBI Ken & Hei','\n\n\nTrained on\nNCBI Ken & Hei\nwith 5-Fold\nCross Validation',
    'Trained on NCBI,\nTested on GRDI','\n\n\nTrained on NCBI\nwith 5-Fold\nCross Validation',
    'Trained on GRDI with\n5-Fold Cross Validation']"""
    bar_titles = ['1','2','3','4']
    for text, label in zip(group._legend.texts, bar_titles):
        text.set_text(label)
    plt.xticks(rotation=-30)
    plt.tight_layout(pad = 5.5)

    if out == 'stdout':
        plt.show()
    else:
        plt.savefig('figures/drug_clusters.png', dpi=600)
        plt.clf()

def label_points(x, y, labels, ax):
    for i, label in enumerate(labels):
        if label in ['CHL', 'NAL']:
            ax.text(x[i]+0.005, y[i]-0.009, label, fontsize = 7)
        else:
            ax.text(x[i]+0.005, y[i]+0.005, label, fontsize = 7)

def simpsons_diversity(results, out): #2
    """
    Shows simpsons diversity of PATRIC/NCBI vs GRDI data
    """
    try:
        simp_df = pd.read_pickle("data/simpsons_diversity.df")
    except:
        print("No diversity index found, run `python acheron/helpers/diversity.py` to generate one")
        raise

    ax = sns.lmplot('salm_amr', 'grdi', data = simp_df, fit_reg = False)
    plt.xlabel('Simpsons Diversity In The NBCI/PATRIC Dataset')
    plt.ylabel("Simpsons Diversity In The GRDI Dataset")
    plt.xlim(0,0.75)
    plt.ylim(0,0.75)
    label_points(simp_df['salm_amr'].values, simp_df['grdi'].values, simp_df.index.values, plt.gca())

    if out == 'stdout':
        plt.show()
    else:
        plt.savefig('figures/simpson_diversity.png', dpi=600)
        plt.clf()


def all_mic_freq_vs_acc(results, out): #3
    """
    Bar chart showing the accuracy of each class
    """
    pass

def all_mic_frequencies(results, out): #4
    """
    Bar chart showing the frequency of each class
    """
    pass

def individual_tests(results, out): #5
    """
    Each drug gets a line graph comparing all models across features
    """
    pass

def group_figures(subset, results, out, fig_num):
    """
    Saves or prints all figures related to whatever subset is declared
    """
    if subset == 'steinkey2021':
        if fig_num == 0:
            acc_vs_features(results, out)
        if fig_num == 1:
            drug_acc_15x7_bars(results, out)
        if fig_num == 2:
            simpsons_diversity(results, out)
        if fig_num == 3:
            all_mic_freq_vs_acc(results, out)
        if fig_num == 4:
            all_mic_frequencies(results, out)
        if fig_num == 5:
            individual_tests(results, out)
        #acc_vs_feats_per_dataset_matchup(results, out)

    else:
        raise Exception("Summaries need to be defined, either make your own or call one in summary.py")
