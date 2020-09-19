# coding=utf-8

#
#

import data_loader as dl
import matplotlib.pyplot as plt
import pandas as pd
import re
import seaborn as sns; sns.set(style="ticks", color_codes=True)
import sys
import numpy as np
import clusterMatrix as cm
import matplotlib.patches as patches
import data_loader as dl

# The function create a correlation matrix to show the relationship between samples
# d : The data object generated by data_loader.py
# type : The matrix
#



def drawCorrMatrixPic(d, type="count", label="des", groupby="tRNA_Families", method="pearson", fig_size= 12, font_size= 12):
    '''
    The function create a correlation matrix to show the relationship between samples
    @param d:  The data object generated by data_loader.py
    @param type: "count"/"pileup" the matrix used for calculated correlation matrix
    @param groupby: group by can be ,"CODE_AA","AA","CODE","tRNA_Family" to sum the value
                    then calculate correlation between them
    @param method: Method to calculate correlation (pearson, kendall, spearman)
    @param font_size: The size of text in matrix
    @return: None
    '''

    df = d["df"]
    s_df = d["s_df"]

    exp_df = d["exp_df"]
    if type=="count":
        exp_df = d["count_exp_df"]
    elif type == "pileup":
        exp_df = d["pileup_exp_df"]

    sample_df = d["s_df"]
    sample_ls = d["sample_ls"]

    #Calulate the expression level for AA
    # col_num = 1
    # fig, axs = plt.subplots(8, col_num, figsize=[25, len(sample_ls)*6/2], sharex=True, squeeze=True,)

    # group = exp_df.groupby(['AA_CODE']).sum()
    # index = 0
    # for s in sample_df["SAMPLE_ID"]:
    #     sample_label = sample_df[sample_df["SAMPLE_ID"]==s]["SAMPLE_LABEL"].values[0]
    #     cur_axs = None
    #     if col_num==1:
    #         cur_axs = axs[int(index / col_num)]
    #     else:
    #         cur_axs = axs[int(index / col_num), index % col_num]
    #     cur_axs.bar(group.index.values, group[s], alpha=0.4)
    #     cur_axs.text(0.02, 0.85, sample_label, transform=cur_axs.transAxes,
    #                                         va='center', fontsize=20, weight='bold')
    #     cur_axs.set_xticklabels(group.index.values, rotation=50, horizontalalignment="right", fontsize=16)
    #     if index % col_num==0 and int(index / col_num)== int(len(sample_ls)/6)+1:
    #         cur_axs.text(-0.3, 1.3, "Normalized Read Number", transform=cur_axs.transAxes,
    #                           va='center', fontsize=20, weight='bold', rotation=90)
    #     index+=1
    # plt.show()

    # Bar charts for all decoder
    # group = exp_df.groupby(['RNA_family']).sum()
    # index = 0
    # for s in sample_df["SAMPLE_ID"]:
    #     sample_label = sample_df[sample_df["SAMPLE_ID"]==s]["SAMPLE_LABEL"].values[0]
    #     cur_axs = None
    #     if col_num==1:
    #         cur_axs = axs[int(index / col_num)]
    #     else:
    #         cur_axs = axs[int(index / col_num), index % col_num]
    #     cur_axs.bar(group.index.values, group[s], alpha=0.4)
    #     cur_axs.text(0.02, 0.85, sample_label, transform=cur_axs.transAxes,
    #                                         va='center', fontsize=20, weight='bold')
    #     cur_axs.set_xticklabels(group.index.values, rotation=50, horizontalalignment="right", fontsize=16)
    #     if index % col_num==0 and int(index / col_num)== int(len(sample_ls)/6)+1:
    #         cur_axs.text(-0.3, 1.3, "Normalized Read Number", transform=cur_axs.transAxes,
    #                           va='center', fontsize=20, weight='bold', rotation=90)
    #     index+=1
    # plt.show()

    code_ls = []
    aa_ls = []
    aa_code_ls = []
    #df = np.log2(exp_df.groupby(['RNA_family']).sum()+1)
    # print(exp_df.columns)
    df = exp_df.groupby([groupby]).sum()
    column_dic = {}
    if label=="des":
        df.rename(columns=d["sample_dic"], inplace=True)
        sample_des_ls = []
        for s in sample_ls:
            sample_des_ls.append(d["sample_dic"][s])
        sample_ls = sample_des_ls
    # Delete column if sum is zero
    for colname in df.columns:
        csum = df[colname].sum()
        if csum == 0:
            del df[colname]
            print("delete the column :" + colname)
    df = df.sort_index(axis=1)
    dl.csv_download_link(df, 'samples_corr.tsv', delete_prompt=False)
    figure = d["report_dir"]+"/cor_matrix.png"
    cm.drawCorrelationClusterMatrix(df, figure, figsize=fig_size, method=method,annot_font_size=font_size) # Correlation matrix

    # for i in df.index:
    #     g = re.search('tRNA-(\w+)-(\w+)', i)
    #     if g:
    #         code_ls.append(g.group(2))
    #         aa_ls.append(g.group(1))
    #         aa_code_ls.append(g.group(1)+"_"+g.group(2))
    #     else:
    #         code_ls.append('None')
    #         aa_ls.append('None')
    #         aa_code_ls.append('None')
    # df['AA'] = aa_ls
    # df["CODE"] = code_ls
    # df["AA_CODE"] = aa_code_ls
    # g = sns.pairplot(df, hue="AA",plot_kws = {'alpha': 0.6, 's': 80, 'edgecolor': 'k'}, height = 2)
    # g.fig.show()
    #g.savefig("test.png")


    # f = plt.figure(figsize=(19, 15))
    # plt.matshow(df.corr(), fignum=f.number)
    # plt.xticks(range(df.shape[1]), df.columns, fontsize=14, rotation=45)
    # plt.yticks(range(df.shape[1]), df.columns, fontsize=14)
    # cb = plt.colorbar()
    # cb.ax.tick_params(labelsize=14)
    # plt.title('Correlation Matrix', fontsize=16)
    # plt.show()

# wdir = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/output"
# sample_tsv = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/sample.tsv"
# trna_anno_bed =  "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/ChIP/hg38_tRNA_60.bed"
# d = dl.LoadWDir(wdir, sample_tsv, trna_anno_bed)
# drawPic(d)