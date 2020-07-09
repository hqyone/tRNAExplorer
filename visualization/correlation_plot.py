# coding=utf-8


#
#  tRNAExplorer v.1.0 (2020.) - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#

import visualization.data_loader as dl
import matplotlib.pyplot as plt
import pandas as pd
import re
import seaborn as sns; sns.set(style="ticks", color_codes=True)
import sys
import numpy as np
import visualization.clusterMatrix as cm
import matplotlib.patches as patches

def drawPic(d):
    # df = pd.read_csv(dl["combined_mean_profile"], sep="\t")
    #
    # s_df = pd.read_csv(dl["sample_source_db"], sep="\t")
    # s_df = s_df[s_df['sample_type']=='total']
    # s_df = s_df.sort_values(by=['cell_type','cell_line_name'])
    #
    # st_df = pd.read_csv(dl["trna_structure_tab"], sep='\t')
    # exp_df = pd.read_csv(dl["mean_trna_sample_matrix"], sep="\t")

    df = d["df"]
    st_df = d["s_df"]
    exp_df= d["exp_df"]
    sample_df = d["s_df"]
    sample_ls = sample_df["#ID"].unique()

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
    print(exp_df.columns)
    df = exp_df.groupby(['tRNA_Families']).sum()
    df = df.sort_index(axis=1)
    cm.drawCorrelationClusterMatrix(df[sample_ls], "matrix.png") # Correlation matrix

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

wdir = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/output"
sample_tsv = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/sample.tsv"
trna_anno_bed =  "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/ChIP/hg38_tRNA_60.bed"
d = dl.LoadWDir(wdir, sample_tsv, trna_anno_bed)
drawPic(d)