# coding=utf-8

#
#  tRNAExplorer v.1.0 (2020.) - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#

import matplotlib.pyplot as plt
import data_loader as dl
import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
import re, os
import math
import numpy as np
import clusterMatrix as cm
import matplotlib.patches as patches

def getDFs(dl):
    df = pd.read_csv(dl["mean_profile"], sep="\t")

    s_df = pd.read_csv(dl["sample_source_db"], sep="\t")
    s_df = s_df[s_df['sample_type']=='total']
    s_df = s_df.sort_values(by=['cell_type','cell_line_name'])

    st_df = pd.read_csv(dl["trna_structure_tab"], sep='\t')
    exp_df = pd.read_csv(dl["mean_trna_sample_matrix"], sep="\t")

    code_ls = []
    aa_ls = []
    code_aa_ls = []
    for i in exp_df["RNA_ID"]:
        g = re.search('tRNA-(\w+)-(\w+)', i)
        if g:
            code_ls.append(g.group(2))
            aa_ls.append(g.group(1))
            code_aa_ls.append(g.group(2)+"_"+g.group(1))
        else:
            code_ls.append('None')
            aa_ls.append('None')
            code_aa_ls.append('None')

    exp_df['AA'] = aa_ls
    exp_df["CODE"] = code_ls
    exp_df["CODE_AA"] = code_aa_ls
    sample_ls = df["#SampleID"].unique()
    sample_df = dl["sample_df"]
    return {"exp_df":exp_df, "sample_df":sample_df,"sample_ls":sample_ls}

# dfs = getDFs(dl)
# exp_df = dfs["exp_df"]
# sample_df = dfs["sample_df"]
# sample_ls = dfs["sample_ls"]
# groupby can be "AA","CODE_AA","CODE"
def drawExpForSamples(exp_df, sample_df, sample_ls, groupby, subplot_col=3, subplot_row=8, figsize_w=10, figsize_h=12,fontsize=12):
    print("Calulate the expression level for "+groupby+" in every samples")
    fig, axs = plt.subplots(subplot_row, subplot_col, figsize=[figsize_w, figsize_h], sharex=True, squeeze=True)
    group = exp_df.groupby([groupby]).sum()
    index = 0
    col_num = subplot_col
    if col_num==1:
        for s in sample_df["SAMPLE_ID"]:
            sample_label = sample_df[sample_df["SAMPLE_ID"]==s]["SAMPLE_LABEL"].values[0]
            axs[int(index / col_num)].bar(group.index.values, group[s], alpha=0.4)
            axs[int(index / col_num)].grid(True)
            axs[int(index / col_num)].text(0.02, 0.85, sample_label, transform=axs[int(index / col_num)].transAxes,
                                                va='center', fontsize=fontsize, weight='bold')
            axs[int(index / col_num)].set_xticklabels(group.index.values, rotation=50, horizontalalignment="right", fontsize=12)
            if index % col_num==0 and int(index / col_num)== int(len(sample_ls)/2)+1:
                axs[int(index / col_num)].text(-0.14, 1.3, "Normalized Read Number", transform=axs[int(index / col_num)].transAxes,
                                  va='center', fontsize=fontsize+4, weight='bold', rotation=90)
            index+=1
    else:
        for s in sample_df["SAMPLE_ID"]:
            sample_label = sample_df[sample_df["SAMPLE_ID"]==s]["SAMPLE_LABEL"].values[0]
            axs[int(index / subplot_col), index % subplot_col].bar(group.index.values, group[s], alpha=0.4)
            axs[int(index / subplot_col), index % subplot_col].grid(True)
            axs[int(index / subplot_col), index % subplot_col].text(0.02, 0.85, sample_label, transform=axs[int(index / subplot_col), index % subplot_col].transAxes,
                                                va='center', fontsize=fontsize, weight='bold')
            axs[int(index / subplot_col), index % subplot_col].set_xticklabels(group.index.values, rotation=50, horizontalalignment="right", fontsize=12)
            if index % subplot_col==0 and int(index / subplot_col)== int(len(sample_ls)/6)+1:
                axs[int(index / subplot_col), index % subplot_col].text(-0.3, 1.3, "Normalized Read Number", transform=axs[int(index / subplot_col), index % subplot_col].transAxes,
                                  va='center', fontsize=fontsize+4, weight='bold', rotation=90)
            index+=1
    plt.show()

# drawExpForAccetpor(exp_df,sample_df,sample_ls, subplot_col=2, subplot_row=12, figsize_w=14)
#
def drawExpForIsotopic(exp_df, sample_df, sample_ls, subplot_col=2, subplot_row=12, figsize_w=10, figsize_h=12,fontsize=12):
    #Calulate the expression level for anti-codon
    group = exp_df.groupby(['CODE_AA']).sum()
    fig, axs = plt.subplots(subplot_row, subplot_col, figsize=[figsize_w, figsize_h], sharex=True, squeeze=True,)
    index = 0
    col_num = subplot_col
    if col_num==1:
        for s in sample_df["SAMPLE_ID"]:
            sample_label = sample_df[sample_df["SAMPLE_ID"]==s]["SAMPLE_LABEL"].values[0]
            axs[int(index / col_num)].bar(group.index.values, group[s], alpha=0.4)
            axs[int(index / col_num)].grid(True)
            axs[int(index / col_num)].text(0.02, 0.85, sample_label, transform=axs[int(index / col_num)].transAxes,
                                                va='center', fontsize=fontsize, weight='bold')
            axs[int(index / col_num)].set_xticklabels(group.index.values, rotation=50, horizontalalignment="right", fontsize=12)
            if index % col_num==0 and int(index / col_num)== int(len(sample_ls)/2)+1:
                axs[int(index / col_num)].text(-0.14, 1.3, "Normalized Read Number", transform=axs[int(index / col_num)].transAxes,
                                  va='center', fontsize=fontsize+4, weight='bold', rotation=90)
            index+=1
    else:
        for s in sample_df["SAMPLE_ID"]:
            sample_label = sample_df[sample_df["SAMPLE_ID"]==s]["SAMPLE_LABEL"].values[0]
            axs[int(index / col_num), index % col_num].bar(group.index.values, group[s], alpha=0.4)
            axs[int(index / col_num), index % col_num].grid(True)
            axs[int(index / col_num), index % col_num].text(0.02, 0.85, sample_label, transform=axs[int(index / col_num), index % col_num].transAxes,
                                                va='center', fontsize=fontsize, weight='bold')
            axs[int(index / col_num), index % col_num].set_xticklabels(group.index.values, rotation=50, horizontalalignment="right", fontsize=12)
            if index % col_num==0 and int(index / col_num)== int(len(sample_ls)/2)+1:
                axs[int(index / col_num), index % col_num].text(-0.14, 1.3, "Normalized Read Number", transform=axs[int(index / col_num), index % col_num].transAxes,
                                  va='center', fontsize=fontsize+4, weight='bold', rotation=90)
            index+=1
    plt.show()

# Test
#drawExpForIsotopic(exp_df, sample_df, subplot_col=1, subplot_row=24,figsize_w=10, figsize_h=18)

# Draw cluster matrix for sample vs anti-codon
# g = sns.clustermap(np.log10(group[sample_df["SAMPLE_ID"]]+1), method='median')
# plt.show()
# plt.savefig('sample_anti_condon_matrix.png')
#
# # Get tRF list
# # tRF_ID SampleNumber, Type,

def normalize(df):
    result = df.copy()
    for feature_name in df.columns:
        max_value = df[feature_name].max()
        min_value = df[feature_name].min()
        result[feature_name] = (df[feature_name] - min_value) / (max_value - min_value)
    return result

def boxplot_sorted(df, by, column, rot=0, ascending=True, fontsize=10,fig_w=12, fig_h =12):
    # use dict comprehension to create new dataframe from the iterable groupby object
    # each group name becomes a column in the new dataframe
    df2 = pd.DataFrame({col:vals[column] for col, vals in df.groupby(by)})
    # find and sort the median values in this new dataframe
    meds = df2.median().sort_values(ascending=ascending)
    # use the columns in the dataframe, ordered sorted by median value
    # return axes so changes can be made outside the function
    return df2[meds.index].boxplot(rot=rot, return_type="axes", fontsize=fontsize,figsize=(fig_w,fig_h))

# groupby can be "tRF_First_Type","exp_samples","tRNA_families","CODE","CODE_AA","AA"
# value can be "max_exp","exp_samples_num","median_exp"
def drawExpBoxPlotForGroup(d,min_exp_cutof=100, sample_ls=[], groupby="tRF_First_Type", value="max_exp",fontsize=14, fig_w=12, fig_h =8):
    '''
    Create a box plot to compare expression level of tRF in multiple levels
    The program will calculate max_exp (max value of expression), exp_samples_num (The number of sample expressed the tRF)
    or median_exp (median expression across samples) and then group by some features (group by)
    to plot a box plot to compare expression between different groups
    @param d: The data object generated by data_loader.py
    @param min_exp_cutof: The minimum of read count for a sample to express the tRFs
    @param sample_ls: A list of IDs of selected sample, if it is empty, it will use all data from all available samples
    @param groupby: feature for group can be one of them ["tRF_First_Type","exp_samples","tRNA_families","CODE","CODE_AA","AA"]
    @param value: The way to calculate the expression ["]median_exp","exp_samples_num","max_exp"]
    @param fontsize: Font size (Default 14)
    @param fig_w: weight of figure (12)
    @param fig_h: height of figure (8)
    @return: None
    '''
    df = d["exp_df"]
    exp_df = df
    sel_s_ls=[]
    if len(sample_ls)>1:
        sel_s_ls, exp_df = FilterExpMartrixWithSampleIDs(df, sample_ls)
        if len(sel_s_ls)<=1:
            exp_df =df
            print("The sample number is too small, sample filtering abort")
        else:
            print("Selected Samples :"+ ",".join(sel_s_ls))

    #tRNA_Families	RNA_IDs	seq	Type
    sample_df = d["s_df"]
    all_sample_ls = sample_df["#SampleID"]
    print("Available Samples :" + ",".join(all_sample_ls))

    # Get tRF tables
    group = exp_df.groupby(['tRF_ID']).sum()
    #group1= exp_df.groupby('tRF_ID').apply(lambda x: ','.join(x.RNA_family))
    group['tRF_Types'] = exp_df.groupby('tRF_ID').apply(lambda x: ','.join(x.Type))
    group['tRF_Unique_Types'] = exp_df.groupby('tRF_ID').apply(lambda x: ",".join(list(set(x.Type))))
    group['tRF_First_Type'] = exp_df.groupby('tRF_ID').apply(lambda x: ",".join(list(set(x.Type))).split(",")[0])
    group['tRNA_families'] = exp_df.groupby('tRF_ID').apply(lambda x: ','.join(x.tRNA_Families))
    group['CODE_AA'] = exp_df.groupby('tRF_ID').apply(lambda x: ",".join(list(set(x.CODE_AA))).split(",")[0])
    group['CODE'] = exp_df.groupby('tRF_ID').apply(lambda x: ",".join(list(set(x.CODE))).split(",")[0])
    group['AA'] = exp_df.groupby('tRF_ID').apply(lambda x: ",".join(list(set(x.AA))).split(",")[0])
    # Filter with min_exp_cutof, above min_exp_cutof means the Gene expressed
    group = group.loc[group.max(axis=1)>min_exp_cutof,]
    group['exp_samples_num'] = group[sel_s_ls].apply(lambda s: (s > 100).sum(), axis=1)
    # Get the max expressing level across samples
    group['max_exp'] = np.log10(group.max(axis=1)+1)
    group['median_exp'] = np.log10(group.median(axis=1) + 1)
    group['tRF_ID'] = group.index
    group['len'], group['str_code'] = group['tRF_ID'].str.split('-', 1).str
    #print(group.groupby('tRF_First_Type').count())
    #group.boxplot(column=['max_exp'], by=['tRF_First_Type'],rot=45, fontsize=10)
    # boxplot_sorted(group,'tRF_First_Type','max_exp',rot=90, ascending=False, fontsize=16)
    ax =boxplot_sorted(group,groupby,value,rot=90, ascending=False, fontsize=fontsize, fig_w=fig_w, fig_h=fig_h)
    ax.set_ylabel(value,fontsize=fontsize+2)
    ax.set_xlabel(groupby, fontsize=fontsize + 2)
    plt.show()
    #group.to_csv('trf.xls', index=False, sep="\t")

def FilterExpMartrixWithSampleIDs(df, sample_ls):
    f_df = df[["tRF_ID","tRNA_Families","RNA_IDs","seq","Type","CODE_AA","AA","CODE"]]
    sel_s_id_ls=[]
    for s in sample_ls:
        if s in df.columns:
            f_df.loc[:,s] = df[s]
            sel_s_id_ls.append(s)
    return [sel_s_id_ls, f_df]



# wdir = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/output"
# sample_tsv = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/sample.tsv"
# trna_anno_bed =  "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/ChIP/hg38_tRNA_60.bed"
# d = dl.LoadWDir(wdir, sample_tsv, trna_anno_bed)
# dfs = getDFs(d)
# exp_df = dfs["exp_df"]
# sample_df = dfs["sample_df"]
# sample_ls = dfs["sample_ls"]
# drawExpBoxPlotForGroup(exp_df,sample_df,groupby="CODE_AA",value="median_exp")

# groupby can be "tRF_First_Type","CODE_AA","AA","CODE","tRNA_Families","CODE_AA","tRNA_ID"
def drawExpMatrixForFamily(d,figure="", sample_ls=[], groupby="tRNA_Families", fontsize=14, fig_w=12, fig_h =18):
    '''
    Create expression matrix of
    @param d: The data object generated by data_loader.py
    @param figure: Absolute path of figure, Default : report_dir+"/exp_matrix.png"
    @param sample_ls: selected sample ID list
    @param groupby: "tRF_First_Type","CODE_AA","AA","CODE","tRNA_Families","CODE_AA","tRNA_ID"
    @param fontsize: Font size (Default 14)
    @param fig_w: weight of figure (12)
    @param fig_h: height of figure (18)
    @return: None
    '''
    df = d["exp_df"]
    exp_df = df
    sel_s_ls = []
    sel_s_des = []
    if len(sample_ls) > 1:
        sel_s_ls, exp_df = FilterExpMartrixWithSampleIDs(df, sample_ls)
        if len(sel_s_ls) <= 1:
            exp_df = df
            print("The sample number is too small, sample filtering abort")
        else:
            print("Selected Samples :" + ",".join(sel_s_ls))

    # tRNA_Families	RNA_IDs	seq	Type
    sample_df = d["s_df"]
    all_sample_ls = sample_df["#SampleID"]
    print("Available Samples :" + ",".join(all_sample_ls))

    for s in sel_s_ls:
        sel_s_des.append(dl.getSampleLabel(s,sample_df))

    # Draw cluster matrix for expression data
    # group = exp_df.groupby(['CODE_AA']).sum()
    # group = exp_df.groupby(['AA']).sum()
    # group = exp_df.groupby(['CODE']).sum()
    # group = exp_df.groupby(['RNA_family']).sum()
    group = exp_df.groupby([groupby]).sum()
    group = np.log10(group+1)
    group = group[sel_s_ls] # Only left sample data
    # Replace columns name with sample Label
    group.columns = sel_s_des
    # df = normalize(group)
    df = group
    # group = exp_df.groupby(['RNA_family']).sum()
    if figure.strip() == "":
        figure = d["report_dir"]+"exp_matrix.png"
    cm.drawExpressionClusterMatrix(df,figure,fig_width=fig_w, fig_height=fig_h)
    # cm.getExpressionDistMatrix(df,'matrix.png',normalized_by=None)
    # cm.testImshow()
    # cm.drawCorrelationClusterMatrix(df, 'matrix.png')

#drawExpMatrixForFamily(exp_df,sample_df)
#print("aa")

# g = sns.clustermap(np.log10(group[sample_df["SAMPLE_ID"]]+1), yticklabels=True, figsize=(35,60))
# g = sns.clustermap(df, yticklabels=True, figsize=(5,10))
# g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 7)
# g.ax_heatmap.set_xlabel('Experiments')
# g.ax_heatmap.set_ylabel('Genes')
# hm = g.ax_heatmap.get_position()
# g.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*0.05, hm.height])
# col = g.ax_col_dendrogram.get_position()
# g.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.15, col.height*0.5])
# row = g.ax_row_dendrogram.get_position()
# g.ax_row_dendrogram.set_position([row.x0, row.y0, row.width*0.15, row.height*2])
# for a in g.ax_row_dendrogram.collections:
#     a.set_linewidth(2)
#
# for a in g.ax_col_dendrogram.collections:
#     a.set_linewidth(2)
#
# #plt.subplots_adjust(left=0.35, right=0.6, top=0.95, bottom=0.15)
# plt.show()


def plot_heatmap(df, row_linkage=None, col_linkage=None, font_size=10, legend_title="", path="", wmult=0.15, hmult=0.95):
    """
    Plots a given df using the provided row_linkage and col_linkage output from
    scipy's linkage function. The other parameters are for adjusting the plot
    formatting
    """
    #cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30)
    #cmap = sns.diverging_palette(220, 20, n=7)
    # ignore genes that have log ratio changes less than 1 for all experiments
    #g = sns.clustermap(df, row_linkage=row_linkage, col_linkage=col_linkage, cmap=cmap, figsize=(20, 20))
    g = sns.clustermap(df, row_linkage=row_linkage, col_linkage=col_linkage,  figsize=(12, 20))
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=30, horizontalalignment="right", fontsize=font_size)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize=font_size)
    plt.setp(g.cax.yaxis.get_majorticklabels(), fontsize=font_size)

    # fix ugly positioning
    hm_pos, rd_pos, cd_pos = (
        g.ax_heatmap.get_position(),
        g.ax_row_dendrogram.get_position(),
        g.ax_col_dendrogram.get_position(),
    )
    # g.ax_heatmap.set_position([hm_pos.x0 - rd_pos.width * (1 - hmult), hm_pos.y0, hm_pos.width * wmult, hm_pos.height])
    # g.ax_row_dendrogram.set_position([rd_pos.x0, rd_pos.y0, rd_pos.width * hmult, rd_pos.height])
    # g.ax_col_dendrogram.set_position(
    #     [cd_pos.x0 - rd_pos.width * (1 - hmult), cd_pos.y0, cd_pos.width * wmult, cd_pos.height / 2]
    # )
    for a in g.ax_row_dendrogram.collections:
        a.set_linewidth(4)

    for a in g.ax_col_dendrogram.collections:
        a.set_linewidth(4)
    cax_pos = g.cax.get_position()
    g.cax.set_position([cax_pos.x0 * 0.75, cax_pos.y0 * 0.25, cax_pos.width * 0.5, cax_pos.height * 0.5])
    g.cax.set_title(legend_title, fontsize=font_size)
    #plt.subplots_adjust(left=0.25, right=0.5, top=0.95, bottom=0.15)
    plt.show()
    #plt.savefig(os.path.join(path, "sig_heatmap.png"), dpi=300, bbox_inches="tight")

#plot_heatmap(df, row_linkage='average', col_linkage='average')

#plt.savefig('sample_tRNA_matrix.png',dpi=300)
# id_ls = []
# id_text_ls = g.ax_heatmap.get_ymajorticklabels()
# score_ls = []
# for i in id_text_ls:
#     trna_id = i._text
#     id_ls.append(trna_id)
#     if len(st_df[st_df['tRNA_name'] == trna_id]) > 0:
#         score_ls.append(st_df[st_df['tRNA_name'] == trna_id]['tRNA_score'].values[0])
#     else:
#         score_ls.append(0)
# plot_data = pd.DataFrame(list(zip(id_ls,score_ls)),
#                            columns=['RNA_ID','Score'])
# plt.figure(figsize=(15, 60))
# ax = sns.barplot(y="RNA_ID", x="Score", data=plot_data, color='grey')
# plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.15)
#plt.show()
#
#
#
# # Draw cluster matrix for expression data
# group = exp_df.groupby(['RNA_ID']).sum()
# group  = exp_df.groupby(['RNA_family']).sum()
# g = sns.clustermap(np.log10(group[sample_df["SAMPLE_ID"]]+1), yticklabels=True, figsize=(35,60))
# g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 5)
# #g.gs.tight_layout(g.fig)
# plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.15)
# plt.show()
# #plt.savefig('sample_tRNA_matrix.png',dpi=300)
# id_ls = []
# id_text_ls = g.ax_heatmap.get_ymajorticklabels()
# score_ls = []
# for i in id_text_ls:
#     trna_id = i._text
#     id_ls.append(trna_id)
#     if len(st_df[st_df['tRNA_name'] == trna_id]) > 0:
#         score_ls.append(st_df[st_df['tRNA_name'] == trna_id]['tRNA_score'].values[0])
#     else:
#         score_ls.append(0)
# plot_data = pd.DataFrame(list(zip(id_ls,score_ls)),
#                            columns=['RNA_ID','Score'])
# plt.figure(figsize=(15, 60))
# axs = sns.barplot(y="RNA_ID", x="Score", data=plot_data, color='grey')
# plt.subplots_adjust(left=0.05, right=0.9, top=0.95, bottom=0.15)
# #plt.show()
# plt.savefig('sample_tRNA_histogram.png',dpi=300)
#
# #Draw the tRNA score vs expression level
# fig, axs = plt.subplots(8, 3, figsize=[10, len(sample_ls)/2], sharex=True, squeeze=True)
# index = 0
# sample_ls.sort()
# for s in sample_ls:
#     sample_label = dl.getSampleLabel(s)
#     trna_ls =[]
#     value_ls = []
#     for i in group.index.values:
#         trna_id = i
#         if len(st_df[st_df['tRNA_name']==trna_id])>0:
#             trna_ls.append(trna_id)
#             score_ls.append(st_df[st_df['tRNA_name'] == trna_id]['tRNA_score'].values[0])
#             value_ls.append(group.loc[trna_id, s])
#     plot_df = pd.DataFrame(list(zip(trna_ls,score_ls,value_ls)),
#                            columns=['RNA_ID','Score','Value'])
#
#     axs[int(index/3),index%3].scatter(plot_df['Score'], plot_df['Value'], alpha=0.2)
#     #ax[int(index/3),index%3].set_xlabel('Score')
#     axs[int(index / 3), index % 3].set_ylim(0,3000)
#     axs[int(index / 3), index % 3].set_xlim(0, 110)
#     axs[int(index / 3), index % 3].text(0.01, 0.85, sample_label, transform=axs[int(index / 3), index % 3].transAxes,
#                                           va='center', fontsize=12, weight='bold')
#     if index % 3==0:
#         axs[int(index/3),index%3].set_ylabel('Read Count')
#     index+=1
# #plt.ylabel("Read Count")
# #plt.show()
# pic_dir = '/Users/hqyone/OneDrive/国内的工作/学术/科研项目/论文课题/tRNA/data_output/rna_seq/pictures'
# plt.savefig(pic_dir + "/exp_score.png", dpi=200)

