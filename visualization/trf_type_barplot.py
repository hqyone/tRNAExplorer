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
import re,numpy as np

plt.rcParams.update({'font.size': 13})

def drawTrfDistForSamples(d):
    st_df = d['st_df']
    df = pd.read_csv(d["mean_sample_trf_type_matrix"], sep="\t",index_col=False)
    s_df = pd.read_csv(d["sample_tsv"], sep='\t',index_col=False)
    st_df = pd.read_csv(d["static_log"], sep='\t', index_col=False)

    df["Description"]= df["#SampleID"]
    df["Total_reads"] = df["#SampleID"]

    sample_ls = d["sample_ls"]
    for s_id in sample_ls:
        description = s_df[s_df["#ID"]==s_id]["Description"].values[0]
        df.loc[(df['#SampleID']==s_id),"Description"]=description

        Total_reads = st_df[st_df["#SampleID"] == s_id]["survived_num"].values[0]/1000
        df.loc[(df['#SampleID'] == s_id), "Total_reads"] = Total_reads

    # #print(s_df)
    # def getSampleIndex(id,df):
    #     return df.index[df['#ID']==id][0]
    #
    # bam_sum_df = df[df.columns[1:len(df.columns)]].sum(1)
    # width=0.8
    # sample_ls = s_df["#ID"]
    # sample_label_ls=[]
    # cell_type_ls = []
    # total_read_counts= []
    # for i in sample_ls:
    #     index = getSampleIndex(i, s_df)
    #     sample_label_ls.append(s_df.at[index, 'Description'])
    #     # cell_type_ls.append(s_df.at[index, 'cell_type'])
    #     # print(static_df.loc[static_df['#BAM_ID']==i,'mapped_read_count'])
    #     total_read_num = st_df[static_df['#SampleID']==i]["survived_num"].values[0]/1000
    #     total_read_counts.append(total_read_num)=total_read_num

    # df["mapped_read_count"] = total_read_counts
    # df["Description"] = sample_label_ls

    width = 0.8
    df = df.sort_values(by=["Description"])
    ratio_df = df.copy()
    for i in range(len(df)):
        #print(df["mapped_read_count"][i])
        #df.iloc[i, 1:len(df.columns)-2]=df.iloc[i, 1:len(df.columns)-2]/df.iloc[i,len(df.columns)-1]
        #df.iloc[i, 1:len(df.columns) - 2] = df.iloc[i, 1:len(df.columns) - 2] / df["mapped_read_count"].tolist()[i]
        df.iloc[i, 1:len(df.columns) - 2] = df.iloc[i, 1:len(df.columns) - 2] / df["Total_reads"].tolist()[i]
        ratio_df.iloc[i, 1:len(df.columns)-2]=df.iloc[i, 1:len(df.columns)-2]/df.iloc[i, 1:len(df.columns)-2].sum()
    sorted_sample_ls = df["#SampleID"]
    #print(df)
    types = ['full_U_tRNA','full_tRNA','5_U_tRNA_halve','5_tRNA_halve','5_U_tRF', '5_tRF','3_U_tRNA_halve' ,'3_tRNA_halve','3_U_tRF','3_tRF' ,'i-tRF','other' ]
    types_colors = ['darkred','red','darkgreen','limegreen', 'lightgreen','greenyellow','navy' ,'blue','dodgerblue','lightblue' ,'gold','grey' ]
    fig, axs = plt.subplots(2,1, figsize=(12,8))

    #https://stackoverflow.com/questions/16006572/plotting-different-colors-in-matplotlib
    #number = len(types)
    #cmap = plt.get_cmap('hsv') #jet gist_rainbow
    #colors = [cmap(i) for i in np.linspace(0, 1, number)]

    buttom1=np.asarray([0]*len(sample_ls))
    buttom2=np.asarray([0]*len(sample_ls))
    index = 0
    for t in types:
        #axs[0].set_xticks(df["#SampleID"])
        axs[0].bar(df["#SampleID"], df[t].values, width, bottom=buttom1, color=types_colors[index])
        axs[1].bar(ratio_df["#SampleID"], ratio_df[t].values, width, bottom=buttom2, color=types_colors[index])
        index+=1
        buttom1 =buttom1+df[t]
        buttom2 =buttom2+ratio_df[t]
    axs[0].get_xaxis().set_ticks([])
    axs[0].set_ylabel("Normalized Reads (K)", fontsize=12)
    axs[1].set_xticklabels(df["Description"], rotation=40, horizontalalignment="right")
    #https://pythonspot.com/matplotlib-legend/
    axs[1].legend(types,loc='upper center', bbox_to_anchor=(1.1, 1.02), ncol=1)
    axs[1].set_ylabel("Ratio", fontsize=12)
    plt.show()
    #plt.savefig(dl.pic_out_dir +"/tRNA_components.png", dpi=100)

def drawTrfDistForAcceptersInSamples(d, isTest=False):
    df = d["df"]
    s_df = d["s_df"]

    st_df = d["st_df"]
    exp_df = d["exp_df"]
    types = d["types"]
    types_colors =d["types_colors"]
    sample_ls = d["sample_ls"]
    sample_df = d["s_df"]

    #Calulate the expression level for AA
    aa_group = exp_df.groupby(['AA']).sum()
    type_aa_group = exp_df.groupby(['Type','AA']).sum()
    index = 0

    for s_id in sample_ls:
        fig, axs = plt.subplots(2,1, figsize=(12,7), sharex=True)
        s_label = dl.getSampleLabel(s_id, sample_df)
        #aa_group = aa_group.sort_values(s_id,ascending=False)
        aa_labels = aa_group.index.tolist()

        type_ls = []
        aa_ls = []
        value_ls=[]
        for t in types:
            for aa in aa_labels:
                type_ls.append(t)
                aa_ls.append(aa)
                if aa in type_aa_group.loc[(t)].index:
                    value_ls.append(type_aa_group.loc[(t, aa)][s_id]/1000)
                else:
                    value_ls.append(0)
        plot_df = pd.DataFrame(list(zip(type_ls, aa_ls, value_ls)), columns=['Type','AA','Value'])

        aa_group = plot_df.groupby(['AA']).sum()
        aa_group = aa_group.sort_values('Value',ascending=False)
        type_group = plot_df.groupby(['Type']).sum()
        plot_type_aa_group = plot_df.groupby(['Type','AA']).sum()

        width = 0.8
        buttom1=np.asarray([0]*len(aa_group.index.tolist()))
        buttom2=np.asarray([0]*len(aa_group.index.tolist()))
        index = 0
        for t in types:
            #axs[0].set_xticks(df["#SampleID"])
            axs[0].bar(aa_group.index.tolist(),
                       plot_type_aa_group.loc[(t)].reindex(aa_group.index.tolist()).Value, width, bottom=buttom1, color=types_colors[index])
            axs[1].bar(aa_group.index.tolist(),
                       plot_type_aa_group.loc[(t)].reindex(aa_group.index.tolist()).Value/aa_group['Value'], width, bottom=buttom2, color=types_colors[index])
            index+=1
            buttom1 =buttom1+plot_type_aa_group.loc[(t)].reindex(aa_group.index.tolist()).Value
            buttom2 =buttom2+plot_type_aa_group.loc[(t)].reindex(aa_group.index.tolist()).Value/aa_group['Value']

        #axs[0].get_xaxis().set_ticks([])
        axs[0].set_ylabel("Normalized Reads (K)", fontsize=15)
        #https://pythonspot.com/matplotlib-legend/
        axs[1].legend(types,loc='upper center', bbox_to_anchor=(1.15, 1.02), ncol=1)
        axs[1].set_ylabel("Ratio", fontsize=14)
        axs[1].set_xlabel("Amino acid", fontsize=15)
        axs[1].set_xticklabels(aa_group.index.tolist(), rotation=40, horizontalalignment="right")
        axs[0].text(0.95, 0.90, s_label, transform=axs[0].transAxes,
                                                  horizontalalignment='right', fontsize=16, weight='bold')
        #plt.show()
        fig.tight_layout()
        if isTest:
            plt.show()
            break
        else:
            plt.savefig(dl.pic_out_dir + "/" + s_id + '_aa.png', dpi=100)
            print("Pictures can be found here : dl.pic_out_dir + "/" + s_id + '_aa.png'")
    print('Finished')

#drawTrfDistForAcceptersInSamples(dl, True)

wdir = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/output"
sample_tsv = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/sample.tsv"
trna_anno_bed =  "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/ChIP/hg38_tRNA_60.bed"
d = dl.LoadWDir(wdir, sample_tsv, trna_anno_bed)
drawTrfDistForSamples(d)
#drawTrfDistForAcceptersInSamples(d, True)