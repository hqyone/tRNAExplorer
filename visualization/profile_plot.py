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
import sys
import numpy as np
import matplotlib.patches as patches


def getSampleIndex(id, df):
    if len(df.index[df['#ID'] == id])>0:
        return df.index[df['#ID'] == id][0]
    else:
        return -1

def getSampleNum(tRNA_ID, df):
    s_num = 0
    sample_ls = df['#SampleID'].unique()
    for s in sample_ls:
        t_df = df[(df['RNA_ID'] == tRNA_ID) & (df['#SampleID'] == s) & (df['type'] == 'total')]
        if len(t_df['profile']) > 0:
            profile = t_df['profile'].values[0].split(",")
            profile_max = -1
            for i in profile:
                if (float(i) > profile_max):
                    profile_max = float(i)
            if (profile_max > 0):
                s_num += 1
    return s_num


def getSampleLabel(sampleID, s_df):
    index = getSampleIndex(sampleID, s_df)
    if index>0:
        return s_df.at[index, 'Description']
    else:
        return sampleID

# Generate rectangle for tRNAs
def getRects(RNA_ID, df):
    rect_ls = []
    sel_df = df[df['name'] == RNA_ID]
    if len(sel_df) > 0:
        rect_ls.append({
            "start": int(sel_df['map_start'].values[0]),
            "end": int(sel_df['map_end'].values[0]),
            "color": "g", "alpha": 0.01
        })
        if sel_df['d_loop_start'].values[0] != -1:
            rect_ls.append({
                "start": int(sel_df['map_start'].values[0]) + int(sel_df['d_loop_start'].values[0]),
                "end": int(sel_df['map_start'].values[0]) + int(sel_df['d_loop_end'].values[0]),
                "color": "b", "alpha": 0.01
            })
        if sel_df['a_loop_start'].values[0] != -1:
            rect_ls.append({
                "start": int(sel_df['map_start'].values[0]) + int(sel_df['a_loop_start'].values[0]),
                "end": int(sel_df['map_start'].values[0]) + int(sel_df['a_loop_end'].values[0]),
                "color": "r", "alpha": 0.01
            })
        if sel_df['t_loop_start'].values[0] != -1:
            rect_ls.append({
                "start": int(sel_df['map_start'].values[0]) + int(sel_df['t_loop_start'].values[0]),
                "end": int(sel_df['map_start'].values[0]) + int(sel_df['t_loop_end'].values[0]),
                "color": "b", "alpha": 0.01
            })
        if sel_df['anticodon_start'].values[0] != -1:
            rect_ls.append({
                "start": int(sel_df['anticodon_start'].values[0]),
                "end": int(sel_df['anticodon_end'].values[0]),
                "color": "r", "alpha": 0.03
            })
    return rect_ls

def drawPic(dl, test=False):
    df = pd.read_csv(dl["combined_mean_profile"], sep="\t")

    s_df = pd.read_csv(dl["sample_tsv"], sep="\t",index_col=False)
    # s_df = s_df[s_df['sample_type']=='total']
    # s_df = s_df.sort_values(by=['cell_type','cell_line_name'])

    st_df = pd.read_csv(dl["trna_anno_bed"], sep='\t',index_col=False)

    #profile_types = ['start_pos', 'end_pos', 'uniq_read', 'multiple_read', 'total']

    profile_types = ['start_pos', 'end_pos', 'uniq_read', 'total']
    profile_types_label = ['start_pos', 'end_pos', 'uniq_read', 'total_read']
    profile_types_color = ['r', 'b', 'pink', 'g']
    profile_types_alpha = [0.4, 0.4, 0.3, 0.2]
    profile_types_line_type = ['-', '-', '--', '-']
    profile_types_line_color = ['red', 'black', 'black', 'green']
    # plt.figure(figsize=(8,6))  # sets the window to 8 x 6 inches

    sample_ls = df["#SampleID"]



    legend_patchs = []
    for i in range(len(profile_types)):
        legend_patchs.append(patches.Patch(facecolor=profile_types_color[i],
                                           alpha=profile_types_alpha[i],
                                           label=profile_types_label[i],
                                           linestyle=profile_types_line_type[i],
                                           edgecolor=profile_types_line_color[i],
                                           linewidth=2))

    sample_ls = df['#SampleID'].unique()
    rna_ls = df['RNA_ID'].unique()

    print(len(rna_ls))

    sorted_sample_matrix = pd.DataFrame(sample_ls, columns=['SampleID'])
    label_ls = []
    for s in sample_ls:
        label_ls.append(getSampleLabel(s, s_df))
    sorted_sample_matrix['label'] = label_ls
    sorted_sample_matrix = sorted_sample_matrix.sort_values(by=["label"])

    for tRNA_ID in rna_ls:
        # tRNA_ID = 'tRNA-Lys-CTT-7-1' #"nm-tRNA-Tyr-GTA-chr17-33"  "tRNA-Lys-CTT-7-1"
        s_num = getSampleNum(tRNA_ID, df)
        fig, axs = plt.subplots(s_num, 1, figsize=[12, s_num+2])
        if s_num>1:
            s_index = 0
            for sample in sorted_sample_matrix['SampleID']:
                t_df = df[(df['RNA_ID'] == tRNA_ID) & (df['#SampleID'] == sample) & (df['type'] == 'total')]
                if len(t_df['profile']) > 0:
                    profile = np.array(t_df['profile'].values[0].split(",")).astype(np.float)
                    profile_max = profile.max()
                    if profile_max > 0:
                        # print(t_df['end_pos'].values[0].split(','))
                        profile = np.array(t_df['profile'].values[0].split(',')).astype(np.float)
                        profile_max = profile.max()
                        if profile_max > 0:
                            p_index = 0
                            for i in profile_types:
                                # print(sample+"_"+i)
                                sourse_df = df[(df['RNA_ID'] == tRNA_ID) & (df['#SampleID'] == sample) & (df['type'] == i)]
                                # print(sourse_df)
                                if len(sourse_df.index) > 0:
                                    profile = np.array(sourse_df['profile'].values[0].split(',')).astype(np.float)
                                    max_value = profile.max()
                                    rects = getRects(tRNA_ID, st_df)
                                    # Create a Rectangle patch
                                    for re in rects:
                                        rect = patches.Rectangle((re['start'], 0), re['end'] - re['start'], profile_max,
                                                                 linewidth=1, edgecolor=re['color'], facecolor=re['color'],
                                                                 alpha=re['alpha'])
                                        axs[s_index].add_patch(rect)
                                    axs[s_index].fill_between(range(1, len(profile) + 1), profile,
                                                              color=profile_types_color[p_index],
                                                              alpha=profile_types_alpha[p_index],
                                                              linestyle=profile_types_line_type[p_index], linewidth=2)
                                    axs[s_index].set_xlim([20, 160])
                                    if (profile_max > 10):
                                        axs[s_index].set_ylabel(str(int(profile_max)))
                                    else:
                                        axs[s_index].set_ylabel(str(float('{0:.1f}'.format(profile_max))))
                                p_index += 1
                            if s_index < s_num - 1:
                                axs[s_index].get_xaxis().set_ticks([])
                            else:
                                axs[s_index].set_xlabel("Location (nt)", fontsize=12)
                                axs[s_index].legend(handles=legend_patchs, loc='lower center', bbox_to_anchor=(0.5, -1.5),
                                                    ncol=4, fontsize=14)
                            if s_index == 0:
                                axs[s_index].text(0.0, 1.3, tRNA_ID, transform=axs[s_index].transAxes,
                                                  va='center', fontsize=16, weight='bold')
                                axs[s_index].text(0.38, 1.2, 'D_Loop', transform=axs[s_index].transAxes,
                                                  va='center', fontsize=10, weight='bold', bbox=dict(facecolor='b', alpha=0.1))
                                axs[s_index].text(0.5, 1.2, 'Anticode', transform=axs[s_index].transAxes,
                                                  va='center', fontsize=10, weight='bold',
                                                  bbox=dict(facecolor='red', alpha=0.1))
                                location = 0.67
                                sel_df = st_df[st_df['name'] == tRNA_ID]
                                if len(sel_df) > 0 and sel_df['t_loop_start'].values[0] != -1:
                                    location = float((int(sel_df['map_start'].values[0]) - 24) +
                                                     (int(sel_df['t_loop_start'].values[0]) + int(
                                                         sel_df['t_loop_end'].values[0])) / 2) / 140
                                axs[s_index].text(location, 1.2, 'T_Loop', transform=axs[s_index].transAxes,
                                                  va='center', fontsize=10, weight='bold', bbox=dict(facecolor='b', alpha=0.1))
                            if s_index == int(s_num/2):
                                axs[s_index].text(-0.06, 1.3, "Pileup deepth", transform=axs[s_index].transAxes,
                                                  va='center', fontsize=16, weight='bold',rotation=90)
                            axs[s_index].get_yaxis().set_ticks([])
                            txt = axs[s_index].text(0.04, 0.8, getSampleLabel(sample, s_df), transform=axs[s_index].transAxes,
                                                    va='center', fontsize=10, weight='bold',
                                                    bbox=dict(facecolor='none', alpha=0.0))
                        s_index += 1
        #plt.show()
        plt.savefig(dl["report_dir"] + "/" + tRNA_ID + '_profile.png', dpi=200)
        if (test):
            plt.show()
            return
    print("The profile pictures are stored in : "+dl.pic_out_dir)

wdir = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/output"
sample_tsv = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/sample.tsv"
trna_anno_bed =  "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/ChIP/hg38_tRNA_60.bed"
d = dl.LoadWDir(wdir, sample_tsv, trna_anno_bed)
drawPic(d, True)