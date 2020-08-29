# coding=utf-8

#
#

import matplotlib.pyplot as plt
import data_loader as dl
import pandas as pd
import sys
import numpy as np
import matplotlib.patches as patches
import os
import re

plt.rcParams['axes.facecolor'] = 'black'

base_color_dic = {
    "A":"r",
    "T":"g",
    "G":"y",
    "C":"b"
}

def getSampleIndex(id, df):
    if len(df.index[df['#SampleID'] == id])>0:
        return df.index[df['#SampleID'] == id][0]
    else:
        return -1

def getSampleTotalCount(s_df, sampleID):
    return s_df[s_df['#SampleID'] == sampleID]["total"].values[0]

def getProfilesStatic(tRNA_ID, p_df,s_df, min_depth, normalize=False):
    s_num = 0
    prof_max = 0

    sample_ls = p_df['#SampleID'].unique()
    for s in sample_ls:
        t_df = p_df[(p_df['tRNA_ID'] == tRNA_ID) & (p_df['#SampleID'] == s) & (p_df['type'] == 'total')]
        if len(t_df['profile']) > 0:
            profile = np.array(t_df['profile'].values[0].split(',')).astype(np.float)
            if normalize:
                sample_total_count = getSampleTotalCount(s_df, s)
                profile = np.divide(profile,sample_total_count/100000.0)
            profile_max = -1
            for i in profile:
                if (float(i) > profile_max):
                    profile_max = float(i)
            if profile_max>prof_max:
                prof_max = profile_max
            if (profile_max >= min_depth):
                s_num += 1
    return [s_num, prof_max]


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
            "color": "g", "alpha": 0.1
        })
        if sel_df['d_loop_start'].values[0] != -1:
            rect_ls.append({
                "start": int(sel_df['map_start'].values[0]) + int(sel_df['d_loop_start'].values[0]),
                "end": int(sel_df['map_start'].values[0]) + int(sel_df['d_loop_end'].values[0]),
                "color": "b", "alpha": 0.1
            })
        if sel_df['a_loop_start'].values[0] != -1:
            rect_ls.append({
                "start": int(sel_df['map_start'].values[0]) + int(sel_df['a_loop_start'].values[0]),
                "end": int(sel_df['map_start'].values[0]) + int(sel_df['a_loop_end'].values[0]),
                "color": "r", "alpha": 0.1
            })
        if sel_df['t_loop_start'].values[0] != -1:
            rect_ls.append({
                "start": int(sel_df['map_start'].values[0]) + int(sel_df['t_loop_start'].values[0]),
                "end": int(sel_df['map_start'].values[0]) + int(sel_df['t_loop_end'].values[0]),
                "color": "b", "alpha": 0.1
            })
        if sel_df['anticodon_start'].values[0] != -1:
            rect_ls.append({
                "start": int(sel_df['anticodon_start'].values[0]),
                "end": int(sel_df['anticodon_end'].values[0]),
                "color": "r", "alpha": 0.3
            })
    return rect_ls

def getBaseRects(sequence):
    rect_ls = []
    index = 0
    for b in sequence:
        color = 'w'
        if b in base_color_dic:
           color = base_color_dic[b]
        rect_ls.append({
            "start": index,
            "end": index+1,
            "color": color,
            "alpha": 0.1
        })
        index+=1
    return rect_ls

def getMutRecs(Mutation_Strs):
    mut_color_dic={
        'I':'w',
        #'M':'b',
        '-':'maroon',
        'A':'r',
        'T':'b',
        'G':'g',
        'C':'y',
        'D':'k'
    }
    rect_ls =[]
    if Mutation_Strs=="-" or Mutation_Strs=="":
        return rect_ls
    mut_ls = Mutation_Strs.split(",")
    y_dic={}
    for mu in mut_ls:
        c = mu.split("=")
        if len(c)>1:
            type, infor, loc, total = c[0].split(":")
            if type=="M":
                type=infor.split('>')[1]
            base = ""
            base_pair =  infor.split(">")
            base = base_pair[0]
            if len(base_pair)>1:
                base = base_pair[1]
            intensity =  round(float(c[1]))
            mut_color = 'k'
            if type in mut_color_dic:
                mut_color = mut_color_dic[type]
            if not loc in y_dic:
                y_dic[loc] = 0
            rect_ls.append({
                "start": int(loc)-1,
                "end": int(loc),
                "color": mut_color,
                "y_loc":y_dic[loc],
                "rect_height":intensity,
                "alpha": 1.0
            })
            y_dic[loc]+=intensity
    return rect_ls

def CombineMutList(Mutation_Str_ls):
    mut_dic = {}
    combined_mut_str = ""
    combine_total_mut_dic ={}
    for MutStrs in Mutation_Str_ls:
        total_dic = {}
        if MutStrs=="-" or MutStrs=="":
            continue
        else:
            mut_ls = MutStrs.split(",")
            for mu in mut_ls:
                c = mu.split("=")
                if len(c)>1:
                    type, infor, loc, total = c[0].split(":")
                    mut_id = type+':'+infor+":"+loc
                    mut_total_id = type + ':' + loc
                    intensity = float(c[1])
                    if not mut_id in mut_dic:
                        mut_dic[mut_id]={
                            'intensity':0,
                            'total':0
                        }
                    mut_dic[mut_id]['intensity']+=intensity
                    total_dic[mut_total_id]=float(total)
        for key in total_dic:
            if not key in combine_total_mut_dic:
                combine_total_mut_dic[key]=0
            combine_total_mut_dic[key]+=total_dic[key]
    for mut_id in mut_dic:
        mut_obj = mut_dic[mut_id]
        mut_ls = mut_id.split(':')
        mut_total_id = mut_ls[0]+':'+mut_ls[2]
        total_number = combine_total_mut_dic[mut_total_id]
        combined_mut_str+=mut_id+':'+str(round(total_number))+'='+str(round(mut_obj['intensity']))+","
    return combined_mut_str.strip(',')


def CombineProfilesWithIntron(sample, cur_profile, tRNA1_ID,tRNA2_ID, d, profile_type):
    p_df = pd.read_csv(d["profiles"], sep="\t")
    st_df = d["st_df"]
    #print(st_df.columns)

    t2_df = p_df[(p_df['tRNA_ID'] == tRNA2_ID) & (p_df['#SampleID'] == sample) & (p_df['type'] == profile_type)]
    profile_1 = cur_profile
    profile_2 = np.array(t2_df['profile'].values[0].split(",")).astype(np.float)

    intron_1_str = st_df[st_df['name'] == tRNA1_ID]['intron_infor'].values[0]
    #print(intron_1_str)
    intron_1_obj = {}
    m = re.search(r'\((\d+)-(\d+)\)', intron_1_str)
    if m:
        intron_1_obj['start'] = int(m.group(1))
        intron_1_obj['end'] = int(m.group(2))
    #print(intron_1_str)
    intron_2_str = st_df[st_df['name'] == tRNA2_ID]['intron_infor'].values[0]
    intron_2_obj = {}
    m = re.search(r'\((\d+)-(\d+)\)', intron_2_str)
    if m:
        intron_2_obj['start'] = int(m.group(1))
        intron_2_obj['end'] = int(m.group(2))

    index1 = 0
    index2 = 0
    for i in range(0, len(profile_1)):
        if index1 == intron_1_obj['start'] - 1:
            index1 = intron_1_obj['end']
        if index2 == intron_2_obj['start'] - 1:
            index2 = intron_2_obj['end']
        if index1 < len(profile_1) and index2 < len(profile_2):
            profile_1[index1] += profile_2[index2]
        index1 += 1
        index2 += 1
    return profile_1

def CombineProfilesByFamily(d, profile_types):
    p_df = pd.read_csv(d["profiles"], sep="\t")
    sample_ls = p_df['#SampleID'].unique()
    rna_family_ls = p_df['tRNA_family'].unique()
    # Combined profile in a tRNA family
    new_p_df = pd.DataFrame(None, columns=list(p_df.columns))
    for s in sample_ls:
        for f in rna_family_ls:
            for p_type in profile_types:
                t_df = p_df[(p_df['tRNA_family'] == f) & (p_df['#SampleID'] == s) & (p_df['type'] == p_type)]
                new_row = {}
                profile = []
                for index, row in t_df.iterrows():
                    new_row['#SampleID'] = s
                    new_row['tRNA_family'] = f
                    new_row['tRNA_ID'] = f.replace('TF_',"")
                    new_row['type'] = p_type
                    cur_profile = np.array(row['profile'].split(',')).astype(np.float)
                    cur_RNA_ID = row['tRNA_ID']
                    if len(profile)==0:
                        profile = cur_profile
                    else:
                        if (len(profile)!=len(cur_profile)):
                            f_tRNA_ID= f.replace('TF_',"")
                            profile = CombineProfilesWithIntron(s, cur_profile, f_tRNA_ID, cur_RNA_ID, d, p_type)
                        else:
                            profile = np.add(profile, cur_profile)
                profile_str_ls = []
                for v in profile:
                    profile_str_ls.append(str(v))
                new_row['profile'] = ','.join(profile_str_ls)
                new_p_df = new_p_df.append(new_row,ignore_index=True )
            mut_df = p_df[(p_df['tRNA_family'] == f) & (p_df['#SampleID'] == s) & (p_df['type'] == 'mutation_dic_str')]
            new_row = {}
            mutation_str_ls = []
            for index, row in mut_df.iterrows():
                new_row['#SampleID'] = s
                new_row['tRNA_family'] = f
                new_row['tRNA_ID'] = f.replace('TF_',"")
                new_row['type'] = 'mutation_dic_str'
                mutation_str_ls.append(row['profile'])
            new_row['profile'] = CombineMutList(mutation_str_ls)
            new_p_df = new_p_df.append(new_row, ignore_index=True )
            seq_df = p_df[(p_df['tRNA_family'] == f) & (p_df['#SampleID'] == s) & (p_df['type'] == 'isequence')]
            new_p_df = new_p_df.append(seq_df, ignore_index=True )
    return new_p_df

def drawProfiles(d, groupby='F', fig_num=1, output_dir ="", min_depth = 20, share_y_lim=True, normalized=True):
    '''
        Draw tRF profiles for each tRNA across samples
        @param d: The data object generated by data_loader.py
        @param groupby: F(Family)/Other(tRNA genes)
        @param num: Int, If it >0 will print profiles for num tRNAs for testing. If num <=0 means prints all profiles
                     or print all profiles in output_dir
        @param output_dir: The output directory. The default is <output_dir>/reports
        @param min_depth: The minimum depth of pileup for a sample to print the profile
        @param share_y_lim: Whether all profiles share the same depth limitation across samples (Default True)
        @param normalized: whether normalize to total tRNA mapped read number (Default True)
        @return:
    '''
    if groupby=="F":
        drawProfiles2(d, fig_num=fig_num, output_dir=output_dir, min_depth=min_depth,share_y_lim=share_y_lim, normalized=normalized)
    else:
        drawProfiles1(d, fig_num=fig_num, output_dir=output_dir, min_depth=min_depth, share_y_lim=share_y_lim,
                      normalized=normalized)

def drawProfiles1(d, fig_num=1, output_dir ="", min_depth = 20, share_y_lim=True, normalized=True):
    '''
    Draw tRF profiles for each tRNA across samples
    @param d: The data object generated by data_loader.py
    @param num: Int, If it >0 will print profiles for num tRNAs for testing. If num <=0 means prints all profiles
                 or print all profiles in output_dir
    @param output_dir: The output directory. The default is <output_dir>/reports
    @param min_depth: The minimum depth of pileup for a sample to print the profile
    @param share_y_lim: Whether all profiles share the same depth limitation across samples (Default True)
    @param normalized: whether normalize to total tRNA mapped read number (Default True)
    @return:
    '''
    outdir = d["report_dir"]
    if os.path.isdir(output_dir):
        outdir = output_dir
    p_df = pd.read_csv(d["profiles"], sep="\t")
    s_df = d["s_df"]
    # s_df = s_df[s_df['sample_type']=='total']
    # s_df = s_df.sort_values(by=['cell_type','cell_line_name'])

    st_df = pd.read_csv(d["trna_anno_bed"], sep='\t', index_col=False)

    #profile_types = ['start_pos', 'end_pos', 'uniq_read', 'multiple_read', 'total']
    profile_types = ['start_pos', 'end_pos', 'total']
    profile_types_label = ['start_pos', 'end_pos', 'total_read']
    profile_types_color = ['r', 'b', 'g']
    profile_types_alpha = [0.4, 0.4,  0.2]
    profile_types_line_type = ['-', '-',  '-']
    profile_types_line_color = ['red', 'black', 'green']
    # plt.figure(figsize=(8,6))  # sets the window to 8 x 6 inches

    sample_ls = p_df["#SampleID"]

    legend_patchs = []
    for i in range(len(profile_types)):
        legend_patchs.append(patches.Patch(facecolor=profile_types_color[i],
                                           alpha=profile_types_alpha[i],
                                           label=profile_types_label[i],
                                           linestyle=profile_types_line_type[i],
                                           edgecolor=profile_types_line_color[i],
                                           linewidth=2))

    sample_ls = p_df['#SampleID'].unique()
    rna_ls = p_df['tRNA_ID'].unique()

    #print(len(rna_ls))

    sorted_sample_matrix = pd.DataFrame(sample_ls, columns=['SampleID'])
    label_ls = []
    for s in sample_ls:
        label_ls.append(getSampleLabel(s, s_df))
    sorted_sample_matrix['label'] = label_ls
    sorted_sample_matrix = sorted_sample_matrix.sort_values(by=["label"])
    trna_index = 0
    for tRNA_ID in rna_ls:
        # tRNA_ID = 'tRNA-Lys-CTT-7-1' #"nm-tRNA-Tyr-GTA-chr17-33"  "tRNA-Lys-CTT-7-1"
        s_num, g_prof_max = getProfilesStatic(tRNA_ID, p_df, s_df, min_depth, normalized)
        if s_num>=1:
            s_index = 0
            fig, axs = plt.subplots(s_num, 1, figsize=[12, s_num+2], squeeze=False)
            for sample in sorted_sample_matrix['SampleID']:
                sample_total_count = s_df[s_df['#SampleID']==sample]["total"].values[0]
                if sample_total_count<=0:
                    sample_total_count = 0
                t_df = p_df[(p_df['tRNA_ID'] == tRNA_ID) & (p_df['#SampleID'] == sample) & (p_df['type'] == 'total')]
                seq_df = p_df[(p_df['tRNA_ID'] == tRNA_ID) & (p_df['#SampleID'] == sample) & (p_df['type'] == 'isequence')]
                mut_df = p_df[(p_df['tRNA_ID'] == tRNA_ID) & (p_df['#SampleID'] == sample) & (p_df['type'] == 'mutation_dic_str')]
                family_id = t_df['tRNA_family']
                if len(t_df['profile']) > 0:
                    sample_des = s_df[s_df['#SampleID']==sample]['Description'].values[0]
                    profile = np.array(t_df['profile'].values[0].split(",")).astype(np.float)
                    if normalized:
                        profile = np.divide(profile, sample_total_count/1000.0)
                    profile_max = profile.max()
                    if profile_max >= min_depth:
                        # Draw background
                        rects = getRects(tRNA_ID, st_df)
                        # Create a Rectangle patch
                        for re in rects:
                            rect_height = profile_max
                            if share_y_lim:
                                rect_height = g_prof_max
                            rect = patches.Rectangle((re['start']-1, 0), re['end'] - re['start']+1, rect_height,
                                                     linewidth=1, edgecolor=re['color'], facecolor=re['color'],
                                                     alpha=re['alpha'])
                            axs[s_index,0].add_patch(rect)
                        # Add the base rects
                        rects = getBaseRects(seq_df['profile'].values[0])
                        for re in rects:
                            rect_height = profile_max
                            if share_y_lim:
                                rect_height = g_prof_max
                            rect = patches.Rectangle((re['start'], 0), re['end'] - re['start'], rect_height/5,
                                                     linewidth=1, edgecolor=re['color'],
                                                     facecolor=re['color'], ec=re['color'],
                                                     alpha=re['alpha'])
                            axs[s_index,0].add_patch(rect)

                        # Add the mutation rects
                        mut_str = mut_df['profile'].values[0]
                        if mut_str!="-" and mut_str!="":
                            print(sample_des+">>"+sample+">>"+tRNA_ID+">>"+mut_str)
                        rects = getMutRecs(mut_str)
                        for re in rects:
                            rect_height = re['rect_height']
                            rect_y = re['y_loc']
                            if normalized:
                                rect_height =  (rect_height*1000)/sample_total_count
                                rect_y =  (rect_y*1000)/sample_total_count
                            rect = patches.Rectangle((re['start'], rect_y), re['end'] - re['start'], rect_height,
                                                     linewidth=0, edgecolor=re['color'],
                                                     facecolor=re['color'], ec=re['color'],
                                                     alpha=re['alpha'])
                            axs[s_index,0].add_patch(rect)
                        p_index = 0
                        for i in profile_types:
                            # print(sample+"_"+i)
                            sourse_df = p_df[(p_df['tRNA_ID'] == tRNA_ID) & (p_df['#SampleID'] == sample) & (p_df['type'] == i)]
                            # print(sourse_df)
                            if len(sourse_df.index) > 0:
                                profile = np.array(sourse_df['profile'].values[0].split(',')).astype(np.float)
                                if normalized:
                                    profile = np.divide(profile, sample_total_count/1000.0)
                                axs[s_index,0].fill_between(range(1, len(profile) + 1), profile,
                                                          color=profile_types_color[p_index],
                                                          alpha=profile_types_alpha[p_index],
                                                          linestyle=profile_types_line_type[p_index],
                                                          linewidth=2)

                                axs[s_index,0].set_xlim([20, 160])
                                if share_y_lim:
                                    axs[s_index,0].set_ylim([0,g_prof_max])
                                if profile_max > 10:
                                    axs[s_index,0].set_ylabel(str(int(profile_max)))
                                else:
                                    axs[s_index,0].set_ylabel(str(float('{0:.1f}'.format(profile_max))))
                            p_index += 1
                        #axs[s_index,0].set_facecolor('xkcd:white')
                        if s_index < s_num - 1:
                            axs[s_index,0].get_xaxis().set_ticks([])
                        else:
                            #axs[s_index,0].get_xaxis().set_ticks([])
                            axs[s_index,0].set_xlabel("Location (nt)", fontsize=12)
                            axs[s_index,0].legend(handles=legend_patchs, loc='lower center', bbox_to_anchor=(0.5, -1.5),
                                                ncol=4, fontsize=14)
                        if s_index == 0:
                            axs[s_index,0].text(0.0, 1.3, family_id, transform=axs[s_index,0].transAxes,
                                              va='center', fontsize=16, weight='bold')
                            axs[s_index,0].text(0.38, 1.1, 'D_Loop', transform=axs[s_index,0].transAxes,
                                              va='center', fontsize=10, weight='bold',
                                              bbox=dict(facecolor='b', alpha=0.2))
                            axs[s_index,0].text(0.5, 1.1, 'Anticodon', transform=axs[s_index,0].transAxes,
                                              va='center', fontsize=10, weight='bold',
                                              bbox=dict(facecolor='red', alpha=0.2))
                            location = 0.67
                            sel_df = st_df[st_df['name'] == tRNA_ID]
                            if len(sel_df) > 0 and sel_df['t_loop_start'].values[0] != -1:
                                location = float((int(sel_df['map_start'].values[0]) - 24) +
                                                 (int(sel_df['t_loop_start'].values[0]) + int(
                                                     sel_df['t_loop_end'].values[0])) / 2) / 140
                            axs[s_index,0].text(location, 1.1, 'T_Loop', transform=axs[s_index,0].transAxes,
                                              va='center', fontsize=10, weight='bold',
                                              bbox=dict(facecolor='b', alpha=0.2))
                        if s_index == int(s_num / 2):
                            y_label = "Pileup depth"
                            if normalized:
                                y_label= "Normalized "+y_label
                            axs[s_index,0].text(-0.06, 1.3, y_label, transform=axs[s_index,0].transAxes,
                                              va='center', fontsize=16, weight='bold', rotation=90)
                        axs[s_index,0].get_yaxis().set_ticks([])
                        txt = axs[s_index,0].text(0.04, 0.8, sample_des,
                                                transform=axs[s_index,0].transAxes,
                                                va='center', fontsize=10, weight='bold',
                                                bbox=dict(facecolor='none', alpha=0.0))
                        s_index += 1
            #if s_index>0:
            #    plt.savefig(outdir + "/" + tRNA_ID + '_profile.png', dpi=200)
            #plt.show()
            trna_index+=1
            if fig_num>0:
                if trna_index<=fig_num:
                    plt.show()
                if fig_num==trna_index:
                    return
            else:
                plt.savefig(outdir + "/" + tRNA_ID + '_profile.png', dpi=200)
        else:
            print("None profile is above the intensity cutoff, skip")
    print("The profile pictures are stored in : " + d.pic_out_dir)

def drawProfiles2(d, fig_num=1, output_dir ="", min_depth = 20, share_y_lim=True, normalized=True):
    '''
    Draw tRF profiles for each tRNA across samples
    @param d: The data object generated by data_loader.py
    @param num: Int, If it >0 will print profiles for num tRNAs for testing. If num <=0 means prints all profiles
                 or print all profiles in output_dir
    @param output_dir: The output directory. The default is <output_dir>/reports
    @param min_depth: The minimum depth of pileup for a sample to print the profile
    @param share_y_lim: Whether all profiles share the same depth limitation across samples (Default True)
    @param normalized: whether normalize to total tRNA mapped read number (Default True)
    @return:
    '''
    outdir = d["report_dir"]
    if os.path.isdir(output_dir):
        outdir = output_dir
    p_df = pd.read_csv(d["profiles"], sep="\t")
    s_df = d["s_df"]
    # s_df = s_df[s_df['sample_type']=='total']
    # s_df = s_df.sort_values(by=['cell_type','cell_line_name'])

    st_df = pd.read_csv(d["trna_anno_bed"], sep='\t', index_col=False)

    #profile_types = ['start_pos', 'end_pos', 'uniq_read', 'multiple_read', 'total']
    profile_types = ['start_pos', 'end_pos', 'total']
    profile_types_label = ['start_pos', 'end_pos', 'total_read']
    profile_types_color = ['r', 'b', 'g']
    profile_types_alpha = [0.4, 0.4,  0.2]
    profile_types_line_type = ['-', '-',  '-']
    profile_types_line_color = ['red', 'black', 'green']
    # plt.figure(figsize=(8,6))  # sets the window to 8 x 6 inches

    sample_ls = p_df["#SampleID"]

    legend_patchs = []
    for i in range(len(profile_types)):
        legend_patchs.append(patches.Patch(facecolor=profile_types_color[i],
                                           alpha=profile_types_alpha[i],
                                           label=profile_types_label[i],
                                           linestyle=profile_types_line_type[i],
                                           edgecolor=profile_types_line_color[i],
                                           linewidth=2))

    sample_ls = p_df['#SampleID'].unique()
    rna_ls = p_df['tRNA_ID'].unique()
    # Combined profiles in a tRNA family
    p_df = CombineProfilesByFamily(d,profile_types)
    #print(len(rna_ls))

    sorted_sample_matrix = pd.DataFrame(sample_ls, columns=['SampleID'])
    label_ls = []
    for s in sample_ls:
        label_ls.append(getSampleLabel(s, s_df))
    sorted_sample_matrix['label'] = label_ls
    sorted_sample_matrix = sorted_sample_matrix.sort_values(by=["label"])
    trna_index = 0
    for tRNA_ID in rna_ls:
        # tRNA_ID = 'tRNA-Lys-CTT-7-1' #"nm-tRNA-Tyr-GTA-chr17-33"  "tRNA-Lys-CTT-7-1"
        s_num, g_prof_max = getProfilesStatic(tRNA_ID, p_df, s_df, min_depth, normalized)
        if s_num>=1:
            s_index = 0
            fig, axs = plt.subplots(s_num, 1, figsize=[12, s_num+2], squeeze=False)
            for sample in sorted_sample_matrix['SampleID']:
                sample_total_count = s_df[s_df['#SampleID']==sample]["total"].values[0]
                if sample_total_count<=0:
                    sample_total_count = 0
                t_df = p_df[(p_df['tRNA_ID'] == tRNA_ID) & (p_df['#SampleID'] == sample) & (p_df['type'] == 'total')]
                seq_df = p_df[(p_df['tRNA_ID'] == tRNA_ID) & (p_df['#SampleID'] == sample) & (p_df['type'] == 'isequence')]
                mut_df = p_df[(p_df['tRNA_ID'] == tRNA_ID) & (p_df['#SampleID'] == sample) & (p_df['type'] == 'mutation_dic_str')]
                if len(t_df['profile']) > 0:
                    sample_des = s_df[s_df['#SampleID']==sample]['Description'].values[0]
                    profile = np.array(t_df['profile'].values[0].split(",")).astype(np.float)
                    if normalized:
                        profile = np.divide(profile, sample_total_count/1000.0)
                    profile_max = profile.max()
                    if profile_max >= min_depth:
                        # Draw background
                        rects = getRects(tRNA_ID, st_df)
                        # Create a Rectangle patch
                        for re in rects:
                            rect_height = profile_max
                            if share_y_lim:
                                rect_height = g_prof_max
                            rect = patches.Rectangle((re['start']-1, 0), re['end'] - re['start']+1, rect_height,
                                                     linewidth=1, edgecolor=re['color'], facecolor=re['color'],
                                                     alpha=re['alpha'])
                            axs[s_index,0].add_patch(rect)
                        # Add the base rects
                        rects = getBaseRects(seq_df['profile'].values[0])
                        for re in rects:
                            rect_height = profile_max
                            if share_y_lim:
                                rect_height = g_prof_max
                            rect = patches.Rectangle((re['start'], 0), re['end'] - re['start'], rect_height/5,
                                                     linewidth=1, edgecolor=re['color'],
                                                     facecolor=re['color'], ec=re['color'],
                                                     alpha=re['alpha'])
                            axs[s_index,0].add_patch(rect)
                        # Add the mutation rects
                        mut_str = mut_df['profile'].values[0]
                        if mut_str != "-" and mut_str != "":
                            print(sample_des + ">>" + sample + ">>" + tRNA_ID + ">>" + mut_str)
                        rects = getMutRecs(mut_str)
                        for re in rects:
                            rect_height = re['rect_height']
                            rect_y = re['y_loc']
                            if normalized:
                                rect_height = (rect_height * 1000) / sample_total_count
                                rect_y = (rect_y * 1000) / sample_total_count
                            rect = patches.Rectangle((re['start'], rect_y), re['end'] - re['start'], rect_height,
                                                     linewidth=0, edgecolor=re['color'],
                                                     facecolor=re['color'], ec=re['color'],
                                                     alpha=re['alpha'])
                            axs[s_index, 0].add_patch(rect)
                        p_index = 0
                        for i in profile_types:
                            # print(sample+"_"+i)
                            sourse_df = p_df[(p_df['tRNA_ID'] == tRNA_ID) & (p_df['#SampleID'] == sample) & (p_df['type'] == i)]
                            # print(sourse_df)
                            if len(sourse_df.index) > 0:
                                profile = np.array(sourse_df['profile'].values[0].split(',')).astype(np.float)
                                if normalized:
                                    profile = np.divide(profile, sample_total_count/1000.0)
                                axs[s_index,0].fill_between(range(1, len(profile) + 1), profile,
                                                          color=profile_types_color[p_index],
                                                          alpha=profile_types_alpha[p_index],
                                                          linestyle=profile_types_line_type[p_index],
                                                          linewidth=2)

                                axs[s_index,0].set_xlim([20, 160])
                                if share_y_lim:
                                    axs[s_index,0].set_ylim([0,g_prof_max])
                                if profile_max > 10:
                                    axs[s_index,0].set_ylabel(str(int(profile_max)))
                                else:
                                    axs[s_index,0].set_ylabel(str(float('{0:.1f}'.format(profile_max))))
                            p_index += 1
                        #axs[s_index,0].set_facecolor('xkcd:white')
                        if s_index < s_num - 1:
                            axs[s_index,0].get_xaxis().set_ticks([])
                        else:
                            #axs[s_index,0].get_xaxis().set_ticks([])
                            axs[s_index,0].set_xlabel("Location (nt)", fontsize=12)
                            axs[s_index,0].legend(handles=legend_patchs, loc='lower center', bbox_to_anchor=(0.5, -1.5),
                                                ncol=4, fontsize=14)
                        if s_index == 0:
                            axs[s_index,0].text(0.0, 1.3, "TF_"+tRNA_ID, transform=axs[s_index,0].transAxes,
                                              va='center', fontsize=16, weight='bold')
                            axs[s_index,0].text(0.38, 1.1, 'D_Loop', transform=axs[s_index,0].transAxes,
                                              va='center', fontsize=10, weight='bold',
                                              bbox=dict(facecolor='b', alpha=0.2))
                            axs[s_index,0].text(0.5, 1.1, 'Anticode', transform=axs[s_index,0].transAxes,
                                              va='center', fontsize=10, weight='bold',
                                              bbox=dict(facecolor='red', alpha=0.2))
                            location = 0.67
                            sel_df = st_df[st_df['name'] == tRNA_ID]
                            if len(sel_df) > 0 and sel_df['t_loop_start'].values[0] != -1:
                                location = float((int(sel_df['map_start'].values[0]) - 24) +
                                                 (int(sel_df['t_loop_start'].values[0]) + int(
                                                     sel_df['t_loop_end'].values[0])) / 2) / 140
                            axs[s_index,0].text(location, 1.1, 'T_Loop', transform=axs[s_index,0].transAxes,
                                              va='center', fontsize=10, weight='bold',
                                              bbox=dict(facecolor='b', alpha=0.2))
                        if s_index == int(s_num / 2):
                            y_label = "Pileup depth"
                            if normalized:
                                y_label= "Normalized "+y_label
                            axs[s_index,0].text(-0.06, 1.3, y_label, transform=axs[s_index,0].transAxes,
                                              va='center', fontsize=16, weight='bold', rotation=90)
                        axs[s_index,0].get_yaxis().set_ticks([])
                        txt = axs[s_index,0].text(0.04, 0.8, sample_des,
                                                transform=axs[s_index,0].transAxes,
                                                va='center', fontsize=10, weight='bold',
                                                bbox=dict(facecolor='none', alpha=0.0))
                        s_index += 1
            #if s_index>0:
            #    plt.savefig(outdir + "/" + tRNA_ID + '_profile.png', dpi=200)
            #plt.show()
            trna_index+=1
            if fig_num>0:
                if trna_index<=fig_num:
                    plt.show()
                if fig_num==trna_index:
                    return
            else:
                plt.savefig(outdir + "/" + tRNA_ID + '_profile.png', dpi=200)

        else:
            print("None profile is above the intensity cutoff, skip!")
    print("The profile pictures are stored in : " + d.pic_out_dir)


# wdir = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/output"
# sample_tsv = "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/sample.tsv"
# trna_anno_bed =  "/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/ChIP/hg38_tRNA_60.bed"
# d = dl.LoadWDir(wdir, sample_tsv, trna_anno_bed)
# drawPic(d, True)