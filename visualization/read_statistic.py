import matplotlib.pyplot as plt
import data_loader as dl
import pandas as pd
import sys
import numpy as np
import matplotlib.patches as patches
import os

plt.rcParams['axes.facecolor'] = 'black'

base_color_dic = {
    "A":"r",
    "T":"g",
    "G":"y",
    "C":"b"
}

def getAddModificationStatisitc(d,proj_name="test",top_num=10,radius = 1.2, fontsize=14, fig_width=14):
    '''
    The function shows the 5' and 3' base addition modifications of tRFs.
    The function will print the top 5' and 3' addition modifications
    And draw four types of pie charts:
        1. A pie chart showing the ratio of tRFs with or without 5'-Addition modification
        2. A pie chart showing the composition of 5'-Addition modifications
        3. A pie chart showing the ratio of tRFs with or without 3'-Addition modification
        4. A pie chart showing the composition of 3'-Addition modifications
    @param d: The data object generated by data_loader.py
    @param proj_name: the name of project, the default is "test", use can check the name of hit.tab file
            it follows the pattern : <sampleID>+"_"+<proj_name>+"_hit.tab"
    @param top_num: The number of top addition modifications to be printed out
    @param radius: the radius of pies
    @param fontsize:
    @param fig_width:
    @return: None
    '''
    sample_dic = d["sample_dic"]
    s_num = len(sample_dic.keys())

    sample_ls = list(sample_dic.keys())
    sample_ls.sort()

    fig, axs = plt.subplots(s_num, 4, figsize=[fig_width, s_num * 3 + 2])

    index = 0
    for s in sample_ls:
        des = sample_dic[s]
        hit_tab = d['wdir'] + "/" + s + "_" + proj_name + "_hit.tab"
        df = pd.read_csv(hit_tab, sep="\t")

        df = dl.add_aa_column(df, trna_id="tRNA_id")

        df['temp5'] = df['read_5_fragment']
        group_a = df.groupby('read_5_fragment').sum()
        axs[index, 1].pie(group_a["mean_number"], labels=group_a.index, autopct='%1.1f%%', shadow=False, radius=radius,
                          labeldistance=None)
        print("\n" + s + "_" + des)
        sorted_df = group_a.sort_values(["mean_number"], ascending=[False]).head(10)
        with pd.option_context('display.max_rows', None, 'display.max_columns',
                               None):  # more options can be specified also
            print(sorted_df['mean_number'])
        group_c = df.groupby('read_3_fragment').sum()
        sorted_df = group_c.sort_values(["mean_number"], ascending=[False]).head(10)
        with pd.option_context('display.max_rows', None, 'display.max_columns',
                               None):  # more options can be specified also
            print(sorted_df['mean_number'])
        axs[index, 3].pie(group_c["mean_number"], labels=group_c.index, autopct='%1.1f%%', shadow=False, radius=radius,
                          labeldistance=None)
        # axs[index,1].text(-0.06, 0.5, des, transform=axs[index,0].transAxes,va='center', fontsize=12, weight='bold', rotation=90)

        df.loc[df['temp5'].isna(), 'read_5_fragment'] = '-'
        df.loc[df['temp5'].notna(), 'read_5_fragment'] = 'Addition'
        group_b = df.groupby('read_5_fragment').sum()
        axs[index, 0].pie(group_b["mean_number"], labels=group_b.index, autopct='%1.1f%%', shadow=False, radius=radius)
        axs[index, 0].text(-0.16, 0.5, des, transform=axs[index, 0].transAxes, va='center', fontsize=fontsize,
                           weight='bold', rotation=90)

        df['temp3'] = df['read_3_fragment']
        df.loc[df['temp3'].isna(), 'read_3_fragment'] = '-'
        df.loc[df['temp3'].notna(), 'read_3_fragment'] = 'Addition'
        group = df.groupby('read_3_fragment').sum()
        axs[index, 2].pie(group["mean_number"], labels=group.index, autopct='%1.1f%%', shadow=False, radius=radius)
        if index == 0:
            axs[index, 0].text(0.2, 1.1, "5'_Addition_Ratio", transform=axs[index, 0].transAxes, va='center',
                               fontsize=fontsize, weight='bold', rotation=0)
            axs[index, 1].text(1.4, 1.1, "5'_Addition_Comp", transform=axs[index, 0].transAxes, va='center',
                               fontsize=fontsize, weight='bold', rotation=0)
            axs[index, 2].text(2.8, 1.1, "3'_Addition_Ratio", transform=axs[index, 0].transAxes, va='center',
                               fontsize=fontsize, weight='bold', rotation=0)
            axs[index, 3].text(4.0, 1.1, "3'_Addition_Comp", transform=axs[index, 0].transAxes, va='center',
                               fontsize=fontsize, weight='bold', rotation=0)
        index += 1

