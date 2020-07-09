# coding=utf-8

#
#  tRNAExplorer v.1.0 (2020.) - tRNA profiles and regulation networks
#  Copyright (C) 2020  Dr. Quanyuan He
#  School of Medicine, Hunan Normal University
#  Email: hqyone@hotmail.com
#  Freely distributed under the GNU General Public License (GPLv3)
#

import numpy as np
import pandas as pd
import os, io, re


#   1) <out_dir>/mean_trna_sample_matrix.tsv  # trna vs sample matrix
#   2) <out_dir>/mean_trna_trf_type_matrix.tsv  # sample:trna vs trf type matrix
#   3) <out_dir>/mean_sample_trf_type_matrix.tsv  # Sample vs trf types (12) matrix
#   4) <out_dir>/combined_mean_profile.tsv  # Profiles tab files (including starpos, endpos, total)
#   5) <out_dir>/trna_sample_readcount_matrix.tsv  # trna vs sample: read_count matrix
#   6) <out_dir>/trna_sample_pileup_matrix.tsv  # trna vs sample: pileup max matrix

def LoadWDir(wdir, sample_tsv, trna_anno_bed, report_dir=""):
    if os.path.isdir(wdir):
        data={}
        #wdir = '/Users/hqyone/OneDrive/国内的工作/学术/科研项目/论文课题/tRNA/data_output/rna_seq'
        ana_dir = wdir
        if report_dir=="":
            report_dir = wdir+"/reports"
        if not os.path.isdir(report_dir):
            try:
                os.makedirs(report_dir)
            except:
                print("An exception occurred during create pictures dir")

        data["mean_trna_sample_matrix"] = ana_dir+"/mean_trna_sample_matrix.tsv"
        data["mean_trna_trf_type_matrix"]=ana_dir+"/mean_trna_trf_type_matrix.tsv"
        data["mean_sample_trf_type_matrix"] = ana_dir+"/mean_sample_trf_type_matrix.tsv"
        data["combined_mean_profile"] = ana_dir+"/combined_mean_profile.tsv"
        data["trna_sample_readcount_matrix"] = ana_dir + "/trna_sample_readcount_matrix.tsv"
        data["trna_sample_pileup_matrix"] = ana_dir + "/trna_sample_pileup_matrix.tsv"
        data["static_log"] = ana_dir + "/static.log"
        data["trna_anno_bed"] = trna_anno_bed
        data["report_dir"] = report_dir

        data["sample_tsv"] = sample_tsv
        data["sample_static_tab"]= ana_dir+"/bam_stat.tsv"

        df = pd.read_csv(data["combined_mean_profile"], sep="\t", index_col=False)
        data["df"] = df

        if os.path.isfile(sample_tsv):
            print(sample_tsv)
            s_df = pd.read_csv(sample_tsv, sep="\t")
            #s_df = s_df[s_df['sample_type']=='total']
            #s_df = s_df.sort_values(by=['cell_type','cell_line_name'])

            data["s_df"] =s_df
            data["sample_ls"] = s_df["#ID"].unique()
            data["types"] = ['full_U_tRNA','full_tRNA','5_U_tRNA_halve','5_tRNA_halve','5_U_tRF', '5_tRF','3_U_tRNA_halve' ,'3_tRNA_halve','3_U_tRF','3_tRF' ,'i-tRF','other' ]
            data["types_colors"] = ['darkred','red','darkgreen','limegreen', 'lightgreen','greenyellow','navy' ,'blue','dodgerblue','lightblue' ,'gold','grey' ]
            # Get EXP_df
            exp_df = pd.read_csv(data["mean_trna_sample_matrix"], sep="\t", index_col=False)
            code_ls = []
            aa_ls = []
            for i in exp_df["RNA_IDs"]:
                g = re.search('tRNA-(\w+)-(\w+)', i)
                if g:
                    code_ls.append(g.group(2))
                    aa_ls.append(g.group(1))
                else:
                    code_ls.append('None')
                    aa_ls.append('None')

            exp_df['AA'] = aa_ls  # Isoacceptor families
            exp_df["CODE"] = code_ls  #Isodecoders
            data["exp_df"]= exp_df

            # Get structure information
            st_df = pd.read_csv(trna_anno_bed, sep='\t')
            data["st_df"] = st_df
            return data
        else:
            print("The sample sheet ("+sample_tsv+") is not exist!")
            print("Please check it!")
            return None
    else:
        print("Work Directory : "+wdir+" is not exist")
        return None

#wdir = '/Users/hqyone/OneDrive/国内的工作/学术/科研项目/论文课题/tRNA/data_output/rna_seq'
#data = LoadRESQDir(wdir)
#print(data)

def getSampleLabel(s_id,sample_df):
    label = "None"
    if len(sample_df[sample_df['#ID']==s_id])>0:
        label = sample_df[sample_df['#ID'] == s_id]['Description'].values[0]
    return label


# wdir="/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/output"
# sample_tsv="/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/sample.tsv"
# LoadWDir(wdir, sample_tsv)