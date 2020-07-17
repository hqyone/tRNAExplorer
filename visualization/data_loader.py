# coding=utf-8

#
#

import numpy as np
import pandas as pd
import os, io, re

def LoadConfig(visual_config, report_dir=""):
    obj = {}
    if os.path.isfile(visual_config):
        V =  open(visual_config,'r')
        for line in V:
            contents = line.strip().split("=")
            if len(contents)>1:
                obj[contents[0].strip()] = contents[1].strip()
        V.close()
        if "out_dir" in obj and "sample_tsv" in obj and "trna_anno_bed" in obj:
            d =  LoadWDir(obj["out_dir"],obj["sample_tsv"], obj["trna_anno_bed"],report_dir=report_dir)
            if os.path.isdir(report_dir):
                d["report_dir"]=report_dir
            return d
    return None


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

        data["trf_sample_matrix"] = ana_dir+"/trf_sample_matrix.tsv"
        data["trna_trf_type_matrix"]=ana_dir+"/trna_trftype_matrix.tsv"
        data["sample_trf_type_matrix"] = ana_dir+"/sample_trftype_matrix.tsv"
        data["profiles"] = ana_dir+"/profiles.tsv"
        data["trna_sample_readcount_matrix"] = ana_dir + "/trna_sample_readcount_matrix.tsv"
        data["trna_sample_pileup_matrix"] = ana_dir + "/trna_sample_pileup_matrix.tsv"
        data["static_log"] = ana_dir + "/static.log"
        data["trna_anno_bed"] = trna_anno_bed
        data["report_dir"] = report_dir
        data["sample_tsv"] = sample_tsv

        df = pd.read_csv(data["profiles"], sep="\t", index_col=False)
        data["df"] = df

        if os.path.isfile(data["static_log"]):
            print(data["static_log"])
            s_df = pd.read_csv(data["static_log"], sep="\t")
            # Only processed samples count
            data["s_df"] =s_df
            data["sample_ls"] = s_df["#SampleID"].unique()
            data["types"] = ['full_U_tRNA','full_tRNA','5_U_tRNA_halve','5_tRNA_halve','5_U_tRF', '5_tRF','3_U_tRNA_halve' ,'3_tRNA_halve','3_U_tRF','3_tRF' ,'i-tRF','other' ]
            data["types_colors"] = ['darkred','red','darkgreen','limegreen', 'lightgreen','greenyellow','navy' ,'blue','dodgerblue','lightblue' ,'gold','grey' ]

            # Get Expression dataframe
            exp_df = pd.read_csv(data["trf_sample_matrix"], sep="\t", index_col=False)
            data["exp_df"] = add_aa_column(exp_df,trna_id="RNA_IDs")

            count_exp_df = pd.read_csv(data["trna_sample_readcount_matrix"], sep="\t", index_col=False)
            data["count_exp_df"] = add_aa_column(count_exp_df, trna_id="tRNA_ID")

            pileup_exp_df = pd.read_csv(data["trna_sample_pileup_matrix"], sep="\t", index_col=False)
            data["pileup_exp_df"] = add_aa_column(pileup_exp_df, trna_id="tRNA_ID")

            data["trf_exp_df"] = pd.read_csv(data["trf_sample_matrix"], sep="\t", index_col=False)

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

def add_aa_column(df, trna_id="RNA_ID"):
    code_ls = []
    aa_ls = []
    code_aa_ls = []
    for i in df[trna_id]:
        g = re.search('tRNA-(\w+)-(\w+)', i)
        if g:
            code_ls.append(g.group(2))
            aa_ls.append(g.group(1))
            code_aa_ls.append(g.group(2) + "_" + g.group(1))
        else:
            code_ls.append('None')
            aa_ls.append('None')
            code_aa_ls.append(g.group(2) + "_" + g.group(1))

    df['AA'] = aa_ls  # Isoacceptor families
    df["CODE"] = code_ls  # Isodecoders
    df["CODE_AA"] = code_aa_ls
    return df

def getTotalIntensity(exp_df, sample_ls):
    sample_dic = {}
    for s in sample_ls:
        if s in exp_df.columns:
            value = exp_df[s].sum()
            sample_dic[s]=value
    return sample_dic



#wdir = '/Users/hqyone/OneDrive/国内的工作/学术/科研项目/论文课题/tRNA/data_output/rna_seq'
#data = LoadRESQDir(wdir)
#print(data)

def getSampleLabel(s_id,sample_df):
    label = "None"
    if len(sample_df[sample_df['#SampleID']==s_id])>0:
        label = sample_df[sample_df['#SampleID'] == s_id]['Description'].values[0]
    return label


# wdir="/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/output"
# sample_tsv="/Users/hqyone/OneDrive/MyProjects/testrepo/new_tools/tRNAExplorer/test_data/RNASeq/sample.tsv"
# LoadWDir(wdir, sample_tsv)