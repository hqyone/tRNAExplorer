# coding=utf-8

#
#

import numpy as np
import pandas as pd
import os, io, re
from IPython.display import display

def LoadConfig(visual_config, report_dir=""):
    obj = {}
    if os.path.isfile(visual_config):
        V =  open(visual_config,'r',encoding='utf-8')
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

def getSampleIndex(id, df):
    if len(df.index[df['#SampleID'] == id])>0:
        return df.index[df['#SampleID'] == id][0]
    else:
        return -1

def getSampleLabel(sampleID, s_df):
    index = getSampleIndex(sampleID, s_df)
    if index>0:
        return s_df.at[index, 'Description']
    else:
        return sampleID


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
        data['wdir']=wdir
        data["trf_sample_matrix"] = ana_dir+"/trf_sample_matrix.tsv"
        data["trna_trf_type_matrix"]=ana_dir+"/trna_trftype_matrix.tsv"
        data["sample_trf_type_matrix"] = ana_dir+"/sample_trftype_matrix.tsv"
        data["profiles"] = ana_dir+"/profiles.tsv"
        data["trna_sample_readcount_matrix"] = ana_dir + "/trna_sample_readcount_matrix.tsv"
        data["trna_sample_pileup_matrix"] = ana_dir + "/trna_sample_pileup_matrix.tsv"
        data["variants"] = ana_dir + "/variants.tsv"
        data["static_log"] = ana_dir + "/static.log"
        data["cleavage_sites"] = ana_dir + "/cleavage_sites.tsv"
        data["trna_anno_bed"] = trna_anno_bed
        data["report_dir"] = report_dir
        data["sample_tsv"] = sample_tsv
        data["output_dir"] = os.path.abspath(__file__)+"/output"

        df = pd.read_csv(data["profiles"], sep="\t", index_col=False)
        data["df"] = df

        if os.path.isfile(data["static_log"]):
            print(data["static_log"])
            s_df = pd.read_csv(data["static_log"], sep="\t")
            # Only processed samples count
            data["s_df"] =s_df
            data["sample_ls"] = s_df["#SampleID"].unique()
            data["types"] = ['full_U_tRNA','full_tRNA','5_U_tRNA_halve','5_tRNA_halve','5_U_tRF', '5_tRF','3_U_tRNA_halve' ,'3_tRNA_halve','3_U_tRF','3_tRF' ,'i-tRF','other' ]
            data["types_colors"] = {
                'full_U_tRNA': 'darkred',
                'full_tRNA': 'red',
                '5_U_tRNA_halve': 'darkgreen',
                '5_tRNA_halve': 'limegreen',
                '5_U_tRF': 'lightgreen',
                '5_tRF': 'greenyellow',
                '3_U_tRNA_halve': 'navy',
                '3_tRNA_halve': 'blue',
                '3_U_tRF': 'dodgerblue',
                '3_tRF': 'lightblue',
                'i-tRF': 'gold',
                'other': 'grey'
            }
            # Get Expression dataframe
            exp_df = pd.read_csv(data["trf_sample_matrix"], sep="\t", index_col=False)
            data["exp_df"] = add_aa_column(exp_df,trna_id="RNA_IDs")

            count_exp_df = pd.read_csv(data["trna_sample_readcount_matrix"], sep="\t", index_col=False)
            data["count_exp_df"] = add_aa_column(count_exp_df, trna_id="tRNA_ID")

            pileup_exp_df = pd.read_csv(data["trna_sample_pileup_matrix"], sep="\t", index_col=False)
            data["pileup_exp_df"] = add_aa_column(pileup_exp_df, trna_id="tRNA_ID")
            
            cleavage_df = pd.read_csv(data["cleavage_sites"], sep="\t", index_col=False)
            data["cleavage_df"] = add_aa_column(cleavage_df, trna_id="tRNA_ID")

            data["trf_exp_df"] = pd.read_csv(data["trf_sample_matrix"], sep="\t", index_col=False)

            data["sample_dic"] = {}
            for s in data["sample_ls"]:
                sample_des = getSampleLabel(s, s_df)
                data["sample_dic"][s] = sample_des
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

def add_AA_statistic(df, trna_id="3_read"):
    AA_len = []
    AA_tail = []
    for i in df[trna_id]:
        g = re.search('^(A+)(\w*)', i)
        if g:
            AA_len.append(len(str(g.group(1))))
            AA_tail.append(str(g.group(2)))
        else:
            AA_len.append(0)
            AA_tail.append("")

    df['AA_length'] = AA_len  # Isoacceptor families
    df["AA_tail"] = AA_tail  # Isodecoders
    return df

def getTotalIntensity(exp_df, sample_ls):
    sample_dic = {}
    for s in sample_ls:
        if s in exp_df.columns:
            value = exp_df[s].sum()
            sample_dic[s]=value
    return sample_dic

def csv_download_link(df, csv_file_name, delete_prompt=True):
    """Display a download link to load a data frame as csv from within a Jupyter notebook"""
    import os
    file_path = os.path.dirname(os.path.abspath(__file__))+"/output"+"/"+csv_file_name
    df.to_csv(file_path, index=True, sep='\t')
    from IPython.display import FileLink
    display(FileLink("./output"+"/"+csv_file_name))
    if delete_prompt:
        a = input('Press enter to delete the file after you have downloaded it.')
        os.remove(file_path)

#df = d["exp_df"]
#target_cols = list(df.columns[df.columns.str.contains('ENC')])

def filterDF(df, num=1, value=100, col_names=[]):
    if num>len(col_names):
        print('Error, check the parameter. Num>len(col_numes)')
        return df
    else:
        df['exp_sample_num'] = df[df[col_names]>=value].count(axis=1)  #Count Specific Values in rows
        f_df = df[df['exp_sample_num']>=num]
        return f_df

#high_confid_df = filterDF(df, num=2, value=1000, col_names=target_cols)
#dl.csv_download_link(high_confid_df, 'high_confid_trfs.tsv', delete_prompt=False)
        
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