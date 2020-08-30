import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import data_loader as dl

def getVariantRatioTabInFamily(d):
    fv = pd.DataFrame(columns=['sample', 'family','loc','RNA_IDs','mem_num','members','ref','muts','mut_reads','total_reads','ratio'])
    v = pd.read_csv(d["variants"], sep="\t")
    # Combine mutations for each tRNA, here we just sum the mutation reads and keep total_reads not change.
    gv = v.groupby(['#SampleID','family','tRNA_ID','loc','ref']).sum()
    gv['mut_reads'] = v.groupby(['#SampleID','family','tRNA_ID','loc','ref'])['mut_reads'].sum()
    gv['total_reads'] = v.groupby(['#SampleID','family','tRNA_ID','loc','ref'])['total_reads'].mean()
    gv['muts'] = v.groupby(['#SampleID','family','tRNA_ID','loc','ref'])['mut'].apply(','.join)
    gv['mut_num'] = v.groupby(['#SampleID','family','tRNA_ID','loc','ref'])['mut_reads'].apply(lambda x: ','.join(x.astype(int).astype(str)))
    gv = gv.reset_index()
    #print(gv)
    # Combine mutations for each tRNA family, here we just sum both the mutation reads and total_reads not change.
    fv = gv.groupby(['#SampleID','family','loc','ref']).sum()
    fv['mut_reads'] = gv.groupby(['#SampleID','family','loc','ref'])['mut_reads'].sum()
    fv['total_reads'] = gv.groupby(['#SampleID','family','loc','ref'])['total_reads'].sum()
    fv['ratio'] = fv['mut_reads']/fv['total_reads']
    fv['muts'] = gv.groupby(['#SampleID','family','loc','ref'])['muts'].apply(','.join)
    fv['mut_num'] = gv.groupby(['#SampleID','family','loc','ref'])['mut_num'].apply(','.join)
    fv['tRNA_IDs'] = gv.groupby(['#SampleID','family','loc','ref'])['tRNA_ID'].apply(','.join)
    fv['meam_num'] = gv.groupby(['#SampleID','family','loc','ref'])['tRNA_ID'].count()
    fv = fv.reset_index()
    
    #Explain for transform https://pbpython.com/pandas_transform.html#:~:text=Understanding%20the%20Transform%20Function%20in%20Pandas%201%20Introduction.,...%204%20Second%20Approach%20-%20Using%20Transform.%20
    fv['ratio_max'] = fv.groupby(['family','loc','ref'])['ratio'].transform('max')
    fv['mut_read_mean'] = fv.groupby(['family','loc','ref'])['mut_reads'].transform('mean')
    
    print("Download tsv here:")
    dl.csv_download_link(fv,'family_mut.tsv', delete_prompt=False) 
    sns.reset_defaults()
    g = sns.FacetGrid(fv, row="#SampleID", height=1.7, aspect=4)
    g.map(sns.distplot, 'ratio',kde=False, bins=10)
    plt.xlim(0, 1)
    plt.figure()
    g = sns.FacetGrid(fv, row="#SampleID", height=1.7, aspect=4)
    g.map(sns.distplot, 'loc',kde=False, bins=75)
    plt.xlim(0, 75)
    plt.show()
    return fv