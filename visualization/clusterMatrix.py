#
#

import scipy
import pylab
import scipy.cluster.hierarchy as sch
import numpy as np
from scipy.spatial.distance import squareform
import pandas as pd
from scipy.spatial import distance_matrix
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import data_loader as dl
plt.rcParams["axes.grid"] = False

def generateDisMatrix(df):
    return pd.DataFrame(distance_matrix(df.values, df.values), index=df.index, columns=df.index)


# Generate distance matrix.
# data = [[5, 7], [7, 3], [8, 1]]
# ctys = ['Boston', 'Phoenix', 'New York']
# df = pd.DataFrame(data, columns=['xcord', 'ycord'], index=ctys)
def drawCorrelationClusterMatrix(df, fig_name, fig_size=12, method="pearson", annot=True, annot_font_size=7):
    #D = generateDisMatrix(df.transpose())
    D = df.corr(method=method)
    condensedD  = D
    #sns.set(font_scale=0.8)
    sns.set(rc={'figure.figsize': (15, 15)})
    g =sns.clustermap(D,figsize=(fig_size, fig_size), cmap = pylab.cm.YlOrRd,annot=annot, annot_kws={"size": annot_font_size}, fmt=".1f")
    for a in g.ax_row_dendrogram.collections:
        a.set_linewidth(2)

    for a in g.ax_col_dendrogram.collections:
        a.set_linewidth(2)
    fig = g.fig
    #fig.subplots_adjust(bottom=0.3, right=0.7)
    fig.show()
    g.savefig(fig_name)
    #condensedD = squareform(D)

    # # Compute and plot first dendrogram
    # fig = pylab.figure(figsize=(10,10))
    # ax1 = fig.add_axes([0.05,0.25,0.1,0.6],frame_on=False)
    # Y = sch.linkage(condensedD, method='single')
    # Z1 = sch.dendrogram(Y, orientation='left',distance_sort='descending')
    # ax1.set_xticks([])
    # ax1.set_yticks([])
    #
    # # Compute and plot second dendrogram.
    # ax2 = fig.add_axes([0.15,0.86,0.6,0.1],frame_on=False)
    # Y = sch.linkage(condensedD, method='single') # centroid
    # Z2 = sch.dendrogram(Y)
    # ax2.set_xticks([])
    # ax2.set_yticks([])
    #
    # # Plot distance matrix.
    #
    # axmatrix = fig.add_axes([0.15,0.25,0.6,0.6])
    #
    # idx1 = Z1['leaves']
    # idx2 = Z2['leaves']
    #
    # #D = D.loci[idx1,:]
    # #D = D.loci(:,idx2)
    #
    # im = axmatrix.matshow(D, aspect='auto', cmap=pylab.cm.YlOrRd)  #pylab.cm.YlGnBu
    # #im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    # # axmatrix.set_xticks([])
    # # axmatrix.set_yticks([])
    # axmatrix.set_xticks(range(len(idx1)))
    # rownames = [sub.replace('.sorted', '') for sub in D.index.values.tolist()]
    # axmatrix.set_xticklabels(rownames, minor=False, fontsize=12)
    # axmatrix.xaxis.set_label_position('bottom')
    # axmatrix.xaxis.tick_bottom()
    #
    # pylab.xticks(rotation=90, fontsize=12)
    #
    # axmatrix.set_yticks(range(len(idx2)))
    # colnames = [sub.replace('.sorted', '') for sub in D.columns.values.tolist()]
    # axmatrix.set_yticklabels(colnames, minor=False,fontsize=12)
    # axmatrix.yaxis.set_label_position('right')
    # axmatrix.yaxis.tick_right()
    #
    # # Plot colorbar.
    # axcolor = fig.add_axes([0.90,0.25,0.02,0.6])
    # pylab.colorbar(im, cax=axcolor)
    # fig.show()
    # fig.savefig(fig_name)

# Draw the clustered expression matrix in different levels
# Levels: 'RNA','RNA_family','CODE', 'AA', 'CODE_AA'
# group = exp_df.groupby(['CODE_AA']).sum()
# group = np.log10(group+1)
# df = normalize(group)
# drawExpressionClusterMatrix(df,'matrix.png',fig_width=6, fig_height=10)
def drawExpressionClusterMatrix(df, fig_name,
                                xl_method = "single",
                                yl_method = "centroid", #single
                                xaxis_font_size=12,
                                yaxis_font_size=12,
                                cmap=pylab.cm.YlGnBu,
                                fig_width =20,
                                fig_height=20,
                                layout=[0.25,0.18,0.55,0.73]):
    # Compute and plot first dendrogram
    matrix_x_r = layout[0]
    matrix_y_r = layout[1]
    matrix_w_r = layout[2]
    matrix_h_r = layout[3]

    # Compute and plot x dendrogram.
    pylab.rcParams["axes.grid"] = False
    fig = pylab.figure(figsize=(fig_width, fig_height))

    ax1 = fig.add_axes([matrix_x_r*(1-0.8), matrix_y_r,matrix_x_r*0.8,matrix_h_r],frame_on=False)
    Y = sch.linkage(df, method=yl_method)
    Z1 = sch.dendrogram(Y, orientation='left') #,distance_sort='descending'
    ax1.set_xticks([])
    ax1.set_yticks([])

    # Compute and plot y dendrogram.
    ax2 = fig.add_axes([matrix_x_r,matrix_y_r+matrix_h_r,matrix_w_r,0.8*(1-matrix_y_r-matrix_h_r)],frame_on=False)
    Y = sch.linkage(df.T, method=xl_method) # centroid, single, average
    Z2 = sch.dendrogram(Y)  #distance_sort='descending'
    ax2.set_xticks([])
    ax2.set_yticks([])

    # Plot distance matrix.
    axmatrix = fig.add_axes([matrix_x_r,matrix_y_r,matrix_w_r,matrix_h_r])

    idx1 = Z1['leaves']
    idx2 = Z2['leaves']

    df = df.iloc[idx1,:]
    df = df.iloc[:,idx2]

    dl.csv_download_link(df, "exp_matrix.tsv", delete_prompt=False)
    im = axmatrix.imshow(df, aspect='auto', origin='lower',cmap=cmap)  #pylab.cm.YlGnBu
    plt.grid(None)
    # #im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([]) # Hide tick lines
    axmatrix.set_yticks([])
    axmatrix.set_xticks(range(len(idx2))) # Show the trick lines
    axmatrix.set_yticks(range(len(idx1))) # Show the trick lines

    colnames = [sub.replace('.sorted', '') for sub in df.columns.values.tolist()]
    axmatrix.set_xticklabels(colnames, minor=False, fontsize=xaxis_font_size)
    #axmatrix.xaxis.set_label_position('bottom')
    #axmatrix.xaxis.tick_bottom()
    #
    pylab.xticks(rotation=90, fontsize=16)
    rownames = [sub.replace('.sorted', '') for sub in df.index.values.tolist()]
    axmatrix.set_yticklabels(rownames, minor=False, fontsize=yaxis_font_size)
    #axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()

    axmatrix.grid(b=None)

    # Plot colorbar.
    axcolor = fig.add_axes([matrix_x_r*(1-0.8),matrix_y_r+matrix_h_r,0.02,0.8*(1-matrix_y_r-matrix_h_r)])
    pylab.colorbar(im, cax=axcolor)

    # Plot annotations
    row_hist_width_r = 0.03
    ls_data = list(range(len(df.index.values)))
    row_hist_ls = [
       #{"name": 'test', "data": ls_data}
    ]
    col_hist_width_r = 0.03
    col_hist_ls = [
        #{"name": 'test', "data": ls_data}
    ]
    axmatrix.tick_params(axis='y', pad=3+20*len(row_hist_ls))
    axmatrix.tick_params(axis='x', pad=3+20*len(col_hist_ls))
    for i in range(len(row_hist_ls)):
        his_obj = row_hist_ls[i]
        a = fig.add_axes([matrix_x_r + matrix_w_r + 0.01, matrix_y_r, row_hist_width_r, matrix_h_r])
        langs = ['C', 'C++', 'Java', 'Python', 'PHP']
        # ls_data = list(range(len(langs)))
        students = [23, 17, 35, 29, 12]

        # a = fig.add_axes([0,0, 100, 100])
        # b = df.iloc[:, [1]]
        # b[b.columns[0]].mask(b[b.columns[0]]>0.5, 5,inplace=True)
        # b[b.columns[0]].mask(b[b.columns[0]] < 0.2, 0, inplace=True)
        # ls = list(map(lambda s: s.split("_")[1], list(df.index)))
        ls = []
        tRNA_ls = list(df.index)
        for i in tRNA_ls:
            c = i.split("_")
            if len(c)>1:
                ls.append(c[1])
            else:
                ls.append(i)
        label_dic = {}
        for i in range(len(ls)):
            if ls[i] not in label_dic:
                label_dic[ls[i]]=i
            ls[i]=label_dic[ls[i]]
        anno_df = pd.DataFrame(ls)
        a.matshow(anno_df, aspect='auto', origin='lower', cmap="tab20")
        #a.matshow(anno_df)
        # a.barh(langs, students)
        a.grid(b=None)
    fig.show()
    fig.savefig(fig_name)

def getExpressionDistMatrix(df, fig_name, cutoff=None, normalized_by="column", exp_level = 24, cols=[]):
    n_df =None
    if len(cols)==0:
        cols = df.columns
    else:
        df=df[cols]
    # Do normalization
    if normalized_by == "column":
        n_df = df.div(df.sum(axis=0), axis=1)
    elif normalized_by == "row":
        n_df = (df.T / df.T.sum()).T
    else:
        n_df = df
    if not cutoff:
        cutoff = n_df.median().median()/2
    #Count the sample number in whihc the tRNA expression is above some cutoff
    n_df['sample_count']=(n_df[cols] > cutoff).sum(1)
    count_max = int(n_df['sample_count'].values.max())
    count_min = int(n_df['sample_count'].values.min())
    # Calulated average expression level for genes in samples which are above the cutoff
    gene_ls = list(df.index.values)
    D = df.T
    g_mean_exp_dic = {}
    vmin = 1000000
    vmax = 0
    for g in gene_ls:
        subset_df = D[D[g] > cutoff]
        mean = subset_df[g].mean()
        if mean==mean:
            g_mean_exp_dic[g] = mean
        else:
            g_mean_exp_dic[g] = 0
        if g_mean_exp_dic[g]>vmax:
            vmax = g_mean_exp_dic[g]
        if g_mean_exp_dic[g]<vmin:
            vmin = g_mean_exp_dic[g]
    a = np.zeros(shape=(exp_level,count_max+1))
    vmax+=cutoff*0.0001
    for g in gene_ls:
        value_index = int(((g_mean_exp_dic[g]-vmin)*exp_level)/(vmax-vmin))
        count = n_df['sample_count'][g]
        a[value_index,count] += 1
    show_df = pd.DataFrame(a)
    fig = pylab.figure(figsize=(10, 10))
    axe = fig.add_axes([0.3,0.2,0.5,0.5])
    axe.imshow(show_df, aspect='auto', origin='lower', cmap="hot")
    axe.yaxis.set_label_position("right")
    #axe.xaxis.set_ticks(np.arange(0, 10))
    axe.xaxis.set_ticks_position('bottom')
    # axe.xaxis.set_ticks(np.arange(0, 9))
    # axe.xaxis.set_ticklabels(list(np.arange(0, 9)))
    axe.set_ylabel("Expression level", fontsize=16)
    axe.yaxis.set_ticks_position('right')
    # axe.yaxis.set_ticks(np.arange(0, 10))
    # axe.yaxis.set_ticklabels(list(np.arange(0, 10)))
    # We want to show all ticks...
    axe.set_xticks(np.arange(len(show_df.columns+1)))
    axe.set_yticks(np.arange(exp_level))

    axe.set_xticks(np.arange(show_df.shape[1] + 1) - .5, minor=True)
    axe.set_yticks(np.arange(show_df.shape[0] + 1) - .5, minor=True)
    axe.grid(which="minor", color="grey", linestyle='-', linewidth=1)
    axe.tick_params(which="minor", bottom=False, left=False)
    # ... and label them with the respective list entries
    #axe.set_xticklabels(list(np.arange(9)))
    #axe.set_yticklabels(list(np.arange(10)))
    axe.set_xlabel("The number of expressed samples", fontsize=16)
    axe.legend()

    # Loop over data dimensions and create text annotations.
    for i in range(len(show_df.index)):
        for j in range(len(df.columns)):
            if int(show_df.iloc[i, j])>0:
                if (int(show_df.iloc[i, j])<30):
                    text = axe.text(j, i, int(show_df.iloc[i, j]),
                                   ha="center", va="center", color="w")
                else:
                    text = axe.text(j, i, int(show_df.iloc[i, j]),
                                    ha="center", va="center", color="b")

    #axe.tick_params(axis='both', which='both', length=0)


    left_hist_pos = {
        "x": 0.195,
        "y": 0.2,
        "w": 0.1,
        "h": 0.5,
    }
    left_axe = fig.add_axes([left_hist_pos["x"],
                             left_hist_pos["y"],
                             left_hist_pos["w"],
                             left_hist_pos["h"]])
    left_hist_data = list(show_df.sum(axis=1).values)
    left_hist_max = max(left_hist_data)
    left_hist_min = min(left_hist_data)
    for i in range(len(left_hist_data)):
        w =  left_hist_data[i]
        y = i/len(left_hist_data)
        x = 0
        h = 1/len(left_hist_data)
        rect = patches.Rectangle((x, y),w, h,
                             linewidth=1, edgecolor='b', facecolor='b',
                             alpha=0.6)
        left_axe.add_patch(rect)
    left_axe.set_yticks([])
    left_axe.set_xlim(left_hist_min, left_hist_max)
    left_axe.invert_xaxis()
    # left_axe.xaxis.set_ticks(np.arange(left_hist_min,left_hist_max, 4))
    # a = list(np.arange(left_hist_min, left_hist_max, 4))
    # a.reverse()
    # left_axe.xaxis.set_ticklabels(a)
    #left_axe.set_xticks([])
    #left_axe.grid(b=None)

    top_hist_pos = {
        "x": 0.3,
        "y": 0.705,
        "w": 0.5,
        "h": 0.1,
    }
    top_axe = fig.add_axes([top_hist_pos["x"],
                             top_hist_pos["y"],
                             top_hist_pos["w"],
                             top_hist_pos["h"]])
    top_hist_data = list(show_df.sum(axis=0).values)
    top_hist_max = max(top_hist_data)
    top_hist_min = min(top_hist_data)
    for i in range(len(top_hist_data)):
        h = top_hist_data[i]
        y = 0
        x = i/len(top_hist_data)
        w = 1 / len(top_hist_data)
        rect = patches.Rectangle((x, y), w, h,
                                 linewidth=1, edgecolor='b', facecolor='b',
                                 alpha=0.6)
        top_axe.add_patch(rect)
    #top_axe.set_yticks([])
    top_axe.set_xticks([])
    top_axe.set_ylim(0, top_hist_max)
    #left_axe.invert_xaxis()
    #top_axe.grid(b=None)

    #axe.set_xticks([])  # Hide tick lines
    #axe.set_yticks([])  # Hide tick lines
    fig.tight_layout()
    fig.show()
    fig.savefig(fig_name)

def testImshow():
    vegetables = ["cucumber", "tomato", "lettuce", "asparagus",
                  "potato", "wheat", "barley"]
    farmers = ["Farmer Joe", "Upland Bros.", "Smith Gardening",
               "Agrifun", "Organiculture", "BioGoods Ltd.", "Cornylee Corp."]

    harvest = np.array([[0.8, 2.4, 2.5, 3.9, 0.0, 4.0, 0.0],
                        [2.4, 0.0, 4.0, 1.0, 2.7, 0.0, 0.0],
                        [1.1, 2.4, 0.8, 4.3, 1.9, 4.4, 0.0],
                        [0.6, 0.0, 0.3, 0.0, 3.1, 0.0, 0.0],
                        [0.7, 1.7, 0.6, 2.6, 2.2, 6.2, 0.0],
                        [1.3, 1.2, 0.0, 0.0, 0.0, 3.2, 5.1],
                        [0.1, 2.0, 0.0, 1.4, 0.0, 1.9, 6.3]])

    fig, ax = plt.subplots()
    im = ax.imshow(harvest)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(farmers)))
    ax.set_yticks(np.arange(len(vegetables)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(farmers)
    ax.set_yticklabels(vegetables)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(vegetables)):
        for j in range(len(farmers)):
            text = ax.text(j, i, harvest[i, j],
                           ha="center", va="center", color="w")

    ax.set_title("Harvest of local farmers (in tons/year)")
    fig.tight_layout()
    plt.show()

