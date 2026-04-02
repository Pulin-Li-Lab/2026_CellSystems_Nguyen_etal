#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
import anndata as ad
import seaborn as sns

# Input files 
region = "regionC"
# processed expression matrix split by regions 
count_file = pd.read_csv(region+"_sct_260226.csv", header = 0, index_col = 0)
cop = count_file
# cell type ident for each cell 
celltype = pd.read_csv(region+"_celltypes_260225.csv", index_col = 0)
# X, Y, Z coordinate of each cell 
meta = pd.read_csv(region+"_metadata_260225.csv", index_col = 0)
cop.index = cop.index.str.replace('.', '-', regex=False)

# find infected cells and associated information
ind =  cop['NP.gRNA'] >= 1
infect = cop.loc[ind]
infect_cell = celltype.loc[ind]
infect_meta = meta.loc[ind]
infect_coord = pd.concat([infect_meta['x'],infect_meta['y']], axis=1)
coord = np.array(infect_coord)

# use DBSCAN to find clusters of infected cells 
# cells within 175 pixel of distance (20 um) are considered neighbors; 
# clusters of at least 3 cells are considered 
dbscan = DBSCAN(eps=175, min_samples=3) 
clusters = dbscan.fit_predict(coord)
clust = np.unique(clusters)
clust = clust[1:len(clust)]

# find the composition of infected cell types for each plaque 
blank = np.zeros((len(clust), 23))
plaque_info = pd.DataFrame(blank, columns = ["ClubEpi",
               "CiliatedEpi",
               "AT2",
               "AT1",
               "AdventitialFib",
               "AlveolarFib",
               "SMC",
               "SCMF",
               "Pericyte",
               "Artery",
               "Vein",
               "Capillary",
               "LymphaticEndo",
               "AlveolarMac",
               "InterstitialMac",
               "InflamMac",
               "PatrolMac",
               "DC",
               "pDC",
               "Neutrophil",
               "Tcell (CD4)",
               "Tcell (CD8)",
               "Bcell"])
for i in clust:
    sub_ind = clusters == i
    sub_infect = infect.loc[sub_ind]
    sub_cell = infect_cell.loc[sub_ind]
    cell_infect = sub_cell['x'].unique()
    num_infect = sub_cell['x'].value_counts(sort=False)
    # print(num_infect)
    for j in range(len(cell_infect)):
        if num_infect.iloc[j] == 0:
            plaque_info.loc[i, cell_infect[j]] = 0
        else:
            plaque_info.loc[i, cell_infect[j]] = num_infect.iloc[j]
plaque_info.to_csv(region+'_plaque_info.csv')

# find the composition of IFN+ cell types for each plaque 
blank = np.zeros((len(clust), 23))
ifn_info = pd.DataFrame(blank, columns = ["ClubEpi",
               "CiliatedEpi",
               "AT2",
               "AT1",
               "AdventitialFib",
               "AlveolarFib",
               "SMC",
               "SCMF",
               "Pericyte",
               "Artery",
               "Vein",
               "Capillary",
               "LymphaticEndo",
               "AlveolarMac",
               "InterstitialMac",
               "InflamMac",
               "PatrolMac",
               "DC",
               "pDC",
               "Neutrophil",
               "Tcell (CD4)",
               "Tcell (CD8)",
               "Bcell"])
ifn = ['Ifnb1', 'Ifnl2', 'Ifna2', 'Ifna4']
for i in clust:
    sub_ind = clusters == i
    sub_infect = infect.loc[sub_ind]
    sub_cell = infect_cell.loc[sub_ind]
    ifn_pos = sub_infect[ifn].sum(axis=1) > 0
    sub_ifn_cell = sub_cell.loc[ifn_pos]
    num_ifn = sub_ifn_cell['x'].value_counts(sort=False)
    cell_ifn = sub_ifn_cell['x'].unique()
    for j in range(len(cell_ifn)):
        ifn_info.loc[i, cell_ifn[j]] = num_ifn.iloc[j]
ifn_info.T.to_csv(region+'_ifn_info.csv')

#%% some visualization 
# look at flu genes -- normalized values for convex hull plotting only, otherwise use expression values 
flu = ['NP.gRNA']
flu_exp = cop[flu].sum(axis=1)
flu_norm = flu_exp / np.max(flu_exp)
np_level = np.linspace(0, 1, 10)
np_thresh_level = np.digitize(flu_norm, np_level)

# look at Ifn genes -- normalized values for convex hull plotting only, otherwise use expression values 
ifn = ['Ifnb1', 'Ifnl2', 'Ifna2', 'Ifna4']
ifn_exp = cop[ifn].sum(axis=1)
ifn_norm = ifn_exp / np.max(ifn_exp)
ifn_level = np.linspace(0, 1, 10)
ifn_thresh_level = np.digitize(ifn_norm, ifn_level)

isg = ['Isg15']
isg_exp = cop[isg].sum(axis=1)
isg_norm = isg_exp / np.max(isg_exp)
isg_level = np.linspace(0, 1, 10)
isg_thresh_level = np.digitize(isg_norm, isg_level)

# create adata object for convex hull plotting 
from scipy import sparse
mtx_count = sparse.csr_matrix(cop) 
adata = ad.AnnData(mtx_count)
adata.obs.index = meta.index
adata.obs['region'] = meta['region'].tolist()
adata.obs['column'] = meta['x'].tolist()
adata.obs['row'] = meta['y'].tolist()
adata.obs['volume'] = meta['volume'].tolist()
adata.obs['celltype'] = celltype['x'].tolist()
adata.obs['nCount_RNA'] = meta['nCount_RNA'].tolist()
adata.obs['flu'] = np_thresh_level
adata.obs['ifn'] = ifn_thresh_level
adata.obs['isg'] = isg_thresh_level

def createConvexHulls(adata, readspath, scale_t = 1,reads_min = 10, reads_max = 100000):
    '''
    Create convex hull borders for cells containing between `reads_min` and `reads_max` reads.
    Returns hulls and barcodes for the corresponding cells.
    @param adata: AnnData object containing spatial data. The `adata.obs.index` should has the $anything-$cell_barcode format.
    @param readspath: path to the reads file
    @param scale_t: scaling factor for the coordinates
    @param reads_min: minimum number of reads for a cell to be included in the convex hull
    @param reads_max: maximum number of reads for a cell to be included in the convex hull
    e.g.
        hulls, barcodes = createConvexHulls(sdata_t,reads_path,scale_t = scale_dict['63'])
        sdata_t2 = sdata_t[np.isin(sdata_t.obs.index.values,barcodes),:].copy()
        plotHulls(sdata_t2, hulls, barcodes, col_t, ctype_color_dict, figsize_scale = 1/200,output_path_t=f'{output_path_t}/{sample_t}_{tech_t}_convex_hull.pdf',scale_bar_t=[200,'200 µm'], legend = True)
    '''
    from scipy.spatial import ConvexHull
    # Get reads information
    reads = pd.read_csv(readspath, index_col=0)
    # Get a list of cell barcodes
    cell_ids = adata.obs.index.values
    cells_unique = np.unique(cell_ids)
    # print(cell_ids)
    Nlabels = cells_unique.shape[0]
    hulls = []
    barcodes = []
    print('Creating cell convex hulls')
    for i in range(len(cells_unique)):
        cell_barcode = cells_unique[i]
        n_reads = adata.obs.loc[cell_barcode, 'nCount_RNA']
        if isinstance(n_reads, pd.Series):
            n_reads = n_reads[0]
            print(f'{cell_barcode} has multiple reads, using the first one')
        if n_reads >= reads_min and n_reads <= reads_max:
            barcode_num = int(cell_barcode.split('-')[1])
            sub_reads = reads.loc[reads['cell_barcode']==barcode_num]
            # print(sub_reads)
            curr_coords = np.array(sub_reads[['column', 'row']]) * scale_t
            barcodes.append(cell_barcode)
            hulls.append(ConvexHull(curr_coords))
    print("Used %d / %d" % (len(barcodes), Nlabels))
    return hulls, barcodes

def plotHulls(adata, hulls, barcodes, color_by, color_dict,linewidth_t = 0.1,figsize_scale = 0.01, legend=None, output_path_t = None,scale_bar_t = None):
    """ Plotting function for hulls
    @param adata: full or subsetted adata (will take longer for more hulls)
    @param hulls: pre-computed scipy.spatial.qhull for each unique cell (barcode)
    @param barcodes: accompanying list of cell barcodes to relate hulls to obs in adata
    @param color_by: name of column in adata.obs (discrete)
    @param color_dict: dictionary of colors to use for each unique value in color_by
    @param linewidth_t: width of hull lines
    @param figsize_scale: scale of figure size [max(column)-min(column),max(row)-min(row)] * figsize_scale
    @param legend: whether to include a legend
    @param output_path_t: path to save figure
    @param scale_bar_t: tuple of (size, units) for scale bar
    e.g.
        hulls, barcodes = createConvexHulls(sdata_t,reads_path,scale_t = scale_dict['63'])
        sdata_t2 = sdata_t[np.isin(sdata_t.obs.index.values,barcodes),:].copy()
        plotHulls(sdata_t2, hulls, barcodes, col_t, ctype_color_dict, figsize_scale = 1/200,output_path_t=f'{output_path_t}/{sample_t}_{tech_t}_convex_hull.pdf',
        scale_bar_t=[200,'200 µm'], legend = True)
    """
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Patch
    # Initialize figure
    plt.rcParams['axes.facecolor'] = 'white'
    fig = plt.figure(figsize=((adata.obs['column'].max()-adata.obs['column'].min())*figsize_scale,
                        (adata.obs['row'].max() - adata.obs['row'].min())*figsize_scale))
    ax = fig.add_subplot(111)
    plt.xlim(adata.obs['column'].min(), adata.obs['column'].max())
    plt.ylim(adata.obs['row'].min(), adata.obs['row'].max())
    # Subset sample hulls accordingly 
    # Barcodes and hulls are matched, so grab index for each barcode in the subdata to get the hulls
    idxs = [barcodes.index(i) for i in adata.obs.index]
    plotted_hulls = [hulls[i] for i in idxs]
    cellcolors = adata.obs[color_by]
    polys = [hull_to_polygon_color(h, color_t = color_dict[cellcolors[i]]) for i,h in enumerate(plotted_hulls)]
    p = PatchCollection(polys, alpha=1, linewidth=linewidth_t, match_original=True, rasterized=False)
    ax.add_collection(p)
    if legend:
        clist = list(set(cellcolors))
        clist = sorted(clist, key=lambda x:list(color_dict.keys()).index(x)) # Reorder legend
        handles = []
        for i in clist:
            handles.append(Patch(color=color_dict[i], label=i))
        ax.legend(handles=handles, title=color_by, loc='upper left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)
    plt.axis('off')
    plt.grid(False)
    if scale_bar_t is not None:
        add_scale_bar(scale_bar_t[0],scale_bar_t[1])
    if output_path_t is not None:
        plt.savefig(f'{output_path_t}', dpi=300, bbox_inches='tight')
    else:
        plt.show()
def hull_to_polygon_color(hull,k =0.9,color_t = None):
    from matplotlib.patches import Polygon
    cent = np.mean(hull.points, 0)
    pts = []
    for pt in hull.points[hull.simplices]:
        pts.append(pt[0].tolist())
        pts.append(pt[1].tolist())
    pts.sort(key=lambda p: np.arctan2(p[1] - cent[1],
                                    p[0] - cent[0]))
    pts = pts[0::2]  # Deleting duplicates
    pts.insert(len(pts), pts[0])
    if color_t is not None:
        poly = Polygon(k*(np.array(pts)- cent) + cent,edgecolor='k', linewidth=1,facecolor = color_t)
    else:
        poly = Polygon(k*(np.array(pts)- cent) + cent,edgecolor='k', linewidth=1)
    poly.set_capstyle('round')
    return poly
def add_scale_bar(length_t,label_t,ax_t = None):
    '''
    Add scale bar to the plot
    If there are multiple axes, then ax_t should be specified.
    If ax_t is None, then the scale bar will be added to the current axis (entire plot).
    '''
    import matplotlib.font_manager as fm
    from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
    fontprops = fm.FontProperties(size=20, family='Arial')
    if ax_t is None:
        bar = AnchoredSizeBar(plt.gca().transData, length_t, label_t, 4, pad=0.1,
                    sep=5, borderpad=0.5, frameon=False,
                    size_vertical=0.1, color='black',fontproperties = fontprops)
        plt.gca().add_artist(bar)
    else:
        bar = AnchoredSizeBar(ax_t.transData, length_t, label_t, 4, pad=0.1,
                    sep=5, borderpad=0.5, frameon=False,
                    size_vertical=0.1, color='black',fontproperties = fontprops)
        ax_t.add_artist(bar)

# change directory to the folder/region of interest 
reads_path = '/Volumes/T9/TranscriptomicsData/STARmap/36hpi_250genes_sum/lung_250genes_processed_sum/Reads_location/I/remain_reads_1.csv'
hulls, barcodes = createConvexHulls(adata,reads_path)
adata_sub = adata[barcodes,:]

# plot cell types by color
cell_color = {
    "ClubEpi": "#1F77B4", 
    "CiliatedEpi": "#AEC7E8",
    "AT2":"#FF7F0E",
    "AT1": "#FFBB78", 
    "AdventitialFib": "#2CA02C",  
    "AlveolarFib":"#98DF8A", 
    "SMC":"#D62728", 
    "SCMF": "#FF9896", 
    "Pericyte":  "#9467BD", 
    "Artery": "#8C564B",  
    "Vein":"#C49C94", 
    "Capillary":  "#E377C2", 
    "LymphaticEndo":"#F7B6D2",  
    "AlveolarMac":  "#7F7F7F",  
    "InterstitialMac": "#C7C7C7",  
    "InflamMac":"#BCBD22", 
    "PatrolMac": "#DBDB8D",  
    "DC": "#17BECF",
    "pDC": "#9EDAE5",  
    "Neutrophil": "#393B79", 
    "Tcell (CD4)": "#637939",
    "Tcell (CD8)": "#B5CF6B",
    "Bcell":   "#843C39"  
}

plotHulls(adata_sub, hulls, barcodes, 'celltype', cell_color, legend = True, figsize_scale = 1/200,scale_bar_t=[200,'200 µm'])

# plot convex hull representing normlized flu level
color = sns.color_palette("rocket_r", 10)
cell_color = {
    1:color[0], 
    2:color[1], 
    3:color[2], 
    4:color[3], 
    5:color[4],
    6:color[5], 
    7:color[6], 
    8:color[7], 
    9:color[8], 
    10:color[9]           
}
plotHulls(adata_sub, hulls, barcodes, 'flu', cell_color, legend = True, figsize_scale = 1/200,scale_bar_t=[200,'200 µm'])


