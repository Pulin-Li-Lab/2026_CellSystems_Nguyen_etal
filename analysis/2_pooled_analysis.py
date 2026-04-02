#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns
from sklearn.cluster import KMeans
import numpy as np
import matplotlib as mpl

# custom color for each cell type 
custom_hex = [
    "#1F77B4",  # ClubEpi
    "#AEC7E8",  # CiliatedEpi
    "#FF7F0E",  # AT2
    "#FFBB78",  # AT1
    "#2CA02C",  # AdventitialFib
    "#98DF8A",  # AlveolarFib
    "#D62728",  # SMC
    "#FF9896",  # SCMF
    "#9467BD",  # Pericyte
    "#8C564B",  # Artery
    "#C49C94",  # Vein
    "#E377C2",  # Capillary
    "#F7B6D2",  # LymphaticEndo
    "#7F7F7F",  # AlveolarMac
    "#C7C7C7",  # InterstitialMac
    "#BCBD22",  # InflamMac
    "#DBDB8D",  # PatrolMac
    "#17BECF",  # DC
    "#9EDAE5",  # pDC
    "#393B79",  # Neutrophil
    "#637939",  # Tcell (CD4)
    "#B5CF6B",  # Tcell (CD8)
    "#843C39"   # Bcell
]

# load information per region 
region = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J']
all_IFN_num = pd.DataFrame([], columns = ["ClubEpi",
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
for i in range(len(region)):
    region_IFN_frac = pd.read_csv("region"+region[i]+"_ifn_info0302.csv", header = 0, index_col = 0)
    all_IFN_num = pd.concat([all_IFN_num, region_IFN_frac.T])

all_IFN_num = pd.DataFrame(all_IFN_num)

all_plaque_comp = pd.DataFrame([], columns = ["ClubEpi",
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
for i in range(len(region)):
    region_comp = pd.read_csv("region"+region[i]+"_plaque_info0302.csv", header = 0, index_col = 0)
    all_plaque_comp = pd.concat([all_plaque_comp, region_comp])
# all_plaque_comp.to_csv('allplaques.csv')
all_plaque_comp = pd.DataFrame(all_plaque_comp)

# Calculate the sum of each row
row_sums = all_plaque_comp.sum(axis=1)

# Divide each row by its corresponding sum
plaque_freq = all_plaque_comp.div(row_sums, axis=0)
plaque_freq.to_csv('AllRegions_Plaque_freq.csv')

# sort foci by cell type composition 
# look at epithelial cell foci first 
airway = plaque_freq[(plaque_freq['CiliatedEpi'] > 0.00) | (plaque_freq['ClubEpi'] > 0.00)]
airway['sum'] = airway['CiliatedEpi'] + airway['ClubEpi']
airway_sort = airway.sort_values(by = 'sum', ascending = False)
airway_sort.drop(columns = ['sum'], inplace = True)
airway_sort.plot(kind='bar', stacked=True, color = custom_hex,align='edge', width=1.0)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# sort IFN data frame by cell type composition
mask = (plaque_freq['CiliatedEpi'] > 0.00) | (plaque_freq['ClubEpi'] > 0) 

airway = plaque_freq[mask].reset_index(drop=True).copy()
airway_ifn = all_IFN_num[mask].reset_index(drop=True).copy()
airway_sort = airway.assign(
    sum=airway['CiliatedEpi'] + airway['ClubEpi']
).sort_values(by='sum', ascending=False)

ordered_indices = airway_sort.index.tolist()
# print(ordered_indices)
airway_ifn_sort = airway_ifn.loc[ordered_indices,:]

# look at alveolar foci
alveolar = plaque_freq[(plaque_freq['CiliatedEpi'] + plaque_freq['ClubEpi']) ==0]
alveolar['sum'] = alveolar['AT1'] + alveolar['AT2']
alveolar_sort = alveolar.sort_values(by = 'sum', ascending = False)
alveolar_sort.drop(columns = ['sum'], inplace = True)
alveolar_sort.plot(kind='bar', stacked=True, color = custom_hex,align='edge', width=1.0)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# look at IFN info
mask = (plaque_freq['CiliatedEpi'] + plaque_freq['ClubEpi']) <= 0.01

alveolar = plaque_freq[mask].reset_index(drop=True).copy()
alveolar_ifn = all_IFN_num[mask].reset_index(drop=True).copy()
alveolar_sort = alveolar.assign(
    sum=alveolar['AT1'] + alveolar['AT2']
).sort_values(by='sum', ascending=False)

ordered_indices = alveolar_sort.index.tolist()
# print(ordered_indices)
alveolar_ifn_sort = alveolar_ifn.loc[ordered_indices,:]

# calculate fraction of IFN+ cells for each cell type per foci
frac = np.float64(all_IFN_num.values)/ np.float64(all_plaque_comp.values)
frac_IFN = pd.DataFrame(frac, index=all_plaque_comp.index, columns=all_plaque_comp.columns)
frac_IFN.replace([np.inf, -np.inf], np.nan, inplace=True) #remove nan due to division by 0
mean_frac = frac_IFN.mean(axis=0)
std_frac = frac_IFN.std(axis = 0)/np.sqrt(100)
plt.errorbar(mean_frac.keys().tolist(), mean_frac, std_frac, linestyle='None', marker='o')
plt.ylim((0,1))



