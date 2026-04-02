input to analysis code is gene expression/meta data/cell type infomation split by region 
0. process/QC cellxgene matrix and identify cell types based on known markers
1. run per_focus.py to generate infected cell type and IFN+ cell type composition per focus per region
2. run pooled_analysis.py to generate probability of IFN+ cell for each cell type per focus 
3. run visualize.m to plot 