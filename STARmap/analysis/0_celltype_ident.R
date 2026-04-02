# cell type identification and validation 
library(Matrix)
library("reticulate")
library(Seurat)
library(dplyr)
library("stringr")


#---- read in gene starression matrix
setwd('~/STARmap/36hpi_250genes_sum/')

star.data <- Read10X(data.dir = "~/STARmap/36hpi_250genes_sum/all_regions", gene.column = 1)
marker <- read.table("marker_long.txt", header = 1)
epi_marker <- read.tbale("epi_marker.txt", header = 1)
cell <- read.table("~/STARmap/36hpi_250genes_sum/all_regions/barcodes.tsv", sep = "\t", row.names = 1)
colnames(cell) <- c( "x", "y", "z", "volume", "np", "np_nuc", "np_cyto", "pctnuc", "pctcyto", "region")

#----pre-processing

star <- CreateSeuratObject(counts = star.data, project = "all")
FeatureScatter(star, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(star, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
star <- AddMetaData(star, cell)
VlnPlot(star, features = c("volume"))
star <- subset(star, subset = nFeature_RNA >= 3 & 
                volume > 5e3 & volume <= 100000 & 
                nCount_RNA > 2 & nCount_RNA <= 500)

#---- normalization 
star <- SCTransform(star)


#---- dimension reduction - only run on marker genes 
star <- RunPCA(star, features = marker$marker)
ElbowPlot(star, ndims = 50)


#---- find neighbors
star <- FindNeighbors(star, dims = 1:30 )
star <- FindClusters(star, resolution = 0.5, algorithm = 4)

star <- RunUMAP(star, dims = 1:30)
star <- RunTSNE(star, dims = 1:30, check_duplicates =FALSE)

DimPlot(star, reduction = "umap", label=TRUE, repel = TRUE)


#---- find marker genes for each cluster for 1st level clustering 
star.markers <- FindAllMarkers(star, features = marker$marker, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top_marker <- star.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)


## for ambiguous clusters, run FindSubCluster for more pure population
# cluster 1
star = FindSubCluster(star, "1", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
star <- SetIdent(star, value = star@meta.data$lung)
star.markers <- FindAllMarkers(star, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- star.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 2
star = FindSubCluster(star, "2", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
star <- SetIdent(star, value = star@meta.data$lung)
star.markers <- FindAllMarkers(star, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- star.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 3
star = FindSubCluster(star, "3", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
star <- SetIdent(star, value = star@meta.data$lung)
star.markers <- FindAllMarkers(star, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- star.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 4
star = FindSubCluster(star, "4", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
star <- SetIdent(star, value = star@meta.data$lung)
star.markers <- FindAllMarkers(star, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- star.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 5
star = FindSubCluster(star, "5", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
star <- SetIdent(star, value = star@meta.data$lung)
star.markers <- FindAllMarkers(star, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- star.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 2
star = FindSubCluster(star, "7", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
star <- SetIdent(star, value = star@meta.data$lung)
star.markers <- FindAllMarkers(star, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- star.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

# cluster 2
star = FindSubCluster(star, "8", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
star <- SetIdent(star, value = star@meta.data$lung)
star.markers <- FindAllMarkers(star, features = marker$marker, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
top <- star.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

star_label <- RenameIdents(object = star,
                           "1_1" = "AT2", 
                           "1_3" = "AT2", 
                           "1_2_1" = "Tcell (CD4)",
                           "1_2_3" = "Tcell (CD4)",
                           "1_4_3" = "AlveolarFib",
                           "1_4_1" = "AlveolarFib",
                           "1_5_2" = "AlveolarFib",
                           "1_5_1" = "LymphaticEndo",
                           "1_5_4" = "pDC",
                           "1_6_2" = "Neutrophil",
                           "2_1" = "Capillary",
                           "2_2" = "InflamMac",
                           "2_5" = "DC",
                           "2_6" = "pDC",
                           "2_4" = "Artery",
                           "2_3" = "Vein",
                           "3_4" = "AdventitialFib", 
                           "3_3" = "AdventitialFib",
                           "3_5" = "SMC",
                           "3_6" = "SCMF",
                           "3_1" = "SMC",
                           "3_2" = "SMC",
                           "5_1" = "Pericyte",
                           "5_3" = "Pericyte",
                           "7_2" = "PatrolMac",
                           "7_4" = "PatrolMac",
                           "7_3" = "Tcell (CD8)",
                           "8_2" = "InterstitialMac",
                           "8_3" = "InterstitialMac",
                           "8_1" = "AlveolarMac",
                           "4_2" = "AirwayEpi",
                           "4_1_1" = "AirwayEpi",
                           "4_1_2" = "AirwayEpi",
                           "4_1_3" = "AirwayEpi",
                           "9" = "Bcell",
                           "6" = "AT1")


#---- subcluster airway epithelial cell cluster
epi_star = subset(star_label, idents = "AirwayEpi")
epi_star <- SCTransform(epi_star)
epi_star <- RunPCA(epi_star, features = epi_marker$marker)

#---- find clusters within airway epithelia
epi_star <- FindNeighbors(epi_star, dims = 1:7 )
epi_star <- FindClusters(epi_star, resolution = 0.3, algorithm = 4)
DoHeatmap(epi_star, epi_marker$V1)
epi_star <- RunUMAP(epi_star, dims = 1:7)
DimPlot(epi_star, reduction = "umap", label=TRUE, repel = TRUE)


#---- determine markers for each cluster
epi_star.markers <- FindAllMarkers(epi_star, features = marker$marker)
FeaturePlot(epi, features = epi_marker$marker)
epi = FindSubCluster(epi, "1", graph.name = "SCT_snn", resolution = 0.2, subcluster.name = "lung", algorithm = 4)
epi <- SetIdent(epi, value = epi@meta.data$lung)
epi_label <- RenameIdents(object = epi,
                          "1_1" = "CLubEpi", 
                          "1_3" = "ClubEpi", 
                          "1_2" = "CiliatedEpi",
                          "2" = "CiliatedEpi",
                          "3" = "Unidentified",
                          "4" = "Unidentified",
                          "5" = "ClubEpi")
#---- map cells back to the main obj 
target_cells <- WhichCells(star_label, idents = 'AirwayEpi')
Idents(star_label, cells = target_cells) = epi_label

my_levels <- c("ClubEpi",
               "CiliatedEpi",
               "AT2",
               "AT1",
               'AdventitialFib',
               'AlveolarFib',
               'SMC',
               'SCMF',
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
               "Bcell")

# Relevel object@ident
star_label@active.ident <- factor(x = star_label@active.ident, levels = my_levels)
DotPlot(star_label, features = c("Scgb1a1", "Cdhr3", "Rsph1", "Rtkn2", "Lama3", 
                                 "Abca3", "Lamp3", "Dcn", "Col14a1", "Itga8", "Col13a1", 
                                "Hhip", "Pdgfrb", "Notch3", "Myh11", "Acta2", 
                                "Efnb2", "Gja5", "Ephb4", "Sema3c", "Plin2", 
                                "Mrc1", "C1qa","Ly6c2", "Spn", "Itgae", "Tcf4", 
                                "Csf3r", "Ccr7", "Rorc", "Cd19"))
# write output files
dat <- star_label
region = subset(dat, subset = region == 'J')
celltype <- Idents(region)
meta <- region@meta.data
write.table(as.matrix(GetAssayData(object = region, slot = "data")), 
            'regionJ_sct_int.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

write.table(celltype, 'regionJ_celltypes.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

write.table(meta, 'regionJ_metadata.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)


#----average cluster average 
sc = readRDS('231006_day0_idents.RDS') # Seurat object for scRNA-seq
sc.average <- Aggregatestarression(sc)
write.table(sc.average$RNA,
            '251111_sc_cluster_average.txt', sep = '\t', row.names = T, col.names = T, quote = F)
cluster.averages <- Aggregatestarression(star)
write.table(cluster.average$RNA,
            '251111_spatial_cluster_average.txt', sep = '\t', row.names = T, col.names = T, quote = F)

#---- label transfer
# find anchor (scRNA-Seq as ref)
sc = readRDS('231006_day0_idents.RDS') # Seurat object for scRNA-seq
common_genes <- intersect(rownames(star_label), rownames(sc))

anchors <- FindTransferAnchors(reference = sc, query = star_label, 
                               dims = 1:30, 
                               features = common_genes)

# label transfer
predictions_cluster <- TransferData(anchorset = anchors, 
                                    refdata = Idents(sc_adult), 
                                    dims = 1:30)
star_label@meta.data$predicted_cluster = predictions_cluster$predicted.id
star_label@meta.data$predicted_cluster_score = predictions_cluster$prediction.score.max



# subset cells - uninfected with no viral genes and ISGs
uninfect = subset(exp_label, NP.gRNA == 0 & NS1.gRNA == 0 & Isg15 == 0 
                  & Ifit1 ==  0 & Ifitm3 ==  0 & Gbp4 ==  0 & Irgm1 ==  0 
                  & Isg20 == 0 & Bst2 == 0)




