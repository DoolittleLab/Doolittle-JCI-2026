##Download supplementary .xlsx files and place in working directory:
#Apoptosis_resistance_geneset.xlsx
#CORDENONSI_YAP_CONSERVED_SIGNATURE.v2024.1_mouse.xlsx

##Install and load all packages
install.packages(c("Seurat", "SCpubr", "ggplot2", "scCustomize", "dplyr", "ggpubr"))
BiocManager::install(c("escape", "GSEABase", "readxl", "DESeq2",
                       "RColorBrewer", "zellkonverter", "clusterProfiler"))
devtools::install_github("cellgeni/schard")
library(Seurat)
library(schard)
library(SCpubr)
library(ggplot2)
library(scCustomize)
library(GSEABase)
library(escape)
library(readxl)
library(DESeq2)
library(RColorBrewer)
library(dplyr)
library(zellkonverter)
library(ggpubr)
library(clusterProfiler)


##Create Seurat Object for Tabula Muris Old samples
options(timeout = 600)
download.file("https://figshare.com/ndownloader/files/23873090", destfile = "tabula_droplet_marrow.h5ad", mode = "wb")
Tabula = schard::h5ad2seurat('tabula_droplet_marrow.h5ad')
# Can also manually download .h5ad file from 
# https://figshare.com/ndownloader/files/23873090 and place in working directory
Tabula = schard::h5ad2seurat('tabula-muris-senis-droplet-processed-official-annotations-Marrow.h5ad')

#Continue here:
Tabula <- SetIdent(Tabula, value = "age")
Tabula_old <- subset(Tabula, idents = c("24m", "30m"))
Tabula_old <- SetIdent(Tabula_old, value = "cell_ontology_class")


##Fig 1a - Tabula Muris DimPlot
Fig1a <- SCpubr::do_DimPlot(Tabula_old, label = TRUE, label.fill = NULL, 
                            label.color = "white") & NoLegend()
ggsave(filename = "Fig1a.svg", plot = Fig1a, width = 7, height = 6)


##Fig 1b - SenMayo gene enrichment
url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-32552-1/MediaObjects/41467_2022_32552_MOESM4_ESM.xlsx"
download.file(url, destfile = "SenMayo.xlsx", mode = "wb")
SenMayo <- read_excel("SenMayo.xlsx", sheet = "mouse")
View(SenMayo)
gene.sets <- list(Sen_Mayo = SenMayo$`Gene(murine)`)
enrichment.score_SenMayo <- escape.matrix(Tabula_old, 
                                       gene.sets = gene.sets, 
                                       groups = 5000, 
                                       min.size = 0)
Tabula_old_SEN <- AddMetaData(Tabula_old, enrichment.score_SenMayo, col.name = "SenMayo")
Fig1b <- FeaturePlot(Tabula_old_SEN, features = "SenMayo", pt.size = 1, order = TRUE, min.cutoff = "q10", 
                     max.cutoff = "q90")+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(filename = "Fig1b.svg", plot = Fig1b, width = 7, height = 6)


##Fig 1c & Extended Data Fig 1a - Myeloid expression plots
pal <- viridis(n=10, option = "F", direction = -1)
Fig1c <- scCustomize::FeaturePlot_scCustom(Tabula_old_SEN, features =c("Itgam", "Cd14", "Cd68", "Elane"), 
                                           colors_use = pal) & NoAxes()
ggsave(filename = "Fig1c.svg", plot = Fig1c, width = 7, height = 6)


##Fig 5i-j - GSEA plots
#Run DESeq2 on raw data
#BMSCs
BMSC_metadata <- read_excel("BMSC_metadata.xlsx")
BMSC_FC_Filtered <- read_excel("BMSC_FC_Filtered.xlsx")
coldata <- as.matrix(BMSC_metadata[, -1])
rownames(coldata) <- BMSC_metadata$id
cts <- as.matrix(BMSC_FC_Filtered[, -1])
rownames(cts) <- BMSC_FC_Filtered$Gene_Name
head(cts, 2)
all(rownames(coldata) == colnames(cts))
#TRUE
BMSC_dds <- DESeqDataSetFromMatrix(countData = cts,
                                   colData = coldata,
                                   design = ~ condition)
BMSC_dds
BMSC_dds$condition <- relevel(BMSC_dds$condition, ref = "Ctrl")
BMSC_dds <- DESeq(BMSC_dds)
BMSC_res <- results(BMSC_dds)
BMSC_res <- BMSC_res[order(BMSC_res$stat),]
head(BMSC_res)
BMSC <- BMSC_res
BMSC <- BMSC[order(-BMSC$stat),]
head(BMSC)
BMSC_geneList <- BMSC$stat
names(BMSC_geneList) <- rownames(BMSC)
BMSC_geneList

#Macrophages
Mac_metadata <- read_excel("Mac_metadata.xlsx")
Mac_FC_Filtered <- read_excel("Macrophage_FC_Filtered.xlsx")
coldata <- as.matrix(Mac_metadata[, -1])
rownames(coldata) <- Mac_metadata$id
cts <- as.matrix(Mac_FC_Filtered[, -1])
rownames(cts) <- Mac_FC_Filtered$Gene_Name
head(cts, 2)
all(rownames(coldata) == colnames(cts))
#TRUE
Mac_dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ condition)
Mac_dds
Mac_dds$condition <- relevel(Mac_dds$condition, ref = "Ctrl")
Mac_dds <- DESeq(Mac_dds)
Mac_res <- results(Mac_dds)
Mac_res <- Mac_res[order(Mac_res$stat),]
head(Mac_res)
Mac <- Mac_res
Mac <- Mac[order(-Mac$stat),]
head(Mac)
Mac_geneList <- Mac$stat
names(Mac_geneList) <- rownames(Mac)
Mac_geneList

#Senmayo GSEA
#BMSCs
url <- "https://www.gsea-msigdb.org/gsea/msigdb/mouse/download_geneset.jsp?geneSetName=SAUL_SEN_MAYO&fileType=gmt"
download.file(url, destfile = "SenMayo.gmt")
sen <- read.gmt("SenMayo.gmt")
BMSC_SenGSEA <- GSEA(geneList = BMSC_geneList, TERM2GENE = sen)
a<- gseaplot(BMSC_SenGSEA, geneSetID = 1, by="runningScore")
ggsave(filename = "BMSC_SenMayo.svg", plot = a, height = 3, width = 5)
#Macrophages
Mac_SenGSEA <- GSEA(geneList = Mac_geneList, TERM2GENE = sen)
#no term enriched under specific pvalueCutoff...
Mac_SenGSEA <- GSEA(geneList = Mac_geneList, TERM2GENE = sen, pvalueCutoff = 1)
a<- gseaplot(Mac_SenGSEA, geneSetID = 1, by="runningScore")
ggsave(filename = "Mac_SenMayo.svg", plot = a, height = 3, width = 5)

#Anti-apoptosis GSEA
#BMSCs
anti_apop_T <- read_xlsx("Apoptosis_resistance_geneset.xlsx")
BMSC_AntiApop <- GSEA(geneList = BMSC_geneList, TERM2GENE = anti_apop_T)
a<- gseaplot(BMSC_AntiApop, geneSetID = 1, by="runningScore")
ggsave(filename = "BMSC_AntiApop.svg", plot = a, height = 3, width = 5)
#Macrophages
Mac_AntiApop <- GSEA(geneList = Mac_geneList, TERM2GENE = anti_apop_T)
#no term enriched under specific pvalueCutoff...
Mac_AntiApop <- GSEA(geneList = Mac_geneList, TERM2GENE = anti_apop_T, pvalueCutoff = 1)
a<- gseaplot(Mac_AntiApop, geneSetID = 1, by="runningScore")
ggsave(filename = "Mac_AntiApop.svg", plot = a, height = 3, width = 5)

#YAP/TAZ GSEA
YAP <- read_xlsx("CORDENONSI_YAP_CONSERVED_SIGNATURE.v2024.1_mouse.xlsx")
BMSC_YAP <- GSEA(geneList = BMSC_geneList, TERM2GENE = YAP)
a<- gseaplot(BMSC_YAP, geneSetID = 1, by="runningScore")
ggsave(filename = "BMSC_YAP.svg", plot = a, height = 3, width = 5)
#Macrophages
Mac_YAP <- GSEA(geneList = Mac_geneList, TERM2GENE = YAP)
#no term enriched under specific pvalueCutoff...
Mac_YAP <- GSEA(geneList = Mac_geneList, TERM2GENE = YAP, pvalueCutoff = 1)
a<- gseaplot(Mac_YAP, geneSetID = 1, by="runningScore")
ggsave(filename = "Mac_YAP.svg", plot = a, height = 3, width = 5)


##Extended Data Fig 1a - Cdkn2a+SenMayo+ plot
EFig1a <- SCpubr::do_NebulosaPlot(Tabula_old_SEN, features = c("Cdkn2a", "SenMayo"), 
                                  joint = TRUE, return_only_joint = TRUE, 
                                  use_viridis = TRUE, viridis.palette = "F", 
                                  viridis.direction = -1)
ggsave(filename = "EFig1a.svg", plot = EFig1a, width = 7, height = 7)


##Extended Data Fig 1b-c - Cyclin-dependent kinase inhibitor Plots
p16 <- SCpubr::do_NebulosaPlot(Tabula_old, features = c("Cdkn2a"), 
                               viridis.palette = "F", 
                               use_viridis = TRUE, viridis.direction = -1)
p21 <- SCpubr::do_NebulosaPlot(Tabula_old, features = c("Cdkn1a"), 
                               viridis.palette = "F", 
                               use_viridis = TRUE, viridis.direction = -1)
p53 <- SCpubr::do_NebulosaPlot(Tabula_old, features = c("Trp53"), 
                               viridis.palette = "F", 
                               use_viridis = TRUE, viridis.direction = -1)
ggsave(filename = "Cdkn2a.svg", plot = p16, width = 7, height = 7)
ggsave(filename = "Cdkn1a.svg", plot = p21, width = 7, height = 7)
ggsave(filename = "Trp53.svg", plot = p53, width = 7, height = 7)


##Extended Data Fig 1d - Cdkn2a+ cells across clusters
p16_pos <- subset(Tabula_old_SEN, subset = Cdkn2a > 0)
p16_pos<- SetIdent(p16_pos, value = "cell_ontology_class")
cell_counts <- p16_pos@meta.data %>%
  group_by(cluster = Idents(p16_pos)) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  mutate(total_count = sum(cell_count), proportion = cell_count / total_count)
set3_palette <- brewer.pal(12, "Set3") 
extended_palette <- colorRampPalette(set3_palette)(18)
cell_counts <- cell_counts %>%
  mutate(cluster = recode(cluster, "megakaryocyte-erythroid progenitor cell" = 'MEP', ))
EFig1d <- ggplot(cell_counts, aes(x = as.factor(cluster), y = proportion, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", width = 1, colour="black", linewidth = 0.25) +
  scale_fill_manual(values = extended_palette, name = "Cluster") +
  labs(title = "Proportion of Cells per Cluster (Cdkn2a > 0)", 
       y = "Proportion of Cells", x = "Cluster") +
  theme_classic() +
  theme(panel.grid = element_blank(), legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=11, colour = "black")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave(filename = "EFig1d.svg", plot = EFig1d, width = 5, height = 5)


##Extended Data Fig 1e - SenePy enrichment analysis
#export  Seurat object as a h5ad file for Python input
sce <- as.SingleCellExperiment(Tabula_old_SEN)
writeH5AD(sce, file = 'PATH/TO/OUTPUT/Tabula_old_SEN.h5ad')

#perform SenePy enrichment in Python and export as 'adata_obs.csv'
obs_data <- read.csv('adata_obs.csv', row.names = 1)
#Check data structure - this should look identical to your seurat object's meta.data 
#Just with a newly added column
head(obs_data)
#Merge the obj_data into your seurat object
#If successful, it should print "Universal score added succesfully"
#We will troubleshoot if that is not the case 
# change <data> to your seurat object
Tabula_old_senepy <- Tabula_old_SEN
if (all(rownames(Tabula_old_senepy@meta.data) %in% rownames(obs_data))) {
  Tabula_old_senepy[['senepy_score']] <- obs_data[rownames(Tabula_old_senepy@meta.data), 'senepy_score']
  print("Senepy score added successfully.")
} else {
  print("Error: Not all cell identifiers from the Seurat object are present in the CSV data.")
}
head(Tabula_old_senepy@meta.data)
EFig1e <- SCpubr::do_FeaturePlot(Tabula_old_senepy, features = "senepy_score", use_viridis = TRUE, viridis.palette = "F", viridis.direction = -1)
ggsave(filename = "EFig1e.svg", plot = EFig1e, width = 7, height = 7)

##Extended Data 1f - Senepy v SenMayo
EFig1f <- VlnPlot_scCustom(Tabula_old_senepy, features = c("senepy_score", "SenMayo"), 
                 stack = TRUE, colors_use = c("orange2", "mediumslateblue")) & 
  NoLegend() + theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 270, vjust = 0.2))
ggsave(filename = "EFig1f.svg", plot = EFig1f, width = 7, height = 7)


##Extended Data Fig 4a
Tabula_mid <- subset(Tabula, idents = "18m")
gene.sets <- list(Sen_Mayo = SenMayo$`Gene(murine)`)
enrichment.score_SenMayo <- escape.matrix(Tabula_mid, 
                                          gene.sets = gene.sets, 
                                          groups = 5000, 
                                          min.size = 0)
Tabula_mid_SEN <- AddMetaData(Tabula_mid, enrichment.score_SenMayo, col.name = "SenMayo")
Tabula_mid_SEN_Cd14 <- subset(Tabula_mid_SEN, Cd14 > 0.01)
a<- VlnPlot(Tabula_mid_SEN_Cd14, features = "SenMayo", split.by = "sex", 
            cols = c("antiquewhite", "indianred")) + stat_compare_means(label = "p.format")
ggsave(filename = "male_female_Cd14_SenMayo.svg", plot = a, height = 3, width = 3)


##Extended Data Fig 4b
Tabula_mid_SEN <- SetIdent(Tabula_mid_SEN, value = "sex")
Tabula_mid_SEN_Cd14 <- SetIdent(Tabula_mid_SEN_Cd14, value = "sex")
a <- VlnPlot(Tabula_mid_SEN, features = c("Esr1", "Esr2", "Ar"), split.by = "sex")
b <- VlnPlot(Tabula_mid_SEN_Cd14, features = c("Esr1", "Esr2", "Ar"), split.by = "sex")
a/b


##Extended Data Fig 4c
#Download .gmt files from https://www.gsea-msigdb.org/gsea/msigdb
WP_url <- "https://www.gsea-msigdb.org/gsea/msigdb/mouse/download_geneset.jsp?geneSetName=WP_ESTROGEN_SIGNALING&fileType=gmt"
download.file(WP_url, destfile = "WP_EST.gmt")
GOBP_url <- "https://www.gsea-msigdb.org/gsea/msigdb/mouse/download_geneset.jsp?geneSetName=GOBP_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY&fileType=gmt"
download.file(GOBP_url, destfile = "GOBP_EST.gmt")
GOBP_POS_url <- "https://www.gsea-msigdb.org/gsea/msigdb/mouse/download_geneset.jsp?geneSetName=GOBP_POSITIVE_REGULATION_OF_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY&fileType=gmt"
download.file(GOBP_POS_url, destfile = "GOBP_POS_EST.gmt")
#WP
WP_EST <- read.gmt("WP_EST.gmt")
WP_EST <- read.gmt("WP_ESTROGEN_SIGNALING.v2024.1.Mm.gmt")
WPgene.set <- list(WP = WP_EST$gene)
enrichment.score_WP <- escape.matrix(Tabula_mid_SEN_Cd14, 
                                     gene.sets = WPgene.set, 
                                     groups = 5000, 
                                     min.size = 0)
Tabula_mid_SEN_Cd14_EST <- AddMetaData(Tabula_mid_SEN_Cd14, enrichment.score_WP, col.name = "WP")
#GOBP
GOBP_EST <- read.gmt("GOBP_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY.v2024.1.Mm.gmt")
GOBPgene.set <- list(GOBP = GOBP_EST$gene)
enrichment.score_GOBP <- escape.matrix(Tabula_mid_SEN_Cd14, 
                                       gene.sets = GOBPgene.set, 
                                       groups = 5000, 
                                       min.size = 0)
Tabula_mid_SEN_Cd14_EST <- AddMetaData(Tabula_mid_SEN_Cd14_EST, enrichment.score_GOBP, col.name = "GOBP")
#GOBP_POS
GOBP_POS_EST <- read.gmt("GOBP_POSITIVE_REGULATION_OF_INTRACELLULAR_ESTROGEN_RECEPTOR_SIGNALING_PATHWAY.v2024.1.Mm.gmt")
GOBP_POSgene.set <- list(GOBP_POS = GOBP_POS_EST$gene)
enrichment.score_GOBP_POS <- escape.matrix(Tabula_mid_SEN_Cd14, 
                                           gene.sets = GOBP_POSgene.set, 
                                           groups = 5000, 
                                           min.size = 0)
Tabula_mid_SEN_Cd14_EST <- AddMetaData(Tabula_mid_SEN_Cd14_EST, enrichment.score_GOBP_POS, col.name = "GOBP_POS")
#Plot
Tabula_mid_SEN_Cd14_EST <- SetIdent(Tabula_mid_SEN_Cd14_EST, value="sex")
a<- VlnPlot(Tabula_mid_SEN_Cd14_EST, features = "WP", cols = c("antiquewhite", "indianred")) + 
  stat_compare_means(label = "p.format") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  labs(y= "Enrichment Score") + coord_cartesian(clip = "off")
b<- VlnPlot(Tabula_mid_SEN_Cd14_EST, features = "GOBP", cols = c("antiquewhite", "indianred")) + 
  stat_compare_means(label = "p.format") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + coord_cartesian(clip = "off")
c<- VlnPlot(Tabula_mid_SEN_Cd14_EST, features = "GOBP_POS", cols = c("antiquewhite", "indianred")) + 
  stat_compare_means(label = "p.format") + 
  theme(axis.title.x = element_blank()) + coord_cartesian(clip = "off")
a|b|c
ggsave(filename = "EFig4c.svg", plot = a|b|c, height = 4, width = 12)


##Extended Data Fig 4d
Tabula_mid_Cd14SENmayo <- subset(Tabula_mid_SEN_Cd14, SenMayo > 500) 
#WP
enrichment.score_WP <- escape.matrix(Tabula_mid_Cd14SENmayo, 
                                     gene.sets = WPgene.set, 
                                     groups = 5000, 
                                     min.size = 0)
Tabula_mid_Cd14SENmayo_EST <- AddMetaData(Tabula_mid_Cd14SENmayo, enrichment.score_WP, col.name = "WP")
#GOBP
enrichment.score_GOBP <- escape.matrix(Tabula_mid_Cd14SENmayo, 
                                       gene.sets = GOBPgene.set, 
                                       groups = 5000, 
                                       min.size = 0)
Tabula_mid_Cd14SENmayo_EST <- AddMetaData(Tabula_mid_Cd14SENmayo_EST, enrichment.score_GOBP, col.name = "GOBP")
#GOBP_POS
enrichment.score_GOBP_POS <- escape.matrix(Tabula_mid_Cd14SENmayo, 
                                           gene.sets = GOBP_POSgene.set, 
                                           groups = 5000, 
                                           min.size = 0)
Tabula_mid_Cd14SENmayo_EST <- AddMetaData(Tabula_mid_Cd14SENmayo_EST, enrichment.score_GOBP_POS, col.name = "GOBP_POS")
#Plot
Tabula_mid_Cd14SENmayo_EST <- SetIdent(Tabula_mid_Cd14SENmayo_EST, value="sex")
a<- VlnPlot(Tabula_mid_Cd14SENmayo_EST, features = "WP", cols = c("antiquewhite", "indianred")) + 
  stat_compare_means(label = "p.format") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  labs(y= "Enrichment Score") + coord_cartesian(clip = "off")
b<- VlnPlot(Tabula_mid_Cd14SENmayo_EST, features = "GOBP", cols = c("antiquewhite", "indianred")) + 
  stat_compare_means(label = "p.format") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + coord_cartesian(clip = "off")
c<- VlnPlot(Tabula_mid_Cd14SENmayo_EST, features = "GOBP_POS", cols = c("antiquewhite", "indianred")) + 
  stat_compare_means(label = "p.format") + 
  theme(axis.title.x = element_blank()) + coord_cartesian(clip = "off")
a|b|c
ggsave(filename = "EFig4d.svg", plot = a|b|c, height = 4, width = 12)


##Extended Data Fig 4e
#WP
enrichment.score_WP <- escape.matrix(Tabula_mid_SEN, 
                                       gene.sets = WPgene.set, 
                                       groups = 5000, 
                                       min.size = 0)
Tabula_mid_SEN_EST <- AddMetaData(Tabula_mid_SEN, enrichment.score_WP, col.name = "WP")
#GOBP
enrichment.score_GOBP <- escape.matrix(Tabula_mid_SEN, 
                                      gene.sets = GOBPgene.set, 
                                      groups = 5000, 
                                      min.size = 0)
Tabula_mid_SEN_EST <- AddMetaData(Tabula_mid_SEN_EST, enrichment.score_GOBP, col.name = "GOBP")
#GOBP_POS
enrichment.score_GOBP_POS <- escape.matrix(Tabula_mid_SEN, 
                                       gene.sets = GOBP_POSgene.set, 
                                       groups = 5000, 
                                       min.size = 0)
Tabula_mid_SEN_EST <- AddMetaData(Tabula_mid_SEN_EST, enrichment.score_GOBP_POS, col.name = "GOBP_POS")
#Plot
Tabula_mid_SEN_EST <- SetIdent(Tabula_mid_SEN_EST, value="cell_ontology_class")
a<- VlnPlot(Tabula_mid_SEN_EST, features = "WP", group.by = "cell_ontology_class", split.by = "sex", 
        idents = c("granulocyte", "granulocytopoietic cell", "macrophage", "monocyte", "promonocyte"), 
        cols = c("antiquewhite", "indianred")) + stat_compare_means(label = "p.format") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  labs(y= "Enrichment Score") + coord_cartesian(clip = "off")
b<- VlnPlot(Tabula_mid_SEN_EST, features = "GOBP", group.by = "cell_ontology_class", split.by = "sex", 
            idents = c("granulocyte", "granulocytopoietic cell", "macrophage", "monocyte", "promonocyte"), 
            cols = c("antiquewhite", "indianred")) + stat_compare_means(label = "p.format") + 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
  labs(y= "Enrichment Score") + coord_cartesian(clip = "off")
c<- VlnPlot(Tabula_mid_SEN_EST, features = "GOBP_POS", group.by = "cell_ontology_class", split.by = "sex", 
            idents = c("granulocyte", "granulocytopoietic cell", "macrophage", "monocyte", "promonocyte"), 
            cols = c("antiquewhite", "indianred")) + stat_compare_means(label = "p.format") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 14)) + 
  labs(y= "Enrichment Score") + coord_cartesian(clip = "off")
a/b/c
ggsave(filename = "EFig4e.svg", plot = a/b/c, height = 10, width = 9)

##Extended Data Fig 8k



##Extended Data Fig 8l
p16_pos <- subset(Tabula, subset = Cdkn2a > 0)
p16_pos <- SetIdent(p16_pos, value = "age")
cell_counts <- p16_pos@meta.data %>%
  group_by(cluster = Idents(p16_pos)) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  mutate(total_count = sum(cell_count), proportion = cell_count / total_count)
set3_palette <- brewer.pal(12, "Set3") 
extended_palette <- colorRampPalette(set3_palette)(18)
EFig8l <- ggplot(cell_counts, aes(x = as.factor(cluster), y = proportion, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", width = 1, colour="black", linewidth = 0.25) +
  scale_fill_manual(values = extended_palette, name = "Cluster") +
  labs(title = "Proportion of Cells per Age (Cdkn2a > 0)", 
       y = "Proportion of Cells", x = "Cluster") +
  theme_classic() +
  theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=11, colour = "black")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave(filename = "EFig8l.svg", plot = EFig8l, width = 5, height = 5)


p16_Cd14_pos <- subset(p16_pos, subset = Cd14 > 0)
p16_Cd14_pos<- SetIdent(p16_Cd14_pos, value = "age")
cell_counts <- p16_Cd14_pos@meta.data %>%
  group_by(cluster = Idents(p16_Cd14_pos)) %>%
  summarise(cell_count = n(), .groups = 'drop') %>%
  mutate(total_count = sum(cell_count), proportion = cell_count / total_count)
set3_palette <- brewer.pal(12, "Set3") 
extended_palette <- colorRampPalette(set3_palette)(18)
EFig8m <- ggplot(cell_counts, aes(x = as.factor(cluster), y = proportion, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", width = 1, colour="black", linewidth = 0.25) +
  scale_fill_manual(values = extended_palette, name = "Cluster") +
  labs(title = "Proportion of Cells per Age (Cdkn2a & Cd14 > 0)", 
       y = "Proportion of Cells", x = "Cluster") +
  theme_classic() +
  theme(panel.grid = element_blank(), legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text = element_text(size=11, colour = "black")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
ggsave(filename = "EFig8m.svg", plot = EFig8m, width = 5, height = 5)

