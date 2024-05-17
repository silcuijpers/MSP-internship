#clean workspace TEST
rm(list=ls())
setwd("~/Documents/GitHub/MSP-internship/Thesis")
########################################################################################
## You can add the code below to install the R-packages if you don't have them yet
## this code is for packages that are originating from Bioconductor

# Install BiocManager from CRAN
install.packages("BiocManager")

# Load BiocManager
library(BiocManager)

# Install Bioconductor packages
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("EnhancedVolcano")
BiocManager::install("VennDiagram")
BiocManager::install("biomaRt")
BiocManager::install("Biostrings")
BiocManager::install("clusterProfiler")
BiocManager::install("GO.db")
BiocManager::install("HDO.db")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("reactome.db")
install.packages("writexl")
BiocManager::install("SPIA")

######################################################################################## 
#load installed packages
library(GEOquery)
library(limma)
library(EnhancedVolcano)
library(VennDiagram)
library(ggfortify)
library(Biostrings)
library(biomaRt)
library(dplyr)
library(clusterProfiler)
library(org.HS.eg.db)
library(enrichplot)
library(ReactomePA)
library(writexl)
library(SPIA)
library(CePa)
library(netgsa)
library(tibble)
library(graphite)


#loading data from GEO
gset <- getGEO("GSE28358",GSEMatrix=TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL571", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#removing sample that misses baseline measurement
#row_to_remove <- which(gset$geo_accession == "GSM701152")
#gset <- gset[-row_to_remove, ]

# Make column names by using the sample names
# In the current dataset this is for example "GSM701076"
fvarLabels(gset) <- make.names(fvarLabels(gset))


# assigning each sample to one of the six groups
gsms <- "10101010101010101032323232323232323232544441054545454554545410324"
sml <- strsplit(gsms, split="")[[1]]

# assigning a groupname to each sample

gset$Group[gset$`intervention:ch1` == "olive oil" & gset$`timepoint:ch1` == "baseline"] <- "baseline"
gset$Group[gset$`intervention:ch1` == "olive oil" & gset$`timepoint:ch1` == "3 months"] <-"olive oil"
gset$Group[gset$`intervention:ch1` == "nuts" & gset$`timepoint:ch1` == "baseline"] <-"nuts baseline"
gset$Group[gset$`intervention:ch1` == "nuts" & gset$`timepoint:ch1` == "3 months"] <-"nuts"
gset$Group[gset$`intervention:ch1` == "low fat" & gset$`timepoint:ch1` == "baseline"] <-"low fat baseline"
gset$Group[gset$`intervention:ch1` == "low fat" & gset$`timepoint:ch1` == "3 months"] <-"low fat"

gs <- factor(sml)
groups <- make.names(c("olive oil","baseline","nuts","nuts baseline","low fat","low fat baseline"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

expres=exprs(gset)
expres
#making boxplot of raw expression data
 dev.new(width=3+ncol(gset)/6, height=5)
 png('downloads/plot/boxplotraw.png')
 ord <- order(gs)  # order samples by group
  palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
            "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
  par(mar=c(7,4,2,1))
  title <- paste ("rawdata")
  boxplot(expres[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
  legend("topright", groups, fill=palette(), bty="n")
  dev.off()

# plotting pca of raw data
texpres=t(expres)
pca_result <- prcomp(texpres, scale = TRUE)

pc1 <- pca_result$x[, 1]
pc2 <- pca_result$x[, 2]
pc3 <- pca_result$x[, 3]

pca.13 <- data.frame(x = pc1, y = pc3)
pca.12<- data.frame(x = pc1, y = pc2)
pca.23<- data.frame(x = pc2, y = pc3)

pca_model_13 <- prcomp(pca.13)
pca_model_12 <- prcomp(pca.12)
pca_model_23 <- prcomp(pca.23)

png('downloads/plot/pca_raw_13.png')
autoplot(pca_model_13, data = as.data.frame(gset), colour = "group", scale = TRUE, label = TRUE, label.size = 3) +
labs(title = "pca raw 13", x = "PC1", y = "PC3")
dev.off()
png('downloads/plot/pca_raw_12.png')
autoplot(pca_model_12, data = as.data.frame(gset), colour = "group", scale = TRUE, label = TRUE, label.size = 3) +
labs(title = "pca raw 12", x = "PC1", y = "PC2")
dev.off()
png('downloads/plot/pca_raw_23.png')
autoplot(pca_model_23, data = as.data.frame(gset), colour = "group", scale = TRUE, label = TRUE, label.size = 3) +
labs(title = "pca raw 23", x = "PC2", y = "PC3")
dev.off()


# clustering the raw expression data
sample_distraw = dist(texpres)
clusters <-hclust(sample_distraw, method = "complete")
png('downloads/plot/clustering_raw.png')
plot(clusters, col=group_colors[group_ids], labels = gset$geo_accession, label.size = 1, cex = 0.5,main = "")
title(main= "Clustering raw expression data")
dev.off()

# replacing values smaller then 0 with NaN and log transform each element
ex<-  exprs(gset)
ex[which(ex <= 0)]<- NaN
exprs(gset) <- log2(ex)

#making boxplot of log2 transformed expression data

dev.new(width=3+ncol(gset)/6, height=5)
png('downloads/plot/boxplotlog2.png')
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("log2 data")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topright", groups, fill=palette(), bty="n")
dev.off()

# plotting pca of log transformed data
t.ex= t(exprs(gset))

pca_result_log <- prcomp(t.ex, scale = TRUE)
pc1.log <- pca_result_log$x[, 1]
pc2.log <- pca_result_log$x[, 2]
pc3.log <- pca_result_log$x[, 3]

pca.13_log <- data.frame(x = pc1.log, y = pc3.log)
pca.12_log<- data.frame(x = pc1.log, y = pc2.log)
pca.23_log<- data.frame(x = pc2.log, y = pc3.log)

pca_model_13log <- prcomp(pca.13_log)
pca_model_12log <- prcomp(pca.12_log)
pca_model_23log <- prcomp(pca.23_log)

png('downloads/plot/pcalog_12.png')
autoplot(pca_model_12log, data = as.data.frame(gset), colour = "group",scale = TRUE, label = TRUE, label.size = 3) +
labs(title = "pca log 12", x = "PC1", y = "PC2")
dev.off()
png('downloads/plot/pcalog_13.png')
autoplot(pca_model_13log, data = as.data.frame(gset), colour = "group",scale = TRUE, label = TRUE, label.size = 3) +
labs(title = "pca log 13", x = "PC1", y = "PC3")
dev.off()
png('downloads/plot/pcalog_23.png')
autoplot(pca_model_23log, data = as.data.frame(gset), colour = "group",scale = TRUE, label = TRUE, label.size = 3) +
labs(title = "pca log 23", x = "PC2", y = "PC3")
dev.off()

#clustering the log transformed data
sample_dist = dist(t.ex)
clusters <-hclust(sample_dist, method = "complete")

png('downloads/plot/clustering_log.png')
plot(clusters, col = group_colors[as.numeric(factor(gset$group))], labels = gset$geo_accession, label.size = 1, cex = 0.5,main = "")
title(main= "Clustering log transformed data")

dev.off()

#calculate the number of rows with Na
rows_with_na <- sum(!complete.cases(ex))
print(rows_with_na)

#removing rows with missing values (NA)
gset <- gset[complete.cases(exprs(gset)), ]

#fitting a linear model to gene expression data
fit <- lmFit(gset, design) 

#making a contrast matrix of group1 olive oil and group2 baseline to calculate the log 2 fold change 
cts.olive <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix.olive<- makeContrasts(contrasts=cts.olive, levels=design)
fit.olive <- contrasts.fit(fit, cont.matrix.olive)
fit.olive <- eBayes(fit.olive, 0.01)# standard deviation of 0.01

cts.nuts <- c(paste(groups[3],"-",groups[4],sep=""))
cont.matrix.nuts<- makeContrasts(contrasts=cts.nuts, levels=design)
fit.nuts <- contrasts.fit(fit, cont.matrix.nuts)
fit.nuts <- eBayes(fit.nuts, 0.01)

cts.lowfat <- c(paste(groups[5],"-",groups[6],sep=""))
cont.matrix.lowfat<- makeContrasts(contrasts=cts.lowfat, levels=design)
fit.lowfat <- contrasts.fit(fit, cont.matrix.lowfat)
fit.lowfat <- eBayes(fit.lowfat, 0.01)

#making top tables for each group using benjamini-hochberg to adjust p-values
tT.olive<- topTable(fit.olive, adjust="BH", sort.by="B", number=Inf)
tT.nuts<-topTable(fit.nuts, adjust="BH", sort.by="B", number = Inf)
tT.lowfat<-topTable(fit.lowfat, adjust="BH", sort.by="B", number= Inf)

# writing top tables of top signigcant genes for each group
tT.olive <- subset(tT.olive, select=c("ID", "Gene.symbol", "Gene.ID", "logFC", "P.Value", "adj.P.Val", "B"))
write.table(tT.olive, file=stdout(), row.names=F, sep="\t")

tT.nuts <- subset(tT.nuts, select=c("ID", "Gene.symbol", "Gene.ID", "logFC", "P.Value", "adj.P.Val", "B"))
write.table(tT.nuts, file=stdout(), row.names=F, sep="\t")


tT.lowfat <- subset(tT.lowfat, select=c("ID", "Gene.symbol", "Gene.ID", "logFC", "P.Value", "adj.P.Val", "B"))
write.table(tT.lowfat, file=stdout(), row.names=F, sep="\t")

#removing gene symbols that are double based on the highest log FC

tT.olive <- tT.olive %>%
  group_by(Gene.symbol) %>%
  filter(logFC == max(logFC))

tT.nuts <- tT.nuts %>%
  group_by(Gene.symbol) %>%
  filter(logFC == max(logFC))

tT.lowfat <- tT.lowfat %>%
  group_by(Gene.symbol) %>%
  filter(logFC == max(logFC))

#plot for adjusted p-value distribution olive oil
png('downloads/plot/olive-adjustedpvalue.png')
hist(tT.olive$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "Olive oil: P-adj value distribution")
dev.off()

# Olive oil p-value distribution
png('downloads/plot/olive-pvalue.png')
hist(tT.olive$P.Value, col = "grey", border = "white", xlab = "P.Value",
     ylab = "Number of genes", main = "Olive oil: P.Value distribution")
dev.off()

# Nuts adjusted p-value
png('downloads/plot/nuts-adjustedpvalue.png')
hist(tT.nuts$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "Nuts: P-adj value distribution")
dev.off()

# Nuts p-value distribution
png('downloads/plot/nuts-pvalue.png')
hist(tT.nuts$P.Value, col = "grey", border = "white", xlab = "P.Value",
     ylab = "Number of genes", main = "Nuts: P.Value distribution")
dev.off()

# Lowfat adjusted p-value
png('downloads/plot/lowfat-adjustedpvalue.png')
hist(tT.lowfat$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "Low fat: P-adj value distribution")
dev.off()

# Lowfat p-value distribution
png('downloads/plot/lowfat-pvalue.png')
hist(tT.lowfat$P.Value, col = "grey", border = "white", xlab = "P.Value",
     ylab = "Number of genes", main = "Low fat: P.Value distribution")
dev.off()



#volcanoplot Olive oil
png('downloads/plot/volcanoplot_olive.png')
EnhancedVolcano(tT.olive, title = "Olive oil", lab = tT.olive$Gene.symbol, 
                labSize = 3, x = 'logFC', xlim = c(-1,1), y = 'P.Value', ylim = c(0,5), pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

#volcano plot Nuts
png('downloads/plot/volcanoplot_nuts.png')
EnhancedVolcano(tT.nuts, title = "Nuts", lab = tT.nuts $Gene.symbol, 
                labSize = 3, x = 'logFC', xlim = c(-1,1), y = 'P.Value', ylim = c(0,5), pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

#volcano plot Low fat
png('downloads/plot/volcanoplot_Lowfat.png')
EnhancedVolcano(tT.olive, title = "Low fat", lab = tT.lowfat$Gene.symbol, 
                labSize = 3, x = 'logFC', xlim = c(-1,1), y = 'P.Value', ylim = c(0,5), pCutoff = 0.05, FCcutoff = 0.26)

dev.off()

### Determine the differentially expressed genes (DEGs) olive oil
DEG_tT.olive <- tT.olive[tT.olive$P.Value< 0.05, c(1:6)]
dim(DEG_tT.olive)
colnames(DEG_tT.olive)

#determening all the DEGs of olive which are downregulated
DEG_tT.olive2.0 <- tT.olive[tT.olive$P.Value< 0.05 & tT.olive$logFC<0, c(1:6)]

#determening all the DEGs of olive which are upregulated
DEG_tT.olive3.0 <- tT.olive[tT.olive$P.Value< 0.05 & tT.olive$logFC>0, c(1:6)]

# comaring DEGs to the DEGs from paper
tT.olive[tT.olive$Gene.symbol == "IL1B", ]
tT.olive[tT.olive$Gene.symbol == "IGFR2", ]
tT.olive[tT.olive$Gene.symbol == "ICAM1", ]
tT.olive[tT.olive$Gene.symbol == "TNF", ] #paper has a log2 ratio around -2
tT.olive[tT.olive$Gene.symbol == "PTGS2", ]
tT.olive[tT.olive$Gene.symbol == "VEGF", ]


#determine DEGs nuts
DEG_tT.nuts <- tT.nuts[tT.nuts$P.Value< 0.05, c(1:6)]
dim(DEG_tT.nuts)

#determening all the DEGs of nuts which are downregulated
DEG_tT.nuts2.0 <- tT.nuts[tT.nuts$P.Value< 0.05 & tT.nuts$logFC<0, c(1:6)]

#determening all the DEGs of nuts which are upregulated
DEG_tT.nuts3.0 <- tT.nuts[tT.nuts$P.Value< 0.05 & tT.nuts$logFC>0, c(1:6)]
tT.nuts[tT.nuts$Gene.symbol == "IL1B", ]
tT.nuts[tT.nuts$Gene.symbol == "IGFR2", ]
tT.nuts[tT.nuts$Gene.symbol == "ICAM1", ]
tT.nuts[tT.nuts$Gene.symbol == "TNF", ]
tT.nuts[tT.nuts$Gene.symbol == "PTGS2", ]
tT.nuts[tT.nuts$Gene.symbol == "VEGF", ]


#determine DEGs lowfat
DEG_tT.lowfat <- tT.lowfat[tT.lowfat$P.Value< 0.05, c(1:6)]
dim(DEG_tT.lowfat)

#determening all the DEGs of low fat which are downregulated
DEG_tT.lowfat2.0 <- tT.lowfat[tT.lowfat$logFC<0 &tT.lowfat$P.Value< 0.05,(1:6)]

#determening all the DEGs of low fat which are upregulated
DEG_tT.lowfat3.0 <- tT.lowfat[tT.lowfat$P.Value< 0.05 & tT.lowfat$logFC>0,(1:6)]

DEG_tT.lowfat3.0
tT.lowfat[tT.lowfat$Gene.symbol == "IL1B", ]
tT.lowfat[tT.lowfat$Gene.symbol == "IGFR2", ]
tT.lowfat[tT.lowfat$Gene.symbol == "ICAM1", ]
tT.lowfat[tT.lowfat$Gene.symbol == "TNF", ]
tT.lowfat[tT.lowfat$Gene.symbol == "PTGS2", ]
tT.lowfat[tT.lowfat$Gene.symbol == "VEGF", ]

### Create a Venn diagram to compare the genes in the after the intervention of olive oil
### nuts and low fat using the R-package: VennDiagram

venn.diagram(x = list(DEG_tT.olive$ID, DEG_tT.nuts$ID, DEG_tT.lowfat$ID),
             category.names = c("Olive","Nuts","Lowfat"),
             output=FALSE,
             filename = 'downloads/plot/venn_diagram_comparisondiets.png',
             col=c("blue","red","yellow"),
             cex = 1.5,
             cat.pos = 4,
             main = "DEG overlap between 3 different diets")

#venndiagram downregulated genes
venn.diagram(x = list(DEG_tT.olive2.0$ID, DEG_tT.nuts2.0$ID, DEG_tT.lowfat2.0$ID),
             category.names = c("Olive","Nuts","Lowfat"),
             output=FALSE,
             filename = 'downloads/plot/venn_diagram_comparisondietsdown.png',
             col=c("blue","red","yellow"),
             cex = 1.5,
             cat.pos = 4,
             main = "down regulated DEG overlap between 3 different diets")

#venndiagram upregulated genes
venn.diagram(x = list(DEG_tT.olive3.0$ID, DEG_tT.nuts3.0$ID, DEG_tT.lowfat3.0$ID),
             category.names = c("Olive","Nuts","Lowfat"),
             output=FALSE,
             filename = 'downloads/plot/venn_diagram_comparisondietsup.png',
             col=c("blue","red","yellow"),
             cex = 1.5,
             cat.pos = 3,
             main = "up regulated DEG overlap between 3 different diets")


# adding Ensemb gene id and entrezgene id to affy gene list

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

idlist=getBM(attributes = c('entrezgene_id', "ensembl_gene_id","external_gene_name"),
                  filters = 'entrezgene_id',
                  values = tT.olive$Gene.ID, 
                  mart = ensembl)

idlist<- idlist %>%
  group_by(entrezgene_id) %>%
  filter(ensembl_gene_id == min(ensembl_gene_id))



# change column name of id_list so that they are similar to DEG_tT
colnames(idlist)[colnames(idlist) == "entrezgene_id"] <- "Gene.ID"

# merging idlist with DEG_tT entrezgene id identifiers  of each group
merged_lowfat <- merge(DEG_tT.lowfat, idlist, by = "Gene.ID")
merged_nuts <- merge(DEG_tT.nuts, idlist, by = "Gene.ID")
merged_olive <- merge(DEG_tT.olive, idlist, by = "Gene.ID")

# merging idlist with DEG_tT entrezgene id identifiers  of each group 
#for later comparison

merged_olivetT <- merge(tT.olive, idlist, by = "Gene.ID")
merged_nutstT <- merge(tT.nuts, idlist, by = "Gene.ID")
merged_lowfatT <- merge(tT.lowfat, idlist, by = "Gene.ID")


# over representation analyses with GO database
ego_olive <- enrichGO(gene          = merged_olive$Gene.ID,
                    universe        = tT.olive$Gene.ID,  # only works when tT and not idlist is used  
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  =  1, # should not be to low otherwise goplot function does not work
                     qvalueCutoff  = 1,
                     readable      = TRUE)

ego_nuts <- enrichGO(gene =   merged_nuts$Gene.ID,
                universe      = tT.nuts$Gene.ID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)

ego_lowfat <- enrichGO(gene   = merged_lowfat$Gene.ID,
                universe      = tT.lowfat$Gene.ID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 1,
                qvalueCutoff  = 1,
                readable      = TRUE)


png("downloads/Plot/ora.dotplot_olive.png")
dotplot(ego_olive, showCategory=15, title = "DOTplot ORA Olive", font.size = 10) 
dev.off()

png("downloads/Plot/ora.dotplot_nuts.png")
dotplot(ego_nuts, showCategory=15, title = "Dotplot ORA Nuts", font.size = 10)
dev.off()

png("downloads/Plot/ora.dotplot_lowfat.png")
dotplot(ego_lowfat, showCategory=15, title = "Dotplot ORA Lowfat", font.size = 10)
dev.off()

#making excel sheets for GO
ego_olive.df <- as.data.frame(ego_olive)
excel_file <- "ora.go_olive.xlsx"
write_xlsx(ego_olive.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

ego_nuts.df <- as.data.frame(ego_nuts)
excel_file <- "ora.go_nuts.xlsx"
write_xlsx(ego_nuts.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

ego_lowfat.df <- as.data.frame(ego_lowfat)
excel_file <- "ora.go_lowfat.xlsx"
write_xlsx(ego_lowfat.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

# Over representation analysis KEGG
ora.kegg_olive <- enrichKEGG(gene         =  merged_olive$Gene.ID,
                  organism     = 'hsa',
                  qvalueCutoff = 1,
                  pvalueCutoff = 1)

ora.kegg_nuts <- enrichKEGG(gene         =  merged_nuts$Gene.ID,
                 organism     = 'hsa',
                 qvalueCutoff = 1,
                 pvalueCutoff = 1)

ora.kegg_lowfat <- enrichKEGG(gene         =  merged_lowfat$Gene.ID,
                   organism     = 'hsa',
                   qvalueCutoff = 0.5,
                   pvalueCutoff = 0.5)

browseKEGG(kk@result, 'hsa04024')
kk@result

#making excel sheets for ora KEGG
ora.kegg_olive.df <- as.data.frame(ora.kegg_olive)
excel_file <- "ora.kegg_olive.xlsx"
write_xlsx(ora.kegg_olive.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

ora.kegg_nuts.df <- as.data.frame(ora.kegg_nuts)
excel_file <- "ora.kegg_nuts.xlsx"
write_xlsx(ora.kegg_nuts.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

ora.kegg_lowfat.df <- as.data.frame(ora.kegg_lowfat)
excel_file <- "ora.kegg_lowfat.xlsx"
write_xlsx(ora.kegg_lowfat.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

#making csv files of the ora KEGG results for Decopath analysis

write.csv( ora.kegg_lowfat ,"ora.kegg_lowfat.csv", row.names = FALSE)

# Over representation analysis with wikipathways as database

ora.wiki_olive <- enrichWP(merged_olive$Gene.ID, organism = "Homo sapiens", pvalueCutoff = 1, qvalueCutoff = 1)
ora.wiki_nuts <- enrichWP(merged_nuts$Gene.ID, organism = "Homo sapiens", pvalueCutoff = 1, qvalueCutoff = 1)
ora.wiki_lowfat <- enrichWP(merged_lowfat$Gene.ID, organism = "Homo sapiens", pvalueCutoff = 1, qvalueCutoff = 1)

#making excel sheets for ora wikipathways

ora.wiki_olive.df <- as.data.frame(ora.wiki_olive)
excel_file <- "ora.wiki_olive.xlsx"
write_xlsx(ora.wiki_olive.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

ora.wiki_nuts.df <- as.data.frame(ora.wiki_nuts)
excel_file <- "ora.wiki_nuts.xlsx"
write_xlsx(ora.wiki_nuts.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

ora.wiki_lowfat.df <- as.data.frame(ora.wiki_lowfat)
excel_file <- "ora.wiki_lowfat.xlsx"
write_xlsx(ora.wiki_lowfat.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

# Over representation analysis with reactome
ora.reactome_olive <- enrichPathway(gene=merged_olive$Gene.ID, pvalueCutoff = 1, readable=TRUE, qvalueCutoff = 1)
ora.reactome_nuts <- enrichPathway(gene=merged_nuts$Gene.ID, pvalueCutoff = 1, readable=TRUE, qvalueCutoff = 1)
ora.reactome_lowfat <- enrichPathway(gene=merged_lowfat$Gene.ID, pvalueCutoff = 1, readable=TRUE, qvalueCutoff = 1)

#making excel sheets for ora reactome
ora.reactome_olive.df <- as.data.frame(ora.reactome_olive)
excel_file <- "ora.reactome_olive.xlsx"
write_xlsx(ora.reactome_olive.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

ora.reactome_nuts.df <- as.data.frame(ora.reactome_nuts)
excel_file <- "ora.reactome_nuts.xlsx"
write_xlsx(ora.reactome_nuts.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

ora.reactome_lowfat.df <- as.data.frame(ora.reactome_lowfat)
excel_file <- "ora.reactome_lowfat.xlsx"
write_xlsx(ora.reactome_lowfat.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

# gene set enrichment analysis (gse)
#making list of ranked genes with expression values for each group
geneList.olive = merged_olive[,4]
names(geneList.olive) = as.character(merged_olive[,1])
geneList.olive = sort(geneList.olive, decreasing = TRUE)

geneList.nuts = merged_nuts[,4]
names(geneList.nuts) = as.character(merged_nuts[,1])
geneList.nuts = sort(geneList.nuts, decreasing = TRUE)

geneList.lowfat = merged_lowfat[,4]
names(geneList.lowfat) = as.character(merged_lowfat[,1])
geneList.lowfat = sort(geneList.lowfat, decreasing = TRUE)


# Gene set enrichment analysis with GO
gsea.GO_olive <- gseGO(geneList = geneList.olive,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       minGSSize    = 50,
                       maxGSSize    = 500,
                       pvalueCutoff = 1,
                       verbose      = FALSE)

gsea.GO_nuts <- gseGO(geneList = geneList.nuts,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       minGSSize    = 50,
                       maxGSSize    = 500,
                       pvalueCutoff = 1,
                       verbose      = FALSE)

gsea.GO_lowfat <- gseGO(geneList = geneList.lowfat,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       minGSSize    = 50,
                       maxGSSize    = 500,
                       pvalueCutoff = 1,
                       verbose      = FALSE)

#making excel sheets for gsea GO
gsea.GO_olive.df <- as.data.frame(gsea.GO_olive)
excel_file <- "gsea.GO_olive.xlsx"
write_xlsx(gsea.GO_olive.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

gsea.GO_nuts.df <- as.data.frame(gsea.GO_nuts)
excel_file <- "gsea.GO_nuts.xlsx"
write_xlsx(gsea.GO_nuts.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

gsea.GO_lowfat.df <- as.data.frame(gsea.GO_lowfat)
excel_file <- "gsea.GO_lowfat.xlsx"
write_xlsx(gsea.GO_lowfat.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

## Gene set enrichment analysis with KEGG with p valueCutoff = 1
gsea.KEGG_olive<- gseKEGG(geneList     = geneList.olive,
               organism     = 'hsa',
               pvalueCutoff = 1,
               verbose      = FALSE)

gsea.KEGG_nuts<- gseKEGG(geneList     = geneList.nuts,
                          organism     = 'hsa',
                          pvalueCutoff = 1,
                          verbose      = FALSE)

gsea.KEGG_lowfat<- gseKEGG(geneList     = geneList.lowfat,
                          organism     = 'hsa',
                          pvalueCutoff = 1, 
                          verbose      = FALSE)

#making excel sheets for gsea KEGG
gsea.KEGG_olive.df <- as.data.frame(gsea.KEGG_olive)
excel_file <- "gsea.KEGG_olive.xlsx"
write_xlsx(gsea.KEGG_olive.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

gsea.KEGG_nuts.df <- as.data.frame(gsea.KEGG_nuts)
excel_file <- "gsea.KEGG_nuts.xlsx"
write_xlsx(gsea.KEGG_nuts.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

gsea.KEGG_lowfat.df <- as.data.frame(gsea.KEGG_lowfat)
excel_file <- "gsea.KEGG_lowfat.xlsx"
write_xlsx(gsea.KEGG_lowfat.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

#writing csv file of gsea KEGG for decopath analysis
write.csv( gsea.KEGG_olive ,"gsea.kegg_olive.csv", row.names = FALSE)

## Gene set enrichment analysis with wikipathways using a pvalueCutoff =1
gsea.WIKI_olive <- gseWP(geneList.olive, organism = "Homo sapiens",pvalueCutoff = 1)
gsea.WIKI_nuts <- gseWP(geneList.nuts, organism = "Homo sapiens",pvalueCutoff = 1)

#making excel sheets for gsea wikipathways
gsea.WIKI_olive.df <- as.data.frame(gsea.WIKI_olive)
excel_file <- "gsea.WIKI_olive.xlsx"
write_xlsx(gsea.WIKI_olive.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

gsea.WIKI_nuts.df <- as.data.frame(gsea.WIKI_nuts)
excel_file <- "gsea.WIKI_nuts.xlsx"
write_xlsx(gsea.WIKI_nuts.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

gsea.WIKI_lowfat.df <- as.data.frame(gsea.WIKI_lowfat)
excel_file <- "gsea.WIKI_lowfat.xlsx"
write_xlsx(gsea.WIKI_lowfat.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

## Gene set enrichment analysis with reactome with p valueCutoff= 1
gsea.reactome_olive <- gsePathway(geneList = geneList.olive, 
                pvalueCutoff = 1,
                pAdjustMethod = "BH", 
                verbose = FALSE)

gsea.reactome_nuts <- gsePathway(geneList = geneList.nuts, 
                                  pvalueCutoff = 1,
                                  pAdjustMethod = "BH", 
                                  verbose = FALSE)

gsea.reactome_lowfat <- gsePathway(geneList = geneList.lowfat, 
                                  pvalueCutoff = 1,
                                  pAdjustMethod = "BH", 
                                  verbose = FALSE)

viewPathway("gene expresseion diet nuts", 
            readable = TRUE, 
            foldChange = geneList.nuts)

#making excel sheets for gsea reactome
gsea.reactome_olive.df <- as.data.frame(gsea.reactome_olive)
excel_file <- "gsea.reactome_olive.xlsx"
write_xlsx(gsea.reactome_olive.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

gsea.reactome_nuts.df <- as.data.frame(gsea.reactome_nuts)
excel_file <- "gsea.reactome_nuts.xlsx"
write_xlsx(gsea.reactome_nuts.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

gsea.reactome_lowfat.df <- as.data.frame(gsea.reactome_lowfat)
excel_file <- "gsea.reactome_lowfat.xlsx"
write_xlsx(gsea.reactome_lowfat.df, excel_file)
cat("Excel file", excel_file, "has been created.\n")

#topology based analysis
#Making named vector for each group containing the logFC and Entrez.ID of significantly expressed genes (used for SPIA)
# Bioconductors SPIA package https://bioconductor.org/packages/release/bioc/vignettes/SPIA/inst/doc/SPIA.pdf

named.vector_olive<- setNames(merged_olive$logFC, merged_olive$Gene.ID)
named.vector_nuts<- setNames(merged_nuts$logFC, merged_nuts$Gene.ID)
named.vector_lowfat <- setNames(merged_lowfat$logFC, merged_lowfat$Gene.ID)

# Signaling Pathway Impact Analysis (SPIA) using fisher method to combine the two p-values with fisher to increase the amount of significant samples
spia_olive <- spia(de=named.vector_olive,all=merged_olivetT$Gene.ID,organism="hsa",data.dir=NULL,pathids=NULL,nB=2000,plots=TRUE,verbose=TRUE,beta = NULL,combine = "fisher")
spia_nuts <- spia(de=named.vector_nuts,all=merged_nutstT$Gene.ID,organism="hsa",data.dir=NULL,pathids=NULL,nB=2000,plots=TRUE,verbose=TRUE,beta = NULL,combine = "fisher")
spia_lowfat <- spia(de=named.vector_lowfat,all=merged_lowfatT$Gene.ID,organism="hsa",data.dir=NULL,pathids=NULL,nB=2000,plots=TRUE,verbose=TRUE,beta = NULL,combine = "fisher")

#plotting spia results
#minimum threshold to have significant results for all thre plots is 0.152, so threshold is set to 0,2
plotP(spia_olive, threshold = 0.2)

plotP(spia_nuts, threshold = 0.2)

plotP(spia_lowfat, threshold = 0.2)
#The pathways at the right of the red oblique line are significant after Bonferroni correction of the global p-values
#The pathways ath the right of the blue oblique line are significant after a FDR correction of global pvalue

# cepa for Over representation analysis
# Centrality based Pathway Enrichment https://cran.r-project.org/web/packages/CePa/vignettes/CePa.pdf

cepa.olive_Reactome <- cepa.all(dif = merged_olive$Gene.symbol, bk = tT.olive$Gene.symbol, pc = PID.db$Reactome)
cepa.nuts_Reactome <- cepa.all(dif = merged_nuts$Gene.symbol, bk = tT.nuts$Gene.symbol, pc = PID.db$Reactome)
cepa.lowfat_Reactome <- cepa.all(dif = merged_lowfat$Gene.symbol, bk = tT.lowfat$Gene.symbol, pc = PID.db$Reactome)

cepa.olive_KEGG <- cepa.all(dif = merged_olive$Gene.symbol, bk = tT.olive$Gene.symbol, pc = PID.db$KEGG)
cepa.nuts_KEGG <- cepa.all(dif = merged_nuts$Gene.symbol, bk = tT.nuts$Gene.symbol, pc = PID.db$KEGG)
cepa.lowfat_KEGG <- cepa.all(dif = merged_lowfat$Gene.symbol, bk = tT.lowfat$Gene.symbol, pc = PID.db$KEGG)

PID.db$Reactome$version

# making plots for cepa ora with different types of data catalogues showing all the pathways
plot(cepa.olive_Reactome)
plot(cepa.nuts_Reactome)
plot(cepa.lowfat_Reactome)

plot(cepa.olive_KEGG)
plot(cepa.nuts_KEGG)
plot(cepa.lowfat_KEGG)

# making ORA Cepa plots with filtered data p-value < 0.05 for KEGG and Reactome
plot(cepa.olive_KEGG, adj.method = "none", only.sig = TRUE, cutoff = 0.05)
plot(cepa.nuts_KEGG, adj.method = "none", only.sig = TRUE, cutoff = 0.05)
plot(cepa.lowfat_KEGG, adj.method = "none", only.sig = TRUE, cutoff = 0.05)

plot(cepa.olive_Reactome, adj.method = "none", only.sig = TRUE, cutoff = 0.05)
plot(cepa.nuts_Reactome, adj.method = "none", only.sig = TRUE, cutoff = 0.05)
plot(cepa.lowfat_Reactome, adj.method = "none", only.sig = TRUE, cutoff = 0.05)

# making expression matrix for GSEA Cepa and Netgsa
ex.top<- data.frame(ex)
ex.top <- rownames_to_column(ex.top, var = "affy.ID")

# Cepa requires external gene name and Netgsa requires entrezgene ID
idlist.ex=getBM(attributes = c('entrezgene_id',"affy_hg_u95av2", "ensembl_gene_id","external_gene_name"),
             filters = 'entrezgene_id',
             values = tT.olive$Gene.ID, 
             mart = ensembl)

# double IDs are filtered out
idlist.ex<- idlist.ex %>%
  group_by(entrezgene_id) %>%
  filter(ensembl_gene_id == min(ensembl_gene_id))

idlist.ex<- idlist.ex[idlist.ex$affy_hg_u95av2 != "", ]

colnames(idlist.ex)[colnames(idlist.ex) == "affy_hg_u95av2"] <- "affy.ID"

matrix.ex <- merge(idlist.ex, ex.top, by = "affy.ID")

matrix.ex<- matrix.ex %>%
  group_by(entrezgene_id) %>%
  filter(ensembl_gene_id == min(ensembl_gene_id))

matrix.ex<- matrix.ex %>%
  group_by(entrezgene_id) %>%
  filter(affy.ID == min(affy.ID))

#making matrix for cepa gsa
matrix.ex.cepa <- matrix.ex[, -c(1, 2, 3)]
matrix.ex.cepa <- as.data.frame(matrix.ex.cepa)
rownames(matrix.ex.cepa) <- matrix.ex.cepa$external_gene_name
rownames(matrix.ex.cepa)[rownames(matrix.ex.cepa) == "DTX2P1-UPK3BP1-PMS2P11"] <- "PMS2P11"
matrix.ex.cepa <- matrix.ex.cepa[, -1]
  
matrix.ex <- matrix.ex[, -1]
matrix.ex <- matrix.ex[, -2]

matrix.ex <- as.data.frame(matrix.ex)
rownames(matrix.ex) <- matrix.ex$entrezgene_id
rownames(matrix.ex) <- paste0("ENTREZID:", rownames(matrix.ex))
matrix.ex <- matrix.ex[, -1]
matrix.ex <- matrix.ex[, -1]

# making a numerical vector where each number stand for a group 1=olive 2=nuts 3=lowfat
group.vector <-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,1,1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,1,1,2,2,3)
unique_groups <- unique(group.vector)
split_matrices <- vector("list", length(unique_groups))

#splitting the matrix into three groups based on their number for matrix with gene symbols (cepa)
for (i in seq_along(unique_groups)) {
  group_id <- unique_groups[i]
  columns_in_group <- which(group.vector == group_id)
  split_matrices[[i]] <- matrix.ex.cepa[, columns_in_group, drop = FALSE]}

for (i in seq_along(split_matrices)) {
  cat("Group", unique_groups[i], ":\n")
  print(split_matrices[[i]])}

olive_matrix.ex.cepa<-split_matrices[[1]]
nuts_matrix.ex<-split_matrices[[2]]
low.fat_matrix.ex<-split_matrices[[3]]

#making vector for each groups to label the control and intervention based on the gsms vector
#gsms = "10101010101010101032323232323232323232544441054545454554545410324"
olive.vector<-c(1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0) # 1 stands for baseline and 0 stands for intervention with olive oil
nuts.vector<-c(3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2) # 3 stands for baseline and 2 stands for intervention wiht nuts 
low.fat.vector<-c(5,4,4,4,4,5,4,5,4,5,4,5,4,5,5,4,5,4,5,4,4) # 4 stands for baseline and 5 stands for intervention with lowfat

label = read.cls("https://mcube.nju.edu.cn/jwang/lab/soft/cepa/P53.cls", 
                 treatment="MUT", control="WT")
n_samples <- length(olive.vector)
treatment_label <- "1"
control_label <- "0"
class_labels <- paste(ifelse(olive.vector == 1, treatment_label, control_label), collapse = " ")
writeLines(paste(n_samples, "2 1", class_labels), "olive.vector.cls")
label.olive = read.cls("olive.vector.cls", treatment=treatment_label, control=control_label)

cepa.all(mat = "olive.output.gct", label = label.olive, pc = PID.db$NCI, nlevel = "tvalue_abs", plevel = "mean")

output_file <- "olive.output.gct"
write.table(olive_matrix.ex, file = output_file, sep = "\t", quote = FALSE, col.names = NA)

#splitting the matrix into three groups based on their number for matrix with entrez ID
for (i in seq_along(unique_groups)) {
  group_id <- unique_groups[i]
  columns_in_group <- which(group.vector == group_id)
  split_matrices[[i]] <- matrix.ex[, columns_in_group, drop = FALSE]}

for (i in seq_along(split_matrices)) {
  cat("Group", unique_groups[i], ":\n")
  print(split_matrices[[i]])}

olive_matrix.ex<-split_matrices[[1]]
nuts_matrix.ex<-split_matrices[[2]]
low.fat_matrix.ex<-split_matrices[[3]]

#Netgsa

paths <- graphite::pathways('hsapiens','kegg')
database_search <- obtainEdgeList(rownames(olive_matrix.ex), c("kegg", "reactome"))

network_info <- prepareAdjMat(olive_matrix.ex,olive.vector, database_search)
?prepareAdjMat()
netgsa<-NetGSA(network_info[["Adj"]], olive_matrix.ex, olive.vector, pathways_mat, 
       lklMethod = "REHE", sampling = TRUE, 
       sample_n = 0.25, sample_p = 0.25)
?netgsa
str(pathways_mat)

