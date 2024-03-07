#clean workspace TEST
rm(list=ls())

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

######################################################################################## 

#load installed packages
library(GEOquery)
library(limma)
library(EnhancedVolcano)
library(VennDiagram)


#loading data from GEO
gset <- getGEO("GSE28358",GSEMatrix=TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL571", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Make column names by using the sample names
# In the current dataset this is for example "GSM701076"
fvarLabels(gset) <- make.names(fvarLabels(gset))


# assigning each sample to one of the six groups
# Did you check whether this maps with the assignement of the groups to the samples? 
gsms <- "10101010101010101032323232323232323232544441054545454554545410324"
sml <- strsplit(gsms, split="")[[1]]

# replacing values smaller then 0 with NaN and log transform each element
ex<-  exprs(gset)
ex[which(ex <= 0)]<- NaN
exprs(gset) <- log2(ex)

#calculate the number of rows with Na
rows_with_na <- sum(!complete.cases(ex))
print(rows_with_na)


# assigning a groupname to each sample
gs <- factor(sml)
groups <- make.names(c("olive oil","baseline","nuts","nuts baseline","low fat","low fat baseline"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

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

#plot for adjusted p-value distribution olive oil
hist(tT.olive$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "Olive oil: P-adj value distribution")

# Olive oil p-value distribution
hist(tT.olive$P.Value, col = "grey", border = "white", xlab = "P.Value",
     ylab = "Number of genes", main = "Olive oil: P.Value distribution")

# Nuts adjusted p-value
hist(tT.nuts$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "Nuts: P-adj value distribution")

# Nuts p-value distribution
hist(tT.nuts$P.Value, col = "grey", border = "white", xlab = "P.Value",
     ylab = "Number of genes", main = "Nuts: P.Value distribution")

# Lowfat adjusted p-value
hist(tT.lowfat$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "Low fat: P-adj value distribution")
#why does is show 20000 genes eventhough I removed some of them
# Nuts p-value distribution
hist(tT.lowfat$P.Value, col = "grey", border = "white", xlab = "P.Value",
     ylab = "Number of genes", main = "Low fat: P.Value distribution")




#volcanoplot Olive oil
png('volcanoplot_olive.png')
EnhancedVolcano(tT.olive, title = "Olive oil", lab = tT.olive$Gene.symbol, 
                labSize = 3, x = 'logFC', xlim = c(-1,1), y = 'P.Value', ylim = c(0,5), pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

#volcano plot Nuts
png('volcanoplot_nuts.png')
EnhancedVolcano(tT.nuts, title = "Nuts", lab = tT.olive$Gene.symbol, 
                labSize = 3, x = 'logFC', xlim = c(-1,1), y = 'P.Value', ylim = c(0,5), pCutoff = 0.05, FCcutoff = 0.26)
dev.off()

#volcano plot Low fat
png('volcanoplot_Lowfat.png')
EnhancedVolcano(tT.olive, title = "Low fat", lab = tT.olive$Gene.symbol, 
                labSize = 3, x = 'logFC', xlim = c(-1,1), y = 'P.Value', ylim = c(0,5), pCutoff = 0.05, FCcutoff = 0.26)

dev.off()

### Determine the differentially expressed genes (DEGs)
DEG_tT.olive <- tT.olive[tT.olive$P.Value< 0.05, c(2:6)]
dim(DEG_tT.olive)
colnames(DEG_tT.olive)

DEG_tT.nuts <- tT.nuts[tT.nuts$P.Value< 0.05, c(2:6)]
dim(DEG_tT.nuts)
colnames(DEG_tT.nuts)

DEG_tT.lowfat <- tT.lowfat[tT.lowfat$P.Value< 0.05, c(2:6)]
dim(DEG_tT.n)
colnames(DEG_tT.nuts)


dT.olive <- decideTests(fit.olive, adjust.method="BH", p.value=0.05, lfc=0)
dT.nuts <- decideTests(fit.nuts, adjust.method="BH", p.value=0.05, lfc=0)
dT.lowfat<-decideTests(fit.lowfat, adjust.method="BH", p.value=0.05, lfc=0)


### Create a Venn diagram to compare the genes in the after the intervention of olive oil
### nuts and low fat using the R-package: VennDiagram



vennDiagram( x = list(
    olive = dT.olive,
    nuts = dT.nuts,
    lowfat = dT.lowfat), category.names = c("olive", "nuts","lowfat"))

