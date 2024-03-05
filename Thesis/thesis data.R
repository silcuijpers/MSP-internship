rm(list=ls())
library(GEOquery)
library(limma)

#loading data from GEO
gset <- getGEO("GSE28358",GSEMatrix=TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL571", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make column names
fvarLabels(gset) <- make.names(fvarLabels(gset))


#assigning eacht sample to one of the six groups
gsms <- "10101010101010101032323232323232323232544441054545454554545410324"
sml <- strsplit(gsms, split="")[[1]]

# replacing values smaller then 0 with NaN and log transfomr each element
ex<-  exprs(gset)
ex[which(ex <= 0)]<- NaN
exprs(gset) <- log2(ex)

# assigning a groupname to each sample
gs <- factor(sml)
groups <- make.names(c("olive oil","baseline","nuts","nuts baseline","low list","low fat baseline"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

#removing rows with missing values (NA)
gset <- gset[complete.cases(exprs(gset)), ]

#fitting a linear model to gene expression data
fit <- lmFit(gset, design) 

#??
cts <- paste(groups, c(tail(groups, -1), head(groups, 1)), sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
show(tT)

fit2 <- eBayes(fit2, 0.01)# standard devitaion of 0.01
tT<- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

