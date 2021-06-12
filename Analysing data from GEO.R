#install packages
install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")

#load library
library(GEOquery)

##load series and platform data from GEO
#change my_id to the dataset that you want
my_id <- "GSE32863"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

#extract data
gse <- gse[[1]]
gse

pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data

## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse))

# log2 transform
ex <- exprs(gse)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

# box-and-whisker plot
dev.new(width=3+ncol(gse)/6, height=5)
par(mar=c(7,4,2,1))
title <- paste ("GSE32863", "/", annotation(gse), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
dev.off()

library(dplyr)
library(limma)
library(umap)

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE32863", "/", annotation(gse), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE32863")

# UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)


#Inspect the clinical variables
sampleInfo <- pData(gse)
sampleInfo

## source_name_ch1 and characteristics_ch1.1 seem to contain factors we might need for the analysis. Let's pick just those columns
sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1.1)

## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo,group = source_name_ch1, patient=characteristics_ch1.1)
sampleInfo

library(pheatmap)
## argument use="c" stops an error if there are any missing data points

corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix)      

#We can incorporate sample information onto the plot to try and understand the clustering.
#We have already created such a data frame previously (sampleInfo).
#However, we need to take care that the rownames of these data match the columns of the correlation matrix.

## Print the rownames/colnames of the sample information and check it matches the correlation matrix
rownames(sampleInfo)
colnames(corMatrix)

pheatmap(corMatrix,
         annotation_col=sampleInfo)    

library(ggplot2)
library(ggrepel)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(exprs(gse)))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste("Patient", patient))) + geom_point() + geom_text_repel()

### lets' say are outliers are samples 1,2 and 3
## replace 1,2,3 with the outliers in your dataset
outlier_samples <- c(1,2,3)

gse <- gse[,-outlier_samples]

library(limma)
gse$characteristics_ch1.1 [gse$characteristics_ch1.1 == "NA"] <- NA

design <- model.matrix(~0+sampleInfo$group)
design

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Adjacent","Adenocarcinoma")

summary(exprs(gse))

## calculate median expression level
cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff

## Identify genes expressed in more than 2 samples
keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse <- gse[keep,]

###SOMETHING HAPPENED HERE AND HAD TO RE-LOAD GSE, no running previous code from cutoff
fit <- lmFit(exprs(gse), design)
head(fit$coefficients)

#define the contrast that we are interested i
contrasts <- makeContrasts(Adenocarcinoma - Adjacent, levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)

#apply the empirical Bayes’ step to get our differential expression statistics and p-values.
fit2 <- eBayes(fit2)
#get a look
topTable(fit2)

#If we want to know how many genes are differentially-expressed overall we can use
#the decideTests function.
decideTests(fit2)
table(decideTests(fit2))

## calculate relative array weights
aw <- arrayWeights(exprs(gse),design)
aw

fit <- lmFit(exprs(gse), design,
             weights = aw)
contrasts <- makeContrasts(Adenocarcinoma - Adjacent, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

anno <- fData(gse)
anno

anno <- select(anno,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
fit2$genes <- anno
topTable(fit2)

full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

## Make sure you have ggplot2 loaded
library(ggplot2)
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()

## change according to your needs
p_cutoff <- 0.05
fc_cutoff <- 1

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 20

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Symbol,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")


## Get the results for particular gene of interest
filter(full_results, Symbol == "SLC22A1")


#Heatmaps of selected genes
## Use to top 20 genes for illustration

topN <- 20
##
ids_of_interest <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(ID)
#extract the corresponding gene symbols.
gene_names <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(Symbol) 

## Get the rows corresponding to ids_of_interest and all columns
gene_matrix <- exprs(gse)[ids_of_interest,]

pheatmap(gene_matrix,
         labels_row = gene_names)
#It is often preferable to scale each row to highlight the differences in each gene across the dataset.
pheatmap(gene_matrix,
         labels_row = gene_names,
         scale="row")


my_genes <- c("SLC22A1","SLC22A2","SLC22A4", "SLC22A5")
ids_of_interest <-  filter(full_results,Symbol %in% my_genes) %>% 
  pull(ID)

gene_names <-  filter(full_results,Symbol %in% my_genes) %>% 
  pull(Symbol)

gene_matrix <- exprs(gse)[ids_of_interest,]
pheatmap(gene_matrix,
         labels_row = gene_names,
         scale="row")