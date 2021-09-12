#load libraries
library(GEOquery)
library(Biobase)
library(dplyr)
library(limma)

##load series and platform data from GEO
#change my_id to the dataset that you want
my_id <- "GSE75037"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

#extract data
gse <- gse[[1]]
gse

## if more than one dataset is present, you can analyse the other dataset by changing the number inside the [[...]]
## e.g. gse <- gse[[2]]

#check the dimension
#view number of features(genes/rows) and samples(columns)
dim(gse)
#in this case, there are 166 samples and 48803 genes

#print information
p <- pData(gse) ## print the sample information
f <-fData(gse) ## print the gene annotation
x <- exprs(gse) ## print the expression data

### CHECK NORMALISATION ####
#inspect the data to discover what scale the data are presented in.
#exprs function can retrieve the expression values as a data frame;
#with one column per-sample and one row per-gene.
#summary function can then be used to print the distributions
summary(exprs(gse))
#A boxplot can also be generated to see if the data have been normalised
#If so, the distributions of each sample should be highly similar.
boxplot(exprs(gse),outline=FALSE)

#visualize the distribution of gene expression levels for each sample
plotDensities(gse, group = pData(gse)[,"source_name_ch1"], legend = "topright")
# Quantile normalize???

### ORGANISING SAMPLE INFORMATION ####

#information about the processing protocol can be extracted by the pData function.
#For your own data, you will have to decide which columns will be useful in
#the analysis. The column giving the main comparison(s)
#of interest and any potential confounding factors. 
sampleInfo <- p
head(sampleInfo)

##Let's pick just those columns that we might need for the analysis
sampleInfo <- select(sampleInfo, 44,characteristics_ch1.1)
## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo,group = 1, patient=characteristics_ch1.1)
sampleInfo

### PRINCIPAL COOMPONENTS ANALYSIS ####
#Plot principal components labeled by group
plotMDS(gse, labels=sampleInfo[,"group"],
        gene.selection="common")
#is there a batch effect?

library(pheatmap)
corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix)  

#Incorporate sample information onto the plot to try and understand the clustering
## Print the rownames of the sample information and check it matches the correlation matrix
rownames(sampleInfo)
colnames(corMatrix)

pheatmap(corMatrix,
         annotation_col=sampleInfo)
#are separated due to normal vs tumour?

#### SAVE CSV ####
### Look at the features data frame and decide the names of the columns you want to keep
#We can export the expression data to a csv for inspection in Excel
library(readr)
features <- select(f,Symbol,Entrez_Gene_ID,RefSeq_ID,Chromosome,Cytoband)
full_output <- cbind(features,x)
write_csv(full_output, path="gse_full_output.csv")

#### DIFFERENTIAL EXPRESSION ANALYSIS #####
library(limma)
design <- model.matrix(~0+sampleInfo$group)
design
#rename
colnames(design) <- c("adenocarcinoma","non_malignant")

# Count the number of samples modeled by each coefficient
colSums(design)

#In order to perform the differential analysis, we have to define the
#contrast that we are interested in. In our case we only have two groups
#and one contrast of interest.
contrasts <- makeContrasts(adenocarcinoma - non_malignant, levels=design)
contrasts

#fit the model
fit <- lmFit(exprs(gse), design)
head(fit$coefficients)

# Fit the contrasts
fit2 <- contrasts.fit(fit,contrasts = contrasts)
#calculate the t-statistics for the contrasts
#apply the empirical Bayes’ step to get our differential
#expression statistics and p-values
efit <- eBayes(fit2)

#### ADD GENES ID ####
anno <- f
anno

anno <- select(anno,Symbol, Entrez_Gene_ID,Chromosome,Cytoband)
efit$genes <- anno
topTable(efit)

#visualise
plotSA(efit, main="Final model: Mean-variance trend")
# Summarize results
#quick look at differential expression levels
topTable(efit)
#If we want to know how many genes are differentially-expressed overall
results <- (decideTests(efit))
summary(results)
 #Significance is defined using an adjusted p-value cutoff that
#is set at 5% by default
#coparison between expression levels in adenocarcinoma and non_mallignant
#8622 genes are found to be down-regulated in adenocarcinoma relative to non_malignant
#10214 are up-regulated un adenocarcinoma relative to non_malignant


# Plot a histogram of the p-values
stats <- topTable(efit, number = nrow(efit), sort.by = "none")
hist(stats[, "P.Value"],col = "grey50", border = "white")

#Some studies require more than an adjusted p-value cut-off
#For a stricter definition on significance, one may require
#log-fold-changes (log-FCs) to be above a minimum value.
tfit <- treat(efit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

adeno_vs_non <- topTreat(tfit, coef=1, n=Inf)
head(adeno_vs_non)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(2,16))

group=as.factor(sampleInfo$group)

##### GLIMMA ########3
#BiocManager::install("Glimma")
library(Glimma)
library(edgeR)

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

d <-glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
         side.main="Entrez_Gene_ID", counts=x, groups=group, launch=TRUE)
d

library(plotly)
htmlwidgets::saveWidget(as_widget(d), "test.html")

###### VISUALIZE INDIVIDUAL GENES#########
# Find the row which contains SLC22A5 expression data
SLC <- which(f[, "Symbol"] == "SLC22A5")

# Plot SLC22A5 expression versus genotype
boxplot(x[SLC, ] ~ p[, "source_name_ch1"],
        main = anno[SLC, ])

