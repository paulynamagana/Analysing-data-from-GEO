#install packages
#install.packages("BiocManager")
#install.packages("forcats")
#install.packages("stringr")
#install.packages("ggplot2")
#install.packages("ggrepel")
#install.packages("readr")
#install.packages("tidyr")
#install.packages("survminer")
#BiocManager::install("GEOquery")
#BiocManager::install("limma")
#BiocManager::install("pheatmap")
#BiocManager::install("org.Hs.eg.db")

#load library
library(GEOquery)
library(Biobase)
library(dplyr)
library(ggplot2)

##load series and platform data from GEO
#change my_id to the dataset that you want
my_id <- "GSE75037"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

#extract data
gse <- gse[[1]]
gse

p <- pData(gse) ## print the sample information/ phenotype data
f <- fData(gse) ## print the gene annotation/ feature data
x <- exprs(gse) ## print the expression matrix

# Create a boxplot of the first gene in the expression matrix
boxplot(x[1,] ~ p[,"source_name_ch1"], main = f[1,"symbol"])

#create expression object
eset <- ExpressionSet(assayData = x,
                      phenoData = AnnotatedDataFrame(p),
                      featureData= AnnotatedDataFrame(f))

#check the dimension
#view number of features(genes/rows) and samples(columns)
dim(eset)

#boxplot
boxplot(exprs(eset)[1,]~pData(eset)[,"source_name_ch1"],
              main=fData(eset)[1,"symbol"])

#visualisation
library(limma)

#plot distribution of samples
plotDensities(eset,legend=FALSE)

#log transform
exprs(eset) <- log(exprs(eset))
plotDensities(eset,legend=FALSE)

# Quantile normalize
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensities(eset, legend = FALSE);abline(v=1.5)

# Determine the genes with mean expression level greater than 1.5
keep <- rowMeans(exprs(eset)) > 1.5
sum(keep)

# Filter the genes
eset2 <- eset[keep]
plotDensities(eset2, legend = FALSE)

#building deisgn matris with model
design <- model.matrix(~0 + source_name_ch1, data=pData(eset))
design

colnames(design) <- c("Adenocarcinoma", "Non_malignant")
design

# Count the number of samples modeled by each coefficient
colSums(design)

# Create a contrasts matrix
cm <- makeContrasts(status = Adenocarcinoma  - Non_malignant,
                    levels = design)
cm

#fir the model
fit <- lmFit(eset,design)

#fit the contrasts
fit2 <- contrasts.fit(fit,contrasts=cm)

#calculate the t-statisics for the contrasts

fit2 <- eBayes(fit2)

#summarise
results <- decideTests(fit2)
summary(results)

#inspecting results
topTable(fit2,number = 5)

stats <- topTable(fit2,number= nrow(fit2), sort.by = "none")
dim(stats)

# Plot a histogram of the p-values
hist(stats[,"P.Value"])
# Create a volcano plot. Highlight the top 5 genes
volcanoplot(fit2,
            highlight=20,
            names=fit2$genes[,"Symbol"])

#ENRICHMENT TESTING
entrez <- fit2$genes[,"Entrez_Gene_ID"]
enrich_kegg <- kegga(fit2,geneid=entrez,species="Hs")
topKEGG(enrich_kegg,number=5)

libbrar




full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

## change according to your needs
p_cutoff <- 0.05
fc_cutoff <- 1

results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()













## source_name_ch1 and characteristics_ch1.1 seem to contain factors we might need for the analysis. Let's pick just those columns
sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1.1)

## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo,disease = source_name_ch1, patient=characteristics_ch1.1)
sampleInfo

#linear model
library(limma)
gse$characteristics_ch1.1 [gse$characteristics_ch1.1 == "NA"] <- NA

design <- model.matrix(~0+sampleInfo$disease)
design
colSums(design)

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Adjacent","Adenocarcinoma")
summary(exprs(gse))

fit <- lmFit(exprs(gse),design)
fit
head(fit$coefficients)

#define the contrast that we are interested in
contrasts <- makeContrasts(status=Adenocarcinoma - Adjacent, levels=design)
contrasts

fit2 <- contrasts.fit(fit, contrasts=contrasts)
head(fit2$coefficients,3)

#calculate the t-statistics
fit2 <- eBayes(fit2)
summary(fit)

#summarise results
#Test for differential expression between 2 groups
results <- decideTests(fit2)
summary(results)

