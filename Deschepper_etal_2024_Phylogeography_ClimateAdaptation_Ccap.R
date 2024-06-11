# INTRODUCTION ----
# This R script goes along with a submitted manuscript titled "Extensive admixture on a global scale and local climatic adaptation of the 
# Mediterranean fruit fly (Diptera, Tephritidae: Ceratitis capitata) revealed by whole genome sequencing by Deschepper et al. (2024).


# PROCEDURE TO MODEL ADMIXTURE USING ADMIXTOOLS 2 ----
# findgraph.R: script to start admixture graph searching based on f2 statistics. This script can be called with findgraph.slurm, 
# in which we specify the number of admixture events in the graph M and the number of graph searching iterations I. The output are
# files containing the best model for each model searching iteration. With testmodels.R, you can combine the outputted models into
# a dataframe and test whether a specific model has significantly better fit than another one. 

# FOLLOWING SCRIPTS ARE ASSOCIATED WITH DETECTING CLIMATE ADAPTATION USING RDA, PCADAPT AND LFMM ----
## PCadapt ----
library(pcadapt)
library(qvalue)
library(SNPRelate)
library(rtracklayer)
setwd("/yourdir")
filename <- read.pcadapt("yourfile.vcf.gz", type = "vcf")
popreg <- read.table("metadate.txt", header = TRUE) # a file containing following info per sample: ID,	Region,	decimalLatitude,	decimalLongitude

# Choosing the number K of Principal Components
# screeplot
x <- pcadapt(input = filename, K = 20, LD.clumping = list(size = 500, thr = 0.1)) 
plot(x, option = "screeplot")
# scoreplot
plot(x, option = "scores", pop = popreg$Region)
plot(x, option = "scores", pop = popreg$Region, i = 1, j = 6)
# Compute test statistics based on optimal number of K
x <- pcadapt(filename, K = 5, method = "mahalanobis")
summary(x)
plot(x , option = "manhattan")
# check the expected uniform distribution of the p-values using a Q-Q plot
plot(x, option = "qqplot")
# Histograms of the test statistic and of the p-values
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x, option = "stat.distribution")

# Choosing a cutoff for outlier detection
# q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qval < alpha)
length(outliers)
# Benjamini-Hochberg procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.01
outliers <- which(padj < alpha)
length(outliers)
# Bonferroni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.01
outliers <- which(padj < alpha)
length(outliers)

# To evaluate if LD might be an issue for your dataset: display loadings
par(mfrow = c(2, 2))
for (i in 1:4) {
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
}

# extract outlier snp positions
bim <- read.table("Ccglobal_DovetailRef_PASS.filtered.renamed.selection.rmvind.rmvindels.bi.maf.05.DP3.miss1.bed.bim", sep = "\t")
outlier.pos <- bim[outliers,]
dim(outlier.pos)
write.table(outlier.pos[,c(1,4,4)], "D:/temp/Cc_Phylogeography/GWAS/PCadapt/pcadapt.outliers.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)


## Set up the genotypes dataset to run RDA and LFMM ----
setwd("/yourdir")
vcf.fn <- paste0(getwd(),"/yourfile.vcf.gz")
snpgdsVCF2GDS(vcf.fn, "Ccglobal_all_miss1maf0.05.gds", method="biallelic.only")
genofile = snpgdsOpen("Ccglobal_all_miss1maf0.05.gds")
coords <- read.table("/yourdir/coords.txt", header = TRUE) # This file contains following columns: ID,	decimalLatitude,	decimalLongitude
#important, check if the coordinates file and genotypes dataset have the same order
read.gdsn(index.gdsn(genofile, "sample.id")) == coords$ID
# get genotypes coded as 0,1,2
genotypes <- snpgdsGetGeno(genofile, snp.id = NULL, snpfirstdim=TRUE) 
# check how many times 0,1,2 and NA occur in snps, no NA allowed
table(c(genotypes), exclude=NULL) 
# make dataframe with all SNPs
df_init <- cbind(read.gdsn(index.gdsn(genofile, "sample.id")),
                 str_extract((read.gdsn(index.gdsn(genofile, "sample.id"))), "[^_]+")) #the region of each sample is pulled from the sample name before the underscore, adapt regex to your needs
df <- data.frame(df_init, t(genotypes))
colnames(df) <- c("ID","Region", 
                  paste0(read.gdsn(index.gdsn(genofile, "snp.chromosome")),"_",read.gdsn(index.gdsn(genofile, "snp.position"))))
dim(df)
df[1:10,1:10] #sanity check
# let's create a test set first with a reduced amount of snps
nr.snps = 50000
testdf <- cbind(df$Region, df[,sample(2:dim(df)[2], nr.snps)]) #sample randomly
testdf <- cbind(df$Region, df[,c(2:nr.snps)]) #sample continuously
df <- testdf #for convenience
colnames(df)[1] <- "Region"
region.vec <- df$Region
# make everything a factor (may take very long if you include all snps)
df[,1:dim(df)[2]] <- lapply(df[,1:dim(df)[2]], as.factor)
gen_df <- df # this is the genotypes dataframe we will use later on for RDA and LFMM, you can also read it in from the file we save below
write.table(x = gen_df, "CcapWorldwide_genotypes_miss1maf0.05.txt", quote = FALSE, row.names = FALSE)

### RDA ----
library(dplyr)
library(readr)
library(stringr)
library(sp)
library(raster)
library(psych)
library(vegan)

# Open the genotypes file created earlier
gen_df <- read_delim("CcapWorldwide_genotypes_miss1maf0.05.txt",
                     col_names = TRUE, num_threads = 10, delim = " ", progress = show_progress())
gen_df <- gen_df %>% arrange(Region)
dim(gen_df)

#### Download climatic data from worldclim and setup dataset ----
spdf <- SpatialPointsDataFrame(coords = coords[, c("decimalLongitude", "decimalLatitude")], data = coords)
# Download the WorldClim data 
worldclim_data <- getData("worldclim", var = c("bio"), res = 2.5)
# Associate worldclim data with coordinates of samples
region <- coords$ID %>% str_extract(., "[^_]+")
env_df <- as.data.frame(cbind(extract(worldclim_data, spdf), as.data.frame(spdf), region))
str(env_df)
# arrange the env dataframe in the same order as the gen dataframe
target <- gen_df$Region
env_df <- env_df %>% arrange(region, levels = target)
popcode <- env_df$ID %>% substr(., 1, nchar(.)-2)


#### checking the dataframes ----
# Confirm that genotypes and environmental data are in the same order
identical(as.vector(gen_df$Region), as.vector(env_df$region))

# assess correlation between predictors, correlation r > 0.7 must be avoided
pairs.panels(env_df[,1:19], scale=T)

# make a selection and review distribution within predictor, select your choice of predictors
pred <- subset(env_df, select=c(bio2, bio5, bio6, bio18, bio19))
pairs.panels(pred, scale=T)


#### Run the RDA ----
cc.rda <- rda(gen_df[,-1] ~ ., data=pred, scale=T)
# partial RDA: factor out the variation explained by geography (highly recommended)
cc.rda <- rda(gen_df[,-1] ~ . + Condition(env_df$decimalLatitude + env_df$decimalLongitude), data = pred, scale=T)
cc.rda
RsquareAdj(cc.rda)
summary(eigenvals(cc.rda, model = "constrained"))
screeplot(cc.rda, main = "RDA")

# check our RDA model for significance using formal tests
# The null hypothesis is that no linear relationship exists between the SNP data and the environmental predictors
signif.full <- anova.cca(cc.rda, parallel = 10) # default is permutation=999
signif.full
# check sign. for each axis (takes a lot of time) 
signif.axis <- anova.cca(cc.rda, by="axis", parallel = 10, permutations = 1)
signif.axis
# checking Variance Inflation Factors for the predictor variables used in the model
# values should be below 10, indicating that collinearity is not a problem
vif.cca(cc.rda)
# basic plotting
plot(cc.rda, scaling = 3)          # default is axes 1 and 2
plot(cc.rda, choices = c(1, 3), scaling = 3)  # axes 1 and 3
# plotting with colors for regions
region <- gen_df$Region %>% as.vector() %>% as.factor()
cols <- c("#66c2a5","#c6d0e5","#FC8D62","#336153","#8da0cb")
cols[region]
# axes 1 & 2
plot(cc.rda, type="n", scaling=3)
points(cc.rda, display="species", pch=20, cex=0.7, col="gray", scaling=3)           # the SNPs
points(cc.rda, display="sites", pch=c(21:25)[region], cex=1.3, col="gray32", scaling=3, bg = cols[region]) # the samples
text(cc.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(region),
       bty="n", col="gray32", pch=21, cex=1, pt.bg=cols)
# axes 1 & 3
plot(cc.rda, type="n", scaling=3, choices = c(1,3))
points(cc.rda, display="species", pch=20, cex=0.7, col="gray", scaling=3)           # the SNPs
points(cc.rda, display="sites", pch=c(21:25)[region], cex=1.3, col="gray32", scaling=3, bg = cols[region]) # the samples
text(cc.rda, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=c("West Africa","East Africa","Mediterranean","Central America","South America"),
       bty="n", col="gray32", pch=21, cex=1.2, pt.bg=bg)

#### Identify candidate SNPs involved in local adaptation ----
load.rda <- scores(cc.rda, choices=c(1:3), display="species") #SNP scores are stored as 'species'
hist(load.rda[,1], main="Loadings on RDA1") # loadings in the tails are likely to be under selection
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
# apply outlier function to significant constraint axes
nr.stdv = 3 # standard deviation cutoff 
cand1 <- outliers(load.rda[,1],nr.stdv) 
cand2 <- outliers(load.rda[,2],nr.stdv)  
cand3 <- outliers(load.rda[,3],nr.stdv)  
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand

# create a data frame summarizing the cand. snps 
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# add in the correlations of each candidate SNP with the eight environmental predictors
foo <- matrix(nrow=(ncand), ncol=ncol(pred))  
colnames(foo) <- colnames(pred)

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen_df[,nam]
  foo[i,] <- apply(pred,2,function(x) cor(x,snp.gen))
}
cand <- cbind.data.frame(cand,foo)  
head(cand)


#### futher investigate candidate snps ----
# looking for duplicate detections
length(cand$snp[duplicated(cand$snp)])
foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # duplicates on axis 1
table(foo[foo[,1]==2,2]) # duplicates on axis 2
table(foo[foo[,1]==3,2]) # duplicates on axis 3

# see which of the predictors each candidate SNP is most strongly correlated with
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,4+(ncol(pred))] <- names(which.max(abs(bar[4:(3+ncol(pred))]))) # gives the variable
  cand[i,4+(ncol(pred)+1)] <- max(abs(bar[4:(3+ncol(pred))]))              # gives the correlation
}

colnames(cand)[4+(ncol(pred))] <- "predictor"
colnames(cand)[4+(ncol(pred)+1)] <- "correlation"
table(cand$predictor)

# save column ids of outliers linked to e.g. bio18 -> use them to filter out outlier loci from whole snp file
ids.bio18 <- match(subset(cand, predictor == "bio18")[,"snp"], colnames(dd)) %>% na.omit()
ids.bio19 <- match(subset(cand, predictor == "bio19")[,"snp"], colnames(dd)) %>% na.omit()
ids.bio2 <- match(subset(cand, predictor == "bio2")[,"snp"], colnames(dd)) %>% na.omit()
ids.bio5 <- match(subset(cand, predictor == "bio5")[,"snp"], colnames(dd)) %>% na.omit()
ids.bio6 <- match(subset(cand, predictor == "bio6")[,"snp"], colnames(dd)) %>% na.omit()

#### plotting the candidate snps ----
sel <- cand$snp
env <- cand$predictor
env[env=="bio2"] <- '#1f78b4'
env[env=="bio5"] <- '#33a02c'
env[env=="bio6"] <- '#6a3d9a'
env[env=="bio18"] <- '#e31a1c'
env[env=="bio19"] <- '#FFEF00'
          
# color by predictor:
col.pred <- rownames(cc.rda$CCA$v) # pull the SNP names
for (i in 1:length(sel)) {           # color code candidate SNPs
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}

col.pred[grep("Scaffold_",col.pred)] <- "#f1eef6" # non-candidate SNPs
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
bg <- c('#1f78b4','#33a02c','#6a3d9a','#e31a1c',"#FFEF00")
# axes 1 & 2
plot(cc.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1))
points(cc.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(cc.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(cc.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("bio2","bio5","bio6","bio18","bio19"), bty="n", col="gray32", pch=21, cex=1.5, pt.bg=bg)

# axes 1 & 3
plot(cc.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices = c(1,3))
points(cc.rda, display="species", pch=21, cex=1, col="gray32", bg=col.pred, scaling=3)
points(cc.rda, display="species", pch=21, cex=1, col=empty.outline, bg=empty, scaling=3)
text(cc.rda, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("bio2","bio4","bio6","bio15"), bty="n", col="gray32", pch=21, cex=1.5, pt.bg=bg)

# save associated snps
chr <- cand$snp %>% str_replace(., "(.+)_\\d+$", "\\1")
pos <- cand$snp %>% str_replace(., ".+_(\\d+)$", "\\1")
write.table(cbind(chr, pos, pos, cand), file = "D:/temp/Cc_Phylogeography/GWAS/RDA/snp.assoc.RDA..txt", 
            quote = FALSE, col.names = c("chr","pos","pos", colnames(cand)),
            row.names = FALSE, sep = "\t")

### LFMM ----
#### Download climatic data from worldclim and setup dataset ----
coords <- read.table("/yourdir/coords.txt", header = TRUE) # This file contains following columns: ID,	decimalLatitude,	decimalLongitude
spdf <- SpatialPointsDataFrame(coords = coords[, c("decimalLongitude", "decimalLatitude")], data = coords)
# Download the WorldClim data 
worldclim_data <- getData("worldclim", var = c("bio"), res = 2.5)
# Associate worldclim data with coordinates of samples
env_df <- as.data.frame(cbind(extract(worldclim_data, spdf), as.data.frame(spdf)))
# Princ. comp. analysis of climatic variables to see variation
env_df2 <- env_df[c("bio2","bio5","bio6","bio18","bio19")]
env_PC <- prcomp(env_df2, center=T, scale=T)
screeplot(env_PC)
g <- ggbiplot(env_PC,
              obs.scale = 1,
              var.scale = 1,
              labels = substr(env_df$ID, 1, 7),
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68)
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
g

biplot(env_PC)

# Or plot every bioClim variables in 1D for all locations to investigate spatial variation
library(ggplot2)
library(ggtext)
library(cowplot)
library(ggrepel)
env_df3 <- cbind(env_df[c("bio2","bio5","bio6","bio18","bio19")], substr(env_df$ID, 1, nchar(env_df$ID)-2), env_df$region) %>% unique()
names(env_df3)[6] <- "pop"
names(env_df3)[7] <- "region"
titles <- c("BIO2 = Mean Diurnal Range","BIO5 = Max Temperature of Warmest Month","BIO6 = Min Temperature of Coldest Month","BIO18 = Precipitation of Warmest Quarter","BIO19 = Precipitation of Coldest Quarter")
colz <- brewer.pal(n = 5, name = 'Dark2')
plot_list = list()
for (i in 5) {
  print(paste("Making plot for", names(env_df3[i])))
  p <- ggplot(env_df3, aes(x = env_df3[,i], y = 0, color = region
  )) +
    ggtitle(titles[i]) +
    annotate("segment", x = min(env_df3[,i]),xend = max(env_df3[,i]), y = 0, yend = 0, size = 2, col = colz[i]) +
    annotate("text", x = min(env_df3[,i]), y = -0.2, label = min(env_df3[,i]), size = 5) +
    annotate("text", x = max(env_df3[,i]), y = -0.2, label = max(env_df3[,i]), size = 5) +
    geom_point(size = 5) + scale_color_manual(values = c("CA" = "#66c2a5",
                                                         "EA" = "#c6d0e5",
                                                         "MED" ="#FC8D62",
                                                         "SA"="#336153",
                                                         "WA"= "#8da0cb")) +
    geom_label_repel(aes(label = pop), col = "black", angle = 90, nudge_y = 0.2,
                     direction = "both", max.overlaps = 20) +
    scale_y_continuous(limits = c(-1,1)) +
    theme(panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20))
  plot_list[[i]] <- p
}
plot_grid(plotlist = plot_list)


#### run the lfmm analysis ----
# set variables for climatic data (env_df) and genotype data (gen_df)
pred <- subset(env_df, select=c(bio2, bio5, bio6, bio18, bio19))
Y <- data.frame(gen_df) #genos
X <- pred
nr.snps <- dim(Y)[2] #you can first test run using a lower amount of snps

# Fit an LFMM with K = 5 factors
mod.lfmm <- lfmm_ridge(Y = Y, X = X, K = 5)


# single core: Perform association testing using the fitted model
X <- pred[,"bio19"] #env variable
mod.lfmm <- lfmm_ridge(Y = Y, X = X, K = 5)
pv <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif")
str(pv)
pvalues.inf <- pv$calibrated.pvalue #calibrated p-values using gif (genome inflation factor)
pos <- cbind(chr = rownames(pv$B), pval.adj = p.adjust(pvalues.inf, method = "BH")) %>% as.data.frame()#correct for mult. testing
colors <- brewer.pal(10,"Set2")[as.factor(substr(pos$chr,1,10))]
plot(-log10(as.numeric(pos$pval)), pch = 19, cex = .5, col = colors, xlab = "SNP", ylab = expression( -log10~italic(Pvalue) ))
legend("topleft", legend=levels(as.factor(substr(pos$chr,1,10))),
       bty="n", col="gray32", pch=21, cex=1, pt.bg=levels(as.factor(colors)))
sig = 0.01
abline(h = -log10(sig), lty = 2, col = "orange")
out.snps <- pos %>% filter(pval.adj < sig)
subset(pos, pval < sig)

# multi core: Perform association testing for multiple predictors in parallel
amp_up_models()
sig = 0.01
out.snps <- data.frame()
out.snps <- foreach(i = 1:length(pred)) %dopar% {
  library(lfmm)
  library(dplyr)
  mod.lfmm <- lfmm_ridge(Y = Y, X = pred[,i], K = 5)
  pv <- lfmm_test(Y = Y, 
                  X = pred[,i], 
                  lfmm = mod.lfmm, 
                  calibrate = "gif")
  pvalues.inf <- pv$calibrated.pvalue #calibrated p-values using gif (genome inflation factor)
  pos <- cbind(chr = rownames(pv$B), pval.adj = p.adjust(pvalues.inf, method = "BH"), predictor = rep(colnames(pred)[i], length(pvalues.inf))) %>% as.data.frame()#correct for mult. testing
  rbind(out.snps, pos %>% filter(pval.adj < sig))
}
out.snps.df <- do.call(rbind.data.frame, out.snps)
out.snps.df
table(out.snps.df$predictor)


# write file with associated loci
chr <- out.snps.df$chr %>% str_replace(., "(.+)_\\d+$", "\\1")
pos <- out.snps.df$chr %>% str_replace(., ".+_(\\d+)$", "\\1")
write.table(cbind(chr, pos, pos, out.snps.df$pval.adj, out.snps.df$predictor),
            "D:/temp/Cc_Phylogeography/GWAS/LFMM/snp.assoc.lfmm.txt",
            col.names = c("chr","pos","pos","pvalue_BH","predictor"), row.names = FALSE, quote = FALSE, sep = "\t")




## Identify candidate SNPs that are present in at least two of the methods applied to assess candidate snps ----
# get row id of intersection candidate snps of LFMM, RDA and PCadapt to extract info from gff
lfmm <- read.table("/yourdir/snp.assoc.lfmm.txt")[-1,1:2]
rda <- read.table("/yourdir/snp.assoc.RDA.txt")[,1:2]
pcadapt <- read.table("/yourdir/pcadapt.outliers.txt")[,1:2]
list_outliers <- list(RDA = paste0(rda[,1],"_", rda[,2]),
                      LFMM = paste0(lfmm[,1],"_", lfmm[,2]),
                      PCadapt = paste0(pcadapt[,1],"_", pcadapt[,2]))

# extract snps that are outliers in at least 2 methods
intLFMM_RDA <- intersect(list_outliers$LFMM, list_outliers$RDA) 
intRDA_PCadapt <- intersect(list_outliers$PCadapt, list_outliers$RDA) 
intLFMM_PCadapt <- intersect(list_outliers$LFMM, list_outliers$PCadapt) 
int <- c(intLFMM_RDA, intRDA_PCadapt, intLFMM_PCadapt)
chr <- int %>% str_replace(., "(.+)_\\d+$", "\\1")
pos <- int %>% str_replace(., ".+_(\\d+)$", "\\1")
write.table(cbind(chr,pos,pos), "/yourdir/2methods.outliers.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
# you can then use bedtools to extract gene annotation for the snps: 
# cat PCadapt.LFMM.outliers.txt | bedtools window -w 1 -a $gff -b stdin | awk '$3=="gene"'
# bedtools intersect -wb -a *.gff -b /mnt/d/temp/Cc_Phylogeography/GWAS/RDA.LFMM.outliers.txt
# venn diagram with RDA, LFMM and PCadapt snps
library(ggvenn)
ggvenn(list_outliers, fill_color = c("#80C673", "#EFC000FF", "#CD534CFF"),
       stroke_color = "transparent", set_name_size = 6, text_size = 5)


