### Working directory and package loading ###
setwd("~/gemma")
library(rutilstimflutre)
library(MM4LMM)

### chosing the chromosome on which the analysis will be performed ###
chr <- 5
genome <- "B"

################################################
#### I . Running the analysis with MM4LMM ####

# Loading the phenotypic file
pheno <- read.table("pheno.csv", header=T, sep=",")
pheno <- pheno[order(pheno$planche, pheno$passage),]
pheno$pair <- paste(pheno$focal, pheno$neighbor, sep="_")
row.names(pheno) <- pheno$pair

# Loading the genotypic file
geno <- fread("geno.csv", data.table = F)
rownames(geno) <- geno$V1
geno <- geno[,-1]
geno <- as.matrix(geno)
geno <- imputeGenosWithMean(geno)

# Loading snp coordinates
anno <- read.table("anno.csv", sep=",", header=T, row.names = 1)
anno$chr <- substr(anno$Chr,1,1)
anno$genome <- substr(anno$Chr,2,2)
anno <- anno[, c("snpId", "chr","genome","pos")]
colnames(anno) <- c("snp","chr","genome","pos")

# Subseting the genotypic file to only keep the desired chromosome
SNPs.to.test<- anno[which(anno$chr==chr & anno$genome==genome),"snp"]
geno.to.test <- geno[,which(colnames(geno)%in%SNPs.to.test)]

# Defining variance-covariance matrices
## Residual matrix = identity 
Ires <- diag(nrow(pheno))

## Genetic matrix = identity matrix (identity between pairs)
load("Kinship.RData")

# Defining incidence matrices
Zpair <- model.matrix(object=~ -1+ pair, data=pheno) 
Zres <- diag(nrow(pheno))

# Estimating variance components 
res <- MMEst(Y=pheno[,"RYT_poids_grains_g_adj"], X = geno.to.test, ZList = list(PAIR = Zpair,  RES=Zres),  VarList =list(PAIR=varcovar, RES=Ires), Method = "Reml")

# Retrieving Wald p-values
pvals <- data.frame(snp = colnames(geno.to.test), pvalue = do.call(rbind,AnovaTest(res))[seq(2,length(AnovaTest(res))*2, by=2),"pval"])
pvals_MM4LMM <- merge(pvals, anno, by="snp")


################################################
#### II . Running the analysis with GEMMA ####

# Modifying the format of the annotation file
row.names(anno)<- anno[,"snp"]
anno$chr <- paste(anno$chr, anno$genome, sep="")
snp.coords <- anno[,c("snp","pos","chr")]
colnames(snp.coords)[2] <- "coord"

# Creating a matrix for the fixed effect (here only intercept)
W <- matrix(rep(1, nrow(pheno)), ncol=1, nrow=nrow(pheno))
row.names(W) <- pheno$pair

# Formating the phenotypic file 
y <-setNames(pheno[,"RYT_poids_grains_g_adj"], pheno[,"pair"])

# Running the analysis
fit <- gemma(model="ulmm", y=y, X=geno.to.test, snp.coords=snp.coords, K.c=varcovar, W=W)
pvals_GEMMA <- fit$tests[,c("rs","p_wald","chr","ps")]

################################################
#### III . Comparing the pvalues obtained with MM4LMM and the pvalues obtained with GEMMA ####
pvals <- merge(pvals_MM4LMM, pvals_GEMMA[,c("rs","p_wald")], by.x="snp", by.y="rs",all.x=T)
X11()
plot(pvalue~p_wald, data=pvals, xlab="GEMMA",ylab="MM4LMM")
