setwd("C:/Users/Nina/Documents/courses/R_summer_school/GWAStutorial/GWAS_Workshop/")
library(GenABEL)
# Load data
genodata <- load.gwaa.data(
  phenofile = "data/phenotype.dat", 
  genofile = "data/genotype.raw", 
  makemap = F, 
  sort = F)
# Get information about number of individuals
nids(genodata)
# Number of SNPs in the data
nsnps(genodata)
# Summary of the first 2 markers
summary(genodata@gtdata)[1:5, ]
chr4.idx <- which(genodata@gtdata@chromosome == 4)
maf <- summary(genodata@gtdata)[chr4.idx[1:500],"Q.2"]
plot(maf, type = 'l', main="MAF on chr 4", 
     xlab = "marker", ylab = "MAF", col = "slateblue")
# Summary of phenotypes
summary(genodata@phdata)
# Make a table with all phenotypes per individual
library(DT)
datatable(genodata@phdata,options=list(pageLength=5))
attach(phdata(genodata))
#Binary trait analysis
tab <- table(pheno, sex)
rownames(tab) <- c("control","case")
colnames(tab) <- c("male","female")
tab
fisher.test(tab)
detach(phdata(genodata))
qc0 <- check.marker(genodata, call = 0.95, perid.call = 0.95, 
                    maf = 1e-08, p.lev = 1e-05, 
                    hweidsubset = genodata@phdata$pheno == 0)
autosomalMarkerNames <- autosomal(genodata.qc0)
# Compute genomic kinship matrix
genodata.qc0.gkin <- ibs(genodata.qc0[, autosomalMarkerNames], weight = "freq")
# Transform it to a distance matrix
genodata.qc0.dist <- as.dist(0.5 - genodata.qc0.gkin)
# Perform multidimensional scaling to display individuals on a 2D
genodata.qc0.pcs  <- cmdscale(genodata.qc0.dist, k=10)
# Label the samples according to case control status
genodata.qc0.pcs <- cbind(genodata.qc0.pcs, genodata.qc0@phdata$pheno)
# Plot the result
{plot(x=genodata.qc0.pcs[, 1], y=genodata.qc0.pcs[, 2], xlab = "PC 1", ylab = "PC 2",
      col=c("grey","slateblue")[as.factor(genodata.qc0.pcs[, 11])], main = "Genomic kinship")
  legend(x=0.06,y=0.09,legend=c("control","case"),
         col=c("grey","slateblue"),pch=19)}
# Subsample 5,000 SNPs that will be used for analysis
genodata.qc0.sub <- genodata.qc0[,sort(sample(1:nsnps(genodata.qc0),5000))]
# Extract allele frequencies
eaf <- summary(gtdata(genodata.qc0.sub))$"Q.2"
# Find relatives
relInfo <- findRelatives(genodata.qc0.sub,q=eaf,gkinCutOff=-1,nmeivec=c(1,2,3,4))
# Plot results
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot")
an.pca <- qtscore(pheno~genodata.qc0.pcs[, 1:10], genodata.qc0, trait="binomial")
plot(an.pca, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot")
estlambda(an.pca[, "P1df"], plot=TRUE)
# Calculate the Bonferroni threshold
bonferroni <- -log10(0.05/nsnps(genodata.qc0))
# Add the line to our plot
plot(an.pca, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot")
abline(h = bonferroni, col = "red")