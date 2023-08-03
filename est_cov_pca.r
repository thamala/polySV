#Script for conducting PCA on mixed ploidy VCF files
#Tuomas Hämälä, April 2023

library(vcfR)
library(ggplot2)
library(ggrepel)

vcf <- read.vcfR("test.vcf") #Filtered and LD-pruned VCF file

#Transform VCF to numeric genotypes
df <- extract.gt(vcf)
df[df == "0|0"] <- 0
df[df == "0|1"] <- 1
df[df == "1|0"] <- 1
df[df == "1|1"] <- 2
df[df == "0/0"] <- 0
df[df == "0/1"] <- 1
df[df == "1/0"] <- 1
df[df == "1/1"] <- 2
df[df == "0/0/0/0"] <- 0
df[df == "0/0/0/1"] <- 1
df[df == "0/0/1/1"] <- 2
df[df == "0/1/1/1"] <- 3
df[df == "1/1/1/1"] <- 4
df[df == "0/0/0/0/0/0"] <- 0
df[df == "0/0/0/0/0/1"] <- 1
df[df == "0/0/0/0/1/1"] <- 2
df[df == "0/0/0/1/1/1"] <- 3
df[df == "0/0/1/1/1/1"] <- 4
df[df == "0/1/1/1/1/1"] <- 5
df[df == "1/1/1/1/1/1"] <- 6
df <- data.frame(apply(df,2,function(x)as.numeric(as.character(x))))

#Remove samples with > 50% missing data
mis <- apply(df,2,function(x)sum(is.na(x))/length(x))
df <- df[,mis <= 0.5]

#Calculate allele frequencies
ploidy <- apply(df,2,max,na.rm=T)
p <- apply(df,1,function(x)sum(x,na.rm=T)/sum(ploidy[!is.na(x)]))

#Removing individuals can change allele frequencies, so we make sure that maf >= 0.05
df <- df[p >= 0.05 & p <= 0.95,]
p <- p[p >= 0.05 & p <= 0.95]

#Estimate a covariance matrix
n <- ncol(df)
cov <- matrix(nrow=n,ncol=n)
for(i in 1:n){
	for(j in 1:i){
		x <- mean(c(ploidy[i],ploidy[j]))
		cov[i,j] <- mean((df[,i]-x*p)*(df[,j]-x*p)/(x*p*(1-p)),na.rm=T)
		cov[j,i] <- cov[i,j]
	}	
}

#Do PCA on the matrix and plot
pc <- prcomp(cov,scale=T)
xlab <- paste0("PC1 (",round(summary(pc)$importance[2]*100),"%)")
ylab <- paste0("PC2 (",round(summary(pc)$importance[5]*100),"%)")
pcs <- data.frame(PC1=pc$x[,1],PC2=pc$x[,2],id=colnames(df),ploidy=ploidy)

ggplot(pcs, aes(PC1, PC2, color=as.factor(ploidy)))+
  geom_point(size=7)+
  labs(x=xlab, y=ylab, color="Ploidy")+
  geom_text_repel(aes(label=id), size=4, force=20, color="black")+
  theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15, "cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12, color="black"),
		plot.title=element_text(size=14, color="black", hjust = 0.5),
		legend.text=element_text(size=11, color="black"),
		legend.title=element_text(size=12, color="black"),
		legend.key=element_blank(),
		aspect.ratio=1)
