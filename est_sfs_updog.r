#Script for estimating SFS and Tajima's D from genotype probabilities.
#Tuomas Hämälä, April 2023

library(updog)

#Functions for folding the SFS and estimating Tajima's D
fold_sfs <- function(sfs){
	out <- NA
	j <- 1
	for(i in 1:(length(sfs-1))){
		if(i <= length(sfs)/2){
			if(i == length(sfs)-i) out[j] <- sfs[i]
			else out[j] <- sfs[i] + sfs[length(sfs)-i]
			j <- j + 1
		}
	}
	out
}

est_var <- function(n, S){
	a1 <- sum(1/1:(n-1))
    a2 <- sum(1/(1:(n-1))^2)
    b1 <- (n+1)/(3.0*(n-1))
    b2 <- 2*((n*n)+n+3.0)/(9*n*(n-1))
    c1 <- b1-1/a1
    c2 <- b2-(n+2)/(a1*n)+(a2/(a1*a1))
    e1 <- c1/a1
    e2 <- c2/((a1*a1)+a2)
    sqrt(e1*S+e2*S*(S-1))
}

est_tajD <- function(sfs){
	n <- length(sfs)
	S <- sum(sfs[-n])
	tW <- S/sum(1/1:(n-1))
	w <- seq(1,n)/n
	tP <- n/(n-1)*2*sum(sfs*w*(1-w))
	(tP-tW)/est_var(n,S)
}

boot_tajD <- function(sfs,n=1000){
	r <- rmultinom(n,sum(sfs),sfs)
	boots <- apply(r,2,est_tajD)
	c(est_tajD(sfs),quantile(boots, c(0.025,0.975))) 
}

#Generate data for Updog
df <- read.table("tetraploids_exon_SV.vcf", header=T) #Filtered VCF file from Sniffles
n <- nrow(df)
c <- ncol(df)
sv_mat <- matrix(nrow=n,ncol=c-9)
size_mat <- matrix(nrow=n,ncol=c-9)
colnames(sv_mat) <- names(df)[10:ncol(df)]
rownames(sv_mat) <- df$ID[1:n]
colnames(size_mat) <- names(df)[10:ncol(df)]
rownames(size_mat) <- df$ID[1:n]
for(i in 1:n){
	if(i==1|i==n|i%%1000==0) cat("Line",i,"/",n,"\n")
	for(j in 10:c){
		temp <- unlist(strsplit(df[i,j],":"))
		sv <- as.numeric(temp[4])
		ref <- as.numeric(temp[3])
		sv_mat[i,j-9] <- sv
		size_mat[i,j-9] <- sv + ref
	}
}

#Use Updog to estimate genotype probabilities
x <- 4
mout <- multidog(refmat=sv_mat, sizemat=size_mat, ploidy=x, model="norm", nc=4)
geno <- mout$inddf

#Estimate SFS assuming a minimum coverage of 10, 20, or 40x. 
#Sites with 20% missing data are excluded and the missing genotypes imputed by drawing them from a binomial distribution
geno$postmean <- ifelse(geno$size==0,NA,geno$postmean)
m <- length(unique(geno$ind))
min10 <- rep(0,m*x)
min20 <- rep(0,m*x)
min40 <- rep(0,m*x)
list <- unique(geno$snp)
n <- length(list)
for(i in 1:n){
	if(i==1|i==n|i%%100==0) cat("Line",i,"/",n,"\n")
	temp <- subset(geno, snp == list[i])

	temp$postmean <- ifelse(temp$size < 10,NA,temp$postmean)
	mis <- sum(is.na(temp$postmean))
	if(mis/nrow(temp) > 0.2) next
	sum <- round(sum(temp$postmean,na.rm=T))
	if(mis > 0) sum <- sum + sum(rbinom(mis,2,sum/((m-mis)*x)))
	if(sum > 0) min10[sum] <- min10[sum] + 1

	temp$postmean <- ifelse(temp$size < 20,NA,temp$postmean)
	mis <- sum(is.na(temp$postmean))
	if(mis/nrow(temp) > 0.2) next
	sum <- round(sum(temp$postmean,na.rm=T))
	if(mis > 0) sum <- sum + sum(rbinom(mis,2,sum/((m-mis)*x)))
	if(sum > 0) min20[sum] <- min20[sum] + 1

	temp$postmean <- ifelse(temp$size < 40,NA,temp$postmean)
	mis <- sum(is.na(temp$postmean))
	if(mis/nrow(temp) > 0.2) next
	sum <- round(sum(temp$postmean,na.rm=T))
	if(mis > 0) sum <- sum + sum(rbinom(mis,2,sum/((m-mis)*x)))
	if(sum > 0) min40[sum] <- min40[sum] + 1
}

#Example SFS
min20 <- c(2349,621,407,338,178,128,71,85,42,35,41,37,21,16,17,18,10,13,12,9,6,10,7,13,5,5,7,25)

#Folded SFS
fold_sfs(min20)

#Tajima's D, including 95% CI
boot_tajD(min20)

