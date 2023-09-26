#Script for plotting the distance between SV and SNP-based adaptive landscapes
#Based on the tutorial by Capblancq and Forester (2021): https://github.com/Capblancq/RDA-landscape-genomics
#Tuomas Hämälä, August 2023

library(vegan)
library(raster)
library(ggplot2)
library(sf)
library(rgdal)
library(geodata)
library(giscoR)

#Function for imputing missing allele frequencies
impute <- function(x){
	x <- as.numeric(x)
	tryCatch({
		u <- mean(x,na.rm=T)
		s <- var(x,na.rm=T)
		a <- (u*(1-u)/s-1)*u
		b <- (1-u)*(u*(1-u)/s-1)
		n <- sum(is.na(x))
		x[is.na(x)] <- rbeta(n,a,b)
	},warning=function(w){x[is.na(x)] <<- median(x,na.rm=T)})
	x
}

#Edges of the Cochlearia range
top <- 70.687081
bottom <- 41.423610
right <- 22.573357
left <- -25.827484
top <- 70.687081
extent <- c(left,right,bottom,top)

#Coordinates of Cochlearia populations and allele frequencies for outlier SVs and SNPs
coords <- read.table("cochearia_pop_coords.txt", header=T)
freqs_sv <- read.table("comb_vg_outl.freq")
freqs_sv <- apply(freqs_sv,1,impute)
freqs_snp <- read.table("comb_snp_outl.freq")
freqs_snp <- apply(freqs_snp,1,impute)

#Get current bioclim data
bc <- worldclim_global(var="bio", res=2.5, path="/Users/tuomashamala/Documents/R/", version="2.1")
names(bc) <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
current <- extract(bc, data.frame(x=coords$lon,y=coords$lat))
current <- scale(current)
scale_env <- attr(current,'scaled:scale') #Mean and SD used to scale both current and future range-wide variables
center_env <- attr(current,'scaled:center')
current <- cbind.data.frame(pop=coords$pop,current)
current_range <- crop(bc,extent)

#Get future bioclim data
cmip <- cmip6_world(model="MRI-ESM2-0",ssp="370",time="2061-2080",var="bioc",res=2.5,path="/Users/tuomashamala/Documents/R/")
names(cmip) <- c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
future_range <- crop(cmip,extent)

#Run RDA on SV and SNP outliers
RDA_sv <- rda(freqs_sv ~ bio2+bio3+bio9+bio15+bio18+bio8+bio10+bio7+bio1+bio5+bio11, current)
RDA_snp <- rda(freqs_snp ~ bio2+bio3+bio9+bio15+bio18+bio8+bio10+bio7+bio1+bio5+bio11, current)

#Do prediction using current climate data and first two RDA axis
var_env_proj <- as.data.frame(scale(current_range[[row.names(RDA_sv$CCA$biplot)]], center_env[row.names(RDA_sv$CCA$biplot)], scale_env[row.names(RDA_sv$CCA$biplot)]))
proj_sv <- list()
proj_snp <- list()
K <- 2
for(i in 1:K){
  ras_sv <- raster(current_range)
  ras_sv[!is.na(ras_sv)] <- as.vector(apply(var_env_proj[,names(RDA_sv$CCA$biplot[,i])], 1, function(x) sum(x*RDA_sv$CCA$biplot[,i])))
  proj_sv[[i]] <- ras_sv
  names(ras_sv) <- paste0("RDA_sv_", as.character(i))
  names(proj_sv)[i] <- paste0("RDA", as.character(i))
  ras_snp <- raster(current_range)
  ras_snp[!is.na(ras_snp)] <- as.vector(apply(var_env_proj[,names(RDA_snp$CCA$biplot[,i])], 1, function(x) sum(x*RDA_snp$CCA$biplot[,i])))
  proj_snp[[i]] <- ras_snp
  names(ras_snp) <- paste0("RDA_snp_", as.character(i))
  names(proj_snp)[i] <- paste0("RDA", as.character(i))
}

weights_sv <- RDA_sv$CCA$eig/sum(RDA_sv$CCA$eig)
weights_snp <- RDA_snp$CCA$eig/sum(RDA_snp$CCA$eig)
  
proj_offset_sv <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(proj_sv[[x]])[,-c(1,2)]))
proj_offset_sv <- as.data.frame(do.call(cbind, lapply(1:K, function(x) proj_offset_sv[,x]*weights_sv[x])))
proj_offset_snp <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(proj_snp[[x]])[,-c(1,2)]))
proj_offset_snp <- as.data.frame(do.call(cbind, lapply(1:K, function(x) proj_offset_snp[,x]*weights_snp[x])))

ras <- proj_offset[[1]]
ras[!is.na(ras)] <- unlist(lapply(1:nrow(proj_offset_sv), function(x) dist(rbind(proj_offset_sv[x,], proj_offset_snp[x,]), method = "euclidean")))
names(ras) <- "offset"
current_offset <- ras

#Do prediction using future climate data and first two RDA axis
var_env_proj <- as.data.frame(scale(future_range[[row.names(RDA_sv$CCA$biplot)]], center_env[row.names(RDA_sv$CCA$biplot)], scale_env[row.names(RDA_sv$CCA$biplot)]))
proj_sv <- list()
proj_snp <- list()
K <- 2
for(i in 1:K){
  ras_sv <- raster(future_range)
  ras_sv[!is.na(ras_sv)] <- as.vector(apply(var_env_proj[,names(RDA_sv$CCA$biplot[,i])], 1, function(x) sum(x*RDA_sv$CCA$biplot[,i])))
  proj_sv[[i]] <- ras_sv
  names(ras_sv) <- paste0("RDA_sv_", as.character(i)) 
  names(proj_sv)[i] <- paste0("RDA", as.character(i))
  ras_snp <- raster(future_range)
  ras_snp[!is.na(ras_snp)] <- as.vector(apply(var_env_proj[,names(RDA_snp$CCA$biplot[,i])], 1, function(x) sum(x*RDA_snp$CCA$biplot[,i])))
  proj_snp[[i]] <- ras_snp
  names(ras_snp) <- paste0("RDA_snp_", as.character(i))
  names(proj_snp)[i] <- paste0("RDA", as.character(i))
}

weights_sv <- RDA_sv$CCA$eig/sum(RDA_sv$CCA$eig)
weights_snp <- RDA_snp$CCA$eig/sum(RDA_snp$CCA$eig)
  
proj_offset_sv <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(proj_sv[[x]])[,-c(1,2)]))
proj_offset_sv <- as.data.frame(do.call(cbind, lapply(1:K, function(x) proj_offset_sv[,x]*weights_sv[x])))
proj_offset_snp <- do.call(cbind, lapply(1:K, function(x) rasterToPoints(proj_snp[[x]])[,-c(1,2)]))
proj_offset_snp <- as.data.frame(do.call(cbind, lapply(1:K, function(x) proj_offset_snp[,x]*weights_snp[x])))

ras <- proj_offset[[1]]
ras[!is.na(ras)] <- unlist(lapply(1:nrow(proj_offset_sv), function(x) dist(rbind(proj_offset_sv[x,], proj_offset_snp[x,]), method = "euclidean")))
names(ras) <- "offset"
future_offset <-  ras

#100 x 100 km grid points across Europe. Downloaded from https://www.eea.europa.eu/data-and-maps/data/eea-reference-grids-2/gis-files/europe-10-km-100-km
shape <- read_sf(dsn="/Users/tuomashamala/Documents/R/Europe/", layer="europe_100km") %>% st_transform(crs = 31370)
#Cleaned occurence data from GBIF
coords <- read.table("cochleria_occurence_coordinates_clean.txt", header=T) 
coch_spatial <- st_as_sf(coords,coords=c("decimalLongitude","decimalLatitude"),crs = 4326) %>% st_transform(crs = 31370)
mask <- st_join(shape, coch_spatial, left = FALSE, largest = TRUE)
current_mask <- mask
future_mask <- mask

#Current grid points
ex <- extract(current_offset, mask)
means <- sapply(ex, mean, na.rm=T)
current_mask$means <- means
current_mask <- na.omit(current_mask)
current_mask$year <- "1970 – 2000"

#Future grid points
ex <- extract(future_offset, mask)
means <- sapply(ex, mean, na.rm=T)
future_mask$means <- means
future_mask <- na.omit(future_mask)
future_mask$year <- "2061 – 2080"

#Combine data and normalise the adaptive distance from 0 to 1
comb_mask <- rbind(current_mask, future_mask)
comb_mask$means <- (comb_mask$means-min(comb_mask$means))/(max(comb_mask$means)-min(comb_mask$means))

#Country borders
svet <- gisco_get_countries(resolution = "20")

#Plot map
ggplot()+
	geom_sf(data=svet,fill=NA,color="grey30")+
	geom_sf(data=comb_mask, aes(fill=means),inherit.aes=FALSE,alpha=0.8)+
	scale_fill_gradient2(low="grey100",high="red2",mid="grey90",midpoint=quantile(comb_mask$means,0.25),name="SNP–SV\ndistance")+
	facet_grid(~year)+
	coord_sf(crs=st_crs(31370),xlim=c(-1285553,1261053),ylim=c(-851506.9,2554251))+
	theme(panel.background=element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_rect(colour="black",fill=NA,size=0.5),
		axis.ticks=element_blank(),
		axis.text=element_blank(),
		axis.title=element_blank(),
		strip.text.x=element_text(size=14,color="black"),
		strip.background =element_rect(fill="white"),
		legend.text=element_text(size=11,color="black"),
		legend.title=element_text(size=12,color="black"))
		
#Biplot of the loadings
load_sv <- as.data.frame(scores(RDA_sv, choices=c(1:2), display="bp"))
load_sv$dat <- "SV"
load_snp <- as.data.frame(scores(RDA_snp, choices=c(1:2), display="bp"))
load_snp$dat <- "SNP"
load <- rbind(load_sv, load_snp)

ord <- factor(load$dat, levels=c("SV","SNP"))
ggplot(load, aes(x=0,y=0,colour=ord))+
	geom_hline(yintercept=0,linetype="dashed",color="grey60",size=0.5) +
 	geom_vline(xintercept=0,linetype="dashed",color="grey60",size=0.5) +
	geom_segment(aes(xend=RDA1,yend=RDA2),size=0.75,linetype=1,arrow=arrow(length=unit(0.02,"npc")))+
	geom_text(aes(x=1.1*RDA1,y=1.1*RDA2,label=c(rownames(load_sv),rownames(load_snp))),size=4)+
	scale_colour_manual(values=c("red2","steelblue2"))+
	labs(x="RDA1",y="RDA2")+
	theme(panel.background = element_blank(),
		panel.grid.major=element_blank(), 
		panel.grid.minor=element_blank(), 
		panel.border=element_blank(),
		axis.line=element_line(color="black",size=0.5),
		axis.text=element_text(size=11,color="black"),
		axis.ticks.length=unit(.15,"cm"),
		axis.ticks=element_line(color="black",size=0.5),
		axis.title=element_text(size=12,color="black"),
		plot.title=element_text(size=14,color="black",hjust = 0.5),
		legend.text=element_text(size=11,color="black"),
		legend.title=element_blank(),
		legend.key=element_blank(),
		aspect.ratio=1)
