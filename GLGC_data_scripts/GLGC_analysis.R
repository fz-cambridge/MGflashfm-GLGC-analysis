#Paper: Leveraging information between multiple population groups and traits improves fine-mapping resolution 
#Authors: F Zhou, O Soremekun, T Chikowore, S Fatumo, I Barroso, AP Morris, JL Asimit

#Data Availability: 
#The GLGC lipids traits GWAS summary statistics from five genetically similar groups 
#were downloaded from http://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/ancestry_specific/. 
#Reference panels for LD and LD scores were generated from the 1000 Genomes data 
#available at https://ctg.cncr.nl/software/MAGMA/ref_data/. 

#Code Availability: 
#Our proposed multi-group fine-mapping methods, MGflashfm and MGfm, 
#are freely available as an R library at our online GitHub. This library also includes updated versions 
#of expanded JAM and flashfm that have dynamic selection of the maximum number of causal variants, 
#as learned from the data. Custom code for the analysis of the GLGC data is also available at our GitHub. 

#Requirements:
#library(MGflashfm)
#QCfunctions.R
#input data relating to GLGC: (1) a csv file contains id of 50 regions; (2) 50 .RData files contain all gwas and ld 

#args=(commandArgs(TRUE))
#regno=as.numeric(args[1])  # region number 

##### LOAD LIBRARY ##### ----
library(MGflashfm)
library(R2BGLiMS)

source("QCfunctions.R")

##### LOAD DATA ##### ----
#Please save the csv file "GLGC_input_50regions_chr_pos_index_100k.csv" in the same directory with this script
#Please save the folder "GLGC_input_50regions_gwas_final" that contains all 50 .RData files in the same directory 
Region_chr_pos_index_unique <- read.csv("GLGC_input_50regions_chr_pos_index_100k.csv")

## Main process----
regno=1  #Run region 1. If run all regions, please use for-loop, i.e. for(regno in 1:50){}

ii = regno
print(ii)
  
id_reg = Region_chr_pos_index_unique$Our_region[ii]
id_chr = Region_chr_pos_index_unique$Chromosome[ii]
id_start = Region_chr_pos_index_unique$Position_start[ii]
id_end   = Region_chr_pos_index_unique$Position_end[ii]
  
load(paste0("./GLGC_input_50regions_gwas_final/Input_reg",id_reg,"_chr",id_chr,"_",id_start,"_",id_end,".RData"))

#######################
# identify traits to finemap across all groups
all_minP <- lapply(gwas_MAF_0_005_raw,function(x) min(x$pvalue))
irm <- grep("EUR_nonFinnish_only",names(all_minP)) #remove EUR with nonFinnish
all_minP <- all_minP[-irm]

hdl_minP <- all_minP[grep("HDL",names(all_minP))]
ldl_minP <- all_minP[grep("LDL",names(all_minP))]
tg_minP <- all_minP[grep("logTG",names(all_minP))]
tc_minP <- all_minP[grep("TC",names(all_minP))]

check <- lapply(list(hdl_minP,ldl_minP,tg_minP,tc_minP),function(x) sum(x<5E-8))
names(check) <- c("HDL","LDL","logTG","TC")
selqt <- which(check>1)

allgwas <- lapply(gwas_MAF_0_005_raw[-irm], basicQC, mafATGC=0.5)  # no filter on maf ATGC	
selgwas <-vector("list",5*length(selqt))
for(j in 1:length(selqt)) {
 nn <- grep(names(selqt[j]),names(allgwas))
 selgwas[(5*(j-1)+1):(5*j)] <- allgwas[nn]
 names(selgwas)[(5*(j-1)+1):(5*j)] <- names(allgwas)[nn]
}

group = c("SAS","HIS","EAS","AFR","EUR")
all_snpinfo <- all_LD <- all_raf <- vector("list",length(group))
names(all_snpinfo) <- names(all_LD) <- names(all_raf) <-  group 
rafgwas <- vector("list",length(group))
names(rafgwas) <- group

##### Reference panels and intersect with gwas and flip if needed
# 1KG HIS, SAS
tmp <- snpinfo_bed_file[grep("MAF_0_001",names(snpinfo_bed_file))]
all_snpinfo[["HIS"]] <- tmp[[grep("HIS",names(tmp))]]
all_snpinfo[["SAS"]] <- tmp[[grep("SAS",names(tmp))]]
all_snpinfo[["EAS"]] <- tmp[[grep("EAS",names(tmp))]]
all_snpinfo[["AFR"]] <- tmp[[grep("AFR",names(tmp))]]
all_snpinfo[["EUR"]] <- tmp[[grep("EUR",names(tmp))]]
rownames(all_snpinfo[["HIS"]]) <- all_snpinfo[["HIS"]]$snpid
rownames(all_snpinfo[["SAS"]]) <- all_snpinfo[["SAS"]]$snpid
rownames(all_snpinfo[["EAS"]]) <- all_snpinfo[["EAS"]]$snpid
rownames(all_snpinfo[["AFR"]]) <- all_snpinfo[["AFR"]]$snpid
rownames(all_snpinfo[["EUR"]]) <- all_snpinfo[["EUR"]]$snpid
rm(tmp)
all_raf[["HIS"]] <- raf_bed_file[["AMR_HIS_1KG_MAF_0_001"]] # wrt allele2
all_raf[["SAS"]] <- raf_bed_file[["SAS_1KG_MAF_0_001"]] # wrt allele2
all_raf[["EAS"]] <- raf_bed_file[["EAS_1KG_MAF_0_001"]] # wrt allele2
all_raf[["AFR"]] <- raf_bed_file[["AFR_1KG_MAF_0_001"]] # wrt allele2
all_raf[["EUR"]] <- raf_bed_file[["EUR_1KG_MAF_0_001"]] # wrt allele2
all_LD[["HIS"]] <- LD_corX_1KG_UKBB[["AMR_HIS_1KG_MAF_0_001"]] # wrt allele2
all_LD[["SAS"]] <- LD_corX_1KG_UKBB[["SAS_1KG_MAF_0_001"]] # wrt allele2
all_LD[["EAS"]] <- LD_corX_1KG_UKBB[["EAS_1KG_MAF_0_001"]] # wrt allele2
all_LD[["AFR"]] <- LD_corX_1KG_UKBB[["AFR_1KG_MAF_0_001"]] # wrt allele2
all_LD[["EUR"]] <- LD_corX_1KG_UKBB[["EUR_1KG_MAF_0_001"]] # wrt allele2

# HIS snps in all traits and in RP----
HISgwas <- selgwas[grep("HIS",names(selgwas))] 
for(i in 1:length(HISgwas)) rownames(HISgwas[[i]]) <- HISgwas[[i]]$rsID
isnp <- Reduce(intersect,lapply(HISgwas,rownames))
for(i in 1:length(HISgwas)) HISgwas[[i]] <- HISgwas[[i]][isnp,]
ksnp <- intersect(isnp, all_snpinfo[["HIS"]]$snpid) 
for(i in 1:length(HISgwas)) HISgwas[[i]] <- HISgwas[[i]][ksnp,]
all_LD[["HIS"]] <- all_LD[["HIS"]][ksnp,ksnp]
all_raf[["HIS"]] <- all_raf[["HIS"]][ksnp]
all_snpinfo[["HIS"]] <- all_snpinfo[["HIS"]][ksnp,]

#flip
rmind <- c()
out <- vector("list",length(HISgwas))
for(i in 1:length(HISgwas)) {
 out[[i]] <- flipgwas(gwas=HISgwas[[i]],RP=all_snpinfo[["HIS"]])
 rmind <- union(rmind,out[[i]][[3]]) # any variant indices that can't be flipped to match RP 
 HISgwas[[i]] <- out[[i]][[1]] # flipped alleles if not matching RP
} 
if(length(rmind)>0) {
 for(i in 1:length(HISgwas)) HISgwas[[i]] <- HISgwas[[i]][-rmind,]
}
sel <- which.max(lapply(HISgwas,function(x) max(x$N)))
rafgwas[["HIS"]] <- HISgwas[[sel]]$POOLED_ALT_AF
names(rafgwas[["HIS"]]) <- HISgwas[[sel]]$rsID
rmind <- which(abs(rafgwas[["HIS"]]-all_raf[["HIS"]])>0.1)
if(length(rmind)>0) {
 for(i in 1:length(HISgwas)) HISgwas[[i]] <- HISgwas[[i]][-rmind,]
}
#snps that have prop > 0.8 for each trait 
HISsnps <- vector("list",length(HISgwas))
for(j in 1:length(HISgwas)) {
 Nprop <- HISgwas[[j]]$N/max(HISgwas[[j]]$N)
 HISsnps[[j]] <- HISgwas[[j]]$rsID[which(Nprop >= 0.8)]
}

# SAS snps in all traits and in RP----
SASgwas <- selgwas[grep("SAS",names(selgwas))] 
for(i in 1:length(SASgwas)) rownames(SASgwas[[i]]) <- SASgwas[[i]]$rsID
isnp <- Reduce(intersect,lapply(SASgwas,rownames))
for(i in 1:length(SASgwas)) SASgwas[[i]] <- SASgwas[[i]][isnp,]
ksnp <- intersect(isnp, all_snpinfo[["SAS"]]$snpid) 
for(i in 1:length(SASgwas)) SASgwas[[i]] <- SASgwas[[i]][ksnp,]
all_LD[["SAS"]] <- all_LD[["SAS"]][ksnp,ksnp]
all_raf[["SAS"]] <- all_raf[["SAS"]][ksnp]
all_snpinfo[["SAS"]] <- all_snpinfo[["SAS"]][ksnp,]

#flip
rmind <- c()
out <- vector("list",length(SASgwas))
for(i in 1:length(SASgwas)) {
 out[[i]] <- flipgwas(gwas=SASgwas[[i]],RP=all_snpinfo[["SAS"]])
 rmind <- union(rmind,out[[i]][[3]]) # any variant indices that can't be flipped to match RP 
 SASgwas[[i]] <- out[[i]][[1]] # flipped alleles if not matching RP
} 
if(length(rmind)>0) {
 for(i in 1:length(SASgwas)) SASgwas[[i]] <- SASgwas[[i]][-rmind,]
}
sel <- which.max(lapply(SASgwas,function(x) max(x$N)))
rafgwas[["SAS"]] <- SASgwas[[sel]]$POOLED_ALT_AF
names(rafgwas[["SAS"]]) <- SASgwas[[sel]]$rsID
rmind <- which(abs(rafgwas[["SAS"]]-all_raf[["SAS"]])>0.1)
if(length(rmind)>0) {
 for(i in 1:length(SASgwas)) SASgwas[[i]] <- SASgwas[[i]][-rmind,]
}
#snps that have prop > 0.8 for each trait 
SASsnps <- vector("list",length(SASgwas))
for(j in 1:length(SASgwas)) {
 Nprop <- SASgwas[[j]]$N/max(SASgwas[[j]]$N)
 SASsnps[[j]] <- SASgwas[[j]]$rsID[which(Nprop >= 0.8)]
}

# EAS snps in all traits and in RP----
EASgwas <- selgwas[grep("EAS",names(selgwas))] 
for(i in 1:length(EASgwas)) rownames(EASgwas[[i]]) <- EASgwas[[i]]$rsID
isnp <- Reduce(intersect,lapply(EASgwas,rownames))
for(i in 1:length(EASgwas)) EASgwas[[i]] <- EASgwas[[i]][isnp,]
ksnp <- intersect(isnp, all_snpinfo[["EAS"]]$snpid) 
for(i in 1:length(EASgwas)) EASgwas[[i]] <- EASgwas[[i]][ksnp,]
all_LD[["EAS"]] <- all_LD[["EAS"]][ksnp,ksnp]
all_raf[["EAS"]] <- all_raf[["EAS"]][ksnp]
all_snpinfo[["EAS"]] <- all_snpinfo[["EAS"]][ksnp,]

#flip
rmind <- c()
out <- vector("list",length(EASgwas))
for(i in 1:length(EASgwas)) {
  out[[i]] <- flipgwas(gwas=EASgwas[[i]],RP=all_snpinfo[["EAS"]])
  rmind <- union(rmind,out[[i]][[3]]) # any variant indices that can't be flipped to match RP 
  EASgwas[[i]] <- out[[i]][[1]] # flipped alleles if not matching RP
} 
if(length(rmind)>0) {
  for(i in 1:length(EASgwas)) EASgwas[[i]] <- EASgwas[[i]][-rmind,]
}
sel <- which.max(lapply(EASgwas,function(x) max(x$N)))
rafgwas[["EAS"]] <- EASgwas[[sel]]$POOLED_ALT_AF
names(rafgwas[["EAS"]]) <- EASgwas[[sel]]$rsID
rmind <- which(abs(rafgwas[["EAS"]]-all_raf[["EAS"]])>0.1)
if(length(rmind)>0) {
  for(i in 1:length(EASgwas)) EASgwas[[i]] <- EASgwas[[i]][-rmind,]
}
#snps that have prop > 0.8 for each trait 
EASsnps <- vector("list",length(EASgwas))
for(j in 1:length(EASgwas)) {
  Nprop <- EASgwas[[j]]$N/max(EASgwas[[j]]$N)
  EASsnps[[j]] <- EASgwas[[j]]$rsID[which(Nprop >= 0.8)]
}

# AFR snps in all traits and in RP----
AFRgwas <- selgwas[grep("AFR",names(selgwas))] 
for(i in 1:length(AFRgwas)) rownames(AFRgwas[[i]]) <- AFRgwas[[i]]$rsID
isnp <- Reduce(intersect,lapply(AFRgwas,rownames))
for(i in 1:length(AFRgwas)) AFRgwas[[i]] <- AFRgwas[[i]][isnp,]
ksnp <- intersect(isnp, all_snpinfo[["AFR"]]$snpid) 
for(i in 1:length(AFRgwas)) AFRgwas[[i]] <- AFRgwas[[i]][ksnp,]
all_LD[["AFR"]] <- all_LD[["AFR"]][ksnp,ksnp]
all_raf[["AFR"]] <- all_raf[["AFR"]][ksnp]
all_snpinfo[["AFR"]] <- all_snpinfo[["AFR"]][ksnp,]

#flip
rmind <- c()
out <- vector("list",length(AFRgwas))
for(i in 1:length(AFRgwas)) {
  out[[i]] <- flipgwas(gwas=AFRgwas[[i]],RP=all_snpinfo[["AFR"]])
  rmind <- union(rmind,out[[i]][[3]]) # any variant indices that can't be flipped to match RP 
  AFRgwas[[i]] <- out[[i]][[1]] # flipped alleles if not matching RP
} 
if(length(rmind)>0) {
  for(i in 1:length(AFRgwas)) AFRgwas[[i]] <- AFRgwas[[i]][-rmind,]
}
sel <- which.max(lapply(AFRgwas,function(x) max(x$N)))
rafgwas[["AFR"]] <- AFRgwas[[sel]]$POOLED_ALT_AF
names(rafgwas[["AFR"]]) <- AFRgwas[[sel]]$rsID
rmind <- which(abs(rafgwas[["AFR"]]-all_raf[["AFR"]])>0.1)
if(length(rmind)>0) {
  for(i in 1:length(AFRgwas)) AFRgwas[[i]] <- AFRgwas[[i]][-rmind,]
}
#snps that have prop > 0.8 for each trait 
AFRsnps <- vector("list",length(AFRgwas))
for(j in 1:length(AFRgwas)) {
  Nprop <- AFRgwas[[j]]$N/max(AFRgwas[[j]]$N)
  AFRsnps[[j]] <- AFRgwas[[j]]$rsID[which(Nprop >= 0.8)]
}

# EUR snps in all traits and in RP----
EURgwas <- selgwas[grep("EUR",names(selgwas))] 
for(i in 1:length(EURgwas)) rownames(EURgwas[[i]]) <- EURgwas[[i]]$rsID
isnp <- Reduce(intersect,lapply(EURgwas,rownames))
for(i in 1:length(EURgwas)) EURgwas[[i]] <- EURgwas[[i]][isnp,]
ksnp <- intersect(isnp, all_snpinfo[["EUR"]]$snpid) 
for(i in 1:length(EURgwas)) EURgwas[[i]] <- EURgwas[[i]][ksnp,]
all_LD[["EUR"]] <- all_LD[["EUR"]][ksnp,ksnp]
all_raf[["EUR"]] <- all_raf[["EUR"]][ksnp]
all_snpinfo[["EUR"]] <- all_snpinfo[["EUR"]][ksnp,]

#flip
rmind <- c()
out <- vector("list",length(EURgwas))
for(i in 1:length(EURgwas)) {
  out[[i]] <- flipgwas(gwas=EURgwas[[i]],RP=all_snpinfo[["EUR"]])
  rmind <- union(rmind,out[[i]][[3]]) # any variant indices that can't be flipped to match RP 
  EURgwas[[i]] <- out[[i]][[1]] # flipped alleles if not matching RP
} 
if(length(rmind)>0) {
  for(i in 1:length(EURgwas)) EURgwas[[i]] <- EURgwas[[i]][-rmind,]
}
sel <- which.max(lapply(EURgwas,function(x) max(x$N)))
rafgwas[["EUR"]] <- EURgwas[[sel]]$POOLED_ALT_AF
names(rafgwas[["EUR"]]) <- EURgwas[[sel]]$rsID
rmind <- which(abs(rafgwas[["EUR"]]-all_raf[["EUR"]])>0.1)
if(length(rmind)>0) {
  for(i in 1:length(EURgwas)) EURgwas[[i]] <- EURgwas[[i]][-rmind,]
}
#snps that have prop > 0.8 for each trait 
EURsnps <- vector("list",length(EURgwas))
for(j in 1:length(EURgwas)) {
  Nprop <- EURgwas[[j]]$N/max(EURgwas[[j]]$N)
  EURsnps[[j]] <- EURgwas[[j]]$rsID[which(Nprop >= 0.8)]
}




##### LOAD INPUT to MSflashfm key functions ##### ----
M <- length(EURgwas)
Neur <- numeric(M); for(i in 1:M) {mm <- which(EURgwas[[i]]$maf>.01);Neur[[i]] <- flashfm::Neff(EURgwas[[i]][mm,"POOLED_ALT_AF"],EURgwas[[i]][mm,"SE"])}
Nafr <- numeric(M); for(i in 1:M) {mm <- which(AFRgwas[[i]]$maf>.01);Nafr[[i]] <- flashfm::Neff(AFRgwas[[i]][mm,"POOLED_ALT_AF"],AFRgwas[[i]][mm,"SE"])}
Neas <- numeric(M); for(i in 1:M) {mm <- which(EASgwas[[i]]$maf>.01);Neas[[i]] <- flashfm::Neff(EASgwas[[i]][mm,"POOLED_ALT_AF"],EASgwas[[i]][mm,"SE"])}
Nsas <- numeric(M); for(i in 1:M) {mm <- which(SASgwas[[i]]$maf>.01);Nsas[[i]] <- flashfm::Neff(SASgwas[[i]][mm,"POOLED_ALT_AF"],SASgwas[[i]][mm,"SE"])}
Nhis <- numeric(M); for(i in 1:M) {mm <- which(HISgwas[[i]]$maf>.01);Nhis[[i]] <- flashfm::Neff(HISgwas[[i]][mm,"POOLED_ALT_AF"],HISgwas[[i]][mm,"SE"])}

afrsnps <- Reduce(intersect,AFRsnps)
hissnps <- Reduce(intersect,HISsnps)
sassnps <- Reduce(intersect,SASsnps)
eassnps <- Reduce(intersect,EASsnps)
eursnps <- Reduce(intersect,EURsnps)
allsnps <- Reduce(union,list(afrsnps,hissnps,sassnps,eassnps,eursnps))

for(i in 1:M) {
  AFRgwas[[i]] <- AFRgwas[[i]][intersect(afrsnps,rownames(AFRgwas[[i]])),]
  EURgwas[[i]] <- EURgwas[[i]][intersect(eursnps,rownames(EURgwas[[i]])),]
  EASgwas[[i]] <- EASgwas[[i]][intersect(eassnps,rownames(EASgwas[[i]])),]
  SASgwas[[i]] <- SASgwas[[i]][intersect(sassnps,rownames(SASgwas[[i]])),]
  HISgwas[[i]] <- HISgwas[[i]][intersect(hissnps,rownames(HISgwas[[i]])),]
}

if (!dir.exists("DIRtmp")){
  dir.create("DIRtmp")
}else{
  print("DIRtmp exists in this directory, please go ahead!")
}
save.path=paste0("DIRtmp/reg-",regno)
dir.create(save.path)


#> names(all_LD)
#[1] "SAS" "HIS" "EAS" "AFR" "EUR"

traits <- sapply(strsplit(names(EURgwas),"_"),"[[",1)

#A <- 5
gwas.list <- list(SAS=SASgwas,HIS=HISgwas,EAS=EASgwas,AFR=AFRgwas,EUR=EURgwas)
A <- length(gwas.list)
for(i in 1:A) {
 for(j in 1:M) {
  gwas.list[[i]][[j]] <- gwas.list[[i]][[j]][,c("rsID","EFFECT_SIZE","POOLED_ALT_AF")]
  names(gwas.list[[i]][[j]]) <- c("rsID","beta","EAF")
 }
 names(gwas.list[[i]]) <- traits
}

Nall <- list(Nsas,Nhis,Neas,Nafr,Neur)
names(Nall) <- names(all_LD)

covY.list <- covY[names(all_LD)]
for(i in 1:A) {
 covY.list[[i]] <- covY.list[[i]][traits,traits]
}
#names(covY.list) <- names(all_LD)

#Run MSflashfm key functions
mgmtCS <- MGFLASHFMwithJAM(gwas.list,
                           LD.list=all_LD,
                           covY.list,
                           Nall,
                           multi=TRUE,
                           TOdds=1,
                           maxcv=1,
                           maxcv_stop = 20,
                           maxcv_autocheck = TRUE,
                           save.path,
                           cpp=0.99,
                           cred=0.99,
                           NCORES=M,    #NCORES=1 if using Windows machine
                           jam.nM.iter=1,
                           flashfmRET=TRUE)

mgCS <- MGflashfmRET(gwas.list,
                     flashfm.list=mgmtCS$flashfm.out,
                     Nall,
                     cred=0.99,
                     multi=FALSE,
                     cpp=0.99,
                     NCORES=M    #NCORES=1 if using Windows machine
)

#unlink("save.path/*")
unlink(paste0(save.path,"/*"))
save(Nall,mgmtCS, mgCS, file = paste0("Output_new_region",id_reg,"_chr",id_chr,"_",id_start,"_",id_end,".RData"))


