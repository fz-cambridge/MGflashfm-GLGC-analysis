basicQC <- function(gwas,mafATGC=0.5) {
# QC for a single gwas
 rmind <- which(nchar(gwas$ALT)>1 | nchar(gwas$REF)>1)
 if(length(rmind)>0) gwas <- gwas[-rmind,]
  
  # rm monomorphic
  rmsnp <- which(gwas$maf==0)
  if(length(rmsnp)>0) gwas <- gwas[-rmsnp,]
  
  #rm non A,T,G,C
  ksnp <- c()
  for(a in c("A","T","G","C")) ksnp <- union(ksnp,grep(a,gwas$REF))
  gwas <- gwas[ksnp,]

# if no rsid, assign chr_pos
ind <- which(is.na(gwas$rsID))
if(length(ind)>0) gwas[ind,"rsID"] <- paste0("chr",gwas[ind,"CHROM"],"_",gwas[ind,"POS_b37"])


 # after removing non-biallelic variants, if there are still duplicate positions in any of the datasets, these are tri-allelic splits and should rm
  tmpqt <- gwas
  dsnp <- unique(tmpqt$rsID[which(duplicated(tmpqt$rsID))]) #  duplicates
  if(length(dsnp)>0){
   for(s in dsnp) {
    ind <- which(tmpqt$rsID==s)
    tmpqt <- tmpqt[-ind,]     
   }
  }
  gwas <- tmpqt #  duplicate positions removed
 
maf <- gwas$maf 
rmsnp <- which((gwas$REF == "A" & gwas$ALT == "T" ) | (gwas$REF == "G" & gwas$ALT == "C" ))
rmsnp <- union(rmsnp,which((gwas$REF == "T" & gwas$ALT == "A" ) | (gwas$REF == "C" & gwas$ALT == "G" )))
if(length(rmsnp)>0) {
   rmind <- which(maf[rmsnp]>mafATGC) 
   if(length(rmind)>0) gwas <- gwas[-rmsnp[rmind],]   
   }

  return(gwas)
 }



flipgwas <- function(gwas,RP) {

flip1 <- which(gwas$ALT != RP$allele1 | gwas$REF != RP$allele2) # check for discrepancies
 check <- c()
 if(length(flip1)>0) check <- which( gwas$ALT[flip1] == RP$allele2[flip1] & gwas$REF[flip1] == RP$allele1[flip1] )  
 if(length(check)>0){
   gwas$REF[flip1[check]] <- RP$allele2[flip1[check]]
   gwas$ALT[flip1[check]] <- RP$allele1[flip1[check]]
   gwas$EFFECT_SIZE[flip1[check]] <- -gwas$EFFECT_SIZE[flip1[check]]
   gwas$POOLED_ALT_AF[flip1[check]] <- 1-gwas$POOLED_ALT_AF[flip1[check]]
 }

 # rm these snps from both datasets since allele codings still don't agree if flipped 
 indrm <-  setdiff(flip1,flip1[check]) 

return(list(gwas,flip1[check],indrm))
}

