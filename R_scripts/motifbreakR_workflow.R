library(motifbreakR)
library(BSgenome)
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library("XtraSNPlocs.Hsapiens.dbSNP144.GRCh37")
library(BSgenome.Hsapiens.UCSC.hg19)
library(MotifDb)
library(BiocParallel)
library(stringr)
library(rtracklayer)

SnpOnMotif <- function(snp_ids, r_scripts, output_dir, motif_dbname="JASPAR",
                       species="Hsapiens", filter_threshold=1e-3,
                       add_motif=F, motif_pwm=""){
  # remove non-rsID SNPS
  snp.ids  <- unique(snp_ids[grepl(snp_ids, pattern="rs")])
  
  all_snps <- SNPlocs.Hsapiens.dbSNP144.GRCh37
  snps.mb <- snps.from.rsid(rsid=snp.ids, dbSNP=SNPlocs.Hsapiens.dbSNP144.GRCh37,
                            search.genome=BSgenome.Hsapiens.UCSC.hg19)
  data(motifbreakR_motif)
  source(paste0(r_scripts, "/", "mitf_motif_break.R"))
  # add extra motifs
  if(add_motif == T) addMitfMotif(motifbreakR_motif, motif_pwm)
  
  motifs <- query(motifbreakR_motif, queryString=motif_dbname, ignore.case=T)
  hsapiens_mitf <- query(motifs, queryString=species)
  
  results <- motifbreakR(snpList=snps.all,
                         pwmList=hsapiens_mitf,
                         filterp=TRUE,
                         threshold=filter_threshold, 
                         method="ic",
                         bkg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                         BPPARAM = BiocParallel::bpparam())
  
  # harangue results GRangess into a dataframe to write out
  res.df <- mcols(results)
  return(res.df)
}