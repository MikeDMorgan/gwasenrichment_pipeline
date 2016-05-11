##########################################################
# What is the effect of RHC SNPs on MITF binding motifs? #
##########################################################

addMitfMotif <- function(MotifDb, pwm.file, ...){
  # need to put just the MITF pwm in as a MotifList object
  mitf.pwm <- as.matrix(read.table(pwm.file, h=F, sep="\t"))
  dimnames(mitf.pwm)[[2]] <- seq_along(1:dim(mitf.pwm)[2])
  dimnames(mitf.pwm)[[1]] <- c("A", "C", "T", "G")
  
  
  mitf_meta <- DataFrame(list(providerName="MA0320.1", providerId="JASPAR_CORE_2014",
                              dataSource="JASPAR_", geneSymbol="MITF", "geneId"="NA",
                              geneIdType="Symbol", proteinId="ENSP00000295600", proteinIdType="ensembl",
                              organism="Hsapiens", sequenceCount="NA", bindingSequence="NA",
                              bindingDomain="NA", tfFamily="bHLH", experimentType="NA",
                              pubmedID="NA"))
  all_mitf <- rbind(attributes(MotifDb)$elementMetadata,
                    mitf_meta)
  attributes(MotifDb)$elementMetadata <- all_mitf
  attributes(MotifDb)$listData[["Hsapiens-JASPAR_CORE_2014-MITF_custom"]] <- mitf.pwm
}
