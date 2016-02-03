##########################################################
# Compare SNP sets with enrichment over epigenetic marks #
# with GAT                                               #
##########################################################
library(reshape2)
library(ggplot2)
library(RColorBrewer)

RA.gat <- read.table("/ifs/projects/proj045/pipeline_enrichment_blond/gat.dir/RA_GWAS-haplotype.gat",
                     h=T, sep="\t", stringsAsFactors=F)

LD.gat <- read.table("/ifs/projects/proj045/pipeline_enrichment_blond/gat.dir/Blond_LDranks-haplotype.gat",
                     h=T, sep="\t", stringsAsFactors=F)

PICS.gat <- read.table("/ifs/projects/proj045/pipeline_enrichment_blond/gat.dir/Blond_PICS-haplotype.gat",
                       h=T, sep="\t", stringsAsFactors=F)

ABF.gat <- read.table("/ifs/projects/proj045/pipeline_enrichment_blond/gat.dir/Blond_ABF-haplotype.gat",
                      h=T, sep="\t", stringsAsFactors=F)

Top.gat <- read.table("/ifs/projects/proj045/pipeline_enrichment_blond/gat.dir/Blond_GwasTop-haplotype.gat",
                      h=T, sep="\t", stringsAsFactors=F)

###########
# LD Rank #
###########

LD.melt <- melt(LD.gat, varnames="annotation", id.vars="annotation")
ld.cell_type = unlist(lapply(strsplit(as.character(LD.melt$variable), fixed=T, split="_"),
                             function(x) {paste0(x[1])}))
LD.melt$cell_type <- ld.cell_type
ld.stat = unlist(lapply(strsplit(as.character(LD.melt$variable), fixed=T, split="_"),
                        function(x) {paste0(x[2])}))
LD.melt$stat <- ld.stat
LD.melt$method <- "LDranks"
LD.melt$annotation <- as.factor(LD.melt$annotation)
LD.melt$value <- as.numeric(LD.melt$value)

n_annots <- length(unique(LD.melt$annotation))
n_cells <- length(unique(LD.melt$cell_type))
annot_cols <- colorRampPalette(brewer.pal(8, "Paired"))

LD_qval <- ggplot(na.omit(LD.melt[LD.melt$stat == "qvalue",]), 
                  aes(x=cell_type, y=-log10(value), 
                      fill=cell_type,
                      group=reorder(annotation, log10(value)))) + 
  geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
  ylim(c(0, 5)) + labs(x="Cell Type", y="-log10 q-value") + 
  geom_hline(yintercept=1.5, linetype="dashed", colour="white") + 
  scale_fill_manual(values=annot_cols(n_cells)) 

LD_fold <- ggplot(na.omit(LD.melt[LD.melt$stat == "fold",]),
                  aes(x=cell_type, y=value, 
                      fill=cell_type,
                      group=reorder(annotation, -value))) + 
  geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
  scale_fill_manual(values=annot_cols(n_cells)) + labs(x="Cell Type", y="Fold Enrichment") +
  geom_hline(yintercept=1.0, linetype="dashed", colour="white")

############################
# Approximate Bayes Factor #
############################

ABF.melt <- melt(ABF.gat, id.vars="annotation")
abf.cell_type = unlist(lapply(strsplit(as.character(ABF.melt$variable), fixed=T, split="_"),
                              function(x) {paste0(x[1])}))
ABF.melt$cell_type <- abf.cell_type
abf.stat = unlist(lapply(strsplit(as.character(ABF.melt$variable), fixed=T, split="_"),
                         function(x) {paste0(x[2])}))
ABF.melt$stat <- abf.stat
ABF.melt$method <- "ABF"
ABF.melt$value <- as.numeric(ABF.melt$value)

ABF_qval <- ggplot(na.omit(ABF.melt[ABF.melt$stat == "qvalue",]),
                   aes(x=cell_type, y=-log10(value), 
                       fill=cell_type,
                       group=reorder(annotation, log10(value)))) + 
  geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
  scale_fill_manual(values=annot_cols(n_cells)) + ylim(c(0, 5)) + 
  labs(x="Cell Type", y="-log10 q-value") + 
  geom_hline(yintercept=1.5, linetype="dashed", colour="white")

ABF_fold <- ggplot(na.omit(ABF.melt[ABF.melt$stat == "fold",]),
                   aes(x=cell_type, y=value, 
                       fill=cell_type,
                       group=reorder(annotation, -value))) + 
  geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
  scale_fill_manual(values=annot_cols(n_cells)) + labs(x="Cell Type", y="Fold Enrichment") +
  geom_hline(yintercept=1.0, linetype="dashed", colour="white")

##########################################
# Probabilistic inference of Causal SNPs #
##########################################

PICS.melt <- melt(PICS.gat, id.vars="annotation")
pics.cell_type = unlist(lapply(strsplit(as.character(PICS.melt$variable), fixed=T, split="_"),
                               function(x) {paste0(x[1])}))
PICS.melt$cell_type <- pics.cell_type
pics.stat = unlist(lapply(strsplit(as.character(PICS.melt$variable), fixed=T, split="_"),
                          function(x) {paste0(x[2])}))
PICS.melt$stat <- pics.stat
PICS.melt$method <- "PICS"
PICS.melt$value <- as.numeric(PICS.melt$value)

PICS_qval <- ggplot(na.omit(PICS.melt[PICS.melt$stat == "qvalue",]),
                    aes(x=cell_type, y=-log10(value), 
                        fill=cell_type,
                        group=reorder(annotation, log10(value)))) + 
  geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
  scale_fill_manual(values=annot_cols(n_cells)) + ylim(c(0, 5)) + 
  labs(x="Cell Type", y="-log10 q-value") + 
  geom_hline(yintercept=1.5, linetype="dashed", colour="white")

PICS_fold <- ggplot(na.omit(PICS.melt[PICS.melt$stat == "fold",]),
                    aes(x=cell_type, y=value, 
                        fill=cell_type,
                        group=reorder(annotation, -value))) + 
  geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
  scale_fill_manual(values=annot_cols(n_cells)) + labs(x="Cell Type", y="Fold Enrichment") + 
  geom_hline(yintercept=1.0, linetype="dashed", colour="white")

####################
# Top 1% GWAS hits #
####################

TOP.melt <- melt(Top.gat, varnames="annotation", id.vars="annotation")
top.cell_type = unlist(lapply(strsplit(as.character(TOP.melt$variable), fixed=T, split="_"),
                              function(x) {paste0(x[1])}))
TOP.melt$cell_type <- top.cell_type
top.stat = unlist(lapply(strsplit(as.character(TOP.melt$variable), fixed=T, split="_"),
                         function(x) {paste0(x[2])}))
TOP.melt$stat <- top.stat
TOP.melt$method <- "TopGWAS"
TOP.melt$value <- as.numeric(TOP.melt$value)

TOP_qval <- ggplot(na.omit(TOP.melt[TOP.melt$stat == "qvalue",]),
                   aes(x=cell_type, y=-log10(value),
                       fill=cell_type,
                       group=reorder(annotation, log10(value)))) + 
  geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
  scale_fill_manual(values=annot_cols(n_cells)) + ylim(c(0, 5)) + labs(x="Cell Type", y="-log10 q-value") + 
  geom_hline(yintercept=1.5, linetype="dashed", colour="white")

TOP_fold <- ggplot(na.omit(TOP.melt[TOP.melt$stat == "fold",]),
                   aes(x=cell_type, y=value, 
                       fill=cell_type,
                       group=reorder(annotation, -value))) + 
  geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
  scale_fill_manual(values=annot_cols(n_cells)) + labs(x="Cell Type", y="Fold Enrichment") +
  geom_hline(yintercept=1.0, linetype="dashed", colour="white")


# ########################
# # ES weighted LD score #
# ########################
# 
# SCORE.melt <- melt(Score.gat, varnames="annotation", id.vars="annotation")
# score.cell_type = unlist(lapply(strsplit(as.character(SCORE.melt$variable), fixed=T, split="_"),
#                                 function(x) {paste0(x[1])}))
# SCORE.melt$cell_type <- score.cell_type
# score.stat = unlist(lapply(strsplit(as.character(SCORE.melt$variable), fixed=T, split="_"),
#                            function(x) {paste0(x[2])}))
# SCORE.melt$stat <- score.stat
# SCORE.melt$method <- "LDscore"
# SCORE.melt$value <- as.numeric(SCORE.melt$value)
# 
# SCORE_qval <- ggplot(na.omit(SCORE.melt[SCORE.melt$stat == "qvalue",]),
#                      aes(x=cell_type, y=-log10(value),
#                          fill=cell_type,
#                          group=reorder(annotation, log10(value)))) + 
#   geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
#   scale_fill_manual(values=annot_cols(n_cells)) + ylim(c(0, 5)) + labs(x="Cell Type", y="-log10 q-value") + 
#   geom_hline(yintercept=1.5, linetype="dashed", colour="white")
# 
# SCORE_fold <- ggplot(na.omit(SCORE.melt[SCORE.melt$stat == "fold",]),
#                      aes(x=cell_type, y=value, 
#                          fill=cell_type,
#                          group=reorder(annotation, -value))) + 
#   geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
#   scale_fill_manual(values=annot_cols(n_cells)) + labs(x="Cell Type", y="Fold Enrichment") +
#   geom_hline(yintercept=1.0, linetype="dashed", colour="white")


################
# RA GWAS SNPs #
################

RA.melt <- melt(RA.gat, varnames="annotation", id.vars="annotation")
ra.cell_type = unlist(lapply(strsplit(as.character(RA.melt$variable), fixed=T, split="_"),
                             function(x) {paste0(x[1])}))
RA.melt$cell_type <- ra.cell_type
ra.stat = unlist(lapply(strsplit(as.character(RA.melt$variable), fixed=T, split="_"),
                        function(x) {paste0(x[2])}))
RA.melt$stat <- ra.stat
RA.melt$method <- "RAgwas"
RA.melt$value <- as.numeric(RA.melt$value)

RA_qval <- ggplot(na.omit(RA.melt[RA.melt$stat == "qvalue",]),
                  aes(x=cell_type, y=-log10(value),
                      fill=cell_type,
                      group=reorder(annotation, log10(value)))) + 
  geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
  scale_fill_manual(values=annot_cols(n_cells)) + ylim(c(0, 5)) + labs(x="Cell Type", y="-log10 q-value") + 
  geom_hline(yintercept=1.5, linetype="dashed", colour="white")

RA_fold <- ggplot(na.omit(RA.melt[RA.melt$stat == "fold",]),
                  aes(x=cell_type, y=value, 
                      fill=cell_type,
                      group=reorder(annotation, -value))) + 
  geom_bar(stat="identity", position="dodge", width=1.0) + theme_dark() + 
  scale_fill_manual(values=annot_cols(n_cells)) + labs(x="Cell Type", y="Fold Enrichment") +
  geom_hline(yintercept=1.0, linetype="dashed", colour="white")


########################
# Merge everything!!!! #
########################

ALL.melt <- rbind.data.frame(ABF.melt, LD.melt, PICS.melt, TOP.melt, RA.melt)
colnames(ALL.melt) <- c("annotation", "variable", "value", "cell_type", "stat", "SNPset")

p_ALL <-ggplot(na.omit(ALL.melt[ALL.melt$stat == "qvalue",]), aes(x=cell_type, y=-log10(value),
                                                                  group=order(annotation, -log10(value)),
                                                                  fill=cell_type)) +
  geom_bar(stat="identity", position="dodge", width=1.0) +
  facet_wrap(~SNPset) + theme_dark() + 
  geom_hline(yintercept=2, linetype="dashed", colour="white") + 
  scale_fill_manual(values=annot_cols(n_cells)) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.7, hjust=0))


#################################################################
# Only plot specific annotations of interest - activation marks #
#################################################################

annots <- c("H3K27ac", "H3K4me1", "H3K4me3", "ATAC", "p300", "DNase", "intergenic",
            "UTR3", "flank", "frameshift", "3flank", "intronic", "5flank", "CDS", 
            "UTR5", "DHS", "MITF")

SELECT.melt <- ALL.melt[ALL.melt$annotation %in% annots,]

p_SELECT <- ggplot(na.omit(SELECT.melt[SELECT.melt$stat == "qvalue",]), aes(x=cell_type, 
                                                                            y=-log10(value),
                                                                            group=order(annotation, -log10(value)),
                                                                            fill=cell_type)) +
  geom_bar(stat="identity", position="dodge", width=1.0) +
  facet_wrap(~SNPset) + theme_dark() + 
  geom_hline(yintercept=2, linetype="dashed", colour="white") + 
  scale_fill_manual(values=annot_cols(n_cells)) + 
  theme(axis.text.x=element_text(angle=90, vjust=0.7, hjust=0)) +
  labs(x="Functional Annotation grouped by Cell Type", y="-log10(q-value)") + 
  theme(text=element_text(size=18))

ggsave(p_SELECT, filename="GAT-compare_SNP_methods-Blond_haplotypes.png", height=14.4, width=18)


##################################################
# Plot the log2 fold Enrichments for qval < 0.05 #
##################################################

select.qval <- na.omit(SELECT.melt[SELECT.melt$stat == "qvalue" & SELECT.melt$value <= 0.01,])
select.merge <- merge(select.qval, SELECT.melt, by=c("annotation", "cell_type", "SNPset"))
# Drop RA Gwas enrichments
select.merge <- select.merge[select.merge$SNPset != "RAgwas",]
new_cols <- colorRampPalette(brewer.pal(8, "Paired"))
n_annot <- length(unique(select.merge$cell_type))

ggplot(na.omit(select.merge[select.merge$stat.y == "l2fold",]), 
       aes(x=annotation, y=value.y, fill=cell_type)) + 
  geom_bar(stat="identity", position="dodge") + 
  theme_dark() + geom_hline(yintercept=1.0, linetype="dashed", colour="white") +
  theme(axis.text.x=element_text(angle=90, vjust=0.7, hjust=0)) +
  labs(x="Functional Annotation grouped by Cell Type", y="log2 Fold Enrichment") +
  theme(text=element_text(size=18)) +  
  scale_fill_manual(values=new_cols(n_annot)) + facet_wrap(~SNPset, scales="free_x")
