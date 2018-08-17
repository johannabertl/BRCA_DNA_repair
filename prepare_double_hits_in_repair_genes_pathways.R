

library("MutationRepair")


load("BRCA_lof_in_repair_genes.Rdata")
load("BRCA_cnv_in_repair_genes.Rdata")
load("BRCA_samples.Rdata")


### annotate double hits ####

BRCA_cnv_in_genes_INV = 2 - as.matrix(BRCA_cnv_in_genes[, -c(1,2)])

BRCA_lof_in_genes_SUB = BRCA_lof_in_genes[BRCA_lof_in_genes$Donor_ID %in% BRCA_samples$Donor_ID[BRCA_samples$CNV], ]
sum(BRCA_lof_in_genes_SUB$Donor_ID == BRCA_cnv_in_genes$Donor_ID) # yes

BRCA_double_prelim = BRCA_cnv_in_genes_INV + as.matrix(BRCA_lof_in_genes_SUB[, -c(1,2)])

BRCA_double_prelim_binary = ifelse(BRCA_double_prelim >=2, 1, 0)
sum(BRCA_double_prelim_binary)

BRCA_double_hits_in_genes = data.frame(BRCA_lof_in_genes_SUB[, c(1,2)], BRCA_double_prelim_binary)


### aggregate over pathways and subpathways ####

load("../PCAWG_analysis/Pearl_overlaps_cosmic.Rdata")
repair_pathways = Pearl_simple[, c("Gene.ID", "Pathway1")]
names(repair_pathways) = c("gene", "set")

BRCA_double_hits_in_pathways = aggregate_over_gene_sets(gene_data = BRCA_double_hits_in_genes, gene_set = repair_pathways, type = "hit", na.rm=T)


### save ####

save(BRCA_double_hits_in_genes, BRCA_double_hits_in_pathways, file="BRCA_double_hits_genes_pathways.Rdata")
