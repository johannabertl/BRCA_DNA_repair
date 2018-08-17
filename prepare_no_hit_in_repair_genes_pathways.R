
library("MutationRepair")

load("BRCA_mutations_in_repair_genes.Rdata")
load("BRCA_cnv_in_repair_genes.Rdata")
load("BRCA_samples.Rdata")

# remove samples without CNV data from SNV data
# (new list BRCA_mutations_in_genes_NEW in the same format as BRCA_mutations_in_genes)

BRCA_mutations_in_genes_NEW = vector("list", length(BRCA_mutations_in_genes))
for(i in 1:length(BRCA_mutations_in_genes)){
  
  BRCA_mutations_in_genes_NEW[[i]] = BRCA_mutations_in_genes[[i]][BRCA_mutations_in_genes[[i]]$Donor_ID %in% BRCA_samples$Donor_ID[BRCA_samples$CNV],]
  
}
names(BRCA_mutations_in_genes_NEW) = names(BRCA_mutations_in_genes)

# The copy number states are annotated in a binary format (0 - neutral or gain (category 2), 1 - loss or deep loss (categories 0 and 1)), then added to the list BRCA_mutations_in_genes_NEW as a new mutation types, and finally with aggregate_over_mutations, I aggregate over all mutation types except 
# - intron_variant
# - downstream_gene_variant
# - upstream_gene_variant

BRCA_cnv_in_genes_BIN = data.frame(BRCA_cnv_in_genes[, c(1,2)], ifelse(BRCA_cnv_in_genes[, -c(1,2)]<2, 1, 0))

BRCA_mutations_in_genes_NEW$CNV = BRCA_cnv_in_genes_BIN

plot(order(BRCA_mutations_in_genes_NEW$`3_prime_UTR_variant`$Donor_ID), order(BRCA_mutations_in_genes_NEW$CNV$Donor_ID))

mut_types = names(BRCA_mutations_in_genes_NEW)
mut_types = mut_types[!(mut_types %in% c("intron_variant", "downstream_gene_variant", "upstream_gene_variant"))]

BRCA_all_hits = aggregate_over_mutations(mut_data = BRCA_mutations_in_genes_NEW, mut_types = mut_types, agg_type="binary")

BRCA_no_hits_in_genes = data.frame(BRCA_all_hits[, c(1,2)], ifelse(BRCA_all_hits[, -c(1,2)]==1, 0, 1))


### aggregate over pathways ####

load("../PCAWG_analysis/Pearl_overlaps_cosmic.Rdata")
repair_pathways = Pearl_simple[, c("Gene.ID", "Pathway1")]
names(repair_pathways) = c("gene", "set")

BRCA_no_hits_in_pathways = aggregate_over_gene_sets(gene_data = BRCA_no_hits_in_genes, gene_set = repair_pathways, type = "nohit", na.rm=T)



### save ####

save(BRCA_no_hits_in_genes, BRCA_no_hits_in_pathways, file="BRCA_no_hits_in_repair_genes_pathways.Rdata")
