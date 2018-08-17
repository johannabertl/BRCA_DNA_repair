
load("BRCA_samples.Rdata")

sum(BRCA_samples$CNV & BRCA_samples$NikZainal560)
samples344 = as.character(BRCA_samples$Donor_ID[BRCA_samples$CNV])


# double hits

load("BRCA_double_hits_genes_pathways.Rdata")

BRCA344_double_hits_in_genes = BRCA_double_hits_in_genes[BRCA_double_hits_in_genes$Donor_ID %in% samples344,]

BRCA344_double_hits_in_pathways = BRCA_double_hits_in_pathways[BRCA_double_hits_in_pathways$Donor_ID %in% samples344,]

save(BRCA344_double_hits_in_genes, BRCA344_double_hits_in_pathways, file="BRCA344_double_hits_genes_pathways.Rdata")

rm(list=c("BRCA344_double_hits_in_genes", "BRCA344_double_hits_in_pathways", "BRCA_double_hits_in_genes", "BRCA_double_hits_in_pathways"))


# mutational signatures

load("BRCA_COSMIC30_mutational_signatures.Rdata")

BRCA344_cosmic30_signatures = vector("list", 5)
names(BRCA344_cosmic30_signatures) = names(BRCA_cosmic30_signatures)
BRCA344_cosmic30_signatures$sig.prop = BRCA_cosmic30_signatures$sig.prop[BRCA_cosmic30_signatures$sig.prop$Donor_ID %in% samples344,]
BRCA344_cosmic30_signatures$sig.count = BRCA_cosmic30_signatures$sig.count[BRCA_cosmic30_signatures$sig.count$Donor_ID %in% samples344,]
BRCA344_cosmic30_signatures$cos_sim_res = BRCA_cosmic30_signatures$cos_sim_res[BRCA_cosmic30_signatures$sig.count$Donor_ID %in% samples344]
BRCA344_cosmic30_signatures$cutoff = BRCA_cosmic30_signatures$cutoff

save(BRCA344_cosmic30_signatures, file="BRCA344_COSMIC30_mutational_signatures.Rdata")

rm(BRCA344_cosmic30_signatures, BRCA_cosmic30_signatures)


# mutation counts

load("BRCA_mutation_counts.Rdata")

BRCA344_mutation_counts = BRCA_mutation_counts[BRCA_mutation_counts$Donor_ID %in% samples344,]
BRCA344_total_mutation_counts = BRCA_total_mutation_count[BRCA_total_mutation_count$Donor_ID %in% samples344,]

save(BRCA344_total_mutation_counts, BRCA344_mutation_counts, file="BRCA344_mutation_counts.Rdata")

rm(BRCA_mutation_counts, BRCA_total_mutation_count, BRCA344_mutation_counts, BRCA344_total_mutation_counts)


# no hits

load("BRCA_no_hits_in_repair_genes_pathways.Rdata")

BRCA344_no_hits_in_genes = BRCA_no_hits_in_genes[BRCA_no_hits_in_genes$Donor_ID %in% samples344,]
BRCA344_no_hits_in_pathways = BRCA_no_hits_in_pathways[BRCA_no_hits_in_pathways$Donor_ID %in% samples344,]

save(BRCA344_no_hits_in_genes, BRCA344_no_hits_in_pathways, file="BRCA344_no_hits_in_repair_genes_pathways.Rdata")

rm(BRCA344_no_hits_in_genes, BRCA344_no_hits_in_pathways, BRCA_no_hits_in_genes, BRCA_no_hits_in_pathways)
