


load("BRCA_mutations_in_repair_genes.Rdata")
names(BRCA_mutations_in_genes)

library("MutationRepair")

BRCA_lof_in_genes = aggregate_over_mutations(BRCA_mutations_in_genes, mut_types = c("disruptive_inframe_deletion", "frameshift_variant", "stop_gained"), agg_type = "sum")

save(BRCA_lof_in_genes, file="BRCA_lof_in_repair_genes.Rdata")
