
library("MutationalPatterns")
library("MutationRepair")
load("BRCA_mutation_counts.Rdata")
data("COSMIC30")

?extract_and_reformat

mutation_matrix = t(as.matrix(BRCA_mutation_counts[, -c(1,2)]))
colnames(mutation_matrix) = BRCA_mutation_counts$Donor_ID


BRCA_cosmic30_signatures = extract_and_reformat(mutations = mutation_matrix, signatures = COSMIC30, samples = BRCA_mutation_counts[, c("Donor_ID", "Project_Code")], cutoff=0.1, min_delta = 0.85)
                                                  

save(BRCA_cosmic30_signatures, file="BRCA_COSMIC30_mutational_signatures.Rdata")
