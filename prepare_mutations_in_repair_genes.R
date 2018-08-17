### Checking the original data ####

# I've subset the original data to the regions where the Pearl repair genes lie and annotated them with the gene names.


load("data/snv_subset.Rdata")

# snv_repair$icgc_mutation_id=NULL

# donor and sample IDs: 

nlevels(snv_repair$icgc_donor_id)
nlevels(snv_repair$icgc_specimen_id)
nlevels(snv_repair$matched_icgc_sample_id)
nlevels(snv_repair$submitted_sample_id)
head(snv_repair$submitted_sample_id)
nlevels(snv_repair$submitted_matched_sample_id)
head(snv_repair$submitted_matched_sample_id)
# the submitted and submitted matched sample_ids seem to be the ones used in Nik-Zainal et al 2016 (supplementary tables)

levels(snv_repair$project_code)

# coordinates: 

levels(snv_repair$chromosome)
hist(snv_repair$chromosome_start)
hist(snv_repair$chromosome_end)

hist(snv_repair$chromosome_end - snv_repair$chromosome_start)
table(snv_repair$chromosome_end - snv_repair$chromosome_start)
table(snv_repair$chromosome_strand) # not informative

# technical information
levels(snv_repair$assembly_version) # fine, ignore
range(snv_repair$gene_build_version) # fine, ignore
table(snv_repair$platform) # 2 different sequencers were used
table(snv_repair$experimental_protocol) # empty
table(snv_repair$sequencing_strategy) # WGS for all samples
table(snv_repair$base_calling_algorithm) # empty
table(snv_repair$alignment_algorithm) # two aligners (BWA and empty), seems to be congruent with the 2 different sequencers
table(snv_repair$variation_calling_algorithm) # 4 levels


# mutations
table(snv_repair$mutation_type)
head(snv_repair$reference_genome_allele)
head(snv_repair$mutated_from_allele)
head(snv_repair$mutated_to_allele)
table(snv_repair$consequence_type) #### MAIN VARIABLE
head(snv_repair$aa_mutation) #?
head(snv_repair$cds_mutation)
head(snv_repair$gene_affected) # gene ID
head(snv_repair$transcript_affected) # transcript ID


# quality
table(snv_repair$quality_score) # empty
table(snv_repair$probability) # empty
hist(snv_repair$total_read_count)
min(snv_repair$total_read_count, na.rm=T)
hist(snv_repair$mutant_allele_read_count)
min(snv_repair$mutant_allele_read_count, na.rm=T)
table(snv_repair$verification_status)
table(snv_repair$verification_platform) # the verified ones were verified by capillary sequencing
table(snv_repair$seq_coverage) # empty


### Checking if the gene names are correct ####

# The dataset also comes with ensemble gene IDs. Here, I check if they are equal to the gene names that I have assigned per region.


snv_repair_red = snv_repair[, c("project_code", "submitted_sample_id", "consequence_type", "gene_affected", "gene_name")]

gene_names_ens = read.table("../PCAWG_analysis/data/Pearl_repair_genes_ensembl.txt", sep="\t", header=T)

snv_repair_red_test = merge(snv_repair_red, gene_names_ens, by.x ="gene_affected", by.y = "Gene.stable.ID", all.x=T, all.y=F)

mean(snv_repair_red_test$gene_name == snv_repair_red_test$Gene.name, na.rm=T)
sum(snv_repair_red_test$gene_name != snv_repair_red_test$Gene.name, na.rm=T)

test = snv_repair_red_test[snv_repair_red_test$gene_name != snv_repair_red_test$Gene.name,]
test = test[!is.na(test$gene_affected),]

# I checked all genes where the gene name I have attached to the region (gene_name) and the gene name I have attached by translating the ensemble gene ID (gene_affected) to a gene name (Gene.name) using the table ../PCAWG_analysis/data/Pearl_repair_genes_ensembl3.txt are not the same. 
# 1) All but one cases are genes that have two names. They can all be found in ../PCAWG_analysis/data/alternative_names_for_GENCODE.txt. That's okay. 
# 2) GEN1 and SMC6 are not the same, but they overlap. 4 intron_variants are assigned ambiguiously. But since they are really considered as having no effect, it does not matter which gene they are assigned to. 

### FINE ###

snv_repair_red = snv_repair_red[, c("project_code", "submitted_sample_id", "consequence_type", "gene_name")]


### save intermediate dataset ####

save(snv_repair_red, file="data/snv_subset_intermediate.Rdata")


### reshape ####

library("reshape2")
library("MutationRepair")

load("data/snv_subset_intermediate.Rdata")
Pearl = unlist(read.table("../PCAWG_analysis/Pearl_repair_gene_symbols.txt", as.is=T))
names(snv_repair_red) = c("Project_Code", "Donor_ID", "Mutation_Class", "Gene_Name")

# remove trailing "a" from Donor_ID (I think it was originally used to distinguish between the cancer and the normal sample ("a" and "b"))
levels(snv_repair_red$Donor_ID) = sub("[ab]\\d?", "", levels(snv_repair_red$Donor_ID))


### create and save samples dataframe ####

BRCA_samples = data.frame(Donor_ID = levels(snv_repair_red$Donor_ID), Project_Code = "BRCA")
save(BRCA_samples, file="BRCA_samples.Rdata")


### create and save mutations in repair genes list ####

# droplevels: to make sure there are not mutation classes with zero observations
snv_repair_red = droplevels(snv_repair_red)

BRCA_mutations_in_genes = reshape_by_types(mutations = snv_repair_red, all_genes = Pearl, all_samples = BRCA_samples)


names(BRCA_mutations_in_genes)
lapply(BRCA_mutations_in_genes, dim)

sumvec = numeric(length(BRCA_mutations_in_genes))
for(i in 1:length(BRCA_mutations_in_genes)) sumvec[i] = sum(BRCA_mutations_in_genes[[i]][, -c(1,2)])

table(snv_repair_red$Mutation_Class)
# same as sumvec -- perfect! 



save(BRCA_mutations_in_genes, file="BRCA_mutations_in_repair_genes.Rdata")

