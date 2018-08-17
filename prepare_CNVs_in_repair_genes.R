### DON't run everything! Note that this script alters the file BRCA_samples.Rdata!!!

### load original data and check ####
# (The data is subset such that it only contains CNV events that affect the Pearl repair genes)

load("data/cnv_subset.Rdata")

names(cnvs)

# are all Pearl genes included?

Pearl = unlist(read.table("../PCAWG_analysis/Pearl_repair_gene_symbols.txt", as.is=T))
Pearl[which(!(Pearl %in% cnvs$Gene.name))]
# they might be genuinely missing, maybe because of too few markers (SNPs on the array) for the region of these genes

# only what is actually contained in this data: 
cnv = droplevels(cnvs)

# genes/locations
nlevels(cnv$gene_affected)
nlevels(cnv$Gene.name)
levels(cnv$chromosome)
hist(cnv$chromosome_start)
hist(cnv$chromosome_end)
hist(cnv$chromosome_end - cnv$chromosome_start)
range(cnv$chromosome_end - cnv$chromosome_start)
sum(is.na(cnv$gene_affected))
sum(is.na(cnv$Gene.name))


# samples
nlevels(cnv$icgc_donor_id)
nlevels(cnv$project_code)
head(cnv$submitted_sample_id)
head(cnv$submitted_matched_sample_id)

# mutations
table(cnv$mutation_type)
# table(cnvs$mutation_type)
hist(cnv$segment_mean)
table(cnv$segment_mean) # => integer values
table(cnv$copy_number) # missing

table(cnv$segment_mean[cnv$mutation_type=="copy neutral"])
table(cnv$segment_mean[cnv$mutation_type=="copy neutral LOH"])
table(cnv$segment_mean[cnv$mutation_type=="gain"])
table(cnv$segment_mean[cnv$mutation_type=="loss"])

table(cnv$copy_number) # missing

# technical details
# a lot of missing variables, incl qualitiy score, probability, etc. (maybe because copy numbers weren't called?)
levels(cnv$assembly_version) # ok
levels(cnv$sequencing_strategy) # ok
levels(cnv$verification_status)
levels(cnv$gene_build_version)
levels(cnv$platform)
levels(cnv$alignment_algorithm)
levels(cnv$variation_calling_algorithm)

# are the copy numbers summarized per gene?
n_samples = nlevels(cnv$submitted_sample_id)
max_vec = numeric(n_samples)
for(i in 1:n_samples) {
  max_vec[i] = max(table(cnv$gene_affected[cnv$submitted_sample_id==levels(cnv$submitted_sample_id)[i]]))
}
max_vec

# Result: the cnvs are not counted by gene, but either by segment or by event. It's not per SNP either, since the affected region (chromomsome_end - chromosome_start) is 800 bp -- 2.5 10^8 bp (!?) long.
# 


### Prepare table for reshaping ####

# reduce information and compute new categories

cnv_red = cnv[, c("submitted_sample_id", "mutation_type", "segment_mean", "Gene.name")]

# new categories: 2 -- two or more alleles
#                 1 -- one allele lost
#                 0 -- deep deletion
# new categories based on mutation_type
# copy neutral -> 2
# copy neutral LOH -> 2
# gain -> 2
# loss and segment mean==0 -> 0
# loss and segment mean==1 -> 1
# loss and segment mean>=2 -> 2

cnv_red$CNV = 2
cnv_red$CNV[cnv_red$mutation_type=="loss" & cnv_red$segment_mean==0] = 0
cnv_red$CNV[cnv_red$mutation_type=="loss" & cnv_red$segment_mean==1] = 1


# checking: 
table(cnv_red$CNV)
table(cnv_red[, c("mutation_type", "segment_mean")])

# preparing table to reshape: 
cnv_red$mutation_type = NULL
cnv_red$segment_mean = NULL
levels(cnv_red$submitted_sample_id) = sub("[ab]\\d?", "", levels(cnv_red$submitted_sample_id))

names(cnv_red) = c("Donor_ID", "Gene_Name", "Mutation_Class")

save(cnv_red, file="data/cnv_subset_intermediate.Rdata")


##########################################################################


# ### prepare samples ####
# load("BRCA_samples.Rdata")
# 
# names(BRCA_samples) = c("Donor_ID", "Project_Code")
# BRCA_samples$CNV = BRCA_samples$Donor_ID %in% levels(cnv_red$Donor_ID)
# BRCA_samples$SNV = T
# 
# save(BRCA_samples, file="BRCA_samples.Rdata")


############################################################################

### reshape ####

load("BRCA_samples.Rdata")
load("data/cnv_subset_intermediate.Rdata")
Pearl = unlist(read.table("../PCAWG_analysis/Pearl_repair_gene_symbols.txt", as.is=T))

library("reshape2")
library("MutationRepair")

cnv_inter = reshape_by_types(cnv_red, all_genes = Pearl, all_samples = BRCA_samples[BRCA_samples$CNV,c("Donor_ID", "Project_Code")])

cnv_final = cnv_inter[[1]]
cnv_final[, -c(1,2)] = 2
cnv_final[cnv_inter[[2]] == 1] = 1
cnv_final[cnv_inter[[1]] == 1] = 0



# test
table(as.matrix(cnv_final[, -c(1,2)]))
table(cnv_red$Mutation_Class)
# the discrepancy stems from multiple events per gene, and the most severe (0) is counted, if there are multiple.

cnv_in_repair_genes = cnv_final

# make GTF2H5, RDM1 and RECQL4 NA (these seem to be genuinely missing measurements, see above)
cnv_in_repair_genes$GTF2H5 = NA
cnv_in_repair_genes$RDM1 = NA
cnv_in_repair_genes$RECQL4 = NA

BRCA_cnv_in_genes = cnv_in_repair_genes


### save ####

save(BRCA_cnv_in_genes, file="BRCA_cnv_in_repair_genes.Rdata")
