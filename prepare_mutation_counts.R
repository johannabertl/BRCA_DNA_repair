

count = read.table("data/BRCA_count_table.txt", header=T, as.is=T)

# remove unused columns

count$genomicSeg=NULL
count$replication_timing=NULL
count$expression=NULL
count$cancer_type=NULL

# remove non-mutated sites

count = count[count$to != count$from,]

# use reverse complements on G and A positions and create the 96 "mutation types"

rev.compl <- function(base_vec) {
  
  out <- base_vec
  out[base_vec=="A"] <- "T"
  out[base_vec=="C"] <- "G"
  out[base_vec=="G"] <- "C"
  out[base_vec=="T"] <- "A"
  
  return(out)
}


count$from_new = ifelse(count$from %in% c("G", "A"), rev.compl(count$from), count$from)
count$to_new = ifelse(count$from %in% c("G", "A"), rev.compl(count$to), count$to)
count$left_new = ifelse(count$from %in% c("G", "A"), rev.compl(count$right), count$left)
count$right_new = ifelse(count$from %in% c("G", "A"), rev.compl(count$left), count$right)

count$mutation = paste0(count$left_new, "[", count$from_new, ">", count$to_new, "]", count$right_new)

count[, c("from", "to", "left", "right", "from_new", "to_new", "left_new", "right_new")]=NULL


# replace sample_id with Donor_ID

load("BRCA_samples.Rdata")

count_merged = merge(BRCA_samples[, c("Donor_ID", "Project_Code", "Guo_ID")], count, by.x = "Guo_ID", by.y = "sample_id", all.y = T)
count_merged$Guo_ID = NULL


# reshape to wide format

library("reshape2")
BRCA_mutation_counts = dcast(count_merged, Donor_ID + Project_Code ~ mutation, value.var = "count", fun.agg = sum)

# reorder colnames so they fit with the signatures, etc. 
library(MutationRepair)
data(COSMIC30)
BRCA_mutation_counts = BRCA_mutation_counts[, c("Donor_ID", "Project_Code", rownames(COSMIC30))]

# total count

BRCA_total_mutation_count = data.frame(BRCA_mutation_counts[, c(1,2)], total_count = rowSums(BRCA_mutation_counts[, -c(1,2)]))

# save

save(BRCA_mutation_counts, BRCA_total_mutation_count, file="BRCA_mutation_counts.Rdata")
