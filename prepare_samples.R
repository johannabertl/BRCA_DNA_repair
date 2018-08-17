### read table ####

sample = read.table("data/sample.tsv", sep="\t", header=T, as.is=T)

# remove unnecessary columns

head(sample)
sample$icgc_specimen_id=NULL
sample$submitted_specimen_id=NULL
sample$submitted_donor_id=NULL
sample$analyzed_sample_interval=NULL
sample$percentage_cellularity=NULL
levels(sample$study)
sample$study=NULL
sample$level_of_cellularity=NULL



### extract donor/sample IDs ####
# that match the ones used by Guo (Juul et al, Bioinformatics, 2018), Nik-Zainal et al, 2016, and ICGC.

# Guo has used the icgc_sample_id, but they are per sample (tumor and normal), i. e. at least 2 per donor. 
# Which one has she used in the count tables? 

count = read.table("data/BRCA_count_table.txt", header=T)

count$genomicSeg=NULL
count$replication_timing=NULL
count$expression=NULL
count$cancer_type=NULL

sl = levels(count$sample_id)

# subset to the sample_ids Guo used: 
sample_subset = sample[sample$icgc_sample_id %in% sl,]

# remove trailing a, b or b2, b3, a2, a3 from the submitted_sample_id to make it fit with the IDs in the supplementary tables of Nik-Zainal et al.
# (this removes the information about the tumor and the blood sample, but I don't need it here)
sample_subset$new_sample_id = sub("[ab]\\d?", "", sample_subset$submitted_sample_id)

# make the new dataframe: 

BRCA_samples = data.frame(Donor_ID = sample_subset$new_sample_id, Project_Code=rep("BRCA", 569), Guo_ID = sample_subset$icgc_sample_id, ICGC_Donor_ID = sample_subset$icgc_donor_id)



### add new columns ####

# add columns about data availability: 

BRCA_samples$SNV = T

load("data/cnv_subset.Rdata")
BRCA_samples$CNV = BRCA_samples$ICGC_Donor_ID %in% levels(cnvs$icgc_donor_id)

# add columns: samples included in Nik-Zainal paper: 

NikZ = read.table("data/Supplementary Table 1 CLINICAL.PATHOLOGY.DATA.FREEZE.ANALYSIS.v4.032015.csv", sep=",", header=T)

BRCA_samples$NikZainal560 = BRCA_samples$Donor_ID %in% NikZ$sample_name



### save ####

save(BRCA_samples, file="BRCA_samples.Rdata")
