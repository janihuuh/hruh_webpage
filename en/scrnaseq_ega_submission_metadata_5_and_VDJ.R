
# Prepare scRNA-seq sample metadata for EGA submission portal
# R (read) and I (index) files need to be in separate experiments (as required by the portal)

# load libraries
library(data.table)
library(dplyr)

# load template
template <- fread("samples-two-fastq-files-paired--20200510-1215.csv", data.table = F)

# load md5 checksums as data frame (includes filenames)
# edit the path to be the path to md5 folder on your computer
path = "/Users/sofie/Downloads/aa_EGA_submission/md5/"
filenames = list.files(path = path, pattern="*md5", recursive = T)

readdata <- function(filename) {
  df <- fread(paste(path, filename, sep = ""), header = F, data.table = F)
  md5 <- df[1,1]
  return(md5)
}

result <- do.call(rbind, lapply(filenames, readdata))
row.names(result) <- filenames

md5 <- result
colnames(md5) <- "md5"
md5 <- data.frame(md5)

# Create metadata annotations:

# sample ID
md5$sample_id <- gsub(".*\\/|_L00.*", "", rownames(md5))

# filename
md5$filename <- rownames(md5)

# add extensions to make sample ids unique
md5 <- md5 %>% 
  mutate(sample_id = ifelse(grepl("L001", filename), paste0(sample_id, "_L001"), paste0(sample_id, "_L002"))) %>%
  mutate(sample_id = ifelse(grepl("_I", filename), paste0(sample_id, "_I"), paste0(sample_id, "_R"))) %>%
  mutate(md5 = "")

# get checksums of encrypted and unencrypted files into separate columns
wrangle <- function(x){
  md5_encrypted <- md5 %>% filter(grepl(x, filename)) %>% filter(!grepl("gpg", filename))
  md5_unencrypted <- md5 %>% filter(grepl(x, filename)) %>% filter(grepl("gpg", filename))
  merge(md5_encrypted, md5_unencrypted, by = "sample_id", suffixes = c("", "_unencrypted"))
}

#create separate data frames for I and R files
md5_i <- merge(wrangle("I1"), wrangle("I1"), by = "sample_id", suffixes = c("_read1", "_read2")) # this adds read2 to I files which will be removed later
md5_r <- merge(wrangle("R1"), wrangle("R2"), by = "sample_id", suffixes = c("_read1", "_read2"))

# add parent directory and correct extension to filenames
md5_r$filename_read1 <- paste0(gsub(".md5", "", md5_r$filename_read1), ".gpg")
md5_r$filename_read2 <- paste0(gsub(".md5", "", md5_r$filename_read2), ".gpg")

md5_i$filename_read1 <- paste0(gsub(".md5", "", md5_i$filename_read1), ".gpg")
md5_i$filename_read2 <- paste0(gsub(".md5", "", md5_i$filename_read2), ".gpg")

#OPTIONAL
#if you want to rename the files in EGA (e.g. AA1 instead of FHRB_XXX, if AA1 is used in the manuscript):
#load a csv file with column "original" for sample IDs used in filenames and "new" for EGA sample IDs
newnames <- fread("manuscript_sample_IDs.csv")
for(i in 1:nrow(newnames)) {
  md5_r$sample_id <- gsub(x = md5_r$sample_id, pattern = newnames$original[i], replacement = newnames$new[i])
  md5_i$sample_id <- gsub(x = md5_i$sample_id, pattern = newnames$original[i], replacement = newnames$new[i])
}

# add column names and order according to template
colnames(md5_r) <- c("Sample alias", "First Checksum", "First Fastq File", "First Unencrypted checksum", "remove", "Second Checksum", "Second Fastq File", "Second Unencrypted checksum", "remove")
md5_r <- md5_r[,colnames(template)]

colnames(md5_i) <- c("Sample alias", "First Checksum", "First Fastq File", "First Unencrypted checksum", "remove", "Second Checksum", "Second Fastq File", "Second Unencrypted checksum", "remove")
md5_i <- md5_i[,colnames(template)]

# bind together metadata of R and I files 
final <- rbind(md5_r, md5_i)

# these are for linking sample ids and files in EGA portal
md5_r_5prime <- md5_r
write.csv(md5_r_5prime, "ega_submission_5prime_files_md5_r.csv", quote = T, row.names = F)

md5_i_5prime <- md5_i %>% select(`Sample alias`, `Fastq file` = `First Fastq File`, `Checksum` = `First Checksum`, `Unencrypted checksum` = `First Unencrypted checksum`)
write.csv(md5_i_5prime, "ega_submission_5prime_files_md5_i.csv", quote = T, row.names = F)

# prepare actual sample metadata

# load template with one line per sample (template can be donwloaded from EGA submission portal if needed)
template_metadata <- fread("ega_aa_scrna_samples.csv", data.table = F)

# add metadata from template for all files
final_simpleid <- final %>% mutate(alias = gsub("_S.*|_TCR.*", "", `Sample alias`)) %>% select(`Sample alias`, alias)
metadata <- merge(template_metadata, final_simpleid) %>% mutate(alias = `Sample alias`, title = `Sample alias`) %>% select(-`Sample alias`)
metadata <- metadata[!duplicated(metadata),]

# write metadata
write.csv(metadata, "ega_aa_submission_metadata.csv", quote = T, row.names = F)

