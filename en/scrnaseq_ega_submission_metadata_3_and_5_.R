
# Prepare scRNA-seq sample metadata for EGA submission portal
# R and I files in separate experiments (as required by the portal)

# load librareis
library(data.table)
library(dplyr)

# load template
template <- fread("samples-two-fastq-files-paired--20200510-1215.csv", data.table = F)

# load md5 checksums as data frame (includes filenames)
path = "/Users/sofie/Downloads/EGA_submission/md5/"
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
md5$sample_id <- gsub("\\/.*", "", rownames(md5))

# filename
md5$filename <- rownames(md5)

# add extensions to make sample ids unique
md5 <- md5 %>% 
  mutate(sample_id = ifelse(grepl("L001", filename), paste0(sample_id, "_L001"), paste0(sample_id, "_L002"))) %>%
  mutate(sample_id = ifelse(grepl("I", filename), paste0(sample_id, "_I"), paste0(sample_id, "_R"))) %>%
  mutate(md5 = "")

# get checksums of encrypted and unencrypted files into separate columns
wrangle <- function(x){
  md5_encrypted <- md5 %>% filter(grepl(x, filename)) %>% filter(!grepl("gpg", filename))
  md5_unencrypted <- md5 %>% filter(grepl(x, filename)) %>% filter(grepl("gpg", filename))
  merge(md5_encrypted, md5_unencrypted, by = "sample_id", suffixes = c("", "_unencrypted"))
}

md5_i <- merge(wrangle("I1"), wrangle("I1"), by = "sample_id", suffixes = c("_read1", "_read2")) # this adds read2 to I files which will be removed later
md5_r <- merge(wrangle("R1"), wrangle("R2"), by = "sample_id", suffixes = c("_read1", "_read2"))

# add parent directory and correct extension to filenames
md5_r$filename_read1 <- paste0("name_of_folder_with_encrypted_files/", gsub(".md5", "", md5_r$filename_read1), ".gpg")
md5_r$filename_read2 <- paste0("name_of_folder_with_encrypted_files/", gsub(".md5", "", md5_r$filename_read2), ".gpg")

md5_i$filename_read1 <- paste0("name_of_folder_with_encrypted_files/", gsub(".md5", "", md5_i$filename_read1), ".gpg")
md5_i$filename_read2 <- paste0("name_of_folder_with_encrypted_files/", gsub(".md5", "", md5_i$filename_read2), ".gpg")

# add columna names and order according to template
colnames(md5_r) <- c("Sample alias", "First Checksum", "First Fastq File", "First Unencrypted checksum", "remove", "Second Checksum", "Second Fastq File", "Second Unencrypted checksum", "remove")
md5_r <- md5_r[,colnames(template)]

colnames(md5_i) <- c("Sample alias", "First Checksum", "First Fastq File", "First Unencrypted checksum", "remove", "Second Checksum", "Second Fastq File", "Second Unencrypted checksum", "remove")
md5_i <- md5_i[,colnames(template)]

# bind together metadata of R and I files 
final <- rbind(md5_r, md5_i)

# write csv files separately for 5' and 3' sequenced samples and R and I files
# these are for linking sample ids and files in EGA portal
md5_r_5prime <- md5_r %>% filter(grepl("3667|5897|6386", `Sample alias`))
write.csv(md5_r_5prime, "ega_submission_5prime_files_md5_r.csv", quote = T, row.names = F)

md5_i_5prime <- md5_i %>% filter(grepl("3667|5897|6386", `Sample alias`)) %>% select(`Sample alias`, `Fastq file` = `First Fastq File`, `Checksum` = `First Checksum`, `Unencrypted checksum` = `First Unencrypted checksum`)
write.csv(md5_i_5prime, "ega_submission_5prime_files_md5_i.csv", quote = T, row.names = F)

md5_r_3prime <- md5_r %>% filter(!grepl("3667|5897|6386", `Sample alias`))
write.csv(md5_r_3prime, "ega_submission_3prime_files_md5_r.csv", quote = T, row.names = F)

md5_i_3prime <- md5_i %>% filter(!grepl("3667|5897|6386", `Sample alias`)) %>% select(`Sample alias`, `Fastq file` = `First Fastq File`, `Checksum` = `First Checksum`, `Unencrypted checksum` = `First Unencrypted checksum`)
write.csv(md5_i_3prime, "ega_submission_3prime_files_md5_i.csv", quote = T, row.names = F)


# prepare actual sample metadata

# load template with one line per sample
template_metadata <- fread("../ega_hemap_scrna_samples.csv", data.table = F) %>% mutate(alias = paste0("FHRB_", alias)) %>% mutate(subjectId = alias)

# add metadata from template for all files
final_simpleid <- final %>% mutate(alias = gsub("_[0-9]_L.*", "", `Sample alias`)) %>% select(`Sample alias`, alias)
metadata <- merge(template_metadata, final_simpleid) %>% mutate(alias = `Sample alias`, title = `Sample alias`) %>% select(-`Sample alias`)

# write metadata
write.csv(metadata, "ega_submission_metadata.csv", quote = T, row.names = F)

# 5'´only
metadata_5 <- metadata %>% filter(grepl("3667|5897|6386", alias))
write.csv(metadata_5, "ega_submission_5prime_metadata.csv", quote = T, row.names = F)

# 3'´only
metadata_3 <- metadata %>% filter(!grepl("3667|5897|6386", alias))
write.csv(metadata_3, "ega_submission_3prime_metadata.csv", quote = T, row.names = F)

