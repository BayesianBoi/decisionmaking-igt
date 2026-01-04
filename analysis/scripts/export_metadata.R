# Export Subject Metadata
# Run this to generate a mapping file (subj_idx -> group)
# Usage: Rscript analysis/scripts/export_metadata.R

source("analysis/utils/load_data.R")

message("Loading all data...")
dat_all <- load_all_igt_data()

# Filter for clinical studies used in fitting
clinical_studies <- c(
    "Ahn2014_HC", "Ahn2014_Amph", "Ahn2014_Hero",
    "Fridberg2010_HC", "Fridberg2010_Cbis"
)
dat_clinical <- dat_all[dat_all$study %in% clinical_studies, ]

# Add index
dat_clinical <- add_subject_index(dat_clinical)

# Create mapping table
metadata <- unique(dat_clinical[, c("subj_idx", "subj_unique", "study")])
metadata$group <- NA
metadata$group[grep("HC", metadata$study)] <- "HC"
metadata$group[grep("Amph", metadata$study)] <- "Amphetamine"
metadata$group[grep("Hero", metadata$study)] <- "Heroin"
metadata$group[grep("Cbis", metadata$study)] <- "Cannabis"

output_file <- "results/subject_metadata.csv"
write.csv(metadata, output_file, row.names = FALSE)

message(sprintf("Metadata saved to %s", output_file))
print(head(metadata))
