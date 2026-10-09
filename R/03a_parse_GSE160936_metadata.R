# ============================================================
# 03a_parse_GSE160936_metadata.R
# ============================================================
# Purpose:
#   Parse GSE160936 sample-level metadata from the GEO series
#   metadata text file and create a clean metadata table with:
#     GSM, donor, region, diagnosis, Braak stage, age, sex, RIN
#
# Fix in this version:
#   The earlier parser accidentally matched "braak tangle stage:"
#   when searching for "age:" because "stage:" contains "age:".
#   This version anchors the search to quoted "age:" values only.
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

setup_path <- file.path(PROJECT_ROOT, "RDS/project_setup.rds")

if (!file.exists(setup_path)) {
  stop("project_setup.rds not found. Run 01_Setup_AD61026.R first.")
}

setup <- readRDS(setup_path)
list2env(setup, envir = .GlobalEnv)

series_path <- file.path(DATA_DIR, "GSE160936_series_metadata_lines.txt")
out_csv     <- file.path(DATA_DIR, "GSE160936_sample_metadata_parsed.csv")

if (!file.exists(series_path)) {
  stop("Missing file: ", series_path)
}

cat("Parsing:", series_path, "\n\n")

lines <- readLines(series_path, warn = FALSE)

# ------------------------------------------------------------
# Helper function to parse one tab-delimited GEO metadata line.
# ------------------------------------------------------------
parse_geo_line <- function(pattern, lines, ignore.case = TRUE) {
  idx <- grep(pattern, lines, ignore.case = ignore.case, perl = TRUE)

  if (length(idx) == 0) {
    stop("Could not find metadata line matching pattern: ", pattern)
  }

  if (length(idx) > 1) {
    cat("NOTE: Multiple lines matched pattern:", pattern, "\n")
    cat("      Using first match.\n")
  }

  line <- lines[idx[1]]
  parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
  values <- parts[-1]
  values <- gsub('^"|"$', "", values)
  return(values)
}

# ------------------------------------------------------------
# Parse metadata lines.
# ------------------------------------------------------------
titles <- parse_geo_line("^!Sample_title", lines)
gsms   <- parse_geo_line("^!Sample_geo_accession", lines)

disease_raw <- parse_geo_line('"disease state:', lines)
braak_raw   <- parse_geo_line('"braak tangle stage:', lines)
region_raw  <- parse_geo_line('"brain region', lines)

# Important: quote anchor prevents accidental match to "stage:".
age_raw     <- parse_geo_line('"age:', lines)
sex_raw     <- parse_geo_line('"Sex:', lines)
rin_raw     <- parse_geo_line('"rin:', lines)

# ------------------------------------------------------------
# Clean field values.
# ------------------------------------------------------------
clean_after_colon <- function(x) {
  trimws(sub("^[^:]+:\\s*", "", x))
}

disease <- clean_after_colon(disease_raw)
braak   <- as.numeric(clean_after_colon(braak_raw))
region  <- clean_after_colon(region_raw)
age     <- as.numeric(clean_after_colon(age_raw))
sex     <- clean_after_colon(sex_raw)
rin     <- as.numeric(clean_after_colon(rin_raw))

# Donor is the title without the final region label.
donor <- sub("\\s+(SSC|EC)$", "", titles)

diagnosis_std <- ifelse(disease == "AD", "AD",
                        ifelse(disease == "Non-disease control", "Control", NA))

# ------------------------------------------------------------
# Link metadata rows to local tar.gz file names by GSM.
# ------------------------------------------------------------
raw_dir <- file.path(DATA_DIR, "GSE160936_RAW")
raw_files <- list.files(raw_dir, pattern = "\\.tar\\.gz$", full.names = FALSE)

file_map <- data.frame(
  file = raw_files,
  gsm_id = sub("^(GSM[0-9]+)_.*$", "\\1", raw_files),
  local_sample_id = sub("_filtered_feature_bc_matrix\\.tar\\.gz$", "",
                        raw_files),
  stringsAsFactors = FALSE
)

meta <- data.frame(
  dataset = "GSE160936",
  gsm_id = gsms,
  title = titles,
  donor_id = donor,
  region = region,
  disease_state = disease,
  diagnosis_std = factor(diagnosis_std, levels = c("Control", "AD")),
  braak_stage = braak,
  age = age,
  sex = sex,
  rin = rin,
  stringsAsFactors = FALSE
)

meta <- merge(meta, file_map, by = "gsm_id", all.x = TRUE, sort = FALSE)

meta <- meta[, c("dataset", "gsm_id", "local_sample_id", "file", "title",
                 "donor_id", "region", "disease_state", "diagnosis_std",
                 "braak_stage", "age", "sex", "rin")]

meta <- meta[match(gsms, meta$gsm_id), ]

# ------------------------------------------------------------
# Sanity checks.
# ------------------------------------------------------------
cat("Rows:", nrow(meta), "\n")
cat("Unique donors:", length(unique(meta$donor_id)), "\n\n")

cat("Sample-level diagnosis table:\n")
print(table(meta$diagnosis_std, useNA = "ifany"))

cat("\nDonor-level diagnosis table:\n")
donor_dx <- meta[!duplicated(meta$donor_id), c("donor_id", "diagnosis_std")]
print(table(donor_dx$diagnosis_std, useNA = "ifany"))

cat("\nRegion table:\n")
print(table(meta$region, useNA = "ifany"))

cat("\nSamples per donor:\n")
print(table(meta$donor_id))

cat("\nAge summary:\n")
print(summary(meta$age))

cat("\nPreview:\n")
print(meta)

if (any(is.na(meta$file))) {
  warning("Some metadata rows did not match a raw tar.gz file.")
  print(meta[is.na(meta$file), ])
}

missing_meta <- setdiff(file_map$gsm_id, meta$gsm_id)
if (length(missing_meta) > 0) {
  warning("Some raw files did not match metadata GSM IDs: ",
          paste(missing_meta, collapse = ", "))
}

write.csv(meta, out_csv, row.names = FALSE)
cat("\nSaved parsed metadata to:\n", out_csv, "\n")
cat("\nDone.\n")
