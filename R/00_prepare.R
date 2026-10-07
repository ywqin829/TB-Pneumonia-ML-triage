## TB/pneumonia revision: source audit and frozen 60/20/20 allocation.
## Run from this directory with Rscript 00_prepare.R.

suppressPackageStartupMessages(library(rsample))

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("models", showWarnings = FALSE)

source_file <- "TB_Pneumonia.csv"
stopifnot(file.exists(source_file))
raw <- read.csv(source_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("Age", "Gender", "PLT", "RBC", "PDW", "MPV", "LYM",
              "NEU", "MONO", "MONO.", "SII", "SIRI", "MLR", "species")
stopifnot(all(required %in% names(raw)))
stopifnot(all(raw$species %in% c("TB", "Pneumonia")))
stopifnot(!anyNA(raw))

is_extra_identical <- duplicated(raw)
clean <- raw[!is_extra_identical, , drop = FALSE]
row.names(clean) <- NULL
clean$row_id <- sprintf("data_row_%04d", seq_len(nrow(clean)))
clean$species <- factor(clean$species, levels = c("TB", "Pneumonia"))
clean$Gender <- factor(clean$Gender, levels = c(1, 2), labels = c("Male", "Female"))
stopifnot(!anyNA(clean$Gender), all(clean$LYM > 0), all(clean$MONO > 0))

## LMR is the reciprocal of the supplied MLR. PIV is reconstructed from
## the four absolute cell counts; IPIV has no defined/source column and is
## not invented. See the method notes for the source of the PIV formula.
clean$LMR <- clean$LYM / clean$MONO
clean$PIV <- clean$NEU * clean$PLT * clean$MONO / clean$LYM
stopifnot(all(is.finite(clean$LMR)), all(is.finite(clean$PIV)))

compare_formula <- function(variable, expected) {
  observed <- as.numeric(clean[[variable]])
  data.frame(variable = variable,
             maximum_absolute_error = max(abs(observed - expected)),
             maximum_relative_error = max(abs(observed - expected) / pmax(abs(expected), 1e-12)),
             stringsAsFactors = FALSE)
}
formula_audit <- do.call(rbind, list(
  compare_formula("SII", clean$PLT * clean$NEU / clean$LYM),
  compare_formula("SIRI", clean$NEU * clean$MONO / clean$LYM),
  compare_formula("MLR", clean$MONO / clean$LYM),
  compare_formula("dNLR", clean$NEU / (clean$WBC - clean$LYM))
))
write.csv(formula_audit, "results/formula_audit.csv", row.names = FALSE)
dn_lower <- clean$WBC - clean$NEU
dn_audit <- data.frame(
  data_row_id = clean$row_id,
  source_dNLR = clean$dNLR,
  standard_dNLR = ifelse(dn_lower > 0, clean$NEU/dn_lower, NA_real_),
  standard_denominator_nonpositive = dn_lower <= 0)
write.csv(dn_audit, "results/dNLR_definition_audit.csv", row.names = FALSE)

## The source is de-identified: these are unique *data rows*, not a verified
## patient-level duplicate audit. Record the original row index separately.
duplicate_audit <- data.frame(source_row = which(is_extra_identical),
                              diagnosis = raw$species[is_extra_identical])
write.csv(duplicate_audit, "results/exact_duplicate_rows.csv", row.names = FALSE)

## Preserve the existing code's seeded, diagnosis-stratified 60/20/20 route.
set.seed(123)
first <- rsample::initial_split(clean, prop = 0.6, strata = species)
train <- rsample::training(first)
remainder <- rsample::testing(first)
second <- rsample::initial_split(remainder, prop = 0.5, strata = species)
validation <- rsample::training(second)
test <- rsample::testing(second)

split_assignment <- rbind(
  data.frame(row_id = train$row_id, split = "training", diagnosis = train$species),
  data.frame(row_id = validation$row_id, split = "validation", diagnosis = validation$species),
  data.frame(row_id = test$row_id, split = "test", diagnosis = test$species)
)
stopifnot(nrow(split_assignment) == nrow(clean), !anyDuplicated(split_assignment$row_id))
write.csv(split_assignment, "results/split_assignment.csv", row.names = FALSE)
split_counts <- as.data.frame(table(split_assignment$split, split_assignment$diagnosis))
names(split_counts) <- c("split", "diagnosis", "n")
write.csv(split_counts, "results/split_counts.csv", row.names = FALSE)

## Retain only columns present in the extract plus explicitly derivable LMR
## and PIV. MLR is omitted from candidate selection because it is the exact
## reciprocal of LMR, but remains available for the descriptive audit.
candidate_features <- setdiff(names(raw), c("Age", "Gender", "MLR", "dNLR", "species"))
candidate_features <- c(candidate_features, "LMR", "PIV")
stopifnot(length(candidate_features) == 25L)
stopifnot(all(vapply(clean[candidate_features], is.numeric, logical(1))))
write.csv(data.frame(feature = candidate_features,
                     source = ifelse(candidate_features %in% c("LMR", "PIV"),
                                     "derived from supplied absolute counts",
                                     "provided CSV")),
          "results/candidate_features.csv", row.names = FALSE)

## A compact source audit, including explicit boundaries of inference.
audit_lines <- c(
  "TB/pneumonia source audit for revision rerun",
  sprintf("Source rows: %d", nrow(raw)),
  sprintf("Excess completely identical rows: %d", sum(is_extra_identical)),
  sprintf("Distinct complete data rows retained: %d", nrow(clean)),
  sprintf("Retained TB rows: %d", sum(clean$species == "TB")),
  sprintf("Retained pneumonia rows: %d", sum(clean$species == "Pneumonia")),
  sprintf("Missing cells in supplied extract: %d", sum(is.na(raw))),
  sprintf("Candidate clinical predictors: %d", length(candidate_features)),
  "Source dNLR equals NEU/(WBC-LYM), not the usual NEU/(WBC-NEU); excluded from model candidates.",
  sprintf("Rows with WBC-NEU <= 0: %d", sum(dn_lower <= 0)),
  "No patient identifier, date, hospital, screening log, or 30-case flag is present.",
  "The 17 identical extra rows cannot be called verified duplicate patients.",
  "Historical incomplete-record exclusions and temporal membership cannot be inferred.",
  "No IPIV or AIS definition/source values are available; neither is fabricated."
)
writeLines(audit_lines, "results/source_audit.txt", useBytes = TRUE)

saveRDS(list(data = clean, train = train, validation = validation, test = test,
             candidate_features = candidate_features, seed = 123L,
             formula_audit = formula_audit),
        "results/prepared_data.rds")

cat(paste(audit_lines, collapse = "\n"), "\n")
print(split_counts)
