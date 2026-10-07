## Delivery smoke check: actual saved predictions, unit conversion, Shiny server.
## Reads inputs/models; writes only results/app_smoke_check.txt.
source("app_predict.R", encoding = "UTF-8")
bundle <- load_tb_bundle(".")
prep <- readRDS("results/prepared_data.rds")
test_predictions <- read.csv("results/test_predictions.csv", stringsAsFactors = FALSE)
reference <- test_predictions[test_predictions$algorithm == "xgbTree", ]
actual <- as.numeric(predict(bundle$model, newdata = prep$test, type = "prob")$TB)
expected <- reference$probability_TB[match(prep$test$row_id, reference$row_id)]
batch_error <- max(abs(actual - expected))
stopifnot(batch_error < 1e-12)

indices <- unique(c(1L, which.min(actual), which.max(actual),
                    which.min(abs(actual - bundle$thresholds[["low"]])),
                    which.min(abs(actual - bundle$thresholds[["high"]])),
                    which.max(replace(actual, actual >= bundle$thresholds[["balanced"]], -Inf)),
                    which.min(replace(actual, actual < bundle$thresholds[["balanced"]], Inf))))
balanced_labels <- ifelse(actual >= bundle$thresholds[["balanced"]], "TB", "Pneumonia")
observed_tb <- as.character(prep$test$species) == "TB"
balanced_confusion <- c(TP = sum(balanced_labels == "TB" & observed_tb),
                        FP = sum(balanced_labels == "TB" & !observed_tb),
                        TN = sum(balanced_labels == "Pneumonia" & !observed_tb),
                        FN = sum(balanced_labels == "Pneumonia" & observed_tb))
stopifnot(identical(unname(balanced_confusion), c(121L, 97L, 238L, 22L)))
single_errors <- vapply(indices, function(i) {
  result <- predict_tb_inputs(bundle, prep$test$MCV[i], prep$test$PDW[i],
                              prep$test$MONO.[i] * 100)
  stopifnot(abs(result$raw_model_input$MONO. - prep$test$MONO.[i]) < 1e-14,
            identical(result$model_class, balanced_labels[i]))
  expected_zone <- if (actual[i] < bundle$thresholds[["low"]]) "low" else {
    if (actual[i] < bundle$thresholds[["high"]]) "intermediate" else "high"
  }
  stopifnot(identical(result$zone, expected_zone))
  abs(result$score - actual[i])
}, numeric(1))
stopifnot(max(single_errors) < 1e-12)
invalid <- tryCatch(predict_tb_inputs(bundle, 90, 12, 101), error = identity)
stopifnot(inherits(invalid, "error"))

app_environment <- new.env(parent = globalenv())
sys.source("app.R", envir = app_environment, keep.source = FALSE)
shiny::testServer(app_environment$server, {
  session$setInputs(mcv = 90, pdw = 12, mono_percent = 6, calculate = 1)
  stopifnot(abs(result()$score - 0.764838218688965) < 1e-7,
            identical(result()$model_class, "TB"), identical(result()$zone, "high"))
  below <- which.max(replace(actual, actual >= bundle$thresholds[["balanced"]], -Inf))
  session$setInputs(mcv = prep$test$MCV[below], pdw = prep$test$PDW[below],
                    mono_percent = prep$test$MONO.[below] * 100, calculate = 2)
  stopifnot(identical(result()$model_class, "Pneumonia"))
})
example <- predict_tb_inputs(bundle, 90, 12, 6)
lines <- c("Local research calculator smoke check: PASS",
           "Delivery date: 2026-10-04",
           sprintf("Batch prediction comparison: %d test rows", length(actual)),
           sprintf("Maximum difference versus saved CSV: %.17g", batch_error),
           sprintf("Single-input checks, including MONO percentage conversion: %d", length(indices)),
           sprintf("Maximum single-input difference: %.17g", max(single_errors)),
           sprintf("Balanced threshold (validation Youden): %.15f", bundle$thresholds[["balanced"]]),
           paste("Balanced test confusion:", paste(names(balanced_confusion), balanced_confusion, collapse = ", ")),
           "Balanced classification above/below threshold: PASS; supplementary zones unchanged: PASS",
           "Invalid MONO% rejected: PASS", "Shiny testServer reactive prediction: PASS",
           sprintf("Example MCV=90, PDW=12, MONO%%=6: score=%.9f, balanced_class=%s, zone=%s", example$score, example$model_class, example$zone),
           "No model refitting, deployment, or source-manuscript modification was performed.")
writeLines(lines, "results/app_smoke_check.txt", useBytes = TRUE)
cat(paste(lines, collapse = "\n"), "\n")
