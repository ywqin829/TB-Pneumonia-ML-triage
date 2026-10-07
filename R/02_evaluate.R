## Locked validation threshold and one-pass test-set evaluation.
## Run after 01_model_pipeline.R from the rerun directory.

suppressPackageStartupMessages({ library(pROC) })
dir.create("results", showWarnings = FALSE)
prep <- readRDS("results/prepared_data.rds")
lock <- readRDS("results/analysis_lock.rds")
val <- read.csv("results/validation_predictions.csv", stringsAsFactors = FALSE)
test <- read.csv("results/test_predictions.csv", stringsAsFactors = FALSE)
winner <- lock$winner
pos <- "TB"
main_val <- val[val$algorithm == winner, ]
main_test <- test[test$algorithm == winner, ]
stopifnot(nrow(main_val) == nrow(prep$validation),
          nrow(main_test) == nrow(prep$test),
          setequal(main_val$row_id, prep$validation$row_id),
          setequal(main_test$row_id, prep$test$row_id))

truth_bin <- function(truth) as.integer(truth == pos)
confusion <- function(truth, prob, threshold) {
  pred <- as.integer(prob >= threshold)
  obs <- truth_bin(truth)
  c(TP = sum(pred == 1 & obs == 1), FP = sum(pred == 1 & obs == 0),
    TN = sum(pred == 0 & obs == 0), FN = sum(pred == 0 & obs == 1))
}
safe_div <- function(x, y) if (y > 0) unname(x/y) else NA_real_
wilson <- function(x, n, z = qnorm(0.975)) {
  if (n == 0) return(c(lower = NA_real_, upper = NA_real_))
  p <- x/n; den <- 1 + z^2/n
  center <- (p + z^2/(2*n))/den
  half <- z * sqrt(p*(1-p)/n + z^2/(4*n^2))/den
  c(lower = max(0, center-half), upper = min(1, center+half))
}
binary_stats <- function(truth, prob, threshold) {
  x <- confusion(truth, prob, threshold)
  n <- sum(x)
  sensitivity <- safe_div(x["TP"], x["TP"] + x["FN"])
  specificity <- safe_div(x["TN"], x["TN"] + x["FP"])
  ppv <- safe_div(x["TP"], x["TP"] + x["FP"])
  npv <- safe_div(x["TN"], x["TN"] + x["FN"])
  precision <- ppv
  recall <- sensitivity
  f1 <- if (is.na(precision) || is.na(recall) || precision + recall == 0)
    NA_real_ else 2*precision*recall/(precision+recall)
  c(x, n=n, sensitivity=sensitivity, specificity=specificity,
    ppv=ppv, npv=npv, accuracy=unname((x["TP"]+x["TN"])/n), f1=f1)
}
auc_ci <- function(truth, prob) {
  roc <- pROC::roc(response = truth, predictor = prob,
                    levels = c("Pneumonia", "TB"), direction = "<", quiet=TRUE)
  c(auc = as.numeric(pROC::auc(roc)),
    lower = unname(pROC::ci.auc(roc, method="delong")[1]),
    upper = unname(pROC::ci.auc(roc, method="delong")[3]))
}
candidate_thresholds <- function(prob) sort(unique(c(0, prob, 1)))
choose_threshold <- function(truth, prob, criterion, target) {
  thresholds <- candidate_thresholds(prob)
  stats <- t(vapply(thresholds, function(t) binary_stats(truth,prob,t),
                    numeric(11)))
  if (criterion == "sensitivity") {
    keep <- which(stats[, "sensitivity"] >= target)
    if (length(keep) == 0) return(NA_real_)
    return(max(thresholds[keep]))
  }
  if (criterion == "specificity") {
    keep <- which(stats[, "specificity"] >= target)
    if (length(keep) == 0) return(NA_real_)
    return(min(thresholds[keep]))
  }
  stop("Unknown criterion")
}
choose_youden <- function(truth, prob) {
  thresholds <- candidate_thresholds(prob)
  stats <- t(vapply(thresholds, function(t) binary_stats(truth,prob,t),
                    numeric(11)))
  score <- stats[, "sensitivity"] + stats[, "specificity"] - 1
  thresholds[which.max(score)]
}

## Sensitivity target is prespecified before test analysis. WHO's TB
## screening TPP lists 90% sensitivity as a minimum for high-sensitivity
## screening; this study's hospital sample is not a TPP qualification trial.
target_sensitivity <- 0.90
target_rulein_specificity <- 0.90
thresholds <- data.frame(
  rule = c("sensitivity_90", "rulein_specificity_90", "youden", "sensitivity_85_reference",
           "historical_0.232", "historical_0.412"),
  threshold = c(
    choose_threshold(main_val$diagnosis, main_val$probability_TB,
                     "sensitivity", target_sensitivity),
    choose_threshold(main_val$diagnosis, main_val$probability_TB,
                     "specificity", target_rulein_specificity),
    choose_youden(main_val$diagnosis, main_val$probability_TB),
    choose_threshold(main_val$diagnosis, main_val$probability_TB,
                     "sensitivity", 0.85),
    0.232, 0.412))
thresholds$source <- c("validation", "validation", "validation", "validation",
                       "historical manuscript only", "historical manuscript only")
write.csv(thresholds, "results/locked_thresholds.csv", row.names=FALSE)
saveRDS(thresholds, "models/locked_thresholds.rds")

threshold_results <- list()
for (i in seq_len(nrow(thresholds))) {
  threshold <- thresholds$threshold[i]
  for (data_name in c("validation", "test")) {
    d <- if (data_name == "validation") main_val else main_test
    s <- binary_stats(d$diagnosis, d$probability_TB, threshold)
    cis <- c(wilson(s["TP"], s["TP"] + s["FN"]),
             wilson(s["TN"], s["TN"] + s["FP"]),
             wilson(s["TP"], s["TP"] + s["FP"]),
             wilson(s["TN"], s["TN"] + s["FN"]),
             wilson(s["TP"] + s["TN"], s["n"]))
    names(cis) <- c("sensitivity_low","sensitivity_high",
                    "specificity_low","specificity_high",
                    "ppv_low","ppv_high","npv_low","npv_high",
                    "accuracy_low","accuracy_high")
    set.seed(17001L)
    f1_boot <- replicate(2000L, {
      idx <- sample.int(nrow(d), nrow(d), replace=TRUE)
      binary_stats(d$diagnosis[idx], d$probability_TB[idx], threshold)["f1"]
    })
    f1_ci <- quantile(f1_boot, c(0.025,0.975), na.rm=TRUE)
    threshold_results[[length(threshold_results)+1L]] <- data.frame(
      dataset=data_name, rule=thresholds$rule[i], threshold=threshold,
      as.list(s), as.list(cis), f1_low=unname(f1_ci[1]),
      f1_high=unname(f1_ci[2]), check.names=FALSE)
  }
}
threshold_results <- do.call(rbind, threshold_results)
write.csv(threshold_results, "results/threshold_performance.csv", row.names=FALSE)

low <- thresholds$threshold[thresholds$rule == "sensitivity_90"]
high <- thresholds$threshold[thresholds$rule == "rulein_specificity_90"]
zone_status <- if (is.finite(low) && is.finite(high) && low < high)
  "valid_three_zone" else "no_ordered_three_zone"
zone_table <- data.frame()
if (zone_status == "valid_three_zone") {
  for (data_name in c("validation", "test")) {
    d <- if (data_name == "validation") main_val else main_test
    zone <- ifelse(d$probability_TB < low, "blue",
                   ifelse(d$probability_TB >= high, "red", "gray"))
    tab <- as.data.frame(table(zone=factor(zone,levels=c("blue","gray","red")),
                               diagnosis=factor(d$diagnosis,levels=c("TB","Pneumonia"))))
    names(tab)[3] <- "n"
    tab$dataset <- data_name
    tab$zone_total <- ave(tab$n, tab$zone, FUN=sum)
    tab$zone_fraction <- tab$zone_total/nrow(d)
    zone_table <- rbind(zone_table, tab)
  }
}
write.csv(zone_table, "results/three_zone_composition.csv", row.names=FALSE)
writeLines(c(paste("zone status:",zone_status),
             paste("lower sensitivity threshold:",low),
             paste("upper specificity threshold:",high)),
           "results/zone_status.txt")

## AUC and confidence interval for all prespecified models on the same
## frozen test cohort. No test result is used to reselect the winner.
model_test_metrics <- do.call(rbind, lapply(unique(test$algorithm), function(nm) {
  d <- test[test$algorithm == nm, ]
  a <- auc_ci(d$diagnosis, d$probability_TB)
  data.frame(algorithm=nm, n=nrow(d), auc=unname(a["auc"]),
             auc_low=unname(a["lower"]), auc_high=unname(a["upper"]),
             brier=mean((d$probability_TB-truth_bin(d$diagnosis))^2))
}))
write.csv(model_test_metrics, "results/model_test_metrics.csv", row.names=FALSE)

main_roc <- pROC::roc(main_test$diagnosis,main_test$probability_TB,
                       levels=c("Pneumonia","TB"),direction="<",quiet=TRUE)
paired_tests <- list()
comparators <- setdiff(unique(test$algorithm), winner)
for (nm in comparators) {
  other <- test[test$algorithm == nm, ]
  other <- other[match(main_test$row_id, other$row_id), ]
  stopifnot(identical(main_test$row_id, other$row_id))
  other_roc <- pROC::roc(other$diagnosis,other$probability_TB,
                         levels=c("Pneumonia","TB"),direction="<",quiet=TRUE)
  comparison <- pROC::roc.test(main_roc, other_roc,
                                method="delong",paired=TRUE)
  difference_ci <- if (!is.null(comparison$conf.int)) unname(comparison$conf.int) else c(NA_real_,NA_real_)
  paired_tests[[nm]] <- data.frame(
    main=winner, comparator=nm,
    auc_main=as.numeric(pROC::auc(main_roc)),
    auc_comparator=as.numeric(pROC::auc(other_roc)),
    auc_difference=as.numeric(pROC::auc(main_roc)-pROC::auc(other_roc)),
    difference_low=difference_ci[1], difference_high=difference_ci[2],
    p_value=as.numeric(comparison$p.value))
}
paired_tests <- do.call(rbind,paired_tests)
paired_tests$p_holm <- p.adjust(paired_tests$p_value,method="holm")
write.csv(paired_tests,"results/paired_auc_comparisons.csv",row.names=FALSE)

## Test-set calibration assessment only: no coefficients are applied back
## to predictions and the test labels never enter model training.
y_test <- truth_bin(main_test$diagnosis)
p_test <- pmin(pmax(main_test$probability_TB,1e-6),1-1e-6)
logit <- qlogis(p_test)
cal_intercept_fit <- glm(y_test ~ offset(logit),family=binomial())
cal_slope_fit <- glm(y_test ~ logit,family=binomial())
calibration <- data.frame(
  n=length(y_test), brier=mean((p_test-y_test)^2),
  intercept=unname(coef(cal_intercept_fit)[1]),
  intercept_low=confint.default(cal_intercept_fit)[1,1],
  intercept_high=confint.default(cal_intercept_fit)[1,2],
  slope=unname(coef(cal_slope_fit)[2]),
  slope_low=confint.default(cal_slope_fit)[2,1],
  slope_high=confint.default(cal_slope_fit)[2,2])
set.seed(16001)
boot_brier <- replicate(2000, {
  idx <- sample.int(length(y_test),length(y_test),replace=TRUE)
  mean((p_test[idx]-y_test[idx])^2)
})
calibration$brier_low <- unname(quantile(boot_brier,0.025))
calibration$brier_high <- unname(quantile(boot_brier,0.975))
write.csv(calibration,"results/test_calibration.csv",row.names=FALSE)

## Net benefit on the same test patients for model, treat all, treat none.
grid <- seq(0.01,0.80,by=0.01)
prevalence <- mean(y_test)
dca <- do.call(rbind,lapply(grid,function(pt) {
  x <- confusion(main_test$diagnosis,main_test$probability_TB,pt)
  weight <- pt/(1-pt)
  data.frame(threshold=pt,
             model=unname(x["TP"]/length(y_test)-x["FP"]/length(y_test)*weight),
             treat_all=prevalence-(1-prevalence)*weight,
             treat_none=0)
}))
write.csv(dca,"results/test_decision_curve.csv",row.names=FALSE)

## Prespecified demographic subgroups; same frozen model and same cutoff.
demographics <- prep$test[,c("row_id","Age","Gender")]
main_test <- merge(main_test,demographics,by="row_id",sort=FALSE)
main_test$age_group <- ifelse(main_test$Age<65,"under_65","65_or_older")
main_test$gender_group <- as.character(main_test$Gender)
subgroups <- list(all=rep(TRUE,nrow(main_test)),
                  under_65=main_test$age_group=="under_65",
                  age_65_or_older=main_test$age_group=="65_or_older",
                  male=main_test$gender_group=="Male",
                  female=main_test$gender_group=="Female")
subgroup_metrics <- list()
for (nm in names(subgroups)) {
  d <- main_test[subgroups[[nm]],]
  s <- binary_stats(d$diagnosis,d$probability_TB,low)
  a <- auc_ci(d$diagnosis,d$probability_TB)
  subgroup_metrics[[nm]] <- data.frame(group=nm,n=nrow(d),
                                       tb_n=sum(d$diagnosis==pos),
                                       pneumonia_n=sum(d$diagnosis!="TB"),
                                       auc=a["auc"],auc_low=a["lower"],auc_high=a["upper"],
                                       sensitivity=s["sensitivity"],
                                       specificity=s["specificity"],
                                       sensitivity_low=wilson(s["TP"],s["TP"]+s["FN"])[1],
                                       sensitivity_high=wilson(s["TP"],s["TP"]+s["FN"])[2],
                                       specificity_low=wilson(s["TN"],s["TN"]+s["FP"])[1],
                                       specificity_high=wilson(s["TN"],s["TN"]+s["FP"])[2])
}
subgroup_metrics <- do.call(rbind,subgroup_metrics)
write.csv(subgroup_metrics,"results/test_subgroup_metrics.csv",row.names=FALSE)

## For the age/sex comparators use each model's validation-derived 90%
## sensitivity cutoff, rather than reusing a probability cutoff from a
## differently calibrated model.
demographic_metrics <- list()
for (nm in c("age_sex","panel_age_sex")) {
  v <- val[val$algorithm==nm,]
  t <- test[test$algorithm==nm,]
  cut <- choose_threshold(v$diagnosis,v$probability_TB,"sensitivity",0.90)
  s <- binary_stats(t$diagnosis,t$probability_TB,cut)
  demographic_metrics[[nm]] <- data.frame(model=nm,validation_cutoff=cut,
                                          as.list(s),check.names=FALSE)
}
write.csv(do.call(rbind,demographic_metrics),
          "results/demographic_comparator_metrics.csv",row.names=FALSE)

writeLines(c(sprintf("Main model: %s",winner),
             sprintf("Selected features: %s",paste(lock$panel,collapse=", ")),
             sprintf("Training/validation/test: %d/%d/%d",lock$training_n,
                     lock$validation_n,lock$test_n),
             sprintf("Validation sensitivity>=90%% cutoff: %.6f",low),
             sprintf("Validation specificity>=90%% rule-in cutoff: %.6f",high),
             sprintf("Three-zone status: %s",zone_status),
             sprintf("Test AUC: %.4f (95%% CI %.4f-%.4f)",
                     model_test_metrics$auc[model_test_metrics$algorithm==winner],
                     model_test_metrics$auc_low[model_test_metrics$algorithm==winner],
                     model_test_metrics$auc_high[model_test_metrics$algorithm==winner])),
           "results/evaluation_summary.txt",useBytes=TRUE)
writeLines(capture.output(sessionInfo()),"results/session_info.txt")
print(threshold_results[threshold_results$dataset=="test",
                        c("rule","threshold","TP","FP","TN","FN",
                          "sensitivity","specificity","ppv","npv")])
print(model_test_metrics)
