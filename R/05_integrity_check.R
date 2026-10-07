## Independent, read-only audit of saved revision outputs. Does not refit prediction models.
## Run from rerun after stages 00-04. Reports are written only under results.
suppressPackageStartupMessages({library(caret);library(recipes);library(themis)
  library(ranger);library(kernlab);library(gbm);library(nnet);library(rpart)
  library(xgboost);library(glmnet)})
checks <- list()
add_check <- function(name, passed, detail="") {
  checks[[length(checks)+1L]] <<- data.frame(check=name,passed=isTRUE(passed),detail=detail)
}
same_num <- function(x,y,tol=1e-9) length(x)==length(y) && all(is.finite(x)) &&
  all(is.finite(y)) && max(abs(x-y))<=tol
read_result <- function(nm) read.csv(file.path("results",nm),stringsAsFactors=FALSE)
auc_parts <- function(y,p) {
  a<-p[y=="TB"]; b<-p[y=="Pneumonia"]
  pair<-outer(a,b,function(x,z)as.numeric(x>z)+.5*as.numeric(x==z))
  list(auc=mean(pair),v10=rowMeans(pair),v01=colMeans(pair),n1=length(a),n0=length(b))
}
auc_stats <- function(y,p) {
  a<-auc_parts(y,p);se<-sqrt(var(a$v10)/a$n1+var(a$v01)/a$n0)
  c(auc=a$auc,low=max(0,a$auc-qnorm(.975)*se),high=min(1,a$auc+qnorm(.975)*se))
}
paired_delong <- function(y,p,q) {
  a<-auc_parts(y,p);b<-auc_parts(y,q)
  d<-a$auc-b$auc
  variance<-var(a$v10-b$v10)/a$n1+var(a$v01-b$v01)/a$n0
  se<-sqrt(variance)
  c(diff=d,low=d-qnorm(.975)*se,high=d+qnorm(.975)*se,
    p=2*pnorm(-abs(d/se)))
}
binary <- function(y,p,cut) {
  z<-p>=cut; pos<-y=="TB"
  tp<-sum(z&pos);fp<-sum(z&!pos);tn<-sum(!z&!pos);fn<-sum(!z&pos)
  c(TP=tp,FP=fp,TN=tn,FN=fn,n=length(y),sensitivity=tp/(tp+fn),
    specificity=tn/(tn+fp),ppv=tp/(tp+fp),npv=tn/(tn+fn),
    accuracy=(tp+tn)/length(y),f1=2*tp/(2*tp+fp+fn))
}
wilson_independent <- function(x,n) {
  z<-qnorm(.975);p<-x/n
  c((p+z*z/(2*n)-z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/(1+z*z/n),
    (p+z*z/(2*n)+z*sqrt(p*(1-p)/n+z*z/(4*n*n)))/(1+z*z/n))
}

raw<-read.csv("TB_Pneumonia.csv",check.names=FALSE)
prep<-readRDS("results/prepared_data.rds");lock<-readRDS("results/analysis_lock.rds")
dedup<-raw[!duplicated(raw),,drop=FALSE];rownames(dedup)<-NULL
add_check("source_rows",nrow(raw)==2399,sprintf("n=%d",nrow(raw)))
add_check("source_no_missing",!anyNA(raw),sprintf("missing=%d",sum(is.na(raw))))
add_check("exact_row_dedup",sum(duplicated(raw))==17 && nrow(dedup)==2382 &&
  sum(raw$species[duplicated(raw)]=="TB")==17,"17 excess identical complete rows, all TB")
add_check("clean_class_counts",sum(dedup$species=="TB")==711 &&
  sum(dedup$species=="Pneumonia")==1671,"TB=711; Pneumonia=1671")
source_columns<-setdiff(names(raw),c("Gender","species"))
add_check("prepared_data_matches_source",identical(dedup[source_columns],prep$data[source_columns]) &&
  identical(as.character(dedup$species),as.character(prep$data$species)),"Original numeric columns preserved after exact-row dedup")
add_check("gender_recode",identical(as.character(prep$data$Gender),ifelse(dedup$Gender==1,"Male","Female")),
          "Matches script recode; source coding still requires author verification")
add_check("derived_LMR_PIV",same_num(prep$data$LMR,dedup$LYM/dedup$MONO) &&
  same_num(prep$data$PIV,dedup$NEU*dedup$PLT*dedup$MONO/dedup$LYM))
formula_candidates<-list(NLR=dedup$NEU/dedup$LYM,PLR=dedup$PLT/dedup$LYM,
  NPR=dedup$NEU/dedup$PLT,NMLR=(dedup$NEU+dedup$MONO)/dedup$LYM,
  MLR=dedup$MONO/dedup$LYM,SIRI=dedup$NEU*dedup$MONO/dedup$LYM,
  SII=dedup$NEU*dedup$PLT/dedup$LYM,dNLR=dedup$NEU/(dedup$WBC-dedup$LYM))
formula_definition<-c(NLR="NEU/LYM",PLR="PLT/LYM",NPR="NEU/PLT",NMLR="(NEU+MONO)/LYM",
  MLR="MONO/LYM",SIRI="NEU*MONO/LYM",SII="NEU*PLT/LYM",dNLR="NEU/(WBC-LYM)")
formula_recheck<-do.call(rbind,lapply(names(formula_candidates),function(nm){
  expected<-formula_candidates[[nm]];observed<-dedup[[nm]]
  residual<-abs(observed-expected);relative<-residual/pmax(abs(expected),1e-12)
  data.frame(variable=nm,formula=formula_definition[nm],n=nrow(dedup),
    maximum_absolute_error=max(residual),maximum_relative_error=max(relative),
    mismatches_relative_gt_1e_7=sum(relative>1e-7),
    mismatches_beyond_rounding=sum(relative>1e-7 & residual>5e-9))
}))
write.csv(formula_recheck,"results/ratio_formula_integrity.csv",row.names=FALSE)
for(i in seq_len(nrow(formula_recheck))){row<-formula_recheck[i,]
  add_check(paste0(row$variable,"_provided_ratio_formula"),row$mismatches_beyond_rounding==0,
    sprintf("%s; maximum relative error %.3g",row$formula,row$maximum_relative_error))}
splits<-list(training=prep$train,validation=prep$validation,test=prep$test)
all_ids<-unlist(lapply(splits,function(x)x$row_id),use.names=FALSE)
add_check("split_disjoint_complete",!anyDuplicated(all_ids) &&
  setequal(all_ids,prep$data$row_id),"1428+476+478=2382, no repeated synthetic row IDs")
expected<-matrix(c(1428,426,1002,476,142,334,478,143,335),nrow=3,byrow=TRUE)
actual<-do.call(rbind,lapply(splits,function(d)c(nrow(d),sum(d$species=="TB"),sum(d$species=="Pneumonia"))))
add_check("split_class_counts",all(unname(actual)==unname(expected)),paste(capture.output(print(actual)),collapse="; "))
assignment<-read_result("split_assignment.csv")
add_check("split_assignment_matches_RDS",all(vapply(names(splits),function(nm)
  setequal(assignment$row_id[assignment$split==nm],splits[[nm]]$row_id),logical(1))))
set.seed(123)
first_recreated<-rsample::initial_split(prep$data,prop=.6,strata=species)
recreated_train<-rsample::training(first_recreated)
second_recreated<-rsample::initial_split(rsample::testing(first_recreated),prop=.5,strata=species)
add_check("seed123_stratified_split_reproduced",identical(recreated_train$row_id,prep$train$row_id) &&
  identical(rsample::training(second_recreated)$row_id,prep$validation$row_id) &&
  identical(rsample::testing(second_recreated)$row_id,prep$test$row_id))
add_check("candidate_definition",length(prep$candidate_features)==25 &&
  !any(c("Age","Gender","MLR","dNLR","species","row_id")%in%prep$candidate_features),
  "25 candidate clinical variables; dNLR/MLR/age/sex excluded from main selector")
add_check("dNLR_source_definition",max(abs(dedup$dNLR-dedup$NEU/(dedup$WBC-dedup$LYM)))<1e-6,
          "Source values match NEU/(WBC-LYM), not standard dNLR")
add_check("standard_dNLR_bad_denominator",sum(dedup$WBC-dedup$NEU<=0)==1,"One row WBC-NEU <=0; row retained, variable excluded")

outer<-read_result("outer_predictions.csv");outer_auc<-read_result("outer_auc.csv")
holdouts<-list()
set.seed(2026);recreated_outer<-caret::createFolds(prep$train$species,k=5,returnTrain=FALSE)
for(i in 1:5){
  ckp<-readRDS(sprintf("models/outer_fold_%d.rds",i))
  sel<-readRDS(sprintf("models/outer_fold_%d_selection.rds",i))
  holdouts[[i]]<-ckp$row_id
  add_check(sprintf("fold_%d_seed2026_holdout_reproduced",i),
    identical(ckp$row_id,prep$train$row_id[recreated_outer[[i]]]))
  add_check(sprintf("fold_%d_selector_checkpoint_consistency",i),identical(ckp$panel,sel$panel) &&
    all(unlist(sel[c("panel","boruta","lasso")],use.names=FALSE)%in%prep$candidate_features) &&
    !"dNLR"%in%unlist(sel[c("panel","boruta","lasso")],use.names=FALSE),paste(sel$panel,collapse=";"))
  add_check(sprintf("fold_%d_truth_and_membership",i),all(ckp$row_id%in%prep$train$row_id) &&
    identical(ckp$diagnosis,as.character(prep$train$species[match(ckp$row_id,prep$train$row_id)])))
  for(method in names(ckp$models)){
    item<-ckp$models[[method]];d<-outer[outer$fold==i&outer$algorithm==method,]
    d<-d[match(ckp$row_id,d$row_id),]
    a<-auc_parts(ckp$diagnosis,item$prob)$auc
    target<-outer_auc$auc[outer_auc$fold==i&outer_auc$algorithm==method]
    add_check(sprintf("fold_%d_%s_predictions_auc",i,method),
      identical(d$row_id,ckp$row_id) && same_num(d$probability_TB,item$prob) &&
      same_num(c(item$auc,target),c(a,a)))
  }
}
add_check("outer_holdouts_disjoint_complete",!anyDuplicated(unlist(holdouts)) &&
  setequal(unlist(holdouts),prep$train$row_id),"Each original training row assessed exactly once per algorithm")
ranking<-read_result("model_ranking.csv")
computed_rank<-aggregate(auc~algorithm,outer_auc,mean)
add_check("model_ranking_and_winner",same_num(ranking$mean_outer_auc,
  computed_rank$auc[match(ranking$algorithm,computed_rank$algorithm)]) &&
  identical(ranking$algorithm[1],lock$winner),paste("winner",lock$winner))
full_sel<-readRDS("models/full_training_selection.rds")
add_check("full_panel_matches_lock",identical(full_sel$panel,lock$panel) &&
  !"dNLR"%in%unlist(full_sel[c("panel","boruta","lasso")],use.names=FALSE),paste(lock$panel,collapse=";"))

val<-read_result("validation_predictions.csv");test<-read_result("test_predictions.csv")
extra_val<-read_result("historical_five_age_sex_validation_predictions.csv")
extra_test<-read_result("historical_five_age_sex_test_predictions.csv")
all_val<-rbind(val,extra_val);all_test<-rbind(test,extra_test)
for(set_nm in c("validation","test")){
  saved<-if(set_nm=="validation")all_val else all_test
  reference<-splits[[set_nm]]
  for(method in unique(saved$algorithm)){
    d<-saved[saved$algorithm==method,]
    add_check(paste(set_nm,method,"case_alignment",sep="_"),nrow(d)==nrow(reference) &&
      !anyDuplicated(d$row_id) && setequal(d$row_id,reference$row_id) &&
      identical(d$diagnosis,as.character(reference$species[match(d$row_id,reference$row_id)])) &&
      all(is.finite(d$probability_TB)&d$probability_TB>=0&d$probability_TB<=1))
    path<-switch(method,age_sex="models/age_sex_xgbTree.rds",
      panel_age_sex="models/panel_age_sex_xgbTree.rds",
      historical_five_xgbTree="models/historical_five_xgbTree.rds",
      historical_five_age_sex="models/historical_five_age_sex_xgbTree.rds",
      paste0("models/final_",method,".rds"))
    fit<-readRDS(path)
    again<-predict(fit,newdata=reference,type="prob")$TB
    add_check(paste(set_nm,method,"saved_model_prediction",sep="_"),
      same_num(again,d$probability_TB[match(reference$row_id,d$row_id)],1e-7))
    if(set_nm=="test"){
      steps<-fit$recipe$steps
      smote<-which(vapply(steps,function(s)inherits(s,"step_smote"),logical(1)))
      add_check(paste(method,"SMOTE_skip_prediction",sep="_"),length(smote)==1 && isTRUE(steps[[smote]]$skip))
      add_check(paste(method,"trainingData_class_n",sep="_"),nrow(fit$trainingData)==1428 &&
        sum(fit$trainingData$species=="TB")==426 && sum(fit$trainingData$species=="Pneumonia")==1002)
      numeric_predictors<-names(fit$trainingData)[vapply(fit$trainingData,is.numeric,logical(1))]
      normalize_step<-which(vapply(steps,function(s)inherits(s,"step_normalize"),logical(1)))
      st<-steps[[normalize_step]]
      add_check(paste(method,"normalization_original_training_only",sep="_"),
        same_num(st$means[numeric_predictors],vapply(prep$train[numeric_predictors],mean,numeric(1))) &&
        same_num(st$sds[numeric_predictors],vapply(prep$train[numeric_predictors],sd,numeric(1))))
      add_check(paste(method,"recipe_no_ID_predictor",sep="_"),
        !any(c("row_id","dNLR","species")%in%fit$recipe$term_info$variable[fit$recipe$term_info$role=="predictor"]))
    }
  }
}
main_v<-val[val$algorithm==lock$winner,];main_t<-test[test$algorithm==lock$winner,]
thresholds<-read_result("locked_thresholds.csv");threshold_table<-read_result("threshold_performance.csv")
for(i in seq_len(nrow(threshold_table))){
  row<-threshold_table[i,];d<-if(row$dataset=="validation")main_v else main_t
  s<-binary(d$diagnosis,d$probability_TB,row$threshold)
  add_check(paste(row$dataset,row$rule,"confusion_metrics",sep="_"),
    same_num(unlist(row[names(s)],use.names=FALSE),s))
  ci<-c(wilson_independent(s["TP"],s["TP"]+s["FN"]),
    wilson_independent(s["TN"],s["TN"]+s["FP"]),
    wilson_independent(s["TP"],s["TP"]+s["FP"]),
    wilson_independent(s["TN"],s["TN"]+s["FN"]),
    wilson_independent(s["TP"]+s["TN"],s["n"]))
  ci_names<-c("sensitivity_low","sensitivity_high","specificity_low","specificity_high",
    "ppv_low","ppv_high","npv_low","npv_high","accuracy_low","accuracy_high")
  add_check(paste(row$dataset,row$rule,"Wilson_CI",sep="_"),same_num(unlist(row[ci_names],use.names=FALSE),ci))
  set.seed(17001L)
  f1_boot<-replicate(2000L,{idx<-sample.int(nrow(d),nrow(d),replace=TRUE)
    binary(d$diagnosis[idx],d$probability_TB[idx],row$threshold)["f1"]})
  add_check(paste(row$dataset,row$rule,"F1_bootstrap_CI",sep="_"),
    same_num(c(row$f1_low,row$f1_high),quantile(f1_boot,c(.025,.975),na.rm=TRUE)))
}
low<-thresholds$threshold[thresholds$rule=="sensitivity_90"]
high<-thresholds$threshold[thresholds$rule=="rulein_specificity_90"]
possible<-sort(unique(c(0,main_v$probability_TB,1)))
sens<-vapply(possible,function(x)binary(main_v$diagnosis,main_v$probability_TB,x)["sensitivity"],numeric(1))
spec<-vapply(possible,function(x)binary(main_v$diagnosis,main_v$probability_TB,x)["specificity"],numeric(1))
add_check("validation_sensitivity90_cutoff_optimal",same_num(low,max(possible[sens>=.90])),
          sprintf("cutoff %.12f; sensitivity %.9f",low,binary(main_v$diagnosis,main_v$probability_TB,low)["sensitivity"]))
add_check("validation_specificity90_cutoff_optimal",same_num(high,min(possible[spec>=.90])),
          sprintf("cutoff %.12f; specificity %.9f",high,binary(main_v$diagnosis,main_v$probability_TB,high)["specificity"]))
add_check("three_zone_threshold_order",low<high)
zones<-read_result("three_zone_composition.csv")
for(set_nm in c("validation","test")){
  d<-if(set_nm=="validation")main_v else main_t
  zone<-ifelse(d$probability_TB<low,"blue",ifelse(d$probability_TB>=high,"red","gray"))
  z<-zones[zones$dataset==set_nm,]
  n<-vapply(seq_len(nrow(z)),function(i)sum(zone==z$zone[i]&d$diagnosis==z$diagnosis[i]),integer(1))
  add_check(paste0(set_nm,"_three_zone_counts"),identical(z$n,n) && sum(z$n)==nrow(d))
}
metrics<-read_result("model_test_metrics.csv")
for(i in seq_len(nrow(metrics))){
  row<-metrics[i,];d<-test[test$algorithm==row$algorithm,]
  a<-auc_stats(d$diagnosis,d$probability_TB)
  add_check(paste0(row$algorithm,"_test_manual_AUC_DeLongCI_Brier"),
    same_num(unlist(row[c("auc","auc_low","auc_high")],use.names=FALSE),a) &&
      same_num(row$brier,mean((d$probability_TB-as.numeric(d$diagnosis=="TB"))^2)))
}
paired<-read_result("paired_auc_comparisons.csv")
for(i in seq_len(nrow(paired))){
  row<-paired[i,];d<-test[test$algorithm==row$comparator,];d<-d[match(main_t$row_id,d$row_id),]
  s<-paired_delong(main_t$diagnosis,main_t$probability_TB,d$probability_TB)
  add_check(paste0(row$comparator,"_paired_manual_DeLong"),
    same_num(unlist(row[c("auc_difference","difference_low","difference_high","p_value")],use.names=FALSE),s))
}
add_check("paired_p_Holm",same_num(paired$p_holm,p.adjust(paired$p_value,method="holm")))
old5<-all_test[all_test$algorithm=="historical_five_xgbTree",]
old5extra<-all_test[all_test$algorithm=="historical_five_age_sex",]
old5extra<-old5extra[match(old5$row_id,old5extra$row_id),]
old5_stats<-paired_delong(old5$diagnosis,old5$probability_TB,old5extra$probability_TB)
old5_saved<-read_result("historical_demographic_paired_test.csv")
add_check("historical5_age_sex_paired_manual_DeLong",same_num(c(old5_saved$auc_difference,old5_saved$p_value),
  old5_stats[c("diff","p")]))
write.csv(data.frame(comparison="original five versus original five plus age/sex",
  auc_original5=auc_parts(old5$diagnosis,old5$probability_TB)$auc,
  auc_original5_age_sex=auc_parts(old5extra$diagnosis,old5extra$probability_TB)$auc,
  auc_difference=unname(old5_stats["diff"]),difference_low=unname(old5_stats["low"]),
  difference_high=unname(old5_stats["high"]),p_value=unname(old5_stats["p"])),
  "results/historical5_demographic_paired_audit.csv",row.names=FALSE)
subgroup<-read_result("test_subgroup_metrics.csv")
demo<-prep$test[match(main_t$row_id,prep$test$row_id),]
for(i in seq_len(nrow(subgroup))){
  row<-subgroup[i,]
  keep<-switch(row$group,all=rep(TRUE,nrow(demo)),under_65=demo$Age<65,
    age_65_or_older=demo$Age>=65,male=demo$Gender=="Male",female=demo$Gender=="Female")
  d<-main_t[keep,];a<-auc_stats(d$diagnosis,d$probability_TB);s<-binary(d$diagnosis,d$probability_TB,low)
  add_check(paste0(row$group,"_subgroup_auc_counts_metrics"),row$n==nrow(d) &&
    row$tb_n==sum(d$diagnosis=="TB") && row$pneumonia_n==sum(d$diagnosis=="Pneumonia") &&
    same_num(unlist(row[c("auc","auc_low","auc_high")],use.names=FALSE),a) &&
    same_num(unlist(row[c("sensitivity","specificity")],use.names=FALSE),s[c("sensitivity","specificity")]))
}
cal<-read_result("test_calibration.csv");y<-as.numeric(main_t$diagnosis=="TB")
p<-pmin(pmax(main_t$probability_TB,1e-6),1-1e-6);logit<-qlogis(p)
ci_fit<-glm(y~offset(logit),family=binomial());slope_fit<-glm(y~logit,family=binomial())
add_check("calibration_coefficients_Wald_CI",same_num(unlist(cal[c("intercept","intercept_low","intercept_high",
  "slope","slope_low","slope_high")],use.names=FALSE),
  c(coef(ci_fit)[1],confint.default(ci_fit)[1,],coef(slope_fit)[2],confint.default(slope_fit)[2,])))
set.seed(16001);bb<-replicate(2000,{ii<-sample.int(length(y),length(y),TRUE);mean((p[ii]-y[ii])^2)})
add_check("Brier_CI_saved_bootstrap",same_num(c(cal$brier,cal$brier_low,cal$brier_high),
  c(mean((p-y)^2),quantile(bb,c(.025,.975)))))
dca<-read_result("test_decision_curve.csv")
expected_nb<-vapply(dca$threshold,function(t){z<-main_t$probability_TB>=t;
  (sum(z&y==1)-sum(z&y==0)*t/(1-t))/length(y)},numeric(1))
add_check("DCA_all_curves",same_num(dca$model,expected_nb) &&
  same_num(dca$treat_all,mean(y)-(1-mean(y))*dca$threshold/(1-dca$threshold)) && all(dca$treat_none==0))
merged_for_cal<-merge(main_t,prep$test[,c("row_id","Age","Gender")],by="row_id",sort=FALSE)
cy<-as.numeric(merged_for_cal$diagnosis=="TB")
cp<-pmin(pmax(merged_for_cal$probability_TB,1e-6),1-1e-6)
set.seed(16002)
calboot<-replicate(1000,{idx<-sample.int(length(cy),length(cy),replace=TRUE)
  yy<-cy[idx];ll<-qlogis(cp[idx]);
  c(unname(coef(suppressWarnings(glm(yy~offset(ll),family=binomial())))[1]),
    unname(coef(suppressWarnings(glm(yy~ll,family=binomial())))[2]))})
calci<-read_result("calibration_bootstrap_intervals.csv")
expected_cal<-t(apply(calboot,1,quantile,c(.025,.975),na.rm=TRUE))
add_check("calibration_bootstrap_CI",same_num(c(calci$lower,calci$upper),
  c(expected_cal[,1],expected_cal[,2])) && all(calci$valid_replicates==rowSums(is.finite(calboot))))

mainfit<-readRDS("models/final_xgbTree.rds")
shap<-readRDS("results/test_shap.rds")
shap_test<-prep$test[match(shap$row_id,prep$test$row_id),]
baked<-recipes::bake(mainfit$recipe,new_data=shap_test)
predictors<-names(baked)[names(baked)!="species"]
matrix_input<-as.matrix(baked[,predictors,drop=FALSE])
contrib<-predict(mainfit$finalModel,newdata=matrix_input,predcontrib=TRUE)
again_prob<-predict(mainfit,newdata=shap_test,type="prob")$TB
add_check("SHAP_stored_contributions_reproduced",same_num(as.numeric(contrib),as.numeric(shap$contribution),1e-7))
add_check("SHAP_additive_prediction_logodds",same_num(plogis(rowSums(contrib)),again_prob,1e-6) &&
  same_num(shap$probability,again_prob,1e-7),"Contributions plus baseline reconstruct TB probabilities within 1e-6")

source_hash<-digest::digest(file="TB_Pneumonia.csv",algo="sha256")
add_check("original_CSV_SHA256_unchanged",identical(source_hash,
  "87afdd4f6c0e9ad2d59b9684886f04e169a21f714dc5c450db7f45ee45d59867"),source_hash)
audit_files<-c("TB_Pneumonia.csv",list.files(pattern="^[0-5][0-9]_.*[.]R$"),
  list.files("models",pattern="[.]rds$",full.names=TRUE))
provenance<-data.frame(file=audit_files,bytes=file.info(audit_files)$size,
  sha256=vapply(audit_files,function(f)digest::digest(file=f,algo="sha256"),character(1)))
write.csv(provenance,"results/integrity_provenance.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"results/integrity_session_info.txt")

all_checks<-do.call(rbind,checks)
write.csv(all_checks,"results/integrity_checks.csv",row.names=FALSE)
writeLines(c(sprintf("Read-only saved-output audit: %d/%d checks passed",sum(all_checks$passed),nrow(all_checks)),
  if(all(all_checks$passed))"No numerical/structural inconsistency detected by these checks." else
    paste("FAILED:",paste(all_checks$check[!all_checks$passed],collapse=", ")),
  "Manual AUC/DeLong were calculated from positive-negative pairwise comparisons, independently of pROC.",
  "Saved model predictions were reproduced without model fitting.",
  "Limits: exact-row dedup does not verify unique patients; no date/site/screening/30-case identifiers exist.",
  "Selectors precede inner tuning and combination CV. Only outer held-out scores evaluate the entire selection route.",
  "The XGBoost combination scorer favors an XGBoost-oriented feature search; algorithm superiority must be qualified.",
  "Historical test-set exposure and reconstructed selection steps make this an internal reanalysis, not new external validation.",
  "Raw stage 01 checkpoint reuse is file-based. The delivered run_all.R verifies an immutable input/code/model manifest; --retrain uses a clean directory."),
  "results/integrity_summary.txt",useBytes=TRUE)
print(all_checks[!all_checks$passed,])
cat(sprintf("Audit complete: %d/%d passed.\n",sum(all_checks$passed),nrow(all_checks)))
if(!all(all_checks$passed))quit(status=1)
