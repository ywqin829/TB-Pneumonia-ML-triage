## Supplemental descriptive, stability, demographic and reproducibility outputs.
## This sources the deterministic evaluation stage; no main model is refitted.
suppressPackageStartupMessages({library(caret);library(recipes);library(themis)})
source("02_evaluate.R")
clean <- prep$data

fmt <- function(x) sprintf("%.3f [%.3f, %.3f]",median(x),quantile(x,.25),quantile(x,.75))
clinical <- setdiff(names(clean),c("row_id","Gender","species"))
baseline <- do.call(rbind,lapply(clinical,function(nm){
  a<-clean[[nm]][clean$species=="TB"];b<-clean[[nm]][clean$species=="Pneumonia"]
  smd<-(mean(a)-mean(b))/sqrt((var(a)+var(b))/2)
  data.frame(variable=nm,TB=fmt(a),Pneumonia=fmt(b),
             p_wilcoxon=suppressWarnings(wilcox.test(a,b,exact=FALSE)$p.value),
             smd=smd,in_model_candidates=nm%in%prep$candidate_features)
}))
write.csv(baseline,"results/baseline_TB_pneumonia.csv",row.names=FALSE)
gender_table<-as.data.frame(table(diagnosis=clean$species,gender=clean$Gender))
names(gender_table)[3]<-"n"
gender_table$fraction<-ave(gender_table$n,gender_table$diagnosis,FUN=function(x)x/sum(x))
gender_table$p_chisquare<-chisq.test(table(clean$species,clean$Gender))$p.value
write.csv(gender_table,"results/baseline_gender.csv",row.names=FALSE)
split_assignment<-read.csv("results/split_assignment.csv")
split_name<-setNames(as.character(split_assignment$split),split_assignment$row_id)
clean$split<-unname(split_name[clean$row_id])
split_baseline<-do.call(rbind,lapply(clinical,function(nm){
  x<-split(clean[[nm]],clean$split)
  data.frame(variable=nm,training=fmt(x$training),validation=fmt(x$validation),test=fmt(x$test),
             p_kruskal=kruskal.test(clean[[nm]],as.factor(clean$split))$p.value)
}))
write.csv(split_baseline,"results/baseline_by_split.csv",row.names=FALSE)
write.csv(as.data.frame(table(split=clean$split,gender=clean$Gender)),
          "results/gender_by_split.csv",row.names=FALSE)
age_summary<-do.call(rbind,lapply(levels(clean$species),function(g){
  x<-clean$Age[clean$species==g]
  data.frame(diagnosis=g,n=length(x),mean=mean(x),sd=sd(x),median=median(x),
             q1=quantile(x,.25),q3=quantile(x,.75),
             skewness=mean((x-mean(x))^3)/(mean((x-mean(x))^2)^1.5))
}))
write.csv(age_summary,"results/age_distribution_summary.csv",row.names=FALSE)

formula_defs<-c(NLR="NEU/LYM",PLR="PLT/LYM",NPR="NEU/PLT",NMLR="(NEU+MONO)/LYM",
                 MLR="MONO/LYM",LMR="LYM/MONO",SIRI="NEU*MONO/LYM",SII="NEU*PLT/LYM",
                 PIV="NEU*PLT*MONO/LYM",dNLR="Provided values equal NEU/(WBC-LYM); excluded")
dictionary<-data.frame(variable=setdiff(names(clean),c("row_id","split")),
                        definition="source CSV field",status="available")
for(nm in names(formula_defs)) dictionary$definition[dictionary$variable==nm]<-formula_defs[nm]
dictionary$status[dictionary$variable=="dNLR"]<-"excluded: definition mismatch"
dictionary$status[dictionary$variable=="MLR"]<-"descriptive only: reciprocal of included LMR"
dictionary$status[dictionary$variable%in%c("LMR","PIV")]<-"reconstructed from supplied absolute counts"
dictionary$definition[dictionary$variable=="Gender"]<-"1=Male, 2=Female; consistent with original Figure 2 percentages"
dictionary$definition[dictionary$variable%in%c("LYM.","NEU.","MONO.")]<-"fraction 0-1; multiply by 100 for percentage display"
dictionary$model_candidate<-dictionary$variable%in%prep$candidate_features
write.csv(dictionary,"results/variable_dictionary.csv",row.names=FALSE)

stability<-data.frame(feature=prep$candidate_features,Boruta=0L,LASSO=0L,final_panel=0L)
for(i in 1:5){
  s<-readRDS(sprintf("models/outer_fold_%d_selection.rds",i))
  stability$Boruta<-stability$Boruta+as.integer(stability$feature%in%s$boruta)
  stability$LASSO<-stability$LASSO+as.integer(stability$feature%in%s$lasso)
  stability$final_panel<-stability$final_panel+as.integer(stability$feature%in%s$panel)
}
stability$final_proportion<-stability$final_panel/5
write.csv(stability,"results/feature_selection_stability.csv",row.names=FALSE)

## Original five-feature panel plus age/sex: preserves the comparator that
## reviewers requested even though the reselected main panel has changed.
legacy_extra_path<-"models/historical_five_age_sex_xgbTree.rds"
if(file.exists(legacy_extra_path)){
  legacy_extra<-readRDS(legacy_extra_path)
} else {
  input<-prep$train[,c("species",lock$legacy_panel,"Age","Gender")]
  r<-recipe(species~.,data=input) |>
    step_dummy(all_nominal_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    step_smote(species,over_ratio=1,neighbors=5)
  set.seed(13002)
  ctrl<-trainControl(method="cv",index=createFolds(input$species,k=3,returnTrain=TRUE),
                      summaryFunction=twoClassSummary,classProbs=TRUE,
                      savePredictions="final",allowParallel=FALSE)
  grid<-expand.grid(nrounds=c(100,250),max_depth=c(2,3),eta=c(.05,.1),gamma=0,
                     colsample_bytree=.8,min_child_weight=1,subsample=.8)
  set.seed(13003)
  legacy_extra<-train(r,data=input,method="xgbTree",metric="ROC",tuneGrid=grid,
                       trControl=ctrl,nthread=1,verbose=0)
  saveRDS(legacy_extra,legacy_extra_path)
}
extra_val<-data.frame(row_id=prep$validation$row_id,diagnosis=prep$validation$species,
                       algorithm="historical_five_age_sex",
                       probability_TB=predict(legacy_extra,prep$validation,type="prob")$TB)
extra_test<-data.frame(row_id=prep$test$row_id,diagnosis=prep$test$species,
                        algorithm="historical_five_age_sex",
                        probability_TB=predict(legacy_extra,prep$test,type="prob")$TB)
write.csv(extra_val,"results/historical_five_age_sex_validation_predictions.csv",row.names=FALSE)
write.csv(extra_test,"results/historical_five_age_sex_test_predictions.csv",row.names=FALSE)
all_v<-rbind(val,extra_val);all_t<-rbind(test,extra_test)
cmp_names<-c(winner,"age_sex","panel_age_sex","historical_five_xgbTree","historical_five_age_sex")
demo<-do.call(rbind,lapply(cmp_names,function(nm){
  v<-all_v[all_v$algorithm==nm,];t<-all_t[all_t$algorithm==nm,]
  cut<-choose_threshold(v$diagnosis,v$probability_TB,"sensitivity",.9)
  s<-binary_stats(t$diagnosis,t$probability_TB,cut);a<-auc_ci(t$diagnosis,t$probability_TB)
  data.frame(model=nm,cutoff=cut,auc=unname(a[1]),auc_low=unname(a[2]),auc_high=unname(a[3]),
             as.list(s),sensitivity_low=wilson(s["TP"],s["TP"]+s["FN"])[1],
             sensitivity_high=wilson(s["TP"],s["TP"]+s["FN"])[2],
             specificity_low=wilson(s["TN"],s["TN"]+s["FP"])[1],
             specificity_high=wilson(s["TN"],s["TN"]+s["FP"])[2])
}))
write.csv(demo,"results/demographic_and_historical_comparison.csv",row.names=FALSE)
legacy_test<-all_t[all_t$algorithm=="historical_five_xgbTree",]
extra_test<-extra_test[match(legacy_test$row_id,extra_test$row_id),]
stopifnot(identical(as.character(legacy_test$diagnosis),as.character(extra_test$diagnosis)))
legacy_roc<-roc(as.character(legacy_test$diagnosis),
                 legacy_test$probability_TB,
                 levels=c("Pneumonia","TB"),direction="<",quiet=TRUE)
extra_roc<-roc(as.character(extra_test$diagnosis),extra_test$probability_TB,
                levels=c("Pneumonia","TB"),direction="<",quiet=TRUE)
comparison<-roc.test(legacy_roc,extra_roc,paired=TRUE,method="delong")
write.csv(data.frame(comparison="historical five versus historical five plus age/sex",
                      auc_original5=as.numeric(auc(legacy_roc)),
                      auc_original5_age_sex=as.numeric(auc(extra_roc)),
                      auc_difference=as.numeric(auc(legacy_roc)-auc(extra_roc)),
                      difference_low=unname(comparison$conf.int[1]),
                      difference_high=unname(comparison$conf.int[2]),
                      p_value=comparison$p.value),
          "results/historical_demographic_paired_test.csv",row.names=FALSE)

## Table S2 replacement includes all actual fitted hyperparameters.
parameters<-list()
for(nm in c("glm","glmnet","rpart","ranger","xgbTree","svmRadial","gbm","nnet","knn")){
  fit<-readRDS(paste0("models/final_",nm,".rds"))
  bt<-fit$bestTune
  parameters[[nm]]<-data.frame(algorithm=nm,parameter=names(bt),value=as.character(bt[1,]))
}
write.csv(do.call(rbind,parameters),"results/selected_hyperparameters.csv",row.names=FALSE)

## Calibration intervals by record-level bootstrap, with failed/singular
## replicates explicitly counted rather than silently accepted.
set.seed(16002)
yy<-truth_bin(main_test$diagnosis);pp<-pmin(pmax(main_test$probability_TB,1e-6),1-1e-6)
cal_boot<-replicate(1000,{
  ii<-sample.int(length(yy),length(yy),replace=TRUE)
  y<-yy[ii];l<-qlogis(pp[ii])
  f1<-suppressWarnings(glm(y~offset(l),family=binomial()))
  f2<-suppressWarnings(glm(y~l,family=binomial()))
  c(intercept=unname(coef(f1)[1]),slope=unname(coef(f2)[2]))
})
boot_summary<-data.frame(parameter=c("intercept","slope"),
                          lower=apply(cal_boot,1,quantile,.025,na.rm=TRUE),
                          upper=apply(cal_boot,1,quantile,.975,na.rm=TRUE),
                          valid_replicates=rowSums(is.finite(cal_boot)))
write.csv(boot_summary,"results/calibration_bootstrap_intervals.csv",row.names=FALSE)
cat("Supplemental tables complete.\n")
