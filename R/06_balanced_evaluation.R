## Author-selected balanced operating strategy (2026-10-04).
## No predictive model is retrained, and earlier 90%-sensitivity results remain.
## Numeric Youden cutoff was computed from validation before this change.
## The author selected the balanced objective after seeing test performance;
## this is a disclosed post-result strategy adjustment, not preregistration.
suppressPackageStartupMessages({library(ggplot2);library(pROC)})
stopifnot(file.exists("results/analysis_lock.rds"))
out <- "balanced_scheme"
dir.create(file.path(out,"results"),showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(out,"figures"),showWarnings=FALSE)
lock <- readRDS("results/analysis_lock.rds")
prep <- readRDS("results/prepared_data.rds")
validation <- read.csv("results/validation_predictions.csv",stringsAsFactors=FALSE)
test <- read.csv("results/test_predictions.csv",stringsAsFactors=FALSE)
main_val <- validation[validation$algorithm==lock$winner,]
main_test <- test[test$algorithm==lock$winner,]
thresholds <- read.csv("results/locked_thresholds.csv",stringsAsFactors=FALSE)
cutoff <- thresholds$threshold[thresholds$rule=="youden"]

confusion <- function(truth,prob,t) {
  positive <- truth=="TB"; predicted <- prob>=t
  c(TP=sum(positive&predicted),FP=sum(!positive&predicted),
    TN=sum(!positive&!predicted),FN=sum(positive&!predicted))
}
stats_at <- function(truth,prob,t) {
  q <- confusion(truth,prob,t)
  sens<-unname(q["TP"]/(q["TP"]+q["FN"]))
  spec<-unname(q["TN"]/(q["TN"]+q["FP"]))
  c(q,n=sum(q),sensitivity=sens,specificity=spec,
    ppv=unname(q["TP"]/(q["TP"]+q["FP"])),
    npv=unname(q["TN"]/(q["TN"]+q["FN"])),
    accuracy=unname((q["TP"]+q["TN"])/sum(q)),
    f1=unname(2*q["TP"]/(2*q["TP"]+q["FP"]+q["FN"])),
    balanced_accuracy=(sens+spec)/2,youden=sens+spec-1)
}
wilson <- function(x,n) {
  z<-qnorm(.975);p<-x/n;den<-1+z*z/n
  center<-(p+z*z/(2*n))/den
  half<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den
  c(low=max(0,center-half),high=min(1,center+half))
}
grid_for_prob <- function(truth,prob) {
  t<-sort(unique(c(0,prob,1)))
  s<-t(vapply(t,function(th)stats_at(truth,prob,th),numeric(13)))
  data.frame(threshold=t,s,check.names=FALSE)
}
val_grid<-grid_for_prob(main_val$diagnosis,main_val$probability_TB)
computed_cutoff<-val_grid$threshold[which.max(val_grid$youden)]
stopifnot(length(cutoff)==1L,abs(cutoff-computed_cutoff)<1e-12)
write.csv(val_grid,file.path(out,"results/validation_threshold_grid.csv"),row.names=FALSE)
write.csv(data.frame(strategy="balanced_youden",threshold=cutoff,
  numerical_threshold_source="validation",
  objective="maximize sensitivity + specificity - 1; equivalently balanced accuracy",
  strategy_timing="author-selected after viewing test results; post-result adjustment"),
  file.path(out,"results/balanced_threshold.csv"),row.names=FALSE)

## Preserve already audited intervals for existing operating points.
previous<-read.csv("results/threshold_performance.csv",stringsAsFactors=FALSE)
comparison<-previous[previous$rule%in%c("sensitivity_90","sensitivity_85_reference","youden"),]
comparison$balanced_accuracy<-(comparison$sensitivity+comparison$specificity)/2
comparison$youden<-comparison$sensitivity+comparison$specificity-1
write.csv(comparison,file.path(out,"results/operating_strategy_comparison.csv"),row.names=FALSE)
balanced<-comparison[comparison$rule=="youden",]
## Stratified record bootstrap: fixed observed TB/pneumonia counts, fixed model
## and threshold. Does not account for threshold selection or unknown clustering.
ba_interval<-function(truth,prob,t,seed){
  set.seed(seed)
  a<-which(truth=="TB");b<-which(truth=="Pneumonia")
  draws<-replicate(2000,{
    ii<-c(sample(a,length(a),replace=TRUE),sample(b,length(b),replace=TRUE))
    unname(stats_at(truth[ii],prob[ii],t)["balanced_accuracy"])
  })
  unname(quantile(draws,c(.025,.975)))
}
balanced$balanced_accuracy_low<-NA_real_;balanced$balanced_accuracy_high<-NA_real_
balanced$auc<-NA_real_;balanced$auc_low<-NA_real_;balanced$auc_high<-NA_real_
for(i in seq_len(nrow(balanced))){
  d<-if(balanced$dataset[i]=="validation")main_val else main_test
  z<-stats_at(d$diagnosis,d$probability_TB,cutoff)
  stopifnot(all(abs(as.numeric(balanced[i,c("TP","FP","TN","FN")])-z[c("TP","FP","TN","FN")])<1e-10))
  ci<-ba_interval(d$diagnosis,d$probability_TB,cutoff,19000L+i)
  balanced$balanced_accuracy_low[i]<-ci[1];balanced$balanced_accuracy_high[i]<-ci[2]
  r<-roc(d$diagnosis,d$probability_TB,levels=c("Pneumonia","TB"),direction="<",quiet=TRUE)
  ac<-ci.auc(r,method="delong")
  balanced$auc[i]<-as.numeric(auc(r));balanced$auc_low[i]<-ac[1];balanced$auc_high[i]<-ac[3]
}
write.csv(balanced,file.path(out,"results/balanced_performance.csv"),row.names=FALSE)

## Compare each already fitted model at its own validation-derived Youden
## cutoff, consistently with the new objective. This does not select a new
## algorithm or panel on test results. AUC comparisons remain unchanged.
all_validation<-validation;all_test<-test
if(file.exists("results/historical_five_age_sex_validation_predictions.csv")){
  all_validation<-rbind(all_validation,read.csv("results/historical_five_age_sex_validation_predictions.csv",stringsAsFactors=FALSE))
  all_test<-rbind(all_test,read.csv("results/historical_five_age_sex_test_predictions.csv",stringsAsFactors=FALSE))
}
val_models<-do.call(rbind,lapply(unique(all_validation$algorithm),function(nm){
  d<-all_validation[all_validation$algorithm==nm,]
  g<-grid_for_prob(d$diagnosis,d$probability_TB)
  q<-g[which.max(g$youden),,drop=FALSE]
  data.frame(algorithm=nm,q,check.names=FALSE)
}))
write.csv(val_models,file.path(out,"results/validation_model_operating_points.csv"),row.names=FALSE)
test_models<-do.call(rbind,lapply(seq_len(nrow(val_models)),function(i){
  nm<-val_models$algorithm[i];t<-val_models$threshold[i]
  d<-all_test[all_test$algorithm==nm,]
  stopifnot(nrow(d)==nrow(main_test),setequal(d$row_id,main_test$row_id))
  z<-stats_at(d$diagnosis,d$probability_TB,t)
  q<-data.frame(algorithm=nm,threshold=t,as.list(z),check.names=FALSE)
  w<-list(sensitivity=wilson(z["TP"],z["TP"]+z["FN"]),
    specificity=wilson(z["TN"],z["TN"]+z["FP"]),
    ppv=wilson(z["TP"],z["TP"]+z["FP"]),npv=wilson(z["TN"],z["TN"]+z["FN"]),
    accuracy=wilson(z["TP"]+z["TN"],z["n"]))
  for(v in names(w)){q[[paste0(v,"_low")]]<-unname(w[[v]][1]);q[[paste0(v,"_high")]]<-unname(w[[v]][2])}
  ci<-ba_interval(d$diagnosis,d$probability_TB,t,19200L+i)
  q$balanced_accuracy_low<-ci[1];q$balanced_accuracy_high<-ci[2]
  r<-roc(d$diagnosis,d$probability_TB,levels=c("Pneumonia","TB"),direction="<",quiet=TRUE)
  ac<-ci.auc(r,method="delong")
  q$auc<-as.numeric(auc(r));q$auc_low<-ac[1];q$auc_high<-ac[3]
  ## Preserve the primary model's previously audited F1 interval exactly;
  ## other models receive deterministic record bootstrap intervals.
  if(nm==lock$winner){
    main<-balanced[balanced$dataset=="test",]
    q$f1_low<-main$f1_low;q$f1_high<-main$f1_high
    q$balanced_accuracy_low<-main$balanced_accuracy_low
    q$balanced_accuracy_high<-main$balanced_accuracy_high
  } else {
    set.seed(19300L+i)
    f<-replicate(2000,{
      ii<-sample.int(nrow(d),nrow(d),replace=TRUE)
      unname(stats_at(d$diagnosis[ii],d$probability_TB[ii],t)["f1"])
    })
    q$f1_low<-unname(quantile(f,.025));q$f1_high<-unname(quantile(f,.975))
  }
  q
}))
write.csv(test_models,file.path(out,"results/balanced_model_comparisons.csv"),row.names=FALSE)

## Subgroups all use the same balanced cutoff and the same trained model.
matched_test<-prep$test[match(main_test$row_id,prep$test$row_id),]
stopifnot(identical(as.character(matched_test$species),main_test$diagnosis))
groups<-list(all=rep(TRUE,nrow(main_test)),under_65=matched_test$Age<65,
  age_65_or_older=matched_test$Age>=65,male=matched_test$Gender=="Male",female=matched_test$Gender=="Female")
subgroups<-do.call(rbind,lapply(names(groups),function(name){
  ii<-groups[[name]];d<-main_test[ii,];z<-stats_at(d$diagnosis,d$probability_TB,cutoff)
  q<-as.data.frame(as.list(z));q$group<-name;q$threshold<-cutoff
  q$tb_n<-sum(d$diagnosis=="TB");q$pneumonia_n<-sum(d$diagnosis=="Pneumonia")
  w<-list(sensitivity=wilson(z["TP"],q$tb_n),specificity=wilson(z["TN"],q$pneumonia_n),
          ppv=wilson(z["TP"],z["TP"]+z["FP"]),npv=wilson(z["TN"],z["TN"]+z["FN"]),
          accuracy=wilson(z["TP"]+z["TN"],z["n"]))
  for(nm in names(w)){q[[paste0(nm,"_low")]]<-unname(w[[nm]][1]);q[[paste0(nm,"_high")]]<-unname(w[[nm]][2])}
  ci<-ba_interval(d$diagnosis,d$probability_TB,cutoff,19100L+match(name,names(groups)))
  q$balanced_accuracy_low<-ci[1];q$balanced_accuracy_high<-ci[2]
  q
}))
write.csv(subgroups,file.path(out,"results/balanced_test_subgroups.csv"),row.names=FALSE)
predictions<-main_test
predictions$balanced_threshold<-cutoff
predictions$balanced_prediction<-ifelse(predictions$probability_TB>=cutoff,"TB","Pneumonia")
write.csv(predictions,file.path(out,"results/balanced_test_predictions.csv"),row.names=FALSE)

red<-"#E64B35";blue<-"#4DBBD5";navy<-"#3C5488"
theme_set(theme_classic(base_size=10,base_family="sans")+
  theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=8),legend.position="bottom",legend.title=element_blank()))
save_panel<-function(g,name,width=6,height=4.5){
  a<-tempfile(fileext=".pdf");b<-tempfile(fileext=".png")
  ggsave(a,g,width=width,height=height,device=cairo_pdf,bg="white")
  ggsave(b,g,width=width,height=height,device=ragg::agg_png,dpi=300,bg="white")
  stopifnot(file.copy(a,file.path(out,"figures",paste0(name,".pdf")),overwrite=TRUE),
            file.copy(b,file.path(out,"figures",paste0(name,".png")),overwrite=TRUE))
  unlink(c(a,b))
}
cm<-as.data.frame(table(actual=factor(main_test$diagnosis,levels=c("Pneumonia","TB")),
  predicted=factor(predictions$balanced_prediction,levels=c("Pneumonia","TB"))))
g_cm<-ggplot(cm,aes(actual,predicted,fill=Freq))+geom_tile(color="white")+
  geom_text(aes(label=Freq),fontface="bold",size=5)+scale_fill_gradient(low="#D7EEF5",high="#43AFC2",guide="none")+
  labs(title="Balanced operating point: test confusion matrix",
    subtitle=sprintf("Validation Youden cutoff = %.6f",cutoff),x="Actual diagnosis",y="Predicted diagnosis")
save_panel(g_cm,"Balanced_test_confusion_matrix")
curve<-rbind(data.frame(threshold=val_grid$threshold,value=val_grid$sensitivity,metric="Sensitivity"),
  data.frame(threshold=val_grid$threshold,value=val_grid$specificity,metric="Specificity"),
  data.frame(threshold=val_grid$threshold,value=val_grid$balanced_accuracy,metric="Balanced accuracy"))
g_tradeoff<-ggplot(curve,aes(threshold,value,color=metric))+geom_step(linewidth=.65)+
  geom_vline(xintercept=cutoff,linetype=2,color="gray40")+
  scale_color_manual(values=c("Sensitivity"=red,"Specificity"=blue,"Balanced accuracy"=navy))+
  scale_y_continuous(labels=scales::percent,limits=c(0,1))+
  labs(title="Operating-point tradeoff on validation data",
    subtitle=sprintf("Youden / balanced accuracy optimum = %.6f",cutoff),x="Score cutoff",y="Validation metric")
save_panel(g_tradeoff,"Balanced_validation_threshold_tradeoff")
display_names<-c(all="All test records",under_65="Age < 65 years",age_65_or_older="Age >= 65 years",male="Male",female="Female")
forest<-rbind(data.frame(group=display_names[subgroups$group],metric="Sensitivity",estimate=subgroups$sensitivity,
  lower=subgroups$sensitivity_low,upper=subgroups$sensitivity_high),
  data.frame(group=display_names[subgroups$group],metric="Specificity",estimate=subgroups$specificity,
    lower=subgroups$specificity_low,upper=subgroups$specificity_high))
forest$group<-factor(forest$group,levels=rev(unname(display_names)))
g_sub<-ggplot(forest,aes(estimate,group,color=metric))+
  geom_errorbarh(aes(xmin=lower,xmax=upper),height=.15,position=position_dodge(width=.4))+
  geom_point(position=position_dodge(width=.4),size=2)+
  scale_color_manual(values=c("Sensitivity"=red,"Specificity"=blue))+
  scale_x_continuous(labels=scales::percent,limits=c(0,1))+
  labs(title="Balanced cutoff: test subgroup performance",subtitle="Same model and cutoff for every subgroup; Wilson 95% CI",x="Sensitivity / specificity",y=NULL)
save_panel(g_sub,"Balanced_test_subgroup_performance",7,4.5)

main<-balanced[balanced$dataset=="test",]
writeLines(c("Author-selected balanced operating strategy: 2026-10-04",
  "The numerical threshold is validation-derived; the strategy choice is a post-result adjustment.",
  sprintf("Model unchanged: %s; %s",lock$winner,paste(lock$panel,collapse=", ")),
  sprintf("Balanced cutoff %.12f",cutoff),
  sprintf("Test: TP %d, FP %d, TN %d, FN %d",main$TP,main$FP,main$TN,main$FN),
  sprintf("Sensitivity %.4f; specificity %.4f; balanced accuracy %.4f",main$sensitivity,main$specificity,main$balanced_accuracy),
  "Earlier sensitivity-oriented results and three-zone thresholds remain available as supplementary analyses.",
  "AUC, scores, calibration, DCA and SHAP do not change when only the binary operating cutoff changes.",
  "No predictive model refitting, source-data editing or manuscript editing was performed."),
  file.path(out,"results/balanced_summary.txt"),useBytes=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"results/session_info.txt"))
print(balanced[,c("dataset","threshold","TP","FP","TN","FN","sensitivity","specificity","balanced_accuracy")])
cat("Balanced evaluation and three independent figure panels complete.\n")
