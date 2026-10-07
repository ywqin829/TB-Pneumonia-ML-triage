## Scientific figures rebuilt from the locked outputs, matching source
## sans-serif fonts and red/blue palette. Each panel is exported separately
## as requested by the author; no combined multi-panel figure is generated.
suppressPackageStartupMessages({library(ggplot2);library(caret);
  library(recipes);library(pROC);library(xgboost)})
dir.create("figures",showWarnings=FALSE)
p<-readRDS("results/prepared_data.rds");lock<-readRDS("results/analysis_lock.rds")
sel<-readRDS("models/full_training_selection.rds")
model<-readRDS(paste0("models/final_",lock$winner,".rds"))
val<-read.csv("results/validation_predictions.csv");test<-read.csv("results/test_predictions.csv")
val<-val[val$algorithm==lock$winner,];test<-test[test$algorithm==lock$winner,]
thr<-read.csv("results/locked_thresholds.csv")
lo<-thr$threshold[thr$rule=="sensitivity_90"];hi<-thr$threshold[thr$rule=="rulein_specificity_90"]
red<-"#E64B35";blue<-"#4DBBD5";navy<-"#3C5488";gray<-"#8491B4"
theme_set(theme_classic(base_size=10,base_family="sans")+
  theme(plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=8),
        legend.position="bottom",legend.title=element_blank(),
        plot.tag=element_text(face="bold",size=12)))
panel_index<-list()
save_panel<-function(g,name,width=6,height=4.5){
  ## Native devices can fail with a non-ASCII Windows output path. Render
  ## into the ASCII runtime temporary directory, then copy verified files.
  tmp_pdf<-tempfile(fileext=".pdf");tmp_png<-tempfile(fileext=".png")
  ggsave(tmp_pdf,g,width=width,height=height,device=cairo_pdf,bg="white")
  ggsave(tmp_png,g,width=width,height=height,device=ragg::agg_png,dpi=300,bg="white")
  stopifnot(file.info(tmp_pdf)$size>0,file.info(tmp_png)$size>0)
  stopifnot(file.copy(tmp_pdf,paste0("figures/",name,".pdf"),overwrite=TRUE),
            file.copy(tmp_png,paste0("figures/",name,".png"),overwrite=TRUE))
  unlink(c(tmp_pdf,tmp_png))
  panel_index[[length(panel_index)+1L]]<<-data.frame(panel=name,width_in=width,height_in=height,dpi=300)
}
label_feature<-function(x){x[x=="MONO."]<-"MONO%";x[x=="LYM."]<-"LYM%";x[x=="NEU."]<-"NEU%";x}

## Figure 1: auditable record flow and resampling architecture.
nodes<-data.frame(x=c(5,5,2,5,8,2,2,5,8,5),y=c(11,9.6,8.1,8.1,8.1,6.4,4.6,4.6,4.6,2.5),
  text=c("Supplied extract: 2,399 data rows\nNo patient ID or date fields",
    "17 excess identical rows removed\n2,382 distinct rows: 711 TB, 1,671 pneumonia\nDiagnosis-stratified random 60/20/20; seed 123",
    "Training: 1,428\n426 TB / 1,002 pneumonia",
    "Validation: 476\n142 TB / 334 pneumonia",
    "Test: 478\n143 TB / 335 pneumonia",
    "Five outer training folds\nBoruta -> LASSO -> direction-aware AUC\n2-8-feature XGBoost combination search",
    "Nine-algorithm comparison\nFold-local preprocessing and SMOTE\nSelect winner, then fit on full training set",
    "Freeze validation thresholds\nSensitivity target 90%\nRule-in specificity target 90%",
    "Test evaluation\nMetrics and 95% CI / paired comparisons\nCalibration / DCA / subgroups",
    paste0("Model evidence: ",lock$winner,"\n",paste(label_feature(lock$panel),collapse=" + "),
           "\nSHAP and local triage dashboard")))
segments<-data.frame(x=c(5,5,5,5,2,2,5,8,3.1,6.1,8),y=c(10.45,9.05,9.05,9.05,7.55,5.85,7.55,7.55,4.6,4.6,4.05),
                      xend=c(5,2,5,8,2,2,5,8,3.9,6.9,5),yend=c(10.2,8.65,8.65,8.65,6.95,5.15,5.15,5.15,4.6,4.6,3.1))
flow<-ggplot()+geom_segment(data=segments,aes(x=x,y=y,xend=xend,yend=yend),
                           arrow=grid::arrow(length=grid::unit(2,"mm")),color="#555555")+
  geom_label(data=nodes,aes(x=x,y=y,label=text),size=2.5,label.size=.3,
             label.padding=grid::unit(.2,"lines"),fill="#F2F5F8")+
  annotate("text",x=5,y=1,label="File-level deduplication cannot establish unique patient identity.\nHistorical missing-record exclusions and the separately collected 30 TB cases cannot be traced in this extract.",
           size=2.6,color="#555555")+
  coord_cartesian(xlim=c(0,10),ylim=c(.5,12),clip="off")+theme_void()+
  labs(title="Reanalysis workflow and cohort accounting")+
  theme(plot.title=element_text(face="bold",hjust=.5,size=13),plot.margin=margin(20,20,20,20))
save_panel(flow,"Fig1_rerun_flow",8.2677,11.6929)

## Figure 2: demographics, inflammatory indices, correlation and Boruta.
d<-p$data
d$age_band<-cut(d$Age,breaks=seq(10,100,10),right=FALSE)
age<-as.data.frame(table(age_band=d$age_band,species=d$species))
age$percent<-ave(age$Freq,age$species,FUN=function(z)100*z/sum(z))
age$signed<-ifelse(age$species=="TB",-age$percent,age$percent)
g_age<-ggplot(age,aes(x=age_band,y=signed,fill=species))+geom_col(width=.85)+coord_flip()+
  scale_fill_manual(values=c(TB=red,Pneumonia=blue))+
  scale_y_continuous(labels=function(x)paste0(abs(x),"%"))+
  labs(title="Age distribution",x="Age group",y="Percentage within diagnosis")+
  theme(axis.text=element_text(size=7),legend.text=element_text(size=7))
sex<-as.data.frame(table(species=d$species,Gender=d$Gender))
sex$proportion<-ave(sex$Freq,sex$species,FUN=function(z)z/sum(z))
g_sex<-ggplot(sex,aes(species,proportion,fill=Gender))+geom_col(position="dodge",color="black",linewidth=.2)+
  coord_flip()+scale_fill_manual(values=c(Male=navy,Female="#EFA18B"))+
  scale_y_continuous(labels=scales::percent)+labs(title="Gender distribution",x=NULL,y="Percentage")+
  theme(legend.text=element_text(size=7))
indices<-c("NLR","PLR","LMR","SII","SIRI","PIV")
radar<-do.call(rbind,lapply(levels(d$species),function(s){
  data.frame(index=indices,species=s,relative_median=sapply(indices,function(nm)
    median(d[d$species==s,nm])/median(d[[nm]])))
}))
radar$index<-factor(radar$index,levels=indices)
radar$angle<-2*pi*(as.integer(radar$index)-1)/length(indices)
radar$x<-radar$relative_median*sin(radar$angle)
radar$y<-radar$relative_median*cos(radar$angle)
radar_max<-ceiling(max(radar$relative_median)*2)/2
radar_levels<-seq(.5,radar_max,.5)
radar_grid<-do.call(rbind,lapply(radar_levels,function(r){
  a<-2*pi*(0:length(indices))/length(indices)
  data.frame(x=r*sin(a),y=r*cos(a),radius=r)
}))
radar_axes<-data.frame(index=indices,angle=2*pi*(seq_along(indices)-1)/length(indices))
radar_axes$x<-radar_max*sin(radar_axes$angle);radar_axes$y<-radar_max*cos(radar_axes$angle)
g_radar<-ggplot()+
  geom_path(data=radar_grid,aes(x,y,group=radius),color="gray80",linewidth=.3)+
  geom_segment(data=radar_axes,aes(x=0,y=0,xend=x,yend=y),color="gray85",linewidth=.3)+
  geom_polygon(data=radar,aes(x,y,group=species,color=species,fill=species),alpha=.12,linewidth=.6)+
  geom_point(data=radar,aes(x,y,color=species),size=1.8)+
  geom_text(data=radar_axes,aes(x=x*1.18,y=y*1.18,label=index),size=3)+
  annotate("text",x=.08,y=radar_levels,label=sprintf("%.1f",radar_levels),size=2.5,hjust=0,color="gray40")+
  scale_color_manual(values=c(TB=red,Pneumonia=blue))+
  scale_fill_manual(values=c(TB=red,Pneumonia=blue))+
  coord_equal(xlim=c(-1,1)*radar_max*1.35,ylim=c(-1,1)*radar_max*1.3,clip="off")+
  labs(title="Inflammatory indices",subtitle="Median / pooled median",x=NULL,y=NULL)+
  theme_void(base_size=10,base_family="sans")+theme(legend.position="bottom",legend.title=element_blank(),
    plot.title=element_text(face="bold",size=10),plot.subtitle=element_text(size=8))
cor_dat<-do.call(rbind,lapply(levels(d$species),function(s){
  m<-cor(d[d$species==s,p$candidate_features],method="spearman")
  z<-as.data.frame(as.table(m));names(z)<-c("feature1","feature2","rho");z$species<-s;z
}))
g_cor<-ggplot(cor_dat,aes(feature1,feature2,fill=rho))+geom_tile()+facet_wrap(~species,nrow=1)+
  scale_fill_gradient2(low="#B2182B",mid="white",high="#2166AC",limits=c(-1,1))+
  scale_x_discrete(labels=label_feature)+scale_y_discrete(labels=label_feature)+
  labs(title="Spearman correlation of candidate predictors",x=NULL,y=NULL)+
  theme(axis.text.x=element_text(angle=90,hjust=1,size=7),axis.text.y=element_text(size=7),
        legend.position="bottom",legend.key.height=grid::unit(2,"mm"),strip.background=element_blank())
ih<-sel$boruta_fit$ImpHistory
boruta_long<-data.frame(feature=rep(colnames(ih),each=nrow(ih)),importance=as.vector(ih))
boruta_long<-boruta_long[is.finite(boruta_long$importance),]
boruta_long$decision<-ifelse(grepl("^shadow",boruta_long$feature),"Shadow reference",
                            as.character(sel$boruta_fit$finalDecision[boruta_long$feature]))
ord<-names(sort(tapply(boruta_long$importance,boruta_long$feature,median)))
boruta_long$feature<-factor(boruta_long$feature,levels=ord)
g_boruta<-ggplot(boruta_long,aes(feature,importance,fill=decision))+geom_boxplot(outlier.size=.3,linewidth=.2)+
  scale_fill_manual(values=c(Confirmed="#B9E68B",Rejected="#EFA18B",Tentative="#F0E442","Shadow reference"="#A9C7E4"))+
  labs(title="Boruta feature importance",subtitle="Full training set only",x=NULL,y="Importance")+
  scale_x_discrete(labels=label_feature)+theme(axis.text.x=element_text(angle=90,hjust=1,size=8))
save_panel(g_age,"Fig2A_age_distribution")
save_panel(g_sex,"Fig2B_gender_distribution")
save_panel(g_radar,"Fig2C_inflammatory_indices")
save_panel(g_cor,"Fig2D_spearman_correlation",9,5.5)
save_panel(g_boruta,"Fig2E_boruta_importance",9,4.5)

## Figure 3: selection and locked-model performance.
cv<-sel$lasso_fit
ld<-data.frame(loglambda=log(cv$lambda),deviance=cv$cvm,low=cv$cvlo,high=cv$cvup)
g_lasso<-ggplot(ld,aes(loglambda,deviance))+geom_errorbar(aes(ymin=low,ymax=high),color="gray65",width=0)+
  geom_point(color=red,size=.8)+geom_vline(xintercept=log(cv$lambda.min),linetype=2)+
  labs(title="LASSO 10-fold cross-validation",x="Log(lambda)",y="Binomial deviance")
au<-head(sel$auc_candidates,9)
ud<-data.frame(feature=factor(label_feature(names(au)),levels=rev(label_feature(names(au)))),auc=as.numeric(au))
g_auc<-ggplot(ud,aes(auc,feature))+geom_segment(aes(x=.5,xend=auc,yend=feature),color="gray60")+
  geom_point(color=red,size=2)+geom_text(aes(label=sprintf("%.3f",auc)),nudge_x=.013,size=2.4,hjust=0)+
  scale_x_continuous(limits=c(.49,max(ud$auc)+.065))+
  labs(title="Single-variable discrimination",subtitle="Direction-aware AUC; training set",x="max(AUC, 1-AUC)",y=NULL)
oa<-read.csv("results/outer_auc.csv")
order_alg<-names(sort(tapply(oa$auc,oa$algorithm,mean)))
oa$algorithm<-factor(oa$algorithm,levels=order_alg)
g_cv<-ggplot(oa,aes(algorithm,auc,fill=algorithm==lock$winner))+geom_boxplot(width=.6,outlier.size=.5)+
  geom_point(position=position_jitter(width=.06,seed=1),size=.6)+coord_flip()+
  scale_fill_manual(values=c("FALSE"="#72AAD3","TRUE"=red),guide="none")+
  labs(title="Five-fold outer CV",x=NULL,y="AUC of the full training pipeline")+
  theme(axis.text=element_text(size=7))
trprob<-predict(model,p$train,type="prob")$TB
roc_df<-function(truth,prob,dataset){r<-roc(truth,prob,levels=c("Pneumonia","TB"),direction="<",quiet=TRUE)
  data.frame(fpr=1-r$specificities,tpr=r$sensitivities,dataset=paste0(dataset," (AUC ",sprintf("%.3f",as.numeric(auc(r))),")"))}
rd<-rbind(roc_df(p$train$species,trprob,"Training apparent"),roc_df(val$diagnosis,val$probability_TB,"Validation"),
           roc_df(test$diagnosis,test$probability_TB,"Test"))
g_roc<-ggplot(rd,aes(fpr,tpr,color=dataset))+geom_line(linewidth=.6)+geom_abline(linetype=3,color="gray60")+
  scale_color_manual(values=c(red,"#879DA6",blue))+labs(title="ROC curves",x="1 - specificity",y="Sensitivity")+
  theme(legend.text=element_text(size=6),legend.position="bottom")
obs<-factor(test$diagnosis,levels=c("Pneumonia","TB"))
pred<-factor(ifelse(test$probability_TB>=lo,"TB","Pneumonia"),levels=c("Pneumonia","TB"))
cm<-as.data.frame(table(actual=obs,predicted=pred))
g_cm<-ggplot(cm,aes(actual,predicted,fill=Freq))+geom_tile(color="white")+
  geom_text(aes(label=Freq),fontface="bold",size=4)+scale_fill_gradient(low="#D7EEF5",high="#43AFC2",guide="none")+
  labs(title="Test confusion matrix",subtitle=sprintf("Validation sensitivity cutoff = %.3f",lo),x="Actual diagnosis",y="Predicted diagnosis")
save_panel(g_lasso,"Fig3A_lasso_CV")
save_panel(g_auc,"Fig3B_direction_aware_AUC")
save_panel(g_cv,"Fig3C_outer_CV")
save_panel(g_roc,"Fig3D_ROC")
save_panel(g_cm,"Fig3E_test_confusion_matrix")

## Figure 4: test calibration and net benefit.
cal<-read.csv("results/test_calibration.csv")
cd<-data.frame(probability=test$probability_TB,observed=as.integer(test$diagnosis=="TB"))
cd$bin<-cut(cd$probability,unique(quantile(cd$probability,seq(0,1,.1))),include.lowest=TRUE)
bins<-aggregate(cbind(probability,observed)~bin,cd,mean)
g_cal<-ggplot(cd,aes(probability,observed))+geom_abline(color="gray70",linewidth=1)+
  geom_smooth(method="loess",formula=y~x,span=1,se=FALSE,color=navy,linewidth=.7)+
  geom_point(data=bins,size=1.8,color=red)+geom_rug(sides="b",alpha=.25)+
  coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
  annotate("text",x=.04,y=.95,hjust=0,vjust=1,size=2.7,
           label=sprintf("Brier = %.3f\nIntercept = %.3f\nSlope = %.3f",cal$brier,cal$intercept,cal$slope))+
  labs(title="Calibration on the test set",x="Predicted TB probability",y="Observed TB fraction")
dca<-read.csv("results/test_decision_curve.csv")
dc<-rbind(data.frame(threshold=dca$threshold,net_benefit=dca$model,strategy="XGBoost"),
           data.frame(threshold=dca$threshold,net_benefit=dca$treat_all,strategy="Treat all"),
           data.frame(threshold=dca$threshold,net_benefit=dca$treat_none,strategy="Treat none"))
g_dca<-ggplot(dc,aes(threshold,net_benefit,color=strategy))+geom_line(linewidth=.7)+
  scale_color_manual(values=c("Treat all"="#879DA6","Treat none"="black",XGBoost=red))+
  coord_cartesian(ylim=c(-.05,.31))+scale_x_continuous(labels=scales::percent)+
  labs(title="Decision curve on the test set",x="Decision threshold",y="Net benefit")
save_panel(g_cal,"Fig4A_test_calibration")
save_panel(g_dca,"Fig4B_test_DCA")

## Figure 5: exact TreeSHAP contributions for the locked XGBoost model.
baked<-recipes::bake(model$recipe,new_data=p$test)
xx<-as.matrix(baked[,setdiff(names(baked),"species"),drop=FALSE])
contrib<-predict(model$finalModel,xx,predcontrib=TRUE)
direct<-predict(model$finalModel,xx)
stopifnot(max(abs(direct-predict(model,p$test,type="prob")$TB))<1e-6,
          max(abs(plogis(rowSums(contrib))-direct))<1e-5)
saveRDS(list(contribution=contrib,row_id=p$test$row_id,probability=direct),"results/test_shap.rds")
features<-setdiff(colnames(contrib),"BIAS")
importance<-sort(colMeans(abs(contrib[,features,drop=FALSE])))
sh<-do.call(rbind,lapply(features,function(nm){
  z<-p$test[[nm]]
  data.frame(feature=nm,shap=contrib[,nm],scaled_value=(z-min(z))/(max(z)-min(z)))
}))
sh$feature<-factor(sh$feature,levels=names(importance))
g_shap<-ggplot(sh,aes(shap,feature,color=scaled_value))+
  geom_point(position=position_jitter(height=.19,width=0,seed=5),alpha=.7,size=1.1)+
  scale_color_gradient(low="#231B43",high="#E83B6B",name="Feature value",breaks=c(0,1),labels=c("Low","High"))+
  scale_y_discrete(labels=label_feature)+labs(title="Global TreeSHAP summary",subtitle="Test records; contributions on log-odds scale",x="SHAP contribution",y=NULL)
g_shap<-g_shap+theme(legend.title=element_text(size=9))
local_plot<-function(index,title){
  z<-contrib[index,features];ord<-order(abs(z),decreasing=FALSE)
  q<-data.frame(feature=factor(label_feature(names(z)[ord]),levels=label_feature(names(z)[ord])),
                contribution=unname(z[ord]))
  ggplot(q,aes(contribution,feature,fill=contribution>0))+geom_col(width=.65)+
    scale_fill_manual(values=c("FALSE"=navy,"TRUE"=red),guide="none")+
    geom_vline(xintercept=0,color="gray50")+
    labs(title=title,subtitle=sprintf("TB probability %.1f%%; illustrative test record",100*direct[index]),
         x="SHAP contribution to log-odds",y=NULL)+theme(axis.text=element_text(size=7))
}
g_high<-local_plot(which.max(direct),"Higher-probability example")
g_low<-local_plot(which.min(direct),"Lower-probability example")
zt<-read.csv("results/three_zone_composition.csv");zt<-zt[zt$dataset=="test",]
zt$zone<-factor(zt$zone,levels=c("blue","gray","red"))
g_zone<-ggplot(zt,aes(zone,n,fill=diagnosis))+geom_col(color="white")+
  geom_text(aes(label=n),position=position_stack(vjust=.5),size=3)+
  scale_fill_manual(values=c(Pneumonia=blue,TB=red))+
  labs(title="Test-set triage zones",subtitle=sprintf("Blue: p < %.3f; gray: %.3f <= p < %.3f; red: p >= %.3f",lo,lo,hi,hi),
       x="Zone",y="Number of records")
save_panel(g_shap,"Fig5A_global_SHAP")
save_panel(g_high,"Fig5B_local_SHAP_high")
save_panel(g_low,"Fig5C_local_SHAP_low")
save_panel(g_zone,"Fig5D_test_triage_zones",7,4.5)

## Supplement: subgroup uncertainty and feature-selection stability.
sg<-read.csv("results/test_subgroup_metrics.csv")
sg$group<-unname(c(all="All test records",under_65="Age < 65 years",age_65_or_older="Age >= 65 years",
                  male="Male",female="Female")[sg$group])
g_sub<-ggplot(sg,aes(auc,reorder(group,auc)))+geom_errorbarh(aes(xmin=auc_low,xmax=auc_high),height=.15)+
  geom_point(color=navy,size=2)+labs(title="Test subgroup AUC with 95% CI",x="AUC",y=NULL)
stab<-read.csv("results/feature_selection_stability.csv")
stab<-stab[stab$final_panel>0,]
stab$feature<-label_feature(stab$feature)
g_stab<-ggplot(stab,aes(reorder(feature,final_panel),final_panel))+geom_col(fill=blue)+coord_flip()+
  scale_y_continuous(breaks=0:5,limits=c(0,5))+labs(title="Feature-selection stability",x=NULL,y="Outer folds selecting the feature")
demo<-read.csv("results/demographic_and_historical_comparison.csv")
demo$model<-unname(c(xgbTree="Reselected 3-feature XGBoost",age_sex="Age + sex",
                    panel_age_sex="Reselected panel + age/sex",historical_five_xgbTree="Original 5-feature XGBoost",
                    historical_five_age_sex="Original 5 features + age/sex")[demo$model])
g_demo<-ggplot(demo,aes(auc,reorder(model,auc)))+
  geom_errorbarh(aes(xmin=auc_low,xmax=auc_high),height=.15)+geom_point(color=red,size=2)+
  labs(title="Prespecified model comparisons on test set",x="AUC with 95% CI",y=NULL)+theme(axis.text.y=element_text(size=7))
save_panel(g_sub,"FigS1A_subgroup_AUC")
save_panel(g_stab,"FigS1B_selection_stability")
save_panel(g_demo,"FigS1C_demographic_comparison",7,4.5)
write.csv(do.call(rbind,panel_index),"figures/panel_index.csv",row.names=FALSE)
cat(length(panel_index),"independent panels exported as PDF and 300-dpi PNG.\n")
