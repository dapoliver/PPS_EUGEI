library(tidyverse)
library(plyr)
library(readxl)
library(gtsummary)
library(gghalves)
library(ggthemes)
library(mice)
library(VIM)
library(glmnet)
library(survminer)
library(survival)
library(survMisc)
library(timeROC)
library(prevalence)
library(yardstick)
library(nestedcv)
library(probably)
library(factoextra)
library(caret)
library(scales)

##### Read in already scored data #####
### Add in polygenic risk score data
PPS_scored <- read.csv("~/Dropbox/Work/PPS/EU-GEI/Databases/PPS_scored_EUGEI_150221.csv")
PRS <- read_excel("~/Dropbox/Work/PPS/EU-GEI/Databases/PRS/EUGEI_WP5_merged.xlsx")
PRS2 <- PRS[c(2,36:45,68,99:101)] # subset to only include ID, PCs1:10, PRS_SCZ in all subjects, PRS_SCZ, BD, MDD in white Europeans only
PRS2 <- PRS2 %>% dplyr::rename(ID=st_subjid) # Rename ID to match PPS_scored
PPS_scored <- join(PPS_scored, PRS2, by="ID") # Add PRS columns to PPS_scored
clinical <- read_excel("~/Dropbox/Work/PPS/EU-GEI/Databases/eugei_data_all.xlsx")
clinical <- clinical %>% subset(select=c(st_subjid, caarms_unu4, caarms_unu5, caarms_non4, caarms_non5, caarms_per4, caarms_per5, caarms_dis4,caarms_dis5, gafex01,gafex02))
clinical[,c(2:11)] <- lapply(clinical[,c(2:11)], FUN=as.numeric)
clinical <- as.data.frame(clinical)
clinical <- clinical %>% mutate(utc = case_when(is.na(caarms_unu5) ~ caarms_unu4*caarms_unu5,
                                                TRUE ~ caarms_unu4),
                                nbi = case_when(is.na(caarms_non5) ~ caarms_non4*caarms_non5,
                                                TRUE ~ caarms_non4),
                                pa = case_when(is.na(caarms_per5) ~ caarms_per4*caarms_per5,
                                                TRUE ~ caarms_per4),
                                ds = case_when(is.na(caarms_dis5) ~ caarms_dis4*caarms_dis5,
                                                TRUE ~ caarms_dis4))
                                
clinical$CAARMS <- rowSums(clinical[,c("utc","nbi", "pa", "ds")], na.rm=TRUE)

clinical <- clinical %>% subset(select=c(st_subjid, CAARMS, gafex01, gafex02))
clinical <- clinical %>% dplyr::rename(ID=st_subjid) # Rename ID to match PPS_scored
PPS_scored <- join(PPS_scored, clinical, by="ID") # Add PRS columns to PPS_scored

### Impute missing data ###
set.seed(123)
mice_imputes = mice(PPS_scored)
PPS_scored <- complete(mice_imputes)


PPS_inc <- read_excel("~/Dropbox/Work/PPS/EU-GEI/Databases/EU-GEI_Dominic_Oliver_single_record_3.8.2020.xlsx")
PPS_inc <- PPS_inc %>% subset(select=c(st_subjid, Subject_status))
PPS_scored <- merge(PPS_scored, PPS_inc, by.x="ID", by.y="st_subjid")
PPS_scored <- PPS_scored %>% mutate(chr = case_when(Subject_status==1 ~ 1,
                                                    Subject_status==3 ~ 0))
PPS_scored_all <- PPS_scored
PPS_scored_hc <- PPS_scored %>% filter(Subject_status==3)
PPS_scored <- PPS_scored %>% filter(Subject_status==1)

##### PLOT SURVIVAL CURVE #####
PPS_scored_uncensored <- PPS_scored
PPS_scored <- PPS_scored %>% mutate(Transition = case_when(day_exit>730 ~ 0,
                                                           TRUE ~ Transition),
                                    day_exit = case_when(day_exit>730 ~ 730,
                                                         TRUE ~ day_exit))
surv <- with(PPS_scored, Surv(day_exit, event = Transition))

#Default setting uses Greenwood CI - default is log transformation
km.surv <- survfit(surv ~ 1, data=PPS_scored)
surv_summ <- data.frame(n=summary(km.surv)$n,
                        time=summary(km.surv)$time,
                        years=summary(km.surv)$time/365,
                        n.risk=summary(km.surv)$n.risk,
                        n.event=summary(km.surv)$n.event,
                        surv=paste0(round(1-(summary(km.surv)$surv),4)*100," (95%CI: ",round(1-(summary(km.surv)$upper),4)*100,
                                    "-",round(1-(summary(km.surv)$lower),4)*100, ", ",summary(km.surv)$n.risk,
                                    " individuals still at risk)"))

png("~/Dropbox/Work/PPS/EU-GEI/Plots/PPS_Surv_140624.png",width = 900, height = 800)
ggsurvplot(km.surv, data=PPS_scored,
           fun="event",
           xscale = 365.25,
           break.time.by = 182.625,
           xlim = c(0,1095.75),
           cumevents = TRUE)
dev.off()

##### Assess proportion of variance explained #####

PPS.lr <- lm(chr ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity + 
                  Immigration  + Paternal_SES + Parental_SMI + 
                  Adult_Life_Events + Tobacco + Cannabis + Childhood_Trauma + 
                  Anhedonia, data = PPS_scored_all)
PPS.lr 
summary(PPS.lr)

PRS.lr <- lm(chr ~ ZRE_SCZ_imp_0.1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PPS_scored_all)
PRS.lr 
summary(PRS.lr)

PPS_PRS.lr <- lm(chr ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity + 
                Immigration  + Paternal_SES + Parental_SMI + 
                Adult_Life_Events + Tobacco + Cannabis + Childhood_Trauma + 
                Anhedonia + ZRE_SCZ_imp_0.1 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = PPS_scored_all)
PPS_PRS.lr 
summary(PPS_PRS.lr)

# Plot
variance_hc <- data.frame(factors=factor(c("PRS","PPS", "PRS + PPS"), levels=c("PRS","PPS", "PRS + PPS")),
                          variance=c(7.2,30.5,34.7))
ggplot(data=variance_hc, aes(x=factors, y=variance, fill=factors)) +
  geom_col() +
  geom_text(aes(label=paste(variance,"%", sep=""),y=variance+5), size=5) +
  ylab("Variance Explained (%)") +
  xlab("") +
  scale_y_continuous(limits=c(0,100)) +
  theme_classic() + 
  scale_fill_discrete_sequential() +
  theme(legend.position = "none")
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Variance_hc_300824.png",width = 24, height = 18, units = "cm")

clinical.cox <- coxph(Surv(day_exit, Transition) ~ CAARMS + gafex01 + gafex02, data = PPS_scored)
clinical.cox
summary(clinical.cox)
survMisc::rsq(clinical.cox)

PPS.cox <- coxph(Surv(day_exit, Transition) ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity + 
                   Immigration  + Paternal_Age + Paternal_SES + Parental_SMI + 
                   Adult_Life_Events + Tobacco + Cannabis + Childhood_Trauma + 
                   Anhedonia, data = PPS_scored)
PPS.cox
summary(PPS.cox)
survMisc::rsq(PPS.cox)

PRS.cox <- coxph(Surv(day_exit, Transition) ~ ZRE_SCZ_imp_0.1, data = PPS_scored)
PRS.cox
summary(PRS.cox)
survMisc::rsq(PRS.cox)

clinical_PRS.cox <- coxph(Surv(day_exit, Transition) ~ CAARMS + gafex01 + gafex02 + ZRE_SCZ_imp_0.1, data = PPS_scored)
clinical_PRS.cox
summary(clinical_PRS.cox)
survMisc::rsq(clinical_PRS.cox)

clinical_PPS.cox <- coxph(Surv(day_exit, Transition) ~ CAARMS + gafex01 + gafex02 + Gender + Handedness + Urbanicity + Pollution + Ethnicity + 
                            Immigration  + Paternal_Age + Paternal_SES + Parental_SMI + 
                            Adult_Life_Events + Tobacco + Cannabis + Childhood_Trauma + 
                            Anhedonia, data = PPS_scored)
clinical_PPS.cox
summary(clinical_PPS.cox)
survMisc::rsq(clinical_PPS.cox)

PPS_PRS.cox <- coxph(Surv(day_exit, Transition) ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity + 
                       Immigration  + Paternal_Age + Paternal_SES + Parental_SMI + 
                       Adult_Life_Events + Tobacco + Cannabis + Childhood_Trauma + 
                       Anhedonia + ZRE_SCZ_imp_0.1 + (PPS*ZRE_SCZ_imp_0.1), data = PPS_scored)
PPS_PRS.cox
summary(PPS_PRS.cox)
survMisc::rsq(PPS_PRS.cox)

clinical_PPS_PRS.cox <- coxph(Surv(day_exit, Transition) ~ CAARMS + gafex01 + gafex02 + Gender + Handedness + Urbanicity + Pollution + Ethnicity + 
                                Immigration  + Paternal_Age + Paternal_SES + Parental_SMI + 
                                Adult_Life_Events + Tobacco + Cannabis + Childhood_Trauma + 
                                Anhedonia + ZRE_SCZ_imp_0.1 + (PPS*ZRE_SCZ_imp_0.1), data = PPS_scored)
clinical_PPS_PRS.cox
summary(clinical_PPS_PRS.cox)
survMisc::rsq(clinical_PPS_PRS.cox)

variance <- data.frame(factors=factor(c("Clinical", "PRS","PPS", "Clinical + PRS", "Clinical + PPS", "PRS + PPS", "Clinical + PRS + PPS"), levels=c("Clinical", "PRS","PPS", "Clinical + PRS", "Clinical + PPS", "PRS + PPS", "Clinical + PRS + PPS")),
                       variance=c(2.6, 5.1, 6.3, 8.1, 9.1, 13, 17))
ggplot(data=variance, aes(x=factors, y=variance, fill=factors)) +
  geom_col() +
  geom_text(aes(label=paste(variance,"%", sep=""),y=variance+5), size=5) +
  ylab("Variance Explained (%)") +
  scale_y_continuous(limits=c(0,100)) +
  xlab("") +
  theme_classic() + 
  scale_fill_discrete_sequential() +
  theme(legend.position = "none")
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Variance_Tx_300824.png",width = 24, height = 18, units = "cm")

##### Detection ####

PPS_scored_all[,c(3:8,10:14,16:17)] <- lapply(PPS_scored_all[,c(3:8,10:14,16:17)], factor)
PPS_scored_mat <- model.matrix(~.-1,PPS_scored_all[,c(3:8,10:14,16:17,23:39)])

PPS_scored_mat <- as.data.frame(PPS_scored_mat) %>% subset(select=c(-ZRE_SCZ_imp_EUR_0.1))
PPS_scored_mat$chr <- PPS_scored_all$chr

table_HC <- tbl_summary(PPS_scored_all,include=c(Gender, Handedness, Urbanicity, Pollution, Ethnicity, Immigration,
                                                 Paternal_Age, Paternal_SES, Parental_SMI, Adult_Life_Events,
                                                 Tobacco, Cannabis, Hearing, Childhood_Trauma, Anhedonia, CAARMS, gafex01, gafex02),
                          by=chr,
                          statistic = all_continuous() ~ c("{mean}", "{sd}"),
                          type = list(CAARMS ~ 'continuous2',
                                      gafex01 ~ 'continuous2',
                                      gafex02 ~ 'continuous2'),
                          digits = list(all_continuous() ~ c(1, 1),
                                        all_categorical() ~ c(0, 1)))  
table_HC %>% as_gt() %>% gt::gtsave("~/Dropbox/Work/PPS/EU-GEI/Table_HC.docx")

results <- data.frame(model=c("PPS","PRS", "PPS+PRS"),
                      C_index=c(1:3),
                      balanced_accuracy=c(1:3),
                      sensitivity=c(1:3),
                      specificity=c(1:3),
                      ppv=c(1:3),
                      npv=c(1:3),
                      calibration_intercept=c(1:3),
                      calibration_slope=c(1:3),
                      Brier=c(1:3),
                      ClinicalUtilityClin=c(1:3),
                      ClinicalUtilityGen=c(1:3)
)
predictors <- list(a=c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                       "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia"),
                   b=c("ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),
                   c=c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                          "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))
plot_data <- data.frame(x1 = seq(0, 1, length.out = 100) , 
                        x2 = seq(0, 1, length.out = 100) , 
                        x3 = seq(0, 1, length.out = 100) , 
                        y1 = c(1:100),
                        y2 = c(1:100),
                        y3 = c(1:100))
plot_pred_norm <- data.frame(z1= c(1:3430),
                             z2 = c(1:3430),
                             z3 = c(1:3430))

calibration_stats <- function(observed, predicted) {
  mean_observed <- mean(observed)
  mean_predicted <- mean(predicted)
  slope <- cov(predicted, observed) / var(predicted)
  intercept <- mean_observed - slope * mean_predicted
  return(list(intercept = intercept, slope = slope))
}

PPS_scored_all[,c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
              "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia")] <- lapply(PPS_scored_all[,c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                                                                                                                                  "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia")], FUN=factor)
###### Fit model and internal validation #######

folds <- repeatfolds(y=factor(PPS_scored_all$chr), repeats = 10, n_outer_folds = 10)
foldids = rep(1,length(folds)*length(folds$Rep1))

for (i in 1:3){
set.seed(123)
temp <- PPS_scored_all %>% select(all_of((predictors[[i]])))
temp_mat <- model.matrix(~.-1,temp)

PRS_multi <- nestcv.glmnet(y=factor(PPS_scored_all$chr), x=temp_mat,family="binomial",
                           n_outer_folds = 10, balance = "randomsample", alpha=1, seed=123) |>
  repeatcv(10, repeat_folds = folds, rep.cores = 2)

PPS_ct <- matrix(0,ncol = 1, nrow = 1)
PPS_ct$FP <- sum(PRS_multi[["output"]]$predy==1 & PRS_multi[["output"]]$testy==0)
PPS_ct$TN <- sum(PRS_multi[["output"]]$predy==0 & PRS_multi[["output"]]$testy==0) 
PPS_ct$TP <- sum(PRS_multi[["output"]]$predy==1 & PRS_multi[["output"]]$testy==1) 
PPS_ct$FN <- sum(PRS_multi[["output"]]$predy==0 & PRS_multi[["output"]]$testy==1)

PPS_ct$total <- sum(PPS_ct$FP, PPS_ct$TN, PPS_ct$TP, PPS_ct$FN)

results$C_index[i] <- paste(round(mean(PRS_multi$result[1:10])*100,1),"% (",
                                      round(mean(PRS_multi$result[1:10]*100,1)-1.96*sd(PRS_multi$result[1:10])*100,1),"%-",
                            round(mean(PRS_multi$result[1:10]*100,1)+1.96*sd(PRS_multi$result[1:10])*100,1), "%)", sep="")
results$balanced_accuracy[i] <- paste(round(sum((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN)),
                                                  (PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN))*PPS_ct$total), 
                                                          n = PPS_ct$total)$lower[2]),
                                                  (propCI(x = ((PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP))*PPS_ct$total), 
                                                          n = PPS_ct$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN))*PPS_ct$total), 
                                                          n = PPS_ct$total)$upper[2]),
                                                  (propCI(x = ((PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP))*PPS_ct$total), 
                                                          n = PPS_ct$total)$upper[2]))*50,1), "%)", sep="")
results$sensitivity[i] <- paste(round((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN))*100,1),"% (",
                                  round(propCI(x = ((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN))*PPS_ct$total), 
                                               n = PPS_ct$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN))*PPS_ct$total), 
                                               n = PPS_ct$total)$upper[2]*100,1), "%)", sep="")
results$specificity[i] <- paste(round((PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP))*100,1),"% (",
                                  round(propCI(x = ((PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP))*PPS_ct$total), 
                                               n = PPS_ct$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP))*PPS_ct$total), 
                                               n = PPS_ct$total)$upper[2]*100,1), "%)", sep="")

PRS_multi[["output"]]$predy <- as.factor(PRS_multi[["output"]]$predy)
results$ppv[i] <- paste0(round(yardstick::ppv(data = PRS_multi[["output"]], truth = testy, estimate = predy)$.estimate*100,1),"%")
results$npv[i] <- paste0(round(yardstick::npv(data = PRS_multi[["output"]], truth = testy, estimate = predy)$.estimate*100,1),"%")

# Calculate calibration slope on the hold-out data (test set)
logistic_calibration <- calibration_stats(observed = as.numeric(PRS_multi[["output"]]$testy)-1, predicted = as.numeric(PRS_multi[["output"]]$predy))
logistic_calibration2 <- predRupdate::pred_val_probs(binary_outcome = as.numeric(PRS_multi[["output"]]$testy)-1, Prob = rescale(PRS_multi[["output"]]$predyp))

results$calibration_intercept[i] <- round(logistic_calibration$intercept[1],2)
results$calibration_slope[i] <- round(logistic_calibration$slope[1],2)

Brier <- data.frame(obs=as.numeric(PRS_multi[["output"]]$testy)-1,
                      pred=scales::rescale((PRS_multi[["output"]]$predyp), to=c(0,1)))
Brier$Brier <- (Brier$obs - Brier$pred)^2
results$Brier[i] <- mean(Brier$Brier)

plot_data[,i+3] <- logistic_calibration$slope[1] * plot_data[,i] + logistic_calibration$intercept[1]
plot_pred_norm[,i] <- (PRS_multi[["output"]]$predyp-min(PRS_multi[["output"]]$predyp))/(max(PRS_multi[["output"]]$predyp)-min(PRS_multi[["output"]]$predyp))

# Decision curve analysis
library(dcurves)
dca_temp <- PRS_multi[["output"]] %>% filter(rep==1)
dca <- data.frame(obs=as.numeric(dca_temp$testy)-1,
                  pred=as.numeric(dca_temp$predy)-1)

dca_assessment <- dca(obs ~ pred,
                      data = dca,
                      prevalence = 0.192,
                      thresholds = seq(0, 0.5, 0.01)
) %>%
  as_tibble() 
results$ClinicalUtilityClin[i] <- dca_assessment %>% 
  filter(label=="pred" & net_benefit>0) %>%
  summarise(net_benefit = mean(net_benefit), .groups = "drop")

  # plot cross validated net benefit values
  ggplot(data=dca_assessment, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2) +
  coord_cartesian(ylim = c(-0.005, 0.2)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,0.35)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  theme_classic()
ggsave(paste0("~/Dropbox/Work/PPS/EU-GEI/Plots/dca_det_clin_", results[i,1],".png"), width=6, height=8, dpi=300)
write_csv(dca_assessment, paste0("~/Dropbox/Work/PPS/EU-GEI/dcaClin_", results[i,1],".csv"))

dca_assessment <- dca(obs ~ pred,
                      data = dca,
                      prevalence = 0.017,
                      thresholds = seq(0, 0.5, 0.001)
) %>%
  as_tibble() 
results$ClinicalUtilityGen[i] <- dca_assessment %>% 
  filter(label=="pred" & net_benefit>0) %>%
  summarise(net_benefit = mean(net_benefit), .groups = "drop")
  # plot cross validated net benefit values
  ggplot(data=dca_assessment,aes( x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2) +
  coord_cartesian(ylim = c(-0.005, 0.02)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,0.1)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  theme_classic()
ggsave(paste0("~/Dropbox/Work/PPS/EU-GEI/Plots/dca_det_gen_", results[i,1],".png"), width=6, height=8, dpi=300)
write_csv(dca_assessment, paste0("~/Dropbox/Work/PPS/EU-GEI/dcaGen_", results[i,1],".csv"))
}

write_csv(results, "~/Dropbox/Work/PPS/EU-GEI/results.csv")

PRS_multi[["output"]]$predyp_norm <- (PRS_multi[["output"]]$predyp-min(PRS_multi[["output"]]$predyp))/(max(PRS_multi[["output"]]$predyp)-min(PRS_multi[["output"]]$predyp))
plot_data <- merge(plot_data, plot_pred_norm)

ggplot() +
geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey") +  # Line segment from (0,0) to (1,1)
geom_smooth(data=plot_data, aes(x1,y1,color = "PPS"), method = "lm", se = FALSE) +
geom_smooth(data=plot_data, aes(x2,y2,color = "PRS"), method = "lm", se = FALSE) +
geom_smooth(data=plot_data, aes(x3,y3,color = "PPS+PRS"), method = "lm", se = FALSE) +
  geom_point(aes(plot_pred_norm$z1, as.numeric(PRS_multi[["output"]]$testy)-1), alpha = 0.5) +
labs(x = "Predicted probabilities", y = "Observed probabilities", color = "Legend") +
scale_color_manual(values = c("PPS" = "#c8526a", "PRS" = "#758EAB", "PPS+PRS" = "#7bae72")) +
theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Calibration_hc_300824_test.png",width = 24, height = 18, units = "cm")

ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey") +  # Line segment from (0,0) to (1,1)
  geom_abline(intercept=0.48, slope=1.03, color = "#599ec4", linewidth=1.2) +
  geom_abline(intercept=0.62, slope=0.51, color = "#ecc363", linewidth=1.2) +
  geom_abline(intercept=0.32, slope=1.10, color = "#7bae72", linewidth=1.2) +
  geom_point(aes(plot_pred_norm$z1, as.numeric(PRS_multi[["output"]]$testy)-1), alpha = 0.5) +
  labs(x = "Predicted probabilities", y = "Observed probabilities", color = "Legend") +
  scale_color_manual(values = c("PPS" = "red", "PRS" = "blue", "PPS+PRS" = "green")) +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Calibration_hc_300824.png",width = 24, height = 18, units = "cm")

calibration <- data.frame(observed=as.numeric(PRS_multi[["output"]]$testy)-1,
                          predicted = as.numeric(PRS_multi[["output"]]$predyp_norm))

calibration %>% cal_plot_logistic(observed, predicted)

library(dcurves)
dca_temp <- PRS_multi[["output"]] %>% filter(rep==1)
dca <- data.frame(obs=as.numeric(dca_temp$testy)-1,
                  pred=as.numeric(dca_temp$predy)-1)

dca_assessment <- dca(obs ~ pred,
                      data = dca,
                      prevalence = 0.192,
                      thresholds = seq(0, 0.5, 0.01)
) %>%
  as_tibble() 
dca_assessment %>% 
  group_by(variable, label, threshold) %>%
  summarise(net_benefit = mean(net_benefit), .groups = "drop") %>%
  # plot cross validated net benefit values
  ggplot(data=dca_assessment, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2) +
  coord_cartesian(ylim = c(-0.005, 0.2)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,0.35)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/dca_det_clin.png", width=6, height=8, dpi=300)

dca_assessment <- dca(obs ~ pred,
    data = dca,
    prevalence = 0.017,
    thresholds = seq(0, 0.5, 0.01)
) %>%
  as_tibble() 
dca_assessment %>% 
  group_by(variable, label, threshold) %>%
  summarise(net_benefit = mean(net_benefit), .groups = "drop") %>%
  # plot cross validated net benefit values
  ggplot(aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2) +
  coord_cartesian(ylim = c(-0.005, 0.02)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,0.1)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/dca_det_gen.png", width=6, height=8, dpi=300)

library(rms)
# Calibration plot for logistic regression
calibration_plot_logistic <- calibrate(PRS_multi, method = "boot", B = 100)
plot(calibration_plot_logistic, main = "Calibration Plot for Logistic Regression")

# Calibration plot for Cox model
calibration_plot_cox <- calibrate(fit_cox, method = "boot", B = 100)
plot(calibration_plot_cox, main = "Calibration Plot for Cox Model")

PRS_multi <- glm(chr ~ ZRE_SCZ_imp_0.1,family=binomial(link='logit'),data=PPS_scored)
summary(PRS_multi)
fmsb::NagelkerkeR2(PRS_multi)

PRS_multi <- glm(chr ~ PPS + ZRE_SCZ_imp_0.1,family=binomial(link='logit'),data=PPS_scored)
summary(PRS_multi)
fmsb::NagelkerkeR2(PRS_multi)

PPS_ind_multi <- glm(chr ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity + 
                       Immigration  + Paternal_Age + Paternal_SES + Parental_SMI + 
                       Adult_Life_Events + Tobacco + Cannabis + Childhood_Trauma + 
                       Anhedonia,family=binomial(link='logit'),data=PPS_scored)
summary(PPS_ind_multi)
fmsb::NagelkerkeR2(PPS_ind_multi)

PPS_ind_PRS_multi <- glm(chr ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity + 
                           Immigration  + Paternal_Age + Paternal_SES + Parental_SMI + 
                           Adult_Life_Events + Tobacco + Cannabis + Childhood_Trauma + 
                           Anhedonia + ZRE_SCZ_imp_0.1,family=binomial(link='logit'),data=PPS_scored)
summary(PPS_ind_PRS_multi)
fmsb::NagelkerkeR2(PPS_ind_PRS_multi)

###### Generate coefficients #######

temp <- PPS_scored_all %>% select(all_of((predictors[[1]])))
temp_mat <- model.matrix(~.-1,temp)
check <- glmnet(y=factor(PPS_scored_all$chr), x=temp_mat,family="binomial",
                n_outer_folds = 10, alpha=1, seed=123)
tmp_coeffs <- coef(check, s = cv.glmnet(as.matrix(temp_mat), PPS_scored_all$chr)$lambda.min)
PPS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored_all %>% select(all_of((predictors[[2]])))
temp_mat <- model.matrix(~.-1,temp)
check <- glmnet(y=factor(PPS_scored_all$chr), x=temp_mat,family="binomial",
                n_outer_folds = 10, alpha=1, seed=123)
tmp_coeffs <- coef(check, s = cv.glmnet(as.matrix(temp_mat), PPS_scored_all$chr)$lambda.min)
PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored_all %>% select(all_of((predictors[[3]])))
temp_mat <- model.matrix(~.-1,temp)
check <- glmnet(y=factor(PPS_scored_all$chr), x=temp_mat,family="binomial",
                n_outer_folds = 10, alpha=1, seed=123)
tmp_coeffs <- coef(check, s = cv.glmnet(as.matrix(temp_mat), PPS_scored_all$chr)$lambda.min)
PPS_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

write_csv(PPS_coef,"~/Dropbox/Work/PPS/EU-GEI/PPS_coef.csv")
write_csv(PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/PRS_coef.csv")
write_csv(PPS_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_coef.csv")

#### Prognosis ####
PPS_scored[,c(3:8,10:14,16:17)] <- lapply(PPS_scored[,c(3:8,10:14,16:17)], factor)

results_cox <- data.frame(model=c("Clinical","PPS","PRS","Clinical+PPS","Clinical+PRS", "PPS+PRS", "All"),
                      C_index=c(1:7),
                      calibration_intercept=c(1:7),
                      calibration_slope=c(1:7),
                      Brier=c(1:7)
)
predictors_cox <- list(a=c("CAARMS", "gafex01", "gafex02"),
                   b=c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                       "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia"),
                   c=c("ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),
                   d=c("CAARMS", "gafex01", "gafex02","Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                       "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia"),
                   e=c("CAARMS", "gafex01", "gafex02","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),
                   f=c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                       "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),
                   g=c("CAARMS", "gafex01", "gafex02","Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                       "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))

calibration_stats <- function(observed, predicted) {
  mean_observed <- mean(observed)
  mean_predicted <- mean(predicted)
  slope <- cov(predicted, observed) / var(predicted)
  intercept <- mean_observed - slope * mean_predicted
  return(list(intercept = intercept, slope = slope))
}

PPS_scored[,c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
              "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia")] <- lapply(PPS_scored[,c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
"Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia")], FUN=factor)

PPS_scored <- PPS_scored %>% mutate(Transition = case_when(day_exit>365.25*3 ~ 0,
                                                           TRUE ~ Transition),
                                    day_exit = case_when(day_exit>365.25*3 ~ 365.25*3,
                                                         TRUE ~ day_exit))
###### Fit model and internal validation #######
folds <- repeatfolds(y=factor(PPS_scored$Transition), repeats = 10, n_outer_folds = 10)

for (i in 1:7){
  set.seed(123)
  temp <- PPS_scored %>% select(all_of((predictors_cox[[i]])))
  temp_mat <- as.data.frame(model.matrix(~.-1,temp))

  set.seed(123)
  cv.PRS_multi <- nestcv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                                n_outer_folds = 10, type.measure = "C", alpha=1, seed=123) |>
    repeatcv(10, repeat_folds = folds, rep.cores = 4)


  results_cox$C_index[i] <- paste(round(mean(cv.PRS_multi$result[1:10]),2)," (",
                                        round(mean(cv.PRS_multi$result[1:10])-1.96*sd(cv.PRS_multi$result[1:10]),2),"-",
                              round(mean(cv.PRS_multi$result[1:10])+1.96*sd(cv.PRS_multi$result[1:10]),2),")", sep="")
    #predicting calibrated scores
  check <- data.frame(observed = as.numeric(cv.PRS_multi[["output"]]$status),
                      #predicted = as.numeric(cv.PRS_multi[["output"]]$predy))
                      predicted = 1 - 0.95402283^exp(as.numeric(cv.PRS_multi[["output"]]$predy)))
  
  logistic_calibration <- calibration_stats(observed = as.numeric(cv.PRS_multi[["output"]]$status), predicted = 1 - 0.95402283^exp(as.numeric(cv.PRS_multi[["output"]]$predy)))
  #logistic_calibration <- calibration_stats(PPS_scored_all$chr, PPS_scored_all$prob)
  results_cox$calibration_intercept[i] <- mean(check$predicted)/mean(check$observed)
  results_cox$calibration_slope[i] <- logistic_calibration$slope[1]
  
  Brier <- data.frame(obs = as.numeric(cv.PRS_multi[["output"]]$status),
                       #predicted = as.numeric(cv.PRS_multi[["output"]]$predy))
                      pred = 1 - 0.95402283^exp(as.numeric(cv.PRS_multi[["output"]]$predy)))
  Brier$Brier <- (Brier$obs * (1 - Brier$pred)^2 + (1 - Brier$obs) * (Brier$pred)^2)
  results_cox$Brier[i] <- mean(Brier$Brier)
}
write_csv(results_cox, "~/Dropbox/Work/PPS/EU-GEI/results_cox_071024.csv")

calibration <- data.frame(observed=as.numeric(cv.PRS_multi[["output"]]$status),
                          predicted = 1 - 0.95402283^exp(as.numeric(cv.PRS_multi[["output"]]$predy)))

calibration %>% cal_plot_logistic(observed, predicted)


ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey") +  # Line segment from (0,0) to (1,1)
  geom_smooth(aes(PPS_scored$prob, PPS_scored$Transition, color = "Original"), method = "loess", se = FALSE) +
  geom_smooth(aes(PPS_scored$prob_recal, PPS_scored$Transition, color = "Re-calibrated"), method = "loess", se = FALSE) +
  geom_point(aes(PPS_scored$prob, PPS_scored$Transition), alpha = 0.5) +
  labs(x = "Predicted probabilities", y = "Observed probabilities", color = "Legend") +
  scale_color_manual(values = c("Original" = "red", "Re-calibrated" = "blue")) +
  theme_classic()

ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Calibration_190324.png",width = 24, height = 18, units = "cm")

###### Generate coefficients #####

temp <- PPS_scored %>% select(all_of((predictors_cox[[1]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[2]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
PPS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[3]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[4]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_PPS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[5]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[6]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
PPS_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[7]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_PPS_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

write_csv(clin_coef,"~/Dropbox/Work/PPS/EU-GEI/clin_coef_surv_200824.csv")
write_csv(PPS_coef,"~/Dropbox/Work/PPS/EU-GEI/PPS_coef_surv_200824.csv")
write_csv(PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/PRS_coef_surv_200824.csv")
write_csv(clin_PPS_coef,"~/Dropbox/Work/PPS/EU-GEI/clinPPS_coef_surv_200824.csv")
write_csv(clin_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/clin_PRS_coef_surv_200824.csv")
write_csv(PPS_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_coef_surv_200824.csv")
write_csv(clin_PPS_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/clin_PPS_PRS_coef_surv_200824.csv")

##### Detection Random Forest ####
    
    # Parameters for nested cross-validation
    outerFolds <- 5
    outerRepeats <- 10
    innerFolds <- 5
    innerRepeats <- 10
  
    # Wrangle data
    PPS_scored_all[,c(3:8,10:14,16:17)] <- lapply(PPS_scored_all[,c(3:8,10:14,16:17)], factor)
    PPS_scored_mat <- as.data.frame(model.matrix(~.-1,PPS_scored_all[,c(3:8,10:14,16:17,23:39)]))
    PPS_scored_mat$chr <- PPS_scored_all$chr
    PPS_hc_rf <- as.data.frame(PPS_scored_mat)
    PPS_hc_rf <- PPS_hc_rf %>% subset(select = c(-ZRE_SCZ_imp_0.1,-ZRE_BD_imp_EUR_0.1:-gafex02))

    PPS_hc_rf[,1:19] <- lapply(PPS_hc_rf[,1:19], factor)
    PPS_hc_rf$chr <- factor(PPS_hc_rf$chr)
    
    PPS_hc_rf$chr <- factor(make.names(levels(PPS_hc_rf$chr))[PPS_hc_rf$chr])

    PPS_hc_rf <- PPS_hc_rf  %>% dplyr::rename(Ethnicity35 = Ethnicity3.5,
                                              Immigration25 = Immigration2.5,
                                              Parental_SMI55 = Parental_SMI5.5,
                                              Adult_Life_Events55 = Adult_Life_Events5.5,
                                              Anhedonia65 = Anhedonia6.5)
    make.names(colnames(PPS_hc_rf), unique=TRUE)

    ## Compute F1-score
    
    F1_score <- function(mat, algoName){
      
      # Remark: left column = prediction // top = real values
      recall <- matrix(1:nrow(mat), ncol = nrow(mat))
      precision <- matrix(1:nrow(mat), ncol = nrow(mat))
      F1_score <- matrix(1:nrow(mat), ncol = nrow(mat))
      
      
      for(i in 1:nrow(mat)){
        recall[i] <- mat[i,i]/rowSums(mat)[i]
        precision[i] <- mat[i,i]/colSums(mat)[i]
      }
      
      for(i in 1:ncol(recall)){
        F1_score[i] <- 2 * ( precision[i] * recall[i] ) / ( precision[i] + recall[i])
      }
      
      # We display the matrix labels
      colnames(F1_score) <- colnames(mat)
      rownames(F1_score) <- algoName
      
      # Display the F1_score for each class
      F1_score
      
      # Display the average F1_score
      mean(F1_score[1,])
    }
    
    # Create a function to perform nested cross-validation with repeats
    nested_cv_with_repeats <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, tuneGrid, seed) {
      
      set.seed(seed)
      seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
      
      best_inner_result_list <- list()
      best_mtry_list <- list()
      best_ntree_list <- list()
      best_nodesize_list <- list()
      all_inner_results <- list()
      c_stat_nested <- data.frame()
      calibration_slopes <- data.frame()
      
      # Outer loop for cross-validation
      for (outer_rep in 1:outerRepeats) {
        cat("Outer Repeat:", outer_rep, "\n")
        set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
        outer_folds <- createFolds(combined_df$chr, k = outerFolds, list = TRUE, returnTrain = TRUE)
        
        for (i in seq_along(outer_folds)) {
          cat("  Outer Fold:", i, "\n")
          
          train_indices <- outer_folds[[i]]
          train <- combined_df[train_indices, ] #  This splits the data into train (inner loop)
          test <- combined_df[-train_indices, ] #  This splits the data into test (outer loop)
          
          inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
          
          # Inner Cross-Validation and Model Training
          for (inner_rep in 1:innerRepeats) {
            cat("    Inner Repeat:", inner_rep, "\n")
            
            # Use a different seed for each inner repeat
            set.seed(seeds[outer_rep] + inner_rep)
            
            # Tune both alpha and lambda
            for (mtry_value in tuneGrid) {
              library(randomForest)
              library(MLmetrics)
              library(caret)
              
              # Define a custom summary function to calculate F1 score
              customSummary <- function(data, lev = NULL, model = NULL) {
                f1 <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[1])
                out <- c(F1 = f1)
                out
              }
              
              # Set up the trainControl with the custom summary function
              control <- trainControl(method = 'cv', 
                                      number = 5, 
                                      #repeats = 10,
                                      classProbs = TRUE,
                                      summaryFunction = customSummary,
                                      search = 'random')
            
              tree <- c(50, 100, 250, 500)
              n.tree <- sample(tree,1)
              nodeSize <- seq(1,(nrow(train)/10), by=1)
              node.size <- sample(nodeSize,1)
              tune_grid_temp <- data.frame(mtry=c(NA,NA,NA,NA,NA))
              tune_grid_temp$mtry <- sample(tuneGrid$mtry,5)
              
              # Train the model using the custom F1 metric
              inner_model <- train(chr ~ .,
                                   data = train,
                                   method = "rf",
                                   metric = "F1",
                                   tuneGrid = tune_grid_temp,
                                   tuneLength=10,
                                   ntree = n.tree,
                                   nodesize=node.size,
                                   trControl = control)
              
              print(inner_model)
              
              if (all(!is.na(inner_model$results$F1))){
              # Store the F1 for each inner repeat
              inner_results[[inner_rep]] <- rbind(inner_results[[inner_rep]], 
                                                    data.frame(min.node.size = node.size,
                                                               tree = n.tree,
                                                               mtry = inner_model$results$mtry, 
                                                               F1 = inner_model$results$F1, 
                                                               repeat_number = inner_rep))
              }
            }
          }
          
          # Combine results from all repeats
          all_inner_results <- do.call(rbind, inner_results)
          print(dim(all_inner_results))
          
          # Calculate the average F1 for each hyperparameter combination
          avg_results <- all_inner_results %>%
            group_by(mtry) %>%
            summarise(avg_F1 = mean(F1, na.rm=TRUE),
                      mtry = mean(mtry, na.rm=TRUE),
                      ntree = mean(tree, na.rm=TRUE),
                      nodesize=mean(min.node.size, na.rm=TRUE),
                      .groups = "keep") 
          
          # Determine the best hyperparameter combination based on the highest average F1
          best_inner_result <- avg_results[which.max(avg_results$avg_F1), ]
          best_min.node.size <- best_inner_result$nodesize
          best_mtry <- best_inner_result$mtry
          best_F1 <- best_inner_result$avg_F1
          best_F1_SE <- best_inner_result$SE_F1
          best_ntree <- best_inner_result$ntree
          
          # Calculate the correct index for the current outer fold and repeat
          index <- (outer_rep - 1) * outerFolds + i
          
          # Store the best hyperparameters for the current outer fold
          best_inner_result_list[[index]] <- best_inner_result
          best_mtry_list[[index]]  <- best_mtry
          best_ntree_list[[index]]  <- best_ntree
          best_nodesize_list[[index]]  <- best_min.node.size
          
          # Train the final model using the best hyperparameters on the inner training set (x_inner_train, y_inner_train)
          control_final <- trainControl(method = 'none')
          repGrid <- data.frame(mtry=best_mtry)  
          
          final_model <- train(chr ~ .,
                         data = train,
                         method = "rf",
                         metric = "F1",
                         ntree = best_ntree,
                         nodesize=best_min.node.size,
                         trControl=control_final,
                         tuneGrid = repGrid)
          
          # Predict the linear predictors (PI) from the Elastic Net model
          test$PI <- predict(final_model, newdata = test, type = "prob")[,2]
          test$chr_pred <- predict(final_model, newdata = test, type = "raw")
            
          cm <- confusionMatrix(data = test$chr_pred, reference = test$chr)
          test <- test %>% mutate(PI=case_when(PI==0 ~ 0.001, PI==1 ~ 0.999, TRUE ~ PI))
          
          # Fit a Cox model on the test set using penalized coefficients
          model_test <- glm(chr ~ PI, data = test, family="binomial")
          
          c_stat_nested <- rbind(c_stat_nested, data.frame(
            C_test = concordance(model_test)$concordance,
            SE_test = concordance(model_test)$cvar,
            Fold = i,
            OuterRepeat = outer_rep,
            n_train = nrow(train),
            events_train = sum(train$chr=="X1"),
            n_test = nrow(test),
            events_test = sum(test$chr=="X1"),
            balanced_accuracy = cm$byClass[11],
            sensitivity = cm$byClass[1],
            specificity = cm$byClass[2],
            ppv = cm$byClass[3],
            npv = cm$byClass[4],
            precision = cm$byClass[5],
            recall = cm$byClass[6],
            f1 = cm$byClass[7]
            
          ))
          
          calibration <- rms::val.prob(p=test$PI, y=as.numeric(test$chr)-1, m=200, pl=F)
          
          # Calculate calibration slope on the hold-out data (test set)
          
          calibration_intercept <- unname(calibration[12])
          calibration_slope <- unname(calibration[13])
          brier <- unname(calibration[11])
          
          
          # Store calibration results
          calibration_slopes <- rbind(calibration_slopes, data.frame(
            fold = i,
            OuterRepeat = outer_rep,
            intercept = calibration_intercept,
            slope = calibration_slope,
            brier = brier
          ))
        }
      }
      
      list(
        best_inner_result_list = best_inner_result_list,
        best_mtry_list = best_mtry_list,
        best_ntree_list = best_ntree_list,
        best_nodesize_list = best_nodesize_list,
        # C_results = C_results,
        # C_SE_results = C_SE_results,
        c_stat_nested = c_stat_nested,
        calibration_slopes = calibration_slopes
      )
    }
    
    # Define tuning grid
    
    tune_grid <- expand.grid(
      mtry=c(3:30)
    )
   
    ###### Run model fitting and internal validation ######
    
    # Perform nested cross-validation with repeats
    PPS_PRS <- PPS_hc_rf %>% subset(select=c(Gender0:Anhedonia65, PC1:PC10, ZRE_SCZ_imp_EUR_0.1, chr))

    results_all <- nested_cv_with_repeats(
      combined_df = PPS_PRS, 
      outerFolds = outerFolds, 
      outerRepeats = outerRepeats, 
      innerFolds = innerFolds, 
      innerRepeats = innerRepeats, 
      tuneGrid = tune_grid,
      seed = 231
    )
    
    PPS_genes <- PPS_hc_rf %>% subset(select=c(PC1:PC10, ZRE_SCZ_imp_EUR_0.1, chr))
    tune_grid_genes <- expand.grid(
      mtry=c(3:11)
    )
    results_genes <- nested_cv_with_repeats(
      combined_df = PPS_genes, 
      outerFolds = outerFolds, 
      outerRepeats = outerRepeats, 
      innerFolds = innerFolds, 
      innerRepeats = innerRepeats, 
      tuneGrid = tune_grid_genes,
      seed = 231
    )
    
    PPS_PPS <- PPS_hc_rf %>% subset(select=c(Gender0:Anhedonia65, chr))
    tune_grid_PPS <- expand.grid(
      mtry=c(3:19)
    )
    results_PPS <- nested_cv_with_repeats(
      combined_df = PPS_PPS, 
      outerFolds = outerFolds, 
      outerRepeats = outerRepeats, 
      innerFolds = innerFolds, 
      innerRepeats = innerRepeats, 
      tuneGrid = tune_grid_PPS,
      seed = 231
    )
    
###### Generate performance metrics summary table######
    
rf_results <- data.frame(model=c("PPS", "PRS", "All"),
                         C=c(paste0(round(mean(results_PPS$c_stat_nested$C),2), " (",
                                    round(mean(results_PPS$c_stat_nested$C) - 1.96*(sd(results_PPS$c_stat_nested$C) / sqrt(nrow(results_PPS$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_PPS$c_stat_nested$C) + 1.96*(sd(results_PPS$c_stat_nested$C) / sqrt(nrow(results_PPS$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_genes$c_stat_nested$C),2), " (",
                                    round(mean(results_genes$c_stat_nested$C) - 1.96*(sd(results_genes$c_stat_nested$C) / sqrt(nrow(results_genes$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_genes$c_stat_nested$C) + 1.96*(sd(results_genes$c_stat_nested$C) / sqrt(nrow(results_genes$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_all$c_stat_nested$C),2), " (",
                                    round(mean(results_all$c_stat_nested$C) - 1.96*(sd(results_all$c_stat_nested$C) / sqrt(nrow(results_all$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_all$c_stat_nested$C) + 1.96*(sd(results_all$c_stat_nested$C) / sqrt(nrow(results_all$c_stat_nested))),2),
                                    ")")),
                         BAC=c(paste0(round(mean(results_PPS$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                      round((mean(results_PPS$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_PPS$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_PPS$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                      round((mean(results_genes$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_genes$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_genes$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                      round((mean(results_all$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_all$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_all$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%)")),
                         sensitivity=c(paste0(round(mean(results_PPS$c_stat_nested$sensitivity)*100,1), "% (",
                                      round((mean(results_PPS$c_stat_nested$sensitivity) - 1.96*(sd(results_PPS$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS$c_stat_nested$sensitivity) + 1.96*(sd(results_PPS$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes$c_stat_nested$sensitivity)*100,1), "% (",
                                      round((mean(results_genes$c_stat_nested$sensitivity) - 1.96*(sd(results_genes$c_stat_nested$sensitivity) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes$c_stat_nested$sensitivity) + 1.96*(sd(results_genes$c_stat_nested$sensitivity) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all$c_stat_nested$sensitivity)*100,1), "% (",
                                      round((mean(results_all$c_stat_nested$sensitivity) - 1.96*(sd(results_all$c_stat_nested$sensitivity) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all$c_stat_nested$sensitivity) + 1.96*(sd(results_all$c_stat_nested$sensitivity) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%)")),
                         specificity=c(paste0(round(mean(results_PPS$c_stat_nested$specificity)*100,1), "% (",
                                              round((mean(results_PPS$c_stat_nested$specificity) - 1.96*(sd(results_PPS$c_stat_nested$specificity) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_PPS$c_stat_nested$specificity) + 1.96*(sd(results_PPS$c_stat_nested$specificity) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                              "%)"),
                                       paste0(round(mean(results_genes$c_stat_nested$specificity)*100,1), "% (",
                                              round((mean(results_genes$c_stat_nested$specificity) - 1.96*(sd(results_genes$c_stat_nested$specificity) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_genes$c_stat_nested$specificity) + 1.96*(sd(results_genes$c_stat_nested$specificity) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                              "%)"),
                                       paste0(round(mean(results_all$c_stat_nested$specificity)*100,1), "% (",
                                              round((mean(results_all$c_stat_nested$specificity) - 1.96*(sd(results_all$c_stat_nested$specificity) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_all$c_stat_nested$specificity) + 1.96*(sd(results_all$c_stat_nested$specificity) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                              "%)")),
                         ppv=c(paste0(round(mean(results_PPS$c_stat_nested$ppv*100),1),"%"),
                               paste0(round(mean(results_genes$c_stat_nested$ppv*100),1),"%"),
                               paste0(round(mean(results_all$c_stat_nested$ppv*100),1),"%")),
                         npv=c(paste0(round(mean(results_PPS$c_stat_nested$npv*100),1),"%"),
                               paste0(round(mean(results_genes$c_stat_nested$npv*100),1),"%"),
                               paste0(round(mean(results_all$c_stat_nested$npv*100),1),"%")),
                         intercept=c(round(mean(results_PPS$calibration_slopes$intercept),2),
                                     round(mean(results_genes$calibration_slopes$intercept),2),
                                     round(mean(results_all$calibration_slopes$intercept),2)),
                         slope=c(round(mean(results_PPS$calibration_slopes$slope),2),
                                 round(mean(results_genes$calibration_slopes$slope),2),
                                 round(mean(results_all$calibration_slopes$slope),2)),
                         Brier=c(round(mean(results_PPS$calibration_slopes$brier),2),
                                 round(mean(results_genes$calibration_slopes$brier),2),
                                 round(mean(results_all$calibration_slopes$brier),2))
                           
                         )
write_csv(rf_results, "~/Dropbox/Work/PPS/EU-GEI/rf_results_181024.csv")

###### Generate variable importance ######
control_final <- trainControl(method = 'none')
check <- as.data.frame(results_all$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_all$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))

final_model <- randomForest(chr ~ .,
                     data = PPS_PRS,
                     method = "rf",
                     metric = "F1",
                     ntree = 500,
                     nodesize=best_min.node.size,
                     trControl=control_final,
                     tuneGrid = repGrid)

# Decision curve analysis
library(dcurves)
dca <- data.frame(obs=as.numeric(predict(final_model, newdata = PPS_PRS, type = "response"))-1,
                  pred=predict(final_model, newdata = PPS_PRS, type = "prob")[,2])

dca_assessment <- dca(obs ~ pred,
                      data = dca,
                      prevalence = 0.192,
                      thresholds = seq(0, 0.5, 0.01)
) %>%
  as_tibble() 

write_csv(dca_assessment, paste0("~/Dropbox/Work/PPS/EU-GEI/dcaClin_rf_", results[i,1],".csv"))

dca_assessment <- dca(obs ~ pred,
                      data = dca,
                      prevalence = 0.017,
                      thresholds = seq(0, 0.5, 0.001)
) %>%
  as_tibble() 

write_csv(dca_assessment, paste0("~/Dropbox/Work/PPS/EU-GEI/dcaGen_rf_", results[i,1],".csv"))

# Function to calculate the depth of a single tree
calculate_tree_depth <- function(tree) {
  # Extract the left and right daughter nodes
  left_daughter <- tree$`left daughter`
  right_daughter <- tree$`right daughter`
  
  # Initialize depth to 0
  depth <- 0
  
  # Traverse through the tree to calculate depth
  while (any(left_daughter != 0 | right_daughter != 0)) {
    depth <- depth + 1
    left_daughter <- left_daughter[left_daughter != 0]
    right_daughter <- right_daughter[right_daughter != 0]
  }
  
  return(depth)
}

# Calculate the depth of each tree
tree_depths <- sapply(1:1, function(i) {
  tree <- getTree(final_model, k = i, labelVar = TRUE)
  calculate_tree_depth(tree)
})

# Calculate the mean depth of all trees
mean_tree_depth <- mean(tree_depths)

# Print the mean tree depth
print(mean_tree_depth)

varImp(final_model)
plot(varImp(final_model), top = 20)

check <- as.data.frame(results_PPS$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_PPS$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))

final_model <- train(chr ~ .,
                     data = PPS_PPS,
                     method = "rf",
                     metric = "F1",
                     #ntree = best_ntree,
                     nodesize=best_min.node.size,
                     trControl=control_final,
                     tuneGrid = repGrid)
varImp(final_model)
vimp_PPS <- varImp(final_model)$importance
vimp_PPS <- tibble::rownames_to_column(vimp_PPS, "Predictor")
vimp_PPS <- vimp_PPS %>% arrange(desc(Overall))
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PPS_VImp.csv")

pdf("~/Dropbox/Work/PPS/EU-GEI/PPS_VImp.pdf")
plot(varImp(final_model), top = 20)
dev.off()

check <- as.data.frame(results_genes$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_genes$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))

final_model <- train(chr ~ .,
                     data = PPS_genes,
                     method = "rf",
                     metric = "F1",
                     #ntree = best_ntree,
                     nodesize=best_min.node.size,
                     trControl=control_final,
                     tuneGrid = repGrid)
varImp(final_model)
vimp_PRS <- varImp(final_model)$importance
vimp_PRS <- tibble::rownames_to_column(vimp_PRS, "Predictor")
vimp_PRS <- vimp_PRS %>% arrange(desc(Overall))
write_csv(vimp_PRS, "~/Dropbox/Work/PPS/EU-GEI/PRS_VImp.csv")

pdf("~/Dropbox/Work/PPS/EU-GEI/PRS_VImp.pdf")
plot(varImp(final_model), top = 20)
dev.off()

check <- as.data.frame(results_all$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_all$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))

final_model <- train(chr ~ .,
                     data = PPS_PRS,
                     method = "rf",
                     metric = "F1",
                     #ntree = best_ntree,
                     nodesize=best_min.node.size,
                     trControl=control_final,
                     tuneGrid = repGrid)
vimp_all <- varImp(final_model)$importance
vimp_all <- tibble::rownames_to_column(vimp_all, "Predictor")
vimp_all <- vimp_all %>% arrange(desc(Overall))

write_csv(vimp_all, "~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_VImp.csv")

pdf("~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_VImp.pdf")
plot(varImp(final_model), top = 20)
dev.off()

##### Prognosis Random Survival Forest ####

# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

# Create a function to perform nested cross-validation with repeats
nested_cv_with_repeats_surv <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, tuneGrid, seed) {
  
  set.seed(seed)
  seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
  
  best_inner_result_list <- list()
  best_mtry_list <- list()
  best_ntree_list <- list()
  best_nodesize_list <- list()
  inner_results <- data.frame()
  all_inner_results <- data.frame()
  c_stat_nested  <- data.frame()
  calibration_slopes <- data.frame()
  
  # Outer loop for cross-validation
  for (outer_rep in 1:outerRepeats) {
    cat("Outer Repeat:", outer_rep, "\n")
    set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
    outer_folds <- createFolds(combined_df$Transition, k = outerFolds, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(outer_folds)) {
      cat("  Outer Fold:", i, "\n")
      
      train_indices <- outer_folds[[i]]
      train <- combined_df[train_indices, ] #  This splits the data into train (inner loop)
      test <- combined_df[-train_indices, ] #  This splits the data into test (outer loop)
      
      inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
      
      # Inner Cross-Validation and Model Training
      for (inner_rep in 1:innerRepeats) {
        cat("    Inner Repeat:", inner_rep, "\n")
        
        # Use a different seed for each inner repeat
        set.seed(seeds[outer_rep] + inner_rep)
        
          library(randomForestSRC)
          library(survcomp)
          
          tree <- c(50, 100, 250, 500)
          nodeSize <- seq(1,(nrow(train)/10), by=1)
          # Train the model using the custom F1 metric
            
            # Store fold-specific C-index scores
            fold_cindex <- numeric(5)
            folds_cv <- createFolds(train$Transition, k = 5, returnTrain = TRUE)
            
          for (fold_idx in 1:5) {
            # Get the training and validation sets for this fold
            
            n.tree <- sample(tree,1)
            node.size <- sample(nodeSize,1)
            mtry <- sample(tuneGrid$mtry,1)
            
            train_indices_cv <- folds_cv[[fold_idx]]

            fold_train <- train[train_indices_cv, ]
            fold_valid <- train[-train_indices_cv, ]
            
            predictor_vars <- setdiff(names(fold_train), c("day_exit", "Transition"))
            
            # Dynamically create the formula for the RSF model
            formula <- as.formula(paste("Surv(day_exit, Transition) ~", paste(predictor_vars, collapse = " + ")))
            
            # Train the RSF model on the training set
            rsf_model <- rfsrc(formula, data = fold_train, ntree = n.tree, mtry = mtry, nodesize = node.size)
            
            # Calculate the C-index on the validation set
            rsf_pred <- predict(rsf_model, fold_valid)
            fold_cindex[fold_idx] <- concordance.index(rsf_pred$survival[,ncol(rsf_pred$survival)], surv.time = fold_valid$day_exit, surv.event = fold_valid$Transition)$c.index
            
          }
          
          # Store the mean C-index for this combination of hyperparameters
          mean_cindex <- mean(fold_cindex)
          inner_results <- rbind(inner_results, data.frame(ntree = n.tree, mtry = mtry, nodesize=node.size, mean_cindex = mean_cindex))
        
        
        # View the results of cross-validation
        print(inner_results)
        
        # Find the best hyperparameter combination based on the highest C-index
        best_params <- inner_results[which.max(inner_results$mean_cindex), ]
        
      }
      
      # Combine results from all repeats
      all_inner_results <- do.call(rbind, inner_results)
      print(dim(all_inner_results))

      # Calculate the average C-index for each hyperparameter combination
      avg_results <- inner_results %>%
        group_by(mtry) %>%
        summarise(avg_C = mean(mean_cindex, na.rm=TRUE),
                  mtry = mean(mtry, na.rm=TRUE),
                  ntree = mean(ntree, na.rm=TRUE),
                  nodesize=mean(nodesize, na.rm=TRUE),
                  .groups = "keep") 
      
      # Determine the best hyperparametercombination based on the highest average C-index
      best_inner_result <- avg_results[which.max(avg_results$avg_C), ]
      best_min.node.size <- best_inner_result$nodesize
      best_mtry <- best_inner_result$mtry
      best_C <- best_inner_result$avg_C
      best_ntree <- best_inner_result$ntree
      
      # Calculate the correct index for the current outer fold and repeat
      index <- (outer_rep - 1) * outerFolds + i
      
      # Store the best hyperparameters for the current outer fold
      best_inner_result_list[[index]] <- best_inner_result
      best_mtry_list[[index]]  <- best_mtry
      best_ntree_list[[index]]  <- best_ntree
      best_nodesize_list[[index]]  <- best_min.node.size
      
      # Train the final model using the best hyperparameters on the inner training set (x_inner_train, y_inner_train)
      control_final <- trainControl(method = 'none')
      repGrid <- data.frame(mtry=best_mtry)  
      
      final_model <- rfsrc(formula, data = train, ntree = best_ntree, mtry = best_mtry, nodesize = best_min.node.size)
      
      # Predict the linear predictors (PI) from the Elastic Net model
      rsf_pred_test <- predict(final_model, test)
      test$PI <- 1-rsf_pred_test$survival[,ncol(rsf_pred_test$survival)]
      test$Tx_pred <- 1-rsf_pred_test$survival[,ncol(rsf_pred_test$survival)]
      
      test <- test %>% mutate(PI=case_when(PI==0 ~ 0.001, PI==1 ~ 0.999, TRUE ~ PI))
      
      # Fit a Cox model on the test set using penalized coefficients
      model_test <- coxph(Surv(day_exit, Transition) ~ PI, data = test)
      
      c_stat_nested <- rbind(c_stat_nested, data.frame(
        C_test = concordance(model_test)$concordance,
        SE_test = concordance(model_test)$cvar,
        Fold = i,
        OuterRepeat = outer_rep,
        n_train = nrow(train),
        events_train = sum(train$Transition==1),
        n_test = nrow(test),
        events_test = sum(test$Transition==1)
        
      ))
      
      calibration <- rms::val.prob(p=test$PI, y=test$Transition, m=200, pl=F)
      
      # Calculate calibration slope on the hold-out data (test set)
      check <- data.frame(observed = as.numeric(test$Transition),
                          #predicted = as.numeric(cv.PRS_multi[["output"]]$predy))
                          predicted = 1 - 0.95402283^exp(as.numeric(test$PI)))
      
      calibration_intercept <- mean(check$predicted)/mean(check$observed)
      calibration_slope <- unname(calibration[13])
      brier <- unname(calibration[11])
      
      # Store calibration results
      calibration_slopes <- rbind(calibration_slopes, data.frame(
        fold = i,
        OuterRepeat = outer_rep,
        intercept = calibration_intercept,
        slope = calibration_slope,
        brier = brier
      ))
    }
  }
  
  list(
    best_inner_result_list = best_inner_result_list,
    best_mtry_list = best_mtry_list,
    best_ntree_list = best_ntree_list,
    best_nodesize_list = best_nodesize_list,
    c_stat_nested = c_stat_nested,
    calibration_slopes = calibration_slopes
  )
}

# Define tuning grid

###### Run model fitting and internal validation ######
predictors_cox <- list(a=c("CAARMS", "gafex01", "gafex02","Transition", "day_exit"),
                       b=c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","Transition", "day_exit"),
                       c=c("ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","Transition", "day_exit"),
                       d=c("CAARMS", "gafex01", "gafex02","Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","Transition", "day_exit"),
                       e=c("CAARMS", "gafex01", "gafex02","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","Transition", "day_exit"),
                       f=c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","Transition", "day_exit"),
                       g=c("CAARMS", "gafex01", "gafex02","Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","Transition", "day_exit"))

PPS_scored[,c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
              "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia")] <- lapply(PPS_scored[,c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                                                                                                                                  "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia")], FUN=factor)

# Clinical
clin_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[1]])))
clin_rf_surv <- as.data.frame(model.matrix(~.-1,clin_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3,3,3)
)

results_clin_surv <- nested_cv_with_repeats_surv(
  combined_df = clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PPS

PPS_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[2]])))
PPS_rf_surv <- as.data.frame(model.matrix(~.-1,PPS_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_rf_surv)-2)
)

results_PPS_surv <- nested_cv_with_repeats_surv(
  combined_df = PPS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PRS
PRS_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[3]])))
PRS_rf_surv <- as.data.frame(model.matrix(~.-1,PRS_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PRS_rf_surv)-2)
)

results_PRS_surv <- nested_cv_with_repeats_surv(
  combined_df = PRS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PPS

PPS_clin_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[4]])))
PPS_clin_rf_surv <- as.data.frame(model.matrix(~.-1,PPS_clin_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_clin_rf_surv)-2)
)

results_PPS_clin_surv <- nested_cv_with_repeats_surv(
  combined_df = PPS_clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PRS

PRS_clin_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[5]])))
PRS_clin_rf_surv <- as.data.frame(model.matrix(~.-1,PRS_clin_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_rf_surv)-2)
)

results_PRS_clin_surv <- nested_cv_with_repeats_surv(
  combined_df = PRS_clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PPS + PRS

PPS_PRS_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[6]])))
PPS_PRS_rf_surv <- as.data.frame(model.matrix(~.-1,PPS_PRS_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_PRS_rf_surv)-2)
)

results_PPS_PRS_surv <- nested_cv_with_repeats_surv(
  combined_df = PPS_PRS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PPS + PRS 

PPS_PRS_clin_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[7]])))
PPS_PRS_clin_rf_surv <- as.data.frame(model.matrix(~.-1,PPS_PRS_clin_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_PRS_clin_rf_surv)-2)
)

results_PPS_PRS_clin_surv <- nested_cv_with_repeats_surv(
  combined_df = PPS_PRS_clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)
###### Generate performance metrics summary table######

rf_surv_results <- data.frame(model=c("Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All"),
                         C=c(paste0(round(mean(results_clin_surv$c_stat_nested$C),2), " (",
                                    round(mean(results_clin_surv$c_stat_nested$C) - 1.96*(sd(results_clin_surv$c_stat_nested$C) / sqrt(nrow(results_clin_surv$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_clin_surv$c_stat_nested$C) + 1.96*(sd(results_clin_surv$c_stat_nested$C) / sqrt(nrow(results_clin_surv$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_PPS_surv$c_stat_nested$C),2), " (",
                                    round(mean(results_PPS_surv$c_stat_nested$C) - 1.96*(sd(results_PPS_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_surv$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_PPS_surv$c_stat_nested$C) + 1.96*(sd(results_PPS_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_surv$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_PRS_surv$c_stat_nested$C),2), " (",
                                    round(mean(results_PRS_surv$c_stat_nested$C) - 1.96*(sd(results_PRS_surv$c_stat_nested$C) / sqrt(nrow(results_PRS_surv$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_PRS_surv$c_stat_nested$C) + 1.96*(sd(results_PRS_surv$c_stat_nested$C) / sqrt(nrow(results_PRS_surv$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_PPS_clin_surv$c_stat_nested$C),2), " (",
                                    round(mean(results_PPS_clin_surv$c_stat_nested$C) - 1.96*(sd(results_PPS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_clin_surv$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_PPS_clin_surv$c_stat_nested$C) + 1.96*(sd(results_PPS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_clin_surv$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_PRS_clin_surv$c_stat_nested$C),2), " (",
                                    round(mean(results_PRS_clin_surv$c_stat_nested$C) - 1.96*(sd(results_PRS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PRS_clin_surv$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_PRS_clin_surv$c_stat_nested$C) + 1.96*(sd(results_PRS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PRS_clin_surv$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_PPS_PRS_clin_surv$c_stat_nested$C),2), " (",
                                    round(mean(results_PPS_PRS_clin_surv$c_stat_nested$C) - 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_PPS_PRS_clin_surv$c_stat_nested$C) + 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_PPS_PRS_clin_surv$c_stat_nested$C),2), " (",
                                    round(mean(results_PPS_PRS_clin_surv$c_stat_nested$C) - 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_PPS_PRS_clin_surv$c_stat_nested$C) + 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))),2),
                                    ")")),
                         intercept=c(round(mean(results_clin_surv$calibration_slopes$intercept),2),
                                     round(mean(results_PPS_surv$calibration_slopes$intercept),2),
                                     round(mean(results_PRS_surv$calibration_slopes$intercept),2),
                                     round(mean(results_PPS_clin_surv$calibration_slopes$intercept),2),
                                     round(mean(results_PRS_clin_surv$calibration_slopes$intercept),2),
                                     round(mean(results_PPS_PRS_surv$calibration_slopes$intercept),2),
                                     round(mean(results_PPS_PRS_clin_surv$calibration_slopes$intercept),2)),
                         slope=c(round(mean(results_clin_surv$calibration_slopes$slope),2),
                                 round(mean(results_PPS_surv$calibration_slopes$slope),2),
                                 round(mean(results_PRS_surv$calibration_slopes$slope),2),
                                 round(mean(results_PPS_clin_surv$calibration_slopes$slope),2),
                                 round(mean(results_PRS_clin_surv$calibration_slopes$slope),2),
                                 round(mean(results_PPS_PRS_surv$calibration_slopes$slope),2),
                                 round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope),2)),
                         Brier=c(round(mean(results_clin_surv$calibration_slopes$brier),2),
                                 round(mean(results_PPS_surv$calibration_slopes$brier),2),
                                 round(mean(results_PRS_surv$calibration_slopes$brier),2),
                                 round(mean(results_PPS_clin_surv$calibration_slopes$brier),2),
                                 round(mean(results_PRS_clin_surv$calibration_slopes$brier),2),
                                 round(mean(results_PPS_PRS_surv$calibration_slopes$brier),2),
                                 round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier),2))
                         
)
write_csv(rf_surv_results, "~/Dropbox/Work/PPS/EU-GEI/rf_surv_results_081024.csv")

###### Generate variable importance ######
# Clinical
control_final <- trainControl(method = 'none')
check <- as.data.frame(results_clin_surv$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_clin_surv$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))
nt <- as.data.frame(results_clin_surv$best_ntree_list[1:50])
nt2 <- as.data.frame(t(nt))
best_ntree <- round(mean(nt2$V1))

predictor_vars <- setdiff(names(clin_rf_surv), c("day_exit", "Transition"))

# Dynamically create the formula for the RSF model
formula <- as.formula(paste("Surv(day_exit, Transition) ~", paste(predictor_vars, collapse = " + ")))

final_model <- rfsrc(formula,
                     data = clin_rf_surv,
                     nodesize=best_min.node.size,
                     ntree=(best_ntree),
                     mtry = repGrid$mtry)

vimp_PPS <- vimp(final_model)$importance
vimp_PPS <- as.data.frame(sort(vimp_PPS, decreasing=TRUE))
vimp_PPS <- tibble::rownames_to_column(vimp_PPS, "Predictor")
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/clin_surv_VImp.csv")

# PPS
control_final <- trainControl(method = 'none')
check <- as.data.frame(results_PPS_surv$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_PPS_surv$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))
nt <- as.data.frame(results_PPS_surv$best_ntree_list[1:50])
nt2 <- as.data.frame(t(nt))
best_ntree <- round(mean(nt2$V1))

predictor_vars <- setdiff(names(PPS_rf_surv), c("day_exit", "Transition"))

# Dynamically create the formula for the RSF model
formula <- as.formula(paste("Surv(day_exit, Transition) ~", paste(predictor_vars, collapse = " + ")))

final_model <- rfsrc(formula,
                     data = PPS_rf_surv,
                     nodesize=best_min.node.size,
                     ntree=(best_ntree),
                     mtry = repGrid$mtry)

vimp_PPS <- vimp(final_model)$importance
vimp_PPS <- as.data.frame(sort(vimp_PPS, decreasing=TRUE))
vimp_PPS <- tibble::rownames_to_column(vimp_PPS, "Predictor")
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PPS_surv_VImp.csv")

# PRS
control_final <- trainControl(method = 'none')
check <- as.data.frame(results_PPS_surv$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_PPS_surv$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))
nt <- as.data.frame(results_PPS_surv$best_ntree_list[1:50])
nt2 <- as.data.frame(t(nt))
best_ntree <- round(mean(nt2$V1))

predictor_vars <- setdiff(names(PRS_rf_surv), c("day_exit", "Transition"))

# Dynamically create the formula for the RSF model
formula <- as.formula(paste("Surv(day_exit, Transition) ~", paste(predictor_vars, collapse = " + ")))

final_model <- rfsrc(formula,
                     data = PRS_rf_surv,
                     nodesize=best_min.node.size,
                     ntree=(best_ntree),
                     mtry = repGrid$mtry)

vimp_PPS <- vimp(final_model)$importance
vimp_PPS <- as.data.frame(sort(vimp_PPS, decreasing=TRUE))
vimp_PPS <- tibble::rownames_to_column(vimp_PPS, "Predictor")
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PRS_surv_VImp.csv")

# PPS+Clinical
control_final <- trainControl(method = 'none')
check <- as.data.frame(results_PPS_clin_surv$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_PPS_clin_surv$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))
nt <- as.data.frame(results_PPS_clin_surv$best_ntree_list[1:50])
nt2 <- as.data.frame(t(nt))
best_ntree <- round(mean(nt2$V1))

predictor_vars <- setdiff(names(PPS_clin_rf_surv), c("day_exit", "Transition"))

# Dynamically create the formula for the RSF model
formula <- as.formula(paste("Surv(day_exit, Transition) ~", paste(predictor_vars, collapse = " + ")))

final_model <- rfsrc(formula,
                     data = PPS_clin_rf_surv,
                     nodesize=best_min.node.size,
                     ntree=(best_ntree),
                     mtry = repGrid$mtry)

vimp_PPS <- vimp(final_model)$importance
vimp_PPS <- as.data.frame(sort(vimp_PPS, decreasing=TRUE))
vimp_PPS <- tibble::rownames_to_column(vimp_PPS, "Predictor")
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PPS_clin_surv_VImp.csv")

# PRS+Clinical
control_final <- trainControl(method = 'none')
check <- as.data.frame(results_PRS_clin_surv$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_PRS_clin_surv$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))
nt <- as.data.frame(results_PRS_clin_surv$best_ntree_list[1:50])
nt2 <- as.data.frame(t(nt))
best_ntree <- round(mean(nt2$V1))

predictor_vars <- setdiff(names(PRS_clin_rf_surv), c("day_exit", "Transition"))

# Dynamically create the formula for the RSF model
formula <- as.formula(paste("Surv(day_exit, Transition) ~", paste(predictor_vars, collapse = " + ")))

final_model <- rfsrc(formula,
                     data = PRS_clin_rf_surv,
                     nodesize=best_min.node.size,
                     ntree=(best_ntree),
                     mtry = repGrid$mtry)

vimp_PPS <- vimp(final_model)$importance
vimp_PPS <- as.data.frame(sort(vimp_PPS, decreasing=TRUE))
vimp_PPS <- tibble::rownames_to_column(vimp_PPS, "Predictor")
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PRS_clin_surv_VImp.csv")

# PPS+PRS
control_final <- trainControl(method = 'none')
check <- as.data.frame(results_PPS_PRS_surv$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_PPS_PRS_surv$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))
nt <- as.data.frame(results_PPS_PRS_surv$best_ntree_list[1:50])
nt2 <- as.data.frame(t(nt))
best_ntree <- round(mean(nt2$V1))

predictor_vars <- setdiff(names(PPS_PRS_rf_surv), c("day_exit", "Transition"))

# Dynamically create the formula for the RSF model
formula <- as.formula(paste("Surv(day_exit, Transition) ~", paste(predictor_vars, collapse = " + ")))

final_model <- rfsrc(formula,
                     data = PPS_PRS_rf_surv,
                     nodesize=best_min.node.size,
                     ntree=(best_ntree),
                     mtry = repGrid$mtry)

vimp_PPS <- vimp(final_model)$importance
vimp_PPS <- as.data.frame(sort(vimp_PPS, decreasing=TRUE))
vimp_PPS <- tibble::rownames_to_column(vimp_PPS, "Predictor")
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_surv_VImp.csv")

# All
control_final <- trainControl(method = 'none')
check <- as.data.frame(results_PPS_PRS_clin_surv$best_mtry_list[1:50])
check2 <- as.data.frame(t(check))
repGrid <- data.frame(mtry=round(mean(check2$V1)))
ns <- as.data.frame(results_PPS_PRS_clin_surv$best_nodesize_list[1:50])
ns2 <- as.data.frame(t(ns))
best_min.node.size <- round(mean(ns2$V1))
nt <- as.data.frame(results_PPS_PRS_clin_surv$best_ntree_list[1:50])
nt2 <- as.data.frame(t(nt))
best_ntree <- round(mean(nt2$V1))

predictor_vars <- setdiff(names(PPS_PRS_clin_rf_surv), c("day_exit", "Transition"))

# Dynamically create the formula for the RSF model
formula <- as.formula(paste("Surv(day_exit, Transition) ~", paste(predictor_vars, collapse = " + ")))

final_model <- rfsrc(formula,
                     data = PPS_PRS_clin_rf_surv,
                     nodesize=best_min.node.size,
                     ntree=(best_ntree),
                     mtry = repGrid$mtry)

vimp_PPS <- vimp(final_model)$importance
vimp_PPS <- as.data.frame(sort(vimp_PPS, decreasing=TRUE))
vimp_PPS <- tibble::rownames_to_column(vimp_PPS, "Predictor")
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_clin_surv_VImp.csv")

##### Summary discrimination figure #####
results_plot <- data.frame(Type=c("LR","LR","LR",
                                  "RF","RF","RF",
                                  "Cox", "Cox","Cox","Cox","Cox","Cox","Cox",
                                  "RSF","RSF","RSF","RSF","RSF","RSF","RSF"),
                              model=c("PPS", "PRS",  "PPS+PRS", 
                                      "PPS", "PRS",  "PPS+PRS",
                                      "Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All",
                                      "Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All"),
                              C=c(0.89, 0.60, 0.90,
                                  0.86, 0.60, 0.87,
                                  0.53,0.44,0.53,0.48,0.52,0.51,0.50,
                                  0.57,0.58,0.55,0.57,0.55,0.54,0.54),
                           lCI=c(0.89,0.56,0.89,
                                 0.85,0.59,0.86,
                                 0.5,0.39,0.52,0.44,0.50,0.46,0.47,
                                 0.56,0.57,0.55,0.56,0.55,0.53,0.53),
                           uCI=c(0.90,0.62,0.91,
                                 0.87,0.62,0.88,
                                 0.57,0.49,0.55,0.52,0.54,0.57,0.52,
                                 0.58,0.59,0.56,0.58,0.56,0.55,0.55),
                           C_se=c(0.0075,0.025,0.01,
                                  0.01,0.015,0.01,
                                  0.035,0.05,0.015,0.04,0.02,0.055,0.025,
                                  0.01,0.01,0.005,0.005,0.01,0.01,0.01),
                           intercept=c(0.48,0.62,0.32,
                                       0.20,1.12,0.13,
                                       0.12,0.27,0.15,0.22,0.16,0.24,0.19,
                                       -1.18,-2.34,-1.19,-1.28,-1.1,-1.83,-1.37
                                       
                           ),
                           slope=c(1.03,0.51,1.10,
                                   0.60,0.29,0.95,
                                   2.18,-1.71,0.90,-0.95,0.95,-1.03,-0.01,
                                   0.23,-0.39,0.22,0.18,0.26,-0.13,0.12
                                   
                           )
                              
)

results_plot$Type <- factor(results_plot$Type,levels = c("LR","RF","Cox","RSF"))
results_plot$model <- factor(results_plot$model,levels = c("Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All"))

ggplot(data=results_plot, aes(x=Type, group=model)) +
  geom_rect(data=results_plot,aes(),xmin=0,xmax=9.5,ymin=0.3,ymax=0.5,
            fill="gray92") +
  geom_rect(data=results_plot,aes(),xmin=0,xmax=9.5,ymin=0.7,ymax=0.8,
            fill="gray92") +
  geom_text(aes(x=2.5, y=0.6,  label="Above chance"), stat="unique", size=12, color="gray80", family="Roboto Condensed") +
  geom_text(aes(x=2.5, y=0.45,  label="Below chance"), stat="unique", size=12, color="gray80", family="Roboto Condensed") +
  geom_text(aes(x=2.5, y=0.75,  label="Acceptable"), stat="unique", size=12, color="gray80", family="Roboto Condensed") +
  geom_text(aes(x=2.5, y=0.85,  label="Excellent"), stat="unique", size=12, color="gray80", family="Roboto Condensed") + 
  geom_text(aes(x=0.66, y=0.95,  label="0.89"), stat="unique", size=8, color="#599ec4", family="Roboto Condensed") +
  geom_text(aes(x=1.00, y=0.95,  label="0.60"), stat="unique", size=8, color="#ecc363", family="Roboto Condensed") +
  geom_text(aes(x=1.33, y=0.95,  label="0.90"), stat="unique", size=8, color="#7bae72", family="Roboto Condensed") +
  geom_text(aes(x=1.66, y=0.95,  label="0.86"), stat="unique", size=8, color="#599ec4", family="Roboto Condensed") +
  geom_text(aes(x=2.00, y=0.95,  label="0.60"), stat="unique", size=8, color="#ecc363", family="Roboto Condensed") +
  geom_text(aes(x=2.33, y=0.95,  label="0.87"), stat="unique", size=8, color="#7bae72", family="Roboto Condensed") +
  
  geom_text(aes(x=2.61, y=0.94,  label="0.53"), stat="unique", size=8, color="#c8526a", family="Roboto Condensed") +
  geom_text(aes(x=2.74, y=0.96,  label="0.44"), stat="unique", size=8, color="#599ec4", family="Roboto Condensed") +
  geom_text(aes(x=2.87, y=0.94,  label="0.53"), stat="unique", size=8, color="#ecc363", family="Roboto Condensed") +
  geom_text(aes(x=3.00, y=0.96,  label="0.48"), stat="unique", size=8, color="#8d828b", family="Roboto Condensed") +
  geom_text(aes(x=3.13, y=0.94,  label="0.52"), stat="unique", size=8, color="#d65f3d", family="Roboto Condensed") +
  geom_text(aes(x=3.26, y=0.96,  label="0.51"), stat="unique", size=8, color="#7bae72", family="Roboto Condensed") +
  geom_text(aes(x=3.39, y=0.94,  label="0.50"), stat="unique", size=8, color="#ab6758", family="Roboto Condensed") +
  
  geom_text(aes(x=3.61, y=0.94,  label="0.57"), stat="unique", size=8, color="#c8526a", family="Roboto Condensed") +
  geom_text(aes(x=3.74, y=0.96,  label="0.58"), stat="unique", size=8, color="#599ec4", family="Roboto Condensed") +
  geom_text(aes(x=3.87, y=0.94,  label="0.55"), stat="unique", size=8, color="#ecc363", family="Roboto Condensed") +
  geom_text(aes(x=4.00, y=0.96,  label="0.57"), stat="unique", size=8, color="#8d828b", family="Roboto Condensed") +
  geom_text(aes(x=4.13, y=0.94,  label="0.55"), stat="unique", size=8, color="#d65f3d", family="Roboto Condensed") +
  geom_text(aes(x=4.26, y=0.96,  label="0.54"), stat="unique", size=8, color="#7bae72", family="Roboto Condensed") +
  geom_text(aes(x=4.39, y=0.94,  label="0.54"), stat="unique", size=8, color="#ab6758", family="Roboto Condensed") +
  geom_pointrange(data=results_plot,mapping=aes(x=Type, y=C, ymin=lCI,ymax=uCI, color=model), size=2, fatten=2, position=position_dodge(width=1)) +
  geom_vline(xintercept=1.5:8.5, linetype=2,color="gray80") +
  scale_color_manual(values=c("#c8526a","#599ec4","#ecc363","#8d828b","#d65f3d","#7bae72","#ab6758")) +
  theme_classic() +
  xlab("Model") +
  ylab("C-index") +
  guides(color = guide_legend(title = "Predictors")) +
  theme(text = element_text(family="Roboto", face="bold", size=21),legend.title = element_text(size=23),legend.text = element_text(size = 23))
ggsave("~/Dropbox/Work/PPS/EU-GEI/Summary_300824_labelled.png",width = 42, height = 32, units = "cm")

##### Calibration summaries ####
lr_cal <- results_plot %>% filter(Type=="LR")
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey") +  # Line segment from (0,0) to (1,1)
  geom_abline(intercept=lr_cal$intercept[1], slope=lr_cal$slope[1], color = "#599ec4", linewidth=1.2) +
  geom_abline(intercept=lr_cal$intercept[2], slope=lr_cal$slope[2], color = "#ecc363", linewidth=1.2) +
  geom_abline(intercept=lr_cal$intercept[3], slope=lr_cal$slope[3], color = "#ab6758", linewidth=1.2) +
  labs(x = "Predicted probabilities", y = "Observed probabilities", color = "Legend") +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Calibration_LR.png",width = 15, height = 15, units = "cm")

rf_cal <- results_plot %>% filter(Type=="RF")
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey") +  # Line segment from (0,0) to (1,1)
  geom_abline(intercept=rf_cal$intercept[1], slope=rf_cal$slope[1], color = "#599ec4", linewidth=1.2) +
  geom_abline(intercept=rf_cal$intercept[2], slope=rf_cal$slope[2], color = "#ecc363", linewidth=1.2) +
  geom_abline(intercept=rf_cal$intercept[3], slope=rf_cal$slope[3], color = "#ab6758", linewidth=1.2) +
  labs(x = "Predicted probabilities", y = "Observed probabilities", color = "Legend") +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Calibration_RF.png",width = 15, height = 15, units = "cm")

cox_cal <- results_plot %>% filter(Type=="Cox")
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey") +  # Line segment from (0,0) to (1,1)
  geom_abline(intercept=cox_cal$intercept[1], slope=cox_cal$slope[1], color = "#c8526a", linewidth=1.2) +
  geom_abline(intercept=cox_cal$intercept[2], slope=cox_cal$slope[2], color = "#599ec4", linewidth=1.2) +
  geom_abline(intercept=cox_cal$intercept[3], slope=cox_cal$slope[3], color = "#ecc363", linewidth=1.2) +
  geom_abline(intercept=cox_cal$intercept[4], slope=cox_cal$slope[4], color = "#8d828b", linewidth=1.2) +
  geom_abline(intercept=cox_cal$intercept[5], slope=cox_cal$slope[5], color = "#d65f3d", linewidth=1.2) +
  geom_abline(intercept=cox_cal$intercept[6], slope=cox_cal$slope[6], color = "#7bae72", linewidth=1.2) +
  geom_abline(intercept=cox_cal$intercept[7], slope=cox_cal$slope[7], color = "#ab6758", linewidth=1.2) +
  labs(x = "Predicted probabilities", y = "Observed probabilities", color = "Legend") +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Calibration_Cox.png",width = 15, height = 15, units = "cm")

rsf_cal <- results_plot %>% filter(Type=="RSF")
ggplot() +
  geom_segment(aes(x = -3, y = -3, xend = 1, yend = 1), linetype = "dashed", color = "grey") +  # Line segment from (0,0) to (1,1)
  geom_abline(intercept=rsf_cal$intercept[1], slope=rsf_cal$slope[1], color = "#c8526a", linewidth=1.2) +
  geom_abline(intercept=rsf_cal$intercept[2], slope=rsf_cal$slope[2], color = "#599ec4", linewidth=1.2) +
  geom_abline(intercept=rsf_cal$intercept[3], slope=rsf_cal$slope[3], color = "#ecc363", linewidth=1.2) +
  geom_abline(intercept=rsf_cal$intercept[4], slope=rsf_cal$slope[4], color = "#8d828b", linewidth=1.2) +
  geom_abline(intercept=rsf_cal$intercept[5], slope=rsf_cal$slope[5], color = "#d65f3d", linewidth=1.2) +
  geom_abline(intercept=rsf_cal$intercept[6], slope=rsf_cal$slope[6], color = "#7bae72", linewidth=1.2) +
  geom_abline(intercept=rsf_cal$intercept[7], slope=rsf_cal$slope[7], color = "#ab6758", linewidth=1.2) +
  labs(x = "Predicted probabilities", y = "Observed probabilities", color = "Legend") +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Calibration_RSF.png",width = 15, height = 15, units = "cm")

##### DCA summary #####
dcaGen_PPS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/dcaGen_PPS.csv")
dcaGen_PRS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/dcaGen_PRS.csv")
dcaGen_PPSPRS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/dcaGen_PPS+PRS.csv")

dcaGen <- dcaGen_PPS %>% filter(label!="pred")
dcaGen_PPS <- dcaGen_PPS %>% filter(label=="pred")
dcaGen_PRS <- dcaGen_PRS %>% filter(label=="pred")
dcaGen_PPSPRS <- dcaGen_PPSPRS %>% filter(label=="pred")

dcaGen_PPS$variable <- "PPS"
dcaGen_PRS$variable <- "PRS"
dcaGen_PPSPRS$variable <- "PPS+PRS"
dcaGen_PPS$label <- "PPS"
dcaGen_PRS$label <- "PRS"
dcaGen_PPSPRS$label <- "PPS+PRS"

dcaGen <- rbind(dcaGen, dcaGen_PPS)
dcaGen <- rbind(dcaGen, dcaGen_PRS)
dcaGen <- rbind(dcaGen, dcaGen_PPSPRS)

dcaGen$label <- factor(dcaGen$label,levels = c("Treat All", "Treat None", "PPS", "PRS", "PPS+PRS"))
ggplot(data=dcaGen, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2) +
  coord_cartesian(ylim = c(-0.005, 0.02)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,0.10)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  scale_color_manual(values=c("gray80","#000000","#599ec4","#ecc363","#7bae72")) +
  theme(text = element_text(family="Roboto", face="bold", size=30),legend.title = element_text(size=23),legend.text = element_text(size = 23)) +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/dca_det_gen_summary.png", width=20, height=15, scale=0.5)

dcaClin_PPS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/dcaClin_rf_a.csv")
dcaClin_PRS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/dcaClin_rf_b.csv")
dcaClin_PPSPRS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/dcaClin_rf_c.csv")

dcaClin <- dcaClin_PPS %>% filter(label!="pred")
dcaClin_PPS <- dcaClin_PPS %>% filter(label=="pred")
dcaClin_PRS <- dcaClin_PRS %>% filter(label=="pred")
dcaClin_PPSPRS <- dcaClin_PPSPRS %>% filter(label=="pred")

dcaClin_PPS$variable <- "PPS"
dcaClin_PRS$variable <- "PRS"
dcaClin_PPSPRS$variable <- "PPS+PRS"
dcaClin_PPS$label <- "PPS"
dcaClin_PRS$label <- "PRS"
dcaClin_PPSPRS$label <- "PPS+PRS"

dcaClin <- rbind(dcaClin, dcaClin_PPS)
dcaClin <- rbind(dcaClin, dcaClin_PRS)
dcaClin <- rbind(dcaClin, dcaClin_PPSPRS)

dcaClin$label <- factor(dcaClin$label,levels = c("Treat All", "Treat None", "PPS", "PRS", "PPS+PRS"))
ggplot(data=dcaClin, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2) +
  coord_cartesian(ylim = c(-0.005, 0.2)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,0.35)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  scale_color_manual(values=c("gray80","#000000","#599ec4","#ecc363","#7bae72")) +
  theme(text = element_text(family="Roboto", face="bold", size=30),legend.title = element_text(size=23),legend.text = element_text(size = 23)) +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/dca_det_clin_summary.png", width=20, height=15, scale=0.5)

#### Remission ####
PPS_scored[,c(3:8,10:14,16:17)] <- lapply(PPS_scored[,c(3:8,10:14,16:17)], factor)

results_rem <- data.frame(model=c("Clinical","PPS","PRS","Clinical+PPS","Clinical+PRS", "PPS+PRS", "All"),
                          C_index=c(1:7),
                          balanced_accuracy=c(1:7),
                          sensitivity=c(1:7),
                          specificity=c(1:7),
                          ppv=c(1:7),
                          npv=c(1:7),
                          calibration_intercept=c(1:7),
                          calibration_slope=c(1:7),
                          Brier=c(1:7),
                          ClinicalUtilityClin=c(1:7)
)
predictors_cox <- list(a=c("CAARMS", "gafex01", "gafex02"),
                       b=c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia"),
                       c=c("ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),
                       d=c("CAARMS", "gafex01", "gafex02","Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia"),
                       e=c("CAARMS", "gafex01", "gafex02","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),
                       f=c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),
                       g=c("CAARMS", "gafex01", "gafex02","Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"))

calibration_stats <- function(observed, predicted) {
  mean_observed <- mean(observed)
  mean_predicted <- mean(predicted)
  slope <- cov(predicted, observed) / var(predicted)
  intercept <- mean_observed - slope * mean_predicted
  return(list(intercept = intercept, slope = slope))
}

PPS_scored[,c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
              "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia")] <- lapply(PPS_scored[,c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                                                                                                                                  "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia")], FUN=factor)

###### Fit model and internal validation #######
folds <- repeatfolds(y=factor(PPS_scored$ARMS_2yr), repeats = 10, n_outer_folds = 10)
foldids = rep(1,length(folds)*length(folds$Rep1))

for (i in 1:7){
  set.seed(123)
  temp <- PPS_scored %>% select(all_of((predictors_cox[[i]])))
  temp_mat <- model.matrix(~.-1,temp)
  
  PRS_multi <- nestcv.glmnet(y=factor(PPS_scored$ARMS_2yr), x=temp_mat,family="binomial",
                             n_outer_folds = 10, balance = "randomsample", alpha=1, seed=123) |>
    repeatcv(10, repeat_folds = folds, rep.cores = 2)
  
  PPS_ct <- matrix(0,ncol = 1, nrow = 1)
  PPS_ct$FP <- sum(PRS_multi[["output"]]$predy==1 & PRS_multi[["output"]]$testy==0)
  PPS_ct$TN <- sum(PRS_multi[["output"]]$predy==0 & PRS_multi[["output"]]$testy==0) 
  PPS_ct$TP <- sum(PRS_multi[["output"]]$predy==1 & PRS_multi[["output"]]$testy==1) 
  PPS_ct$FN <- sum(PRS_multi[["output"]]$predy==0 & PRS_multi[["output"]]$testy==1)
  
  PPS_ct$total <- sum(PPS_ct$FP, PPS_ct$TN, PPS_ct$TP, PPS_ct$FN)
  
  results_rem$C_index[i] <- paste(round(mean(PRS_multi$result[1:10])*100,1),"% (",
                              round(mean(PRS_multi$result[1:10]*100,1)-1.96*sd(PRS_multi$result[1:10])*100,1),"%-",
                              round(mean(PRS_multi$result[1:10]*100,1)+1.96*sd(PRS_multi$result[1:10])*100,1), "%)", sep="")
  results_rem$balanced_accuracy[i] <- paste(round(sum((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN)),
                                                  (PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP)))*50,1),"% (",
                                        round(sum((propCI(x = ((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN))*PPS_ct$total), 
                                                          n = PPS_ct$total)$lower[2]),
                                                  (propCI(x = ((PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP))*PPS_ct$total), 
                                                          n = PPS_ct$total)$lower[2]))*50,1),"%-",
                                        round(sum((propCI(x = ((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN))*PPS_ct$total), 
                                                          n = PPS_ct$total)$upper[2]),
                                                  (propCI(x = ((PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP))*PPS_ct$total), 
                                                          n = PPS_ct$total)$upper[2]))*50,1), "%)", sep="")
  results_rem$sensitivity[i] <- paste(round((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN))*100,1),"% (",
                                  round(propCI(x = ((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN))*PPS_ct$total), 
                                               n = PPS_ct$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((PPS_ct$TP/(PPS_ct$TP+PPS_ct$FN))*PPS_ct$total), 
                                               n = PPS_ct$total)$upper[2]*100,1), "%)", sep="")
  results_rem$specificity[i] <- paste(round((PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP))*100,1),"% (",
                                  round(propCI(x = ((PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP))*PPS_ct$total), 
                                               n = PPS_ct$total)$lower[2]*100,1),"%-",                                   
                                  round(propCI(x = ((PPS_ct$TN/(PPS_ct$TN+PPS_ct$FP))*PPS_ct$total), 
                                               n = PPS_ct$total)$upper[2]*100,1), "%)", sep="")
  
  PRS_multi[["output"]]$predy <- as.factor(PRS_multi[["output"]]$predy)
  results_rem$ppv[i] <- paste0(round(yardstick::ppv(data = PRS_multi[["output"]], truth = testy, estimate = predy)$.estimate*100,1),"%")
  results_rem$npv[i] <- paste0(round(yardstick::npv(data = PRS_multi[["output"]], truth = testy, estimate = predy)$.estimate*100,1),"%")
  
  # Calculate calibration slope on the hold-out data (test set)
  logistic_calibration <- calibration_stats(observed = as.numeric(PRS_multi[["output"]]$testy)-1, predicted = as.numeric(PRS_multi[["output"]]$predy))
  logistic_calibration2 <- predRupdate::pred_val_probs(binary_outcome = as.numeric(PRS_multi[["output"]]$testy)-1, Prob = rescale(PRS_multi[["output"]]$predyp))
  
  results_rem$calibration_intercept[i] <- round(logistic_calibration$intercept[1],2)
  results_rem$calibration_slope[i] <- round(logistic_calibration$slope[1],2)
  
  Brier <- data.frame(obs=as.numeric(PRS_multi[["output"]]$testy)-1,
                      pred=scales::rescale((PRS_multi[["output"]]$predyp), to=c(0,1)))
  Brier$Brier <- (Brier$obs - Brier$pred)^2
  results_rem$Brier[i] <- mean(Brier$Brier)
  
  # plot_data[,i+3] <- logistic_calibration$slope[1] * plot_data[,i] + logistic_calibration$intercept[1]
  # plot_pred_norm[,i] <- (PRS_multi[["output"]]$predyp-min(PRS_multi[["output"]]$predyp))/(max(PRS_multi[["output"]]$predyp)-min(PRS_multi[["output"]]$predyp))
  
  # Decision curve analysis
  library(dcurves)
  dca_temp <- PRS_multi[["output"]] %>% filter(rep==1)
  dca <- data.frame(obs=as.numeric(dca_temp$testy)-1,
                    pred=as.numeric(dca_temp$predy)-1)
  
  dca_assessment <- dca(obs ~ pred,
                        data = dca,
                        prevalence = 0.192,
                        thresholds = seq(0, 0.5, 0.01)
  ) %>%
    as_tibble() 
  results_rem$ClinicalUtilityClin[i] <- dca_assessment %>% 
    filter(label=="pred" & net_benefit>0) %>%
    summarise(net_benefit = mean(net_benefit), .groups = "drop")
  
  # plot cross validated net benefit values
  ggplot(data=dca_assessment, aes(x = threshold, y = net_benefit, color = label)) +
    stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2) +
    coord_cartesian(ylim = c(-0.005, 0.2)) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,0.35)) +
    labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
    theme_classic()
  ggsave(paste0("~/Dropbox/Work/PPS/EU-GEI/Plots/dca_rem_", results_rem[i,1],".png"), width=6, height=8, dpi=300)
  write_csv(dca_assessment, paste0("~/Dropbox/Work/PPS/EU-GEI/dcaRem_", results_rem[i,1],".csv"))
  
  
}

write_csv(results_rem, "~/Dropbox/Work/PPS/EU-GEI/results_remission.csv")

PRS_multi[["output"]]$predyp_norm <- (PRS_multi[["output"]]$predyp-min(PRS_multi[["output"]]$predyp))/(max(PRS_multi[["output"]]$predyp)-min(PRS_multi[["output"]]$predyp))
plot_data <- merge(plot_data, plot_pred_norm)

ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed", color = "grey") +  # Line segment from (0,0) to (1,1)
  geom_smooth(data=plot_data, aes(x1,y1,color = "PPS"), method = "lm", se = FALSE) +
  geom_smooth(data=plot_data, aes(x2,y2,color = "PRS"), method = "lm", se = FALSE) +
  geom_smooth(data=plot_data, aes(x3,y3,color = "PPS+PRS"), method = "lm", se = FALSE) +
  geom_point(aes(plot_pred_norm$z1, as.numeric(PRS_multi[["output"]]$testy)-1), alpha = 0.5) +
  labs(x = "Predicted probabilities", y = "Observed probabilities", color = "Legend") +
  scale_color_manual(values = c("PPS" = "#c8526a", "PRS" = "#758EAB", "PPS+PRS" = "#7bae72")) +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Calibration_remission.png",width = 24, height = 18, units = "cm")

calibration <- data.frame(observed=as.numeric(PRS_multi[["output"]]$testy)-1,
                          predicted = as.numeric(PRS_multi[["output"]]$predyp_norm))

calibration %>% cal_plot_logistic(observed, predicted)


library(rms)
# Calibration plot for logistic regression
calibration_plot_logistic <- calibrate(PRS_multi, method = "boot", B = 100)
plot(calibration_plot_logistic, main = "Calibration Plot for Logistic Regression")

# Calibration plot for Cox model
calibration_plot_cox <- calibrate(fit_cox, method = "boot", B = 100)
plot(calibration_plot_cox, main = "Calibration Plot for Cox Model")

PRS_multi <- glm(chr ~ ZRE_SCZ_imp_0.1,family=binomial(link='logit'),data=PPS_scored)
summary(PRS_multi)
fmsb::NagelkerkeR2(PRS_multi)

PRS_multi <- glm(chr ~ PPS + ZRE_SCZ_imp_0.1,family=binomial(link='logit'),data=PPS_scored)
summary(PRS_multi)
fmsb::NagelkerkeR2(PRS_multi)

PPS_ind_multi <- glm(chr ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity + 
                       Immigration  + Paternal_Age + Paternal_SES + Parental_SMI + 
                       Adult_Life_Events + Tobacco + Cannabis + Childhood_Trauma + 
                       Anhedonia,family=binomial(link='logit'),data=PPS_scored)
summary(PPS_ind_multi)
fmsb::NagelkerkeR2(PPS_ind_multi)

PPS_ind_PRS_multi <- glm(chr ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity + 
                           Immigration  + Paternal_Age + Paternal_SES + Parental_SMI + 
                           Adult_Life_Events + Tobacco + Cannabis + Childhood_Trauma + 
                           Anhedonia + ZRE_SCZ_imp_0.1,family=binomial(link='logit'),data=PPS_scored)
summary(PPS_ind_PRS_multi)
fmsb::NagelkerkeR2(PPS_ind_PRS_multi)

###### Generate coefficients #####

temp <- PPS_scored %>% select(all_of((predictors_cox[[1]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[2]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
PPS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[3]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[4]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_PPS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[5]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[6]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
PPS_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp <- PPS_scored %>% select(all_of((predictors_cox[[7]])))
temp_mat <- model.matrix(~.-1,temp)
check <- cv.glmnet(y=Surv(PPS_scored$day_exit, PPS_scored$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_PPS_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

write_csv(clin_coef,"~/Dropbox/Work/PPS/EU-GEI/clin_coef_surv_remission.csv")
write_csv(PPS_coef,"~/Dropbox/Work/PPS/EU-GEI/PPS_coef_surv_remission.csv")
write_csv(PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/PRS_coef_surv_remission.csv")
write_csv(clin_PPS_coef,"~/Dropbox/Work/PPS/EU-GEI/clinPPS_coef_surv_remission.csv")
write_csv(clin_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/clin_PRS_coef_surv_remission.csv")
write_csv(PPS_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_coef_surv_remission.csv")
write_csv(clin_PPS_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/clin_PPS_PRS_coef_surv_remission.csv")

##### RF #####
# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

# Wrangle data
PPS_scored[,c(3:8,10:14,16:17)] <- lapply(PPS_scored[,c(3:8,10:14,16:17)], factor)
PPS_scored_mat <- as.data.frame(model.matrix(~.-1,PPS_scored[,c(3:8,10:14,16:17,23:39)]))
PPS_hc_rf <- as.data.frame(PPS_scored_mat)
PPS_hc_rf <- PPS_hc_rf %>% subset(select = c(-ZRE_SCZ_imp_0.1,-ZRE_BD_imp_EUR_0.1:-ZRE_MDD_imp_EUR_0.1))
PPS_hc_rf$ARMS_2yr <- PPS_scored$ARMS_2yr

PPS_hc_rf[,1:19] <- lapply(PPS_hc_rf[,1:19], factor)
PPS_hc_rf$ARMS_2yr <- factor(PPS_hc_rf$ARMS_2yrX1)

PPS_hc_rf$ARMS_2yr <- factor(make.names(levels(PPS_hc_rf$ARMS_2yr))[PPS_hc_rf$ARMS_2yr])

PPS_hc_rf <- PPS_hc_rf  %>% dplyr::rename(Ethnicity35 = Ethnicity3.5,
                                          Immigration25 = Immigration2.5,
                                          Parental_SMI55 = Parental_SMI5.5,
                                          Adult_Life_Events55 = Adult_Life_Events5.5,
                                          Anhedonia65 = Anhedonia6.5)
make.names(colnames(PPS_hc_rf), unique=TRUE)

## Compute F1-score

F1_score <- function(mat, algoName){
  
  # Remark: left column = prediction // top = real values
  recall <- matrix(1:nrow(mat), ncol = nrow(mat))
  precision <- matrix(1:nrow(mat), ncol = nrow(mat))
  F1_score <- matrix(1:nrow(mat), ncol = nrow(mat))
  
  
  for(i in 1:nrow(mat)){
    recall[i] <- mat[i,i]/rowSums(mat)[i]
    precision[i] <- mat[i,i]/colSums(mat)[i]
  }
  
  for(i in 1:ncol(recall)){
    F1_score[i] <- 2 * ( precision[i] * recall[i] ) / ( precision[i] + recall[i])
  }
  
  # We display the matrix labels
  colnames(F1_score) <- colnames(mat)
  rownames(F1_score) <- algoName
  
  # Display the F1_score for each class
  F1_score
  
  # Display the average F1_score
  mean(F1_score[1,])
}

# Create a function to perform nested cross-validation with repeats
nested_cv_with_repeats_rem <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, tuneGrid, seed) {
  
  set.seed(seed)
  seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
  combined_df$ARMS_2yr <- factor(combined_df$ARMS_2yr)
  combined_df$ARMS_2yr <- factor(make.names(levels(combined_df$ARMS_2yr))[combined_df$ARMS_2yr])

  best_inner_result_list <- list()
  best_mtry_list <- list()
  best_ntree_list <- list()
  best_nodesize_list <- list()
  all_inner_results <- list()
  c_stat_nested <- data.frame()
  calibration_slopes <- data.frame()
  
  # Outer loop for cross-validation
  for (outer_rep in 1:outerRepeats) {
    cat("Outer Repeat:", outer_rep, "\n")
    set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
    outer_folds <- createFolds(combined_df$ARMS_2yr, k = outerFolds, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(outer_folds)) {
      cat("  Outer Fold:", i, "\n")
      
      train_indices <- outer_folds[[i]]
      train <- combined_df[train_indices, ] #  This splits the data into train (inner loop)
      test <- combined_df[-train_indices, ] #  This splits the data into test (outer loop)
      
      inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
      
      # Inner Cross-Validation and Model Training
      for (inner_rep in 1:innerRepeats) {
        cat("    Inner Repeat:", inner_rep, "\n")
        
        # Use a different seed for each inner repeat
        set.seed(seeds[outer_rep] + inner_rep)
        
        # Tune both alpha and lambda
        for (mtry_value in tuneGrid) {
          library(randomForest)
          library(MLmetrics)
          library(caret)
          
          # Define a custom summary function to calculate F1 score
          customSummary <- function(data, lev = NULL, model = NULL) {
            f1 <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[1])
            out <- c(F1 = f1)
            out
          }
          
          # Set up the trainControl with the custom summary function
          control <- trainControl(method = 'cv', 
                                  number = 5, 
                                  #repeats = 10,
                                  classProbs = TRUE,
                                  summaryFunction = customSummary,
                                  search = 'random')
          
          tree <- c(50, 100, 250, 500)
          n.tree <- sample(tree,1)
          nodeSize <- seq(1,(nrow(train)/10), by=1)
          node.size <- sample(nodeSize,1)
          tune_grid_temp <- data.frame(mtry=c(NA,NA,NA))
          tune_grid_temp$mtry <- sample(tuneGrid$mtry,3)
          
          # Train the model using the custom F1 metric
          inner_model <- train(ARMS_2yr ~ .,
                               data = train,
                               method = "rf",
                               metric = "F1",
                               tuneGrid = tune_grid_temp,
                               tuneLength=10,
                               ntree = n.tree,
                               nodesize=node.size,
                               trControl = control)
          
          print(inner_model)
          
          if (all(!is.na(inner_model$results$F1))){
            # Store the F1 for each inner repeat
            inner_results[[inner_rep]] <- rbind(inner_results[[inner_rep]], 
                                                data.frame(min.node.size = node.size,
                                                           tree = n.tree,
                                                           mtry = inner_model$results$mtry, 
                                                           F1 = inner_model$results$F1, 
                                                           repeat_number = inner_rep))
          }
        }
      }
      
      # Combine results from all repeats
      all_inner_results <- do.call(rbind, inner_results)
      print(dim(all_inner_results))
      
      # Calculate the average F1 for each hyperparameter combination
      avg_results <- all_inner_results %>%
        group_by(mtry) %>%
        summarise(avg_F1 = mean(F1, na.rm=TRUE),
                  mtry = mean(mtry, na.rm=TRUE),
                  ntree = mean(tree, na.rm=TRUE),
                  nodesize=mean(min.node.size, na.rm=TRUE),
                  .groups = "keep") 
      
      # Determine the best hyperparameter combination based on the highest average F1
      best_inner_result <- avg_results[which.max(avg_results$avg_F1), ]
      best_min.node.size <- best_inner_result$nodesize
      best_mtry <- best_inner_result$mtry
      best_F1 <- best_inner_result$avg_F1
      best_F1_SE <- best_inner_result$SE_F1
      best_ntree <- best_inner_result$ntree
      
      # Calculate the correct index for the current outer fold and repeat
      index <- (outer_rep - 1) * outerFolds + i
      
      # Store the best hyperparameters for the current outer fold
      best_inner_result_list[[index]] <- best_inner_result
      best_mtry_list[[index]]  <- best_mtry
      best_ntree_list[[index]]  <- best_ntree
      best_nodesize_list[[index]]  <- best_min.node.size
      
      # Train the final model using the best hyperparameters on the inner training set (x_inner_train, y_inner_train)
      control_final <- trainControl(method = 'none')
      repGrid <- data.frame(mtry=best_mtry)  
      
      final_model <- train(ARMS_2yr ~ .,
                           data = train,
                           method = "rf",
                           metric = "F1",
                           ntree = best_ntree,
                           nodesize=best_min.node.size,
                           trControl=control_final,
                           tuneGrid = repGrid)
      
      # Predict the linear predictors (PI) from the Elastic Net model
      test$PI <- predict(final_model, newdata = test, type = "prob")[,2]
      test$chr_pred <- predict(final_model, newdata = test, type = "raw")
      
      cm <- confusionMatrix(data = test$chr_pred, reference = test$ARMS_2yr)
      test <- test %>% mutate(PI=case_when(PI==0 ~ 0.001, PI==1 ~ 0.999, TRUE ~ PI))
      
      # Fit a Cox model on the test set using penalized coefficients
      model_test <- glm(ARMS_2yr ~ PI, data = test, family="binomial")
      
      c_stat_nested <- rbind(c_stat_nested, data.frame(
        C_test = concordance(model_test)$concordance,
        SE_test = concordance(model_test)$cvar,
        Fold = i,
        OuterRepeat = outer_rep,
        n_train = nrow(train),
        events_train = sum(train$ARMS_2yr=="X1"),
        n_test = nrow(test),
        events_test = sum(test$ARMS_2yr=="X1"),
        balanced_accuracy = cm$byClass[11],
        sensitivity = cm$byClass[1],
        specificity = cm$byClass[2],
        ppv = cm$byClass[3],
        npv = cm$byClass[4],
        precision = cm$byClass[5],
        recall = cm$byClass[6],
        f1 = cm$byClass[7]
        
      ))
      
      calibration <- rms::val.prob(p=test$PI, y=as.numeric(test$ARMS_2yr)-1, m=200, pl=F)
      
      # Calculate calibration slope on the hold-out data (test set)
      
      calibration_intercept <- unname(calibration[12])
      calibration_slope <- unname(calibration[13])
      brier <- unname(calibration[11])
      
      
      # Store calibration results
      calibration_slopes <- rbind(calibration_slopes, data.frame(
        fold = i,
        OuterRepeat = outer_rep,
        intercept = calibration_intercept,
        slope = calibration_slope,
        brier = brier
      ))
    }
  }
  
  list(
    best_inner_result_list = best_inner_result_list,
    best_mtry_list = best_mtry_list,
    best_ntree_list = best_ntree_list,
    best_nodesize_list = best_nodesize_list,
    # C_results = C_results,
    # C_SE_results = C_SE_results,
    c_stat_nested = c_stat_nested,
    calibration_slopes = calibration_slopes
  )
}

# Define tuning grid

tune_grid <- expand.grid(
  mtry=c(3:30)
)

###### Run model fitting and internal validation ######
predictors_cox <- list(a=c("CAARMS", "gafex01", "gafex02","ARMS_2yr"),
                       b=c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ARMS_2yr"),
                       c=c("ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","ARMS_2yr"),
                       d=c("CAARMS", "gafex01", "gafex02","Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ARMS_2yr"),
                       e=c("CAARMS", "gafex01", "gafex02","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","ARMS_2yr"),
                       f=c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","ARMS_2yr"),
                       g=c("CAARMS", "gafex01", "gafex02","Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                           "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia","ZRE_SCZ_imp_0.1","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","ARMS_2yr"))

PPS_scored[,c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
              "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia")] <- lapply(PPS_scored[,c("Gender","Handedness","Urbanicity","Pollution","Ethnicity","Immigration","Paternal_SES",
                                                                                                                                  "Parental_SMI", "Adult_Life_Events","Tobacco", "Cannabis", "Childhood_Trauma","Anhedonia")], FUN=factor)

# Clinical
clin_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[1]])))
clin_rf_surv <- as.data.frame(model.matrix(~.-1,clin_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3,3,3)
)

results_clin_surv <- nested_cv_with_repeats_rem(
  combined_df = clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PPS

PPS_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[2]])))
PPS_rf_surv <- as.data.frame(model.matrix(~.-1,PPS_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_rf_surv)-2)
)

results_PPS_surv <- nested_cv_with_repeats_rem(
  combined_df = PPS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PRS
PRS_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[3]])))
PRS_rf_surv <- as.data.frame(model.matrix(~.-1,PRS_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PRS_rf_surv)-2)
)

results_PRS_surv <- nested_cv_with_repeats_rem(
  combined_df = PRS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PPS

PPS_clin_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[4]])))
PPS_clin_rf_surv <- as.data.frame(model.matrix(~.-1,PPS_clin_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_clin_rf_surv)-2)
)

results_PPS_clin_surv <- nested_cv_with_repeats_rem(
  combined_df = PPS_clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PRS

PRS_clin_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[5]])))
PRS_clin_rf_surv <- as.data.frame(model.matrix(~.-1,PRS_clin_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_rf_surv)-2)
)

results_PRS_clin_surv <- nested_cv_with_repeats_rem(
  combined_df = PRS_clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PPS + PRS

PPS_PRS_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[6]])))
PPS_PRS_rf_surv <- as.data.frame(model.matrix(~.-1,PPS_PRS_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_PRS_rf_surv)-2)
)

results_PPS_PRS_surv <- nested_cv_with_repeats_rem(
  combined_df = PPS_PRS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PPS + PRS 

PPS_PRS_clin_rf_surv <- PPS_scored %>% select(all_of((predictors_cox[[7]])))
PPS_PRS_clin_rf_surv <- as.data.frame(model.matrix(~.-1,PPS_PRS_clin_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_PRS_clin_rf_surv)-2)
)

results_PPS_PRS_clin_surv <- nested_cv_with_repeats_rem(
  combined_df = PPS_PRS_clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)
###### Generate performance metrics summary table######

rf_rem_results <- data.frame(model=c("Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All"),
                              C=c(paste0(round(mean(results_clin_surv$c_stat_nested$C),2), " (",
                                         round(mean(results_clin_surv$c_stat_nested$C) - 1.96*(sd(results_clin_surv$c_stat_nested$C) / sqrt(nrow(results_clin_surv$c_stat_nested))),2),
                                         "-",
                                         round(mean(results_clin_surv$c_stat_nested$C) + 1.96*(sd(results_clin_surv$c_stat_nested$C) / sqrt(nrow(results_clin_surv$c_stat_nested))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_surv$c_stat_nested$C),2), " (",
                                         round(mean(results_PPS_surv$c_stat_nested$C) - 1.96*(sd(results_PPS_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_surv$c_stat_nested))),2),
                                         "-",
                                         round(mean(results_PPS_surv$c_stat_nested$C) + 1.96*(sd(results_PPS_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_surv$c_stat_nested))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_surv$c_stat_nested$C),2), " (",
                                         round(mean(results_PRS_surv$c_stat_nested$C) - 1.96*(sd(results_PRS_surv$c_stat_nested$C) / sqrt(nrow(results_PRS_surv$c_stat_nested))),2),
                                         "-",
                                         round(mean(results_PRS_surv$c_stat_nested$C) + 1.96*(sd(results_PRS_surv$c_stat_nested$C) / sqrt(nrow(results_PRS_surv$c_stat_nested))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_clin_surv$c_stat_nested$C),2), " (",
                                         round(mean(results_PPS_clin_surv$c_stat_nested$C) - 1.96*(sd(results_PPS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_clin_surv$c_stat_nested))),2),
                                         "-",
                                         round(mean(results_PPS_clin_surv$c_stat_nested$C) + 1.96*(sd(results_PPS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_clin_surv$c_stat_nested))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_clin_surv$c_stat_nested$C),2), " (",
                                         round(mean(results_PRS_clin_surv$c_stat_nested$C) - 1.96*(sd(results_PRS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PRS_clin_surv$c_stat_nested))),2),
                                         "-",
                                         round(mean(results_PRS_clin_surv$c_stat_nested$C) + 1.96*(sd(results_PRS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PRS_clin_surv$c_stat_nested))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_surv$c_stat_nested$C),2), " (",
                                         round(mean(results_PPS_PRS_surv$c_stat_nested$C) - 1.96*(sd(results_PPS_PRS_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_surv$c_stat_nested))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_surv$c_stat_nested$C) + 1.96*(sd(results_PPS_PRS_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_surv$c_stat_nested))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_clin_surv$c_stat_nested$C),2), " (",
                                         round(mean(results_PPS_PRS_clin_surv$c_stat_nested$C) - 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_clin_surv$c_stat_nested$C) + 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))),2),
                                         ")")),
                              BAC=c(paste0(round(mean(results_clin_surv$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                           round((mean(results_clin_surv$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_clin_surv$results_clin_surv$balanced_accuracy) / sqrt(nrow(results_clin_surv$c_stat_nested))))*100,1),
                                           "%-",
                                           round((mean(results_clin_surv$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_clin_surv$results_clin_surv$balanced_accuracy) / sqrt(nrow(results_clin_surv$c_stat_nested))))*100,1),
                                           "%)"),
                                    paste0(round(mean(results_PPS_surv$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                           round((mean(results_PPS_surv$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_PPS_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_surv$c_stat_nested))))*100,1),
                                           "%-",
                                           round((mean(results_PPS_surv$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_PPS_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_surv$c_stat_nested))))*100,1),
                                           "%)"),
                                    paste0(round(mean(results_PRS_surv$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                           round((mean(results_PRS_surv$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_PRS_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PRS_surv$c_stat_nested))))*100,1),
                                           "%-",
                                           round((mean(results_PRS_surv$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_PRS_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PRS_surv$c_stat_nested))))*100,1),
                                           "%)"),
                                    paste0(round(mean(results_PPS_clin_surv$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                           round((mean(results_PPS_clin_surv$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_PPS_clin_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_clin_surv$c_stat_nested))))*100,1),
                                           "%-",
                                           round((mean(results_PPS_clin_surv$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_PPS_clin_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_clin_surv$c_stat_nested))))*100,1),
                                           "%)"),
                                    paste0(round(mean(results_PRS_clin_surv$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                           round((mean(results_PRS_clin_surv$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_PRS_clin_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PRS_clin_surv$c_stat_nested))))*100,1),
                                           "%-",
                                           round((mean(results_PRS_clin_surv$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_PRS_clin_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PRS_clin_surv$c_stat_nested))))*100,1),
                                           "%)"),
                                    paste0(round(mean(results_PPS_PRS_surv$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                           round((mean(results_PPS_PRS_surv$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_PPS_PRS_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_PRS_surv$c_stat_nested))))*100,1),
                                           "%-",
                                           round((mean(results_PPS_PRS_surv$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_PPS_PRS_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_PRS_surv$c_stat_nested))))*100,1),
                                           "%)"),
                                    paste0(round(mean(results_PPS_PRS_clin_surv$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                           round((mean(results_PPS_PRS_clin_surv$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))))*100,1),
                                           "%-",
                                           round((mean(results_PPS_PRS_clin_surv$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))))*100,1),
                                           "%)")),
                              sensitivity=c(paste0(round(mean(results_clin_surv$c_stat_nested$sensitivity)*100,1), "% (",
                                                   round((mean(results_clin_surv$c_stat_nested$sensitivity) - 1.96*(sd(results_clin_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_clin_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_clin_surv$c_stat_nested$sensitivity) + 1.96*(sd(results_clin_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_clin_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PPS_surv$c_stat_nested$sensitivity)*100,1), "% (",
                                                   round((mean(results_PPS_surv$c_stat_nested$sensitivity) - 1.96*(sd(results_PPS_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PPS_surv$c_stat_nested$sensitivity) + 1.96*(sd(results_PPS_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PRS_surv$c_stat_nested$sensitivity)*100,1), "% (",
                                                   round((mean(results_PRS_surv$c_stat_nested$sensitivity) - 1.96*(sd(results_PRS_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PRS_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PRS_surv$c_stat_nested$sensitivity) + 1.96*(sd(results_PRS_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PRS_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PPS_clin_surv$c_stat_nested$sensitivity)*100,1), "% (",
                                                   round((mean(results_PPS_clin_surv$c_stat_nested$sensitivity) - 1.96*(sd(results_PPS_clin_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_clin_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PPS_clin_surv$c_stat_nested$sensitivity) + 1.96*(sd(results_PPS_clin_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_clin_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PRS_clin_surv$c_stat_nested$sensitivity)*100,1), "% (",
                                                   round((mean(results_PRS_clin_surv$c_stat_nested$sensitivity) - 1.96*(sd(results_PRS_clin_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PRS_clin_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PRS_clin_surv$c_stat_nested$sensitivity) + 1.96*(sd(results_PRS_clin_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PRS_clin_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PPS_PRS_surv$c_stat_nested$sensitivity)*100,1), "% (",
                                                   round((mean(results_PPS_PRS_surv$c_stat_nested$sensitivity) - 1.96*(sd(results_PPS_PRS_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_PRS_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PPS_PRS_surv$c_stat_nested$sensitivity) + 1.96*(sd(results_PPS_PRS_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_PRS_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PPS_PRS_clin_surv$c_stat_nested$sensitivity)*100,1), "% (",
                                                   round((mean(results_PPS_PRS_clin_surv$c_stat_nested$sensitivity) - 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PPS_PRS_clin_surv$c_stat_nested$sensitivity) + 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))))*100,1),
                                                   "%)")),
                              specificity=c(paste0(round(mean(results_clin_surv$c_stat_nested$specificity)*100,1), "% (",
                                                   round((mean(results_clin_surv$c_stat_nested$specificity) - 1.96*(sd(results_clin_surv$c_stat_nested$specificity) / sqrt(nrow(results_clin_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_clin_surv$c_stat_nested$specificity) + 1.96*(sd(results_clin_surv$c_stat_nested$specificity) / sqrt(nrow(results_clin_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PPS_surv$c_stat_nested$specificity)*100,1), "% (",
                                                   round((mean(results_PPS_surv$c_stat_nested$specificity) - 1.96*(sd(results_PPS_surv$c_stat_nested$specificity) / sqrt(nrow(results_PPS_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PPS_surv$c_stat_nested$specificity) + 1.96*(sd(results_PPS_surv$c_stat_nested$specificity) / sqrt(nrow(results_PPS_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PRS_surv$c_stat_nested$specificity)*100,1), "% (",
                                                   round((mean(results_PRS_surv$c_stat_nested$specificity) - 1.96*(sd(results_PRS_surv$c_stat_nested$specificity) / sqrt(nrow(results_PRS_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PRS_surv$c_stat_nested$specificity) + 1.96*(sd(results_PRS_surv$c_stat_nested$specificity) / sqrt(nrow(results_PRS_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PPS_clin_surv$c_stat_nested$specificity)*100,1), "% (",
                                                   round((mean(results_PPS_clin_surv$c_stat_nested$specificity) - 1.96*(sd(results_PPS_clin_surv$c_stat_nested$specificity) / sqrt(nrow(results_PPS_clin_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PPS_clin_surv$c_stat_nested$specificity) + 1.96*(sd(results_PPS_clin_surv$c_stat_nested$specificity) / sqrt(nrow(results_PPS_clin_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PRS_clin_surv$c_stat_nested$specificity)*100,1), "% (",
                                                   round((mean(results_PRS_clin_surv$c_stat_nested$specificity) - 1.96*(sd(results_PRS_clin_surv$c_stat_nested$specificity) / sqrt(nrow(results_PRS_clin_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PRS_clin_surv$c_stat_nested$specificity) + 1.96*(sd(results_PRS_clin_surv$c_stat_nested$specificity) / sqrt(nrow(results_PRS_clin_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PPS_PRS_surv$c_stat_nested$specificity)*100,1), "% (",
                                                   round((mean(results_PPS_PRS_surv$c_stat_nested$specificity) - 1.96*(sd(results_PPS_PRS_surv$c_stat_nested$specificity) / sqrt(nrow(results_PPS_PRS_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PPS_PRS_surv$c_stat_nested$specificity) + 1.96*(sd(results_PPS_PRS_surv$c_stat_nested$specificity) / sqrt(nrow(results_PPS_PRS_surv$c_stat_nested))))*100,1),
                                                   "%)"),
                                            paste0(round(mean(results_PPS_PRS_clin_surv$c_stat_nested$specificity)*100,1), "% (",
                                                   round((mean(results_PPS_PRS_clin_surv$c_stat_nested$specificity) - 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$specificity) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))))*100,1),
                                                   "%-",
                                                   round((mean(results_PPS_PRS_clin_surv$c_stat_nested$specificity) + 1.96*(sd(results_PPS_PRS_clin_surv$c_stat_nested$specificity) / sqrt(nrow(results_PPS_PRS_clin_surv$c_stat_nested))))*100,1),
                                                   "%)")),
                              ppv=c(paste0(round(mean(results_clin_surv$c_stat_nested$ppv*100),1),"%"),
                                    paste0(round(mean(results_PPS_surv$c_stat_nested$ppv*100),1),"%"),
                                    paste0(round(mean(results_PRS_surv$c_stat_nested$ppv*100),1),"%"),
                                    paste0(round(mean(results_PPS_clin_surv$c_stat_nested$ppv*100),1),"%"),
                                    paste0(round(mean(results_PRS_clin_surv$c_stat_nested$ppv*100),1),"%"),
                                    paste0(round(mean(results_PPS_PRS_surv$c_stat_nested$ppv*100),1),"%"),
                                    paste0(round(mean(results_PPS_PRS_clin_surv$c_stat_nested$ppv*100),1),"%")),
                              npv=c(paste0(round(mean(results_clin_surv$c_stat_nested$npv*100),1),"%"),
                                    paste0(round(mean(results_PPS_surv$c_stat_nested$npv*100),1),"%"),
                                    paste0(round(mean(results_PRS_surv$c_stat_nested$npv*100),1),"%"),
                                    paste0(round(mean(results_PPS_clin_surv$c_stat_nested$npv*100),1),"%"),
                                    paste0(round(mean(results_PRS_clin_surv$c_stat_nested$npv*100),1),"%"),
                                    paste0(round(mean(results_PPS_PRS_surv$c_stat_nested$npv*100),1),"%"),
                                    paste0(round(mean(results_PPS_PRS_clin_surv$c_stat_nested$npv*100),1),"%")),
                              intercept=c(round(mean(results_clin_surv$calibration_slopes$intercept),2),
                                          round(mean(results_PPS_surv$calibration_slopes$intercept),2),
                                          round(mean(results_PRS_surv$calibration_slopes$intercept),2),
                                          round(mean(results_PPS_clin_surv$calibration_slopes$intercept),2),
                                          round(mean(results_PRS_clin_surv$calibration_slopes$intercept),2),
                                          round(mean(results_PPS_PRS_surv$calibration_slopes$intercept),2),
                                          round(mean(results_PPS_PRS_clin_surv$calibration_slopes$intercept),2)),
                              slope=c(round(mean(results_clin_surv$calibration_slopes$slope),2),
                                      round(mean(results_PPS_surv$calibration_slopes$slope),2),
                                      round(mean(results_PRS_surv$calibration_slopes$slope),2),
                                      round(mean(results_PPS_clin_surv$calibration_slopes$slope),2),
                                      round(mean(results_PRS_clin_surv$calibration_slopes$slope),2),
                                      round(mean(results_PPS_PRS_surv$calibration_slopes$slope),2),
                                      round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope),2)),
                              Brier=c(round(mean(results_clin_surv$calibration_slopes$brier),2),
                                      round(mean(results_PPS_surv$calibration_slopes$brier),2),
                                      round(mean(results_PRS_surv$calibration_slopes$brier),2),
                                      round(mean(results_PPS_clin_surv$calibration_slopes$brier),2),
                                      round(mean(results_PRS_clin_surv$calibration_slopes$brier),2),
                                      round(mean(results_PPS_PRS_surv$calibration_slopes$brier),2),
                                      round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier),2))
                              
)
                              
write_csv(rf_rem_results, "~/Dropbox/Work/PPS/EU-GEI/rf_results_remission.csv")