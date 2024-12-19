library(tidyverse) # For data wrangling
library(plyr) # For data wrangling
library(readxl) # For reading in data
library(missForest) # For random forest imputation
library(glmnet) # For logistic regression and Cox
library(survminer) # For survival data
library(survival) # For survival data
library(survMisc) # For survival data
library(yardstick) # To calculate metrics
library(probably) # To calculate metrics
library(caret) # For model fitting and internal validation
library(scales) # For rescaling values

##### Read in already processed data #####
PPS_scored <- read.csv("~/Dropbox/Work/PPS/EU-GEI/Databases/PPS_processed.csv")

PPS_scored_all <- PPS_scored
PPS_scored <- PPS_scored %>% filter(chr==1)

##### Censor survival data at 3 years #####
PPS_scored_uncensored <- PPS_scored
PPS_scored <- PPS_scored %>% mutate(Transition = case_when(day_exit>1095.75 ~ 0,
                                                           TRUE ~ Transition),
                                    day_exit = case_when(day_exit>1095.75 ~ 1095.75,
                                                         TRUE ~ day_exit))

##### Assess proportion of variance explained #####
###### Detection ######
PPS_scored_all$Ethnicity_ED <- as.factor(PPS_scored_all$Ethnicity_ED)

PPS_variance <- PPS_scored_all %>% subset(select=c(chr, Gender, Handedness, Urbanicity, Pollution, Ethnicity_ED, 
                                                     FirstGenImmigrantNA, FirstGenImmigrant, SecondGenImmigrantNA, SecondGenImmigrant,  
                                                     PaternalSES, Parental_SMI, AdultLifeEvents, 
                                                     Tobacco, Cannabis, ChildhoodTrauma, Anhedonia, PRS_resid))

# Impute missing data with random forest imputation
imputed_data <- missForest(PPS_variance)
PPS_scored_all_imp <- as.data.frame(imputed_data$ximp)

PPS.lr <- lm(chr ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity_ED + 
             FirstGenImmigrantNA + FirstGenImmigrant + SecondGenImmigrantNA + SecondGenImmigrant +  
             PaternalSES + Parental_SMI + AdultLifeEvents + 
             Tobacco + Cannabis + ChildhoodTrauma + Anhedonia, data = PPS_scored_all_imp)
PPS.lr 
summary(PPS.lr)

PRS.lr <- lm(chr ~ PRS_resid, data = PPS_scored_all_imp)
PRS.lr 
summary(PRS.lr)

PPS_PRS.lr <- lm(chr ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity_ED + 
                   FirstGenImmigrantNA + FirstGenImmigrant + SecondGenImmigrantNA + SecondGenImmigrant +  
                   PaternalSES + Parental_SMI + AdultLifeEvents + 
                   Tobacco + Cannabis + ChildhoodTrauma + Anhedonia + PRS_resid, data = PPS_scored_all_imp)
PPS_PRS.lr 
summary(PPS_PRS.lr)

# Plot
variance_hc <- data.frame(factors=factor(c("PPS","PRS", "PRS + PPS"), levels=c("PPS", "PRS","PRS + PPS")),
                          variance=c(30.3,1.6,31.3))
ggplot(data=variance_hc, aes(x=factors, y=variance, fill=factors)) +
  geom_col() +
  geom_text(aes(label=paste(variance,"%", sep=""),y=variance+5), size=5) +
  ylab("Variance Explained (%)") +
  xlab("") +
  scale_y_continuous(limits=c(0,100), breaks = seq(0, 100, by = 10)) +
  theme_classic() + 
  scale_fill_manual(values=c("#599ec4","#ecc363","#7bae72")) +
  theme(legend.position = "none")
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Variance_hc_161124.png",width = 24, height = 18, units = "cm")

###### Prognosis ######
PPS_scored$Ethnicity <- as.factor(PPS_scored$Ethnicity)
PPS_scored$Ethnicity_ED <- as.factor(PPS_scored$Ethnicity_ED)

PPS_variance_chr <- PPS_scored %>% subset(select=c(Gender, Handedness, Urbanicity, Pollution, Ethnicity_ED, 
                                                   FirstGenImmigrantNA, FirstGenImmigrant, SecondGenImmigrantNA, SecondGenImmigrant,  
                                                   PaternalSES, Parental_SMI, AdultLifeEvents, 
                                                   Tobacco, Cannabis, ChildhoodTrauma, Anhedonia, PRS_resid, CAARMS, gafex01, gafex02,Transition, day_exit))

# Impute missing data with random forest imputation
imputed_data <- missForest(PPS_variance_chr)
PPS_scored_imp <- as.data.frame(imputed_data$ximp)
PPS_scored_imp <- PPS_scored_imp %>% mutate(Transition = case_when(Transition<0.5 ~ 0,
                                                                   Transition>=0.5 ~ 1))
clinical.cox <- coxph(Surv(day_exit, Transition) ~ CAARMS + gafex01 + gafex02, data = PPS_scored_imp)
clinical.cox
summary(clinical.cox)
survMisc::rsq(clinical.cox)

PPS.cox <- coxph(Surv(day_exit, Transition) ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity_ED + 
                   FirstGenImmigrantNA + FirstGenImmigrant + SecondGenImmigrantNA + SecondGenImmigrant +  
                   PaternalSES + Parental_SMI + AdultLifeEvents + 
                   Tobacco + Cannabis + ChildhoodTrauma + Anhedonia, data = PPS_scored_imp)
PPS.cox
summary(PPS.cox)
survMisc::rsq(PPS.cox)

PRS.cox <- coxph(Surv(day_exit, Transition) ~ PRS_resid, data = PPS_scored_imp)
PRS.cox
summary(PRS.cox)
survMisc::rsq(PRS.cox)

clinical_PRS.cox <- coxph(Surv(day_exit, Transition) ~ CAARMS + gafex01 + gafex02 + PRS_resid, data = PPS_scored_imp)
clinical_PRS.cox
summary(clinical_PRS.cox)
survMisc::rsq(clinical_PRS.cox)

clinical_PPS.cox <- coxph(Surv(day_exit, Transition) ~ CAARMS + gafex01 + gafex02 + Gender + Handedness + Urbanicity + Pollution + Ethnicity_ED + 
                            FirstGenImmigrantNA + FirstGenImmigrant + SecondGenImmigrantNA + SecondGenImmigrant +  
                            PaternalSES + Parental_SMI + AdultLifeEvents + 
                            Tobacco + Cannabis + ChildhoodTrauma + Anhedonia, data = PPS_scored_imp)
clinical_PPS.cox
summary(clinical_PPS.cox)
survMisc::rsq(clinical_PPS.cox)

PPS_PRS.cox <- coxph(Surv(day_exit, Transition) ~ Gender + Handedness + Urbanicity + Pollution + Ethnicity_ED + 
                       FirstGenImmigrantNA + FirstGenImmigrant + SecondGenImmigrantNA + SecondGenImmigrant +  
                       PaternalSES + Parental_SMI + AdultLifeEvents + 
                       Tobacco + Cannabis + ChildhoodTrauma + Anhedonia + PRS_resid, data = PPS_scored_imp)
PPS_PRS.cox
summary(PPS_PRS.cox)
survMisc::rsq(PPS_PRS.cox)

clinical_PPS_PRS.cox <- coxph(Surv(day_exit, Transition) ~ CAARMS + gafex01 + gafex02 + Gender + Handedness + Urbanicity + Pollution + Ethnicity_ED + 
                                FirstGenImmigrantNA + FirstGenImmigrant + SecondGenImmigrantNA + SecondGenImmigrant +  
                                PaternalSES + Parental_SMI + AdultLifeEvents + 
                                Tobacco + Cannabis + ChildhoodTrauma + Anhedonia + PRS_resid, data = PPS_scored_imp)
clinical_PPS_PRS.cox
summary(clinical_PPS_PRS.cox)
survMisc::rsq(clinical_PPS_PRS.cox)

variance <- data.frame(factors=factor(c("Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All"), levels=c("Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All")),
                       variance=c(6.9, 30.0, 3.3, 36.0, 10.0, 35.0, 40.0))
ggplot(data=variance, aes(x=factors, y=variance, fill=factors)) +
  geom_col() +
  geom_text(aes(label=paste(variance,"%", sep=""),y=variance+5), size=5) +
  ylab("Variance Explained (%)") +
  scale_y_continuous(limits=c(0,100), breaks = seq(0, 100, by = 10)) +
  xlab("") +
  theme_classic() + 
  scale_fill_manual(values=c("#c8526a","#599ec4","#ecc363","#8d828b","#d65f3d","#7bae72","#ab6758")) +
  theme(legend.position = "none")
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/Variance_Tx_021224.png",width = 24, height = 18, units = "cm")

##### Plot survival curve #####
surv <- with(PPS_scored_imp, Surv(day_exit, event = Transition))

#Default setting uses Greenwood CI - default is log transformation
km.surv <- survfit(surv ~ 1, data=PPS_scored_imp)
surv_summ <- data.frame(n=summary(km.surv)$n,
                        time=summary(km.surv)$time,
                        years=summary(km.surv)$time/365,
                        n.risk=summary(km.surv)$n.risk,
                        n.event=summary(km.surv)$n.event,
                        surv=paste0(round(1-(summary(km.surv)$surv),4)*100," (95%CI: ",round(1-(summary(km.surv)$upper),4)*100,
                                    "-",round(1-(summary(km.surv)$lower),4)*100, ", ",summary(km.surv)$n.risk,
                                    " individuals still at risk)"))

png("~/Dropbox/Work/PPS/EU-GEI/Plots/PPS_Surv_181224.png",width = 900, height = 800)
ggsurvplot(km.surv, data=PPS_scored_imp,
           fun="event",
           xscale = 365.25,
           break.time.by = 182.625,
           xlim = c(0,1095.75),
           cumevents = TRUE)
dev.off()

##### Detection LR ####

# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

###### Run model fitting and internal validation ######

# Run PPS model
PPS_PPS <- PPS_scored_all %>% subset(select=c(Gender:Tobacco, chr))

results_PPS <- LR_repeated_nested_cv(
  combined_df = PPS_PPS, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PRS model
PPS_genes <- PPS_scored_all %>% subset(select=c(PRS_resid, chr))

results_genes <- LR_repeated_nested_cv(
  combined_df = PPS_genes, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run combined PPS and PRS model
PPS_PRS <- PPS_scored_all %>% subset(select=c(Gender:Tobacco, PRS_resid, chr))

results_all <- LR_repeated_nested_cv(
  combined_df = PPS_PRS, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)
###### Generate performance metrics summary table######

lr_results <- data.frame(model=c("PPS", "PRS", "All"),
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
                         ppv=c(paste0(round(mean(results_PPS$c_stat_nested$ppv, na.rm=TRUE)*100,1), "% (",
                                      round((mean(results_PPS$c_stat_nested$ppv, na.rm=TRUE) - 1.96*(sd(results_PPS$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS$c_stat_nested$ppv, na.rm=TRUE) + 1.96*(sd(results_PPS$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes$c_stat_nested$ppv, na.rm=TRUE)*100,1), "% (",
                                      round((mean(results_genes$c_stat_nested$ppv, na.rm=TRUE) - 1.96*(sd(results_genes$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes$c_stat_nested$ppv, na.rm=TRUE) + 1.96*(sd(results_genes$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all$c_stat_nested$ppv, na.rm=TRUE)*100,1), "% (",
                                      round((mean(results_all$c_stat_nested$ppv, na.rm=TRUE) - 1.96*(sd(results_all$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all$c_stat_nested$ppv, na.rm=TRUE) + 1.96*(sd(results_all$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%)")),
                         npv=c(paste0(round(mean(results_PPS$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_PPS$c_stat_nested$npv) - 1.96*(sd(results_PPS$c_stat_nested$npv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS$c_stat_nested$npv) + 1.96*(sd(results_PPS$c_stat_nested$npv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_genes$c_stat_nested$npv) - 1.96*(sd(results_genes$c_stat_nested$npv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes$c_stat_nested$npv) + 1.96*(sd(results_genes$c_stat_nested$npv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_all$c_stat_nested$npv) - 1.96*(sd(results_all$c_stat_nested$npv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all$c_stat_nested$npv) + 1.96*(sd(results_all$c_stat_nested$npv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%)")),
                         intercept=c(paste0(round(mean(results_PPS$calibration_slopes$intercept),2), " (",
                                            round((mean(results_PPS$calibration_slopes$intercept) - 1.96*(sd(results_PPS$calibration_slopes$intercept) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_PPS$calibration_slopes$intercept) + 1.96*(sd(results_PPS$calibration_slopes$intercept) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_genes$calibration_slopes$intercept),2), " (",
                                            round((mean(results_genes$calibration_slopes$intercept) - 1.96*(sd(results_genes$calibration_slopes$intercept) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_genes$calibration_slopes$intercept) + 1.96*(sd(results_genes$calibration_slopes$intercept) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_all$calibration_slopes$intercept),2), " (",
                                            round((mean(results_all$calibration_slopes$intercept) - 1.96*(sd(results_all$calibration_slopes$intercept) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_all$calibration_slopes$intercept) + 1.96*(sd(results_all$calibration_slopes$intercept) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                            ")")),
                         slope=c(paste0(round(mean(results_PPS$calibration_slopes$slope),2), " (",
                                        round((mean(results_PPS$calibration_slopes$slope) - 1.96*(sd(results_PPS$calibration_slopes$slope) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS$calibration_slopes$slope) + 1.96*(sd(results_PPS$calibration_slopes$slope) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes$calibration_slopes$slope),2), " (",
                                        round((mean(results_genes$calibration_slopes$slope) - 1.96*(sd(results_genes$calibration_slopes$slope) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes$calibration_slopes$slope) + 1.96*(sd(results_genes$calibration_slopes$slope) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all$calibration_slopes$slope),2), " (",
                                        round((mean(results_all$calibration_slopes$slope) - 1.96*(sd(results_all$calibration_slopes$slope) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all$calibration_slopes$slope) + 1.96*(sd(results_all$calibration_slopes$slope) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        ")")),
                         Brier=c(paste0(round(mean(results_PPS$calibration_slopes$brier),2), " (",
                                        round((mean(results_PPS$calibration_slopes$brier) - 1.96*(sd(results_PPS$calibration_slopes$brier) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS$calibration_slopes$brier) + 1.96*(sd(results_PPS$calibration_slopes$brier) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes$calibration_slopes$brier),2), " (",
                                        round((mean(results_genes$calibration_slopes$brier) - 1.96*(sd(results_genes$calibration_slopes$brier) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes$calibration_slopes$brier) + 1.96*(sd(results_genes$calibration_slopes$brier) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all$calibration_slopes$brier),2), " (",
                                        round((mean(results_all$calibration_slopes$brier) - 1.96*(sd(results_all$calibration_slopes$brier) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all$calibration_slopes$brier) + 1.96*(sd(results_all$calibration_slopes$brier) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        ")")),
                         ICI=c(paste0(round(mean(results_PPS$calibration_slopes$ICI),2), " (",
                                      round((mean(results_PPS$calibration_slopes$ICI) - 1.96*(sd(results_PPS$calibration_slopes$ICI) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_PPS$calibration_slopes$ICI) + 1.96*(sd(results_PPS$calibration_slopes$ICI) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_genes$calibration_slopes$ICI),2), " (",
                                      round((mean(results_genes$calibration_slopes$ICI) - 1.96*(sd(results_genes$calibration_slopes$ICI) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_genes$calibration_slopes$ICI) + 1.96*(sd(results_genes$calibration_slopes$ICI) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_all$calibration_slopes$ICI),2), " (",
                                      round((mean(results_all$calibration_slopes$ICI) - 1.96*(sd(results_all$calibration_slopes$ICI) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_all$calibration_slopes$ICI) + 1.96*(sd(results_all$calibration_slopes$ICI) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                      ")"))
                         
)
write_csv(lr_results, "~/Dropbox/Work/PPS/EU-GEI/lr_results_171224.csv")
write_csv(results_all$dca_general, "~/Dropbox/Work/PPS/EU-GEI/lr_all_dcacgen_171224.csv")
write_csv(results_all$dca_clinical, "~/Dropbox/Work/PPS/EU-GEI/lr_all_dcaclin_171224.csv")
write_csv(results_PPS$dca_general, "~/Dropbox/Work/PPS/EU-GEI/lr_PPS_dcacgen_171224.csv")
write_csv(results_PPS$dca_clinical, "~/Dropbox/Work/PPS/EU-GEI/lr_PPS_dcaclin_171224.csv")
write_csv(results_genes$dca_general, "~/Dropbox/Work/PPS/EU-GEI/lr_PRS_dcacgen_171224.csv")
write_csv(results_genes$dca_clinical, "~/Dropbox/Work/PPS/EU-GEI/lr_PRS_dcaclin_171224.csv")

###### Generate coefficients #######
PPS_PPS[,c(1:16)] <- lapply(PPS_PPS[,c(1:16)], factor)
options(na.action='na.pass')
temp_mat <- model.matrix(~.-1,PPS_PPS[,c(1:16)])
temp_mat <- temp_mat %>% subset(select=c(-Gender0))
imputed_data <- missForest(temp_mat)
temp_mat <- imputed_data$ximp

data <- data.frame(temp_mat, chr = factor(PPS_scored_all$chr))

recipe <- recipe(chr ~ ., data = data) %>%
  step_smote(chr, over_ratio = 1)

# Apply the recipe
processed_data <- prep(recipe) %>% bake(new_data = NULL)

# Convert to matrix for glmnet
smote_x <- as.matrix(processed_data[, -ncol(processed_data)])
smote_y <- processed_data$chr

# Train the glmnet model with SMOTE-transformed data
check <- cv.glmnet(
  y = smote_y,
  x = smote_x,
  family = "binomial",
  nfolds = 10,   # Number of folds (adjust if needed)
  alpha = 1,
  seed = 123
)
tmp_coeffs <- coef(check, s = cv.glmnet(as.matrix(temp_mat), PPS_scored_all$chr)$lambda.min)
PPS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

imputed_data <- mice(PPS_genes)
temp_mat <- complete(imputed_data)
temp_mat <- temp_mat[,1]
data <- data.frame(temp_mat, chr = factor(PPS_scored_all$chr))

recipe <- recipe(chr ~ ., data = data) %>%
  step_smote(chr, over_ratio = 1)

# Apply the recipe
processed_data <- prep(recipe) %>% bake(new_data = NULL)

# Convert to matrix for glmnet
smote_x <- as.matrix(processed_data[, -ncol(processed_data)])
smote_y <- processed_data$chr

# Train the glmnet model with SMOTE-transformed data
check <- glm(smote_y ~ smote_x,
             family="binomial")
tmp_coeffs <- coef(check, s = cv.glmnet(as.matrix(temp_mat), PPS_scored_all$chr)$lambda.min)
PRS_coef <- data.frame(name = "PRS", coefficient = check$coefficients[2])

PPS_PRS[,c(1:16)] <- lapply(PPS_PRS[,c(1:16)], factor)
temp_mat <- model.matrix(~.-1,PPS_PRS[,c(1:17)])
temp_mat <- temp_mat %>% subset(select=c(-Gender0))
imputed_data <- missForest(temp_mat)
temp_mat <- imputed_data$ximp
data <- data.frame(temp_mat, chr = factor(PPS_scored_all$chr))

recipe <- recipe(chr ~ ., data = data) %>%
  step_smote(chr, over_ratio = 1)

# Apply the recipe
processed_data <- prep(recipe) %>% bake(new_data = NULL)

# Convert to matrix for glmnet
smote_x <- as.matrix(processed_data[, -ncol(processed_data)])
smote_y <- processed_data$chr

# Train the glmnet model with SMOTE-transformed data
check <- cv.glmnet(
  y = smote_y,
  x = smote_x,
  family = "binomial",
  nfolds = 10,   # Number of folds (adjust if needed)
  alpha = 1,
  seed = 123
)
tmp_coeffs <- coef(check, s = cv.glmnet(as.matrix(temp_mat), PPS_scored_all$chr)$lambda.min)
PPS_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

write_csv(PPS_coef,"~/Dropbox/Work/PPS/EU-GEI/PPS_coef_191224.csv")
write_csv(PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/PRS_coef_191224.csv")
write_csv(PPS_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_coef_191224.csv")

#### Prognosis Cox ####

# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

###### Run model fitting and internal validation ######

# Run clinical model
clinical <- PPS_scored %>% subset(select=c(CAARMS, gafex01, gafex02, Transition, day_exit))

results_clin_surv <- Cox_repeated_nested_cv(
  combined_df = clinical, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS model
PPS_PPS <- PPS_scored %>% subset(select=c(Gender:Tobacco, Transition, day_exit))

results_PPS_surv <- Cox_repeated_nested_cv(
  combined_df = PPS_PPS, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PRS model
PPS_genes <- PPS_scored %>% subset(select=c(PRS_resid, Transition, day_exit))

results_PRS_surv <- Cox_repeated_nested_cv(
  combined_df = PPS_genes, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS+clinical model
PPS_clin <- PPS_scored %>% subset(select=c(Gender:Tobacco, CAARMS, gafex01, gafex02, Transition, day_exit))

results_PPS_clin_surv <- Cox_repeated_nested_cv(
  combined_df = PPS_clin, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PRS+clinical model
PRS_clin <- PPS_scored %>% subset(select=c(CAARMS, gafex01, gafex02, PRS_resid, Transition, day_exit))

results_PRS_clin_surv <- Cox_repeated_nested_cv(
  combined_df = PRS_clin, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS+PRS model
PPS_PRS <- PPS_scored %>% subset(select=c(Gender:Tobacco, PRS_resid, Transition, day_exit))

results_PPS_PRS_surv <- Cox_repeated_nested_cv(
  combined_df = PPS_PRS, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS+PRS+clinical model
PPS_PRS_clin <- PPS_scored %>% subset(select=c(Gender:Tobacco, CAARMS, gafex01, gafex02, PRS_resid, Transition, day_exit))

results_PPS_PRS_clin_surv <- Cox_repeated_nested_cv(
  combined_df = PPS_PRS_clin, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

###### Generate performance metrics summary table######

cox_results <- data.frame(model=c("Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All"),
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
                              calibration_in_the_large=c(paste0(round(mean(results_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PPS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PPS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PPS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PRS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PRS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PRS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                                                ")")),
                              slope=c(paste0(round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                             ")")),
                              Brier=c(paste0(round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                             ")")),
                          ICI=c(paste0(round(mean(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                         round(mean(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                         ")"))
                              
)
write_csv(cox_results, "~/Dropbox/Work/PPS/EU-GEI/cox_results_031224.csv")

###### Generate coefficients #####
imputed_data <- missForest(PPS_variance_chr)
PPS_scored_imp <- as.data.frame(imputed_data$ximp)

clinical <- PPS_scored_imp %>% subset(select=c(CAARMS, gafex01, gafex02))
PPS_PPS <- PPS_scored_imp %>% subset(select=c(Gender:Anhedonia))
PPS_genes <- PPS_scored_imp %>% subset(select=c(PRS_resid))
PPS_clin <- PPS_scored_imp %>% subset(select=c(Gender:Anhedonia,CAARMS, gafex01, gafex02))
PRS_clin <- PPS_scored_imp %>% subset(select=c(PRS_resid,CAARMS, gafex01, gafex02))
PPS_PRS <- PPS_scored_imp %>% subset(select=c(Gender:Anhedonia, PRS_resid))
PPS_PRS_clin <- PPS_scored_imp %>% subset(select=c(Gender:Anhedonia, PRS_resid, CAARMS, gafex01, gafex02))


temp_mat <- model.matrix(~.-1,clinical)
check <- cv.glmnet(y=Surv(PPS_scored_imp$day_exit, PPS_scored_imp$Transition), x=temp_mat,family="cox",
                type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp_mat <- model.matrix(~.-1,PPS_PPS)
check <- cv.glmnet(y=Surv(PPS_scored_imp$day_exit, PPS_scored_imp$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
PPS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp_mat <- model.matrix(~.-1,PPS_genes)
check <- coxph(Surv(time = day_exit, event = Transition) ~ PRS_resid, data=PPS_scored_imp)
PRS_coef <- as.data.frame(coef(check))

temp_mat <- model.matrix(~.-1,PPS_clin)
check <- cv.glmnet(y=Surv(PPS_scored_imp$day_exit, PPS_scored_imp$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_PPS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp_mat <- model.matrix(~.-1,PRS_clin)
check <- cv.glmnet(y=Surv(PPS_scored_imp$day_exit, PPS_scored_imp$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp_mat <- model.matrix(~.-1, PPS_PRS)
check <- cv.glmnet(y=Surv(PPS_scored_imp$day_exit, PPS_scored_imp$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
PPS_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

temp_mat <- model.matrix(~.-1, PPS_PRS_clin)
check <- cv.glmnet(y=Surv(PPS_scored_imp$day_exit, PPS_scored_imp$Transition), x=temp_mat,family="cox",
                   type.measure = "C", alpha=1, seed=123)
tmp_coeffs <- coef(check, s = check$lambda.min)
clin_PPS_PRS_coef <- data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)

write_csv(clin_coef,"~/Dropbox/Work/PPS/EU-GEI/clin_coef_surv_131224.csv")
write_csv(PPS_coef,"~/Dropbox/Work/PPS/EU-GEI/PPS_coef_surv_131224.csv")
write_csv(PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/PRS_coef_surv_131224.csv")
write_csv(clin_PPS_coef,"~/Dropbox/Work/PPS/EU-GEI/clinPPS_coef_surv_131224.csv")
write_csv(clin_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/clin_PRS_coef_surv_131224.csv")
write_csv(PPS_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_coef_surv_131224.csv")
write_csv(clin_PPS_PRS_coef,"~/Dropbox/Work/PPS/EU-GEI/clin_PPS_PRS_coef_surv_131224.csv")

##### Detection Random Forest ####
    # Parameters for nested cross-validation
    outerFolds <- 5
    outerRepeats <- 10
    innerFolds <- 5
    innerRepeats <- 10
   
    ###### Run model fitting and internal validation ######
    # Perform nested cross-validation with repeats
    PPS_PRS <- PPS_scored_all %>% subset(select=c(Gender:Tobacco, PRS_resid, chr))
    
    # Define tuning grid
    tune_grid <- expand.grid(
      mtry=c(3:(ncol(PPS_PRS)-2))
    )
    
    results_all <- RF_repeated_nested_cv(
      combined_df = PPS_PRS, 
      outerFolds = outerFolds, 
      outerRepeats = outerRepeats, 
      innerFolds = innerFolds, 
      innerRepeats = innerRepeats, 
      tuneGrid = tune_grid,
      seed = 231
    )
    
    PPS_genes <- PPS_scored_all %>% subset(select=c(PRS_resid, chr))
    tune_grid_genes <- expand.grid(
      mtry=1
    )
    results_genes <- RF_repeated_nested_cv(
      combined_df = PPS_genes, 
      outerFolds = outerFolds, 
      outerRepeats = outerRepeats, 
      innerFolds = innerFolds, 
      innerRepeats = innerRepeats, 
      tuneGrid = tune_grid_genes,
      seed = 231
    )
    
    PPS_PPS <- PPS_scored_all %>% subset(select=c(Gender:Tobacco, chr))
    tune_grid_PPS <- expand.grid(
      mtry=c(3:(ncol(PPS_PPS)-2))
    )
    results_PPS <- RF_repeated_nested_cv(
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
                         ppv=c(paste0(round(mean(results_PPS$c_stat_nested$ppv)*100,1), "% (",
                                      round((mean(results_PPS$c_stat_nested$ppv) - 1.96*(sd(results_PPS$c_stat_nested$ppv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS$c_stat_nested$ppv) + 1.96*(sd(results_PPS$c_stat_nested$ppv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes$c_stat_nested$ppv)*100,1), "% (",
                                      round((mean(results_genes$c_stat_nested$ppv) - 1.96*(sd(results_genes$c_stat_nested$ppv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes$c_stat_nested$ppv) + 1.96*(sd(results_genes$c_stat_nested$ppv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all$c_stat_nested$ppv)*100,1), "% (",
                                      round((mean(results_all$c_stat_nested$ppv) - 1.96*(sd(results_all$c_stat_nested$ppv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all$c_stat_nested$ppv) + 1.96*(sd(results_all$c_stat_nested$ppv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%)")),
                         npv=c(paste0(round(mean(results_PPS$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_PPS$c_stat_nested$npv) - 1.96*(sd(results_PPS$c_stat_nested$npv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS$c_stat_nested$npv) + 1.96*(sd(results_PPS$c_stat_nested$npv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_genes$c_stat_nested$npv) - 1.96*(sd(results_genes$c_stat_nested$npv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes$c_stat_nested$npv) + 1.96*(sd(results_genes$c_stat_nested$npv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_all$c_stat_nested$npv) - 1.96*(sd(results_all$c_stat_nested$npv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all$c_stat_nested$npv) + 1.96*(sd(results_all$c_stat_nested$npv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%)")),
                         intercept=c(paste0(round(mean(results_PPS$calibration_slopes$intercept),2), " (",
                                            round((mean(results_PPS$calibration_slopes$intercept) - 1.96*(sd(results_PPS$calibration_slopes$intercept) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_PPS$calibration_slopes$intercept) + 1.96*(sd(results_PPS$calibration_slopes$intercept) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_genes$calibration_slopes$intercept),2), " (",
                                            round((mean(results_genes$calibration_slopes$intercept) - 1.96*(sd(results_genes$calibration_slopes$intercept) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_genes$calibration_slopes$intercept) + 1.96*(sd(results_genes$calibration_slopes$intercept) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_all$calibration_slopes$intercept),2), " (",
                                            round((mean(results_all$calibration_slopes$intercept) - 1.96*(sd(results_all$calibration_slopes$intercept) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_all$calibration_slopes$intercept) + 1.96*(sd(results_all$calibration_slopes$intercept) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                            ")")),
                         slope=c(paste0(round(mean(results_PPS$calibration_slopes$slope),2), " (",
                                        round((mean(results_PPS$calibration_slopes$slope) - 1.96*(sd(results_PPS$calibration_slopes$slope) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS$calibration_slopes$slope) + 1.96*(sd(results_PPS$calibration_slopes$slope) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes$calibration_slopes$slope),2), " (",
                                        round((mean(results_genes$calibration_slopes$slope) - 1.96*(sd(results_genes$calibration_slopes$slope) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes$calibration_slopes$slope) + 1.96*(sd(results_genes$calibration_slopes$slope) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all$calibration_slopes$slope),2), " (",
                                        round((mean(results_all$calibration_slopes$slope) - 1.96*(sd(results_all$calibration_slopes$slope) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all$calibration_slopes$slope) + 1.96*(sd(results_all$calibration_slopes$slope) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        ")")),
                         Brier=c(paste0(round(mean(results_PPS$calibration_slopes$brier),2), " (",
                                        round((mean(results_PPS$calibration_slopes$brier) - 1.96*(sd(results_PPS$calibration_slopes$brier) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS$calibration_slopes$brier) + 1.96*(sd(results_PPS$calibration_slopes$brier) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes$calibration_slopes$brier),2), " (",
                                        round((mean(results_genes$calibration_slopes$brier) - 1.96*(sd(results_genes$calibration_slopes$brier) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes$calibration_slopes$brier) + 1.96*(sd(results_genes$calibration_slopes$brier) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all$calibration_slopes$brier),2), " (",
                                        round((mean(results_all$calibration_slopes$brier) - 1.96*(sd(results_all$calibration_slopes$brier) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all$calibration_slopes$brier) + 1.96*(sd(results_all$calibration_slopes$brier) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        ")")),
                         ICI=c(paste0(round(mean(results_PPS$calibration_slopes$ICI),2), " (",
                                      round((mean(results_PPS$calibration_slopes$ICI) - 1.96*(sd(results_PPS$calibration_slopes$ICI) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_PPS$calibration_slopes$ICI) + 1.96*(sd(results_PPS$calibration_slopes$ICI) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_genes$calibration_slopes$ICI),2), " (",
                                      round((mean(results_genes$calibration_slopes$ICI) - 1.96*(sd(results_genes$calibration_slopes$ICI) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_genes$calibration_slopes$ICI) + 1.96*(sd(results_genes$calibration_slopes$ICI) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_all$calibration_slopes$ICI),2), " (",
                                      round((mean(results_all$calibration_slopes$ICI) - 1.96*(sd(results_all$calibration_slopes$ICI) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_all$calibration_slopes$ICI) + 1.96*(sd(results_all$calibration_slopes$ICI) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                      ")"))
                         )
write_csv(rf_results, "~/Dropbox/Work/PPS/EU-GEI/rf_results_081224.csv")

###### Generate variable importance ######

imputed_data <- missForest(PPS_variance)
PPS_scored_all_imp <- as.data.frame(imputed_data$ximp)

PPS_PPS <- PPS_scored_all_imp %>% subset(select=c(Gender:Anhedonia, chr))
PPS_PPS$chr <- factor(PPS_PPS$chr)
PPS_PPS$chr <- factor(make.names(levels(PPS_PPS$chr))[PPS_PPS$chr])

PPS_PRS <- PPS_scored_all_imp %>% subset(select=c(Gender:Anhedonia, PRS_resid, chr))
PPS_PRS$chr <- factor(PPS_PRS$chr)
PPS_PRS$chr <- factor(make.names(levels(PPS_PRS$chr))[PPS_PRS$chr])

PPS_genes <- PPS_scored_all_imp %>% subset(select=c(PRS_resid, chr))
PPS_genes$chr <- factor(PPS_genes$chr)
PPS_genes$chr <- factor(make.names(levels(PPS_genes$chr))[PPS_genes$chr])

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
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PPS_VImp_121224.csv")

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
write_csv(vimp_PRS, "~/Dropbox/Work/PPS/EU-GEI/PRS_VImp_121224.csv")

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

write_csv(vimp_all, "~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_VImp_121224.csv")

pdf("~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_VImp.pdf")
plot(varImp(final_model), top = 20)
dev.off()

##### Prognosis Random Survival Forest ####
# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

###### Run model fitting and internal validation ######
# Clinical
clin_rf_surv <- PPS_scored %>% subset(select=c(CAARMS, gafex01, gafex02, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3,3,3)
)

results_clin_surv <- RSF_repeated_nested_cv(
  combined_df = clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PPS

PPS_rf_surv <- PPS_scored %>% subset(select=c(Gender:Tobacco, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_rf_surv)-2)
)

results_PPS_surv <- RSF_repeated_nested_cv(
  combined_df = PPS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PRS
PRS_rf_surv <- PPS_scored %>% subset(select=c(PRS_resid, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PRS_rf_surv)-2)
)

results_PRS_surv <- RSF_repeated_nested_cv(
  combined_df = PRS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PPS

PPS_clin_rf_surv <- PPS_scored %>% subset(select=c(Gender:Tobacco, CAARMS, gafex01, gafex02, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_clin_rf_surv)-2)
)

results_PPS_clin_surv <- RSF_repeated_nested_cv(
  combined_df = PPS_clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PRS

PRS_clin_rf_surv <- PPS_scored %>% subset(select=c(PRS_resid, CAARMS, gafex01, gafex02, Transition, day_exit))
PRS_clin_rf_surv <- as.data.frame(model.matrix(~.-1,PRS_clin_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_rf_surv)-2)
)

results_PRS_clin_surv <- RSF_repeated_nested_cv(
  combined_df = PRS_clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PPS + PRS

PPS_PRS_rf_surv <- PPS_scored %>% subset(select=c(Gender:Tobacco, PRS_resid, Transition, day_exit))
PPS_PRS_rf_surv <- as.data.frame(model.matrix(~.-1,PPS_PRS_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_PRS_rf_surv)-2)
)

results_PPS_PRS_surv <- RSF_repeated_nested_cv(
  combined_df = PPS_PRS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PPS + PRS 

PPS_PRS_clin_rf_surv <- PPS_scored %>% subset(select=c(Gender:Tobacco, PRS_resid, CAARMS, gafex01, gafex02, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_PRS_clin_rf_surv)-2)
)

results_PPS_PRS_clin_surv <- RSF_repeated_nested_cv(
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
                          calibration_in_the_large=c(paste0(round(mean(results_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PRS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PRS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PRS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                                            ")")),
                          slope=c(paste0(round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                         ")")),
                          Brier=c(paste0(round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                         ")")),
                          ICI=c(paste0(round(mean(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                       ")"))
                         
)
write_csv(rf_surv_results, "~/Dropbox/Work/PPS/EU-GEI/rf_surv_results_081224.csv")

###### Generate variable importance ######

PPS_variance_chr2 <- as.data.frame(model.matrix(~.-1, data=PPS_variance_chr))
imputed_data <- missForest(PPS_variance_chr2)
PPS_scored_imp <- as.data.frame(imputed_data$ximp)

clinical <- PPS_scored_imp %>% subset(select=c(CAARMS, gafex01, gafex02))
PPS_PPS <- PPS_scored_imp %>% subset(select=c(Gender:Anhedonia))
PPS_genes <- PPS_scored_imp %>% subset(select=c(PRS_resid))
PPS_clin <- PPS_scored_imp %>% subset(select=c(Gender:Anhedonia,CAARMS, gafex01, gafex02))
PRS_clin <- PPS_scored_imp %>% subset(select=c(PRS_resid,CAARMS, gafex01, gafex02))
PPS_PRS <- PPS_scored_imp %>% subset(select=c(Gender:Anhedonia, PRS_resid))
PPS_PRS_clin <- PPS_scored_imp %>% subset(select=c(Gender:Anhedonia, PRS_resid, CAARMS, gafex01, gafex02))

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
                     data = PPS_scored_imp,
                     nodesize=best_min.node.size,
                     ntree=(best_ntree),
                     mtry = repGrid$mtry)
final_model <- rfsrc(formula,
                     data = PPS_scored_imp,
                     nodesize=4,
                     ntree=(500),
                     mtry = 6)

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

predictor_vars <- setdiff(names(PPS_scored_imp), c("day_exit", "Transition", "CAARMS", "gafex01", 'gafex02', "PRS_resid"))

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
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PPS_surv_VImp_131224.csv")

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

vimp_PPS<- vimp(final_model)$importance
vimp_PPS <- as.data.frame(sort(vimp_PPS, decreasing=TRUE))
vimp_PPS <- tibble::rownames_to_column(vimp_PPS, "Predictor")
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PRS_surv_VImp_131224.csv")

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

predictor_vars <- setdiff(names(PPS_scored_imp), c("day_exit", "Transition",'PRS_resid'))

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
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PPS_clin_surv_VImp_131224.csv")

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
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PRS_clin_surv_VImp_131224.csv")

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

predictor_vars <- setdiff(names(PPS_scored_imp), c("day_exit", "Transition", "CAARMS", "gafex01", "gafex02"))

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
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_surv_VImp_131224.csv")

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

predictor_vars <- setdiff(names(PPS_scored_imp), c("day_exit", "Transition"))

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
write_csv(vimp_PPS, "~/Dropbox/Work/PPS/EU-GEI/PPS_PRS_clin_surv_VImp_131224.csv")

##### Summary discrimination figure #####
results_plot <- data.frame(Type=c("LR","LR","LR",
                                  "RF","RF","RF",
                                  "Cox", "Cox","Cox","Cox","Cox","Cox","Cox",
                                  "RSF","RSF","RSF","RSF","RSF","RSF","RSF"),
                              model=c("PPS", "PRS",  "PPS+PRS", 
                                      "PPS", "PRS",  "PPS+PRS",
                                      "Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All",
                                      "Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All"),
                              C=c(0.88, 0.62, 0.88,
                                  0.83, 0.55, 0.82,
                                  0.61,0.54,0.58,0.55,0.57,0.55,0.55,
                                  0.58,0.58,0.59,0.57,0.58,0.57,0.57),
                           lCI=c(0.86,0.60,0.87,
                                 0.81,0.52,0.79,
                                 0.59,0.53,0.56,0.53,0.55,0.53,0.53,
                                 0.51,0.57,0.57,0.55,0.57,0.56,0.56),
                           uCI=c(0.89,0.64,0.89,
                                 0.86,0.59,0.86,
                                 0.63,0.56,0.59,0.57,0.59,0.53,0.53,
                                 0.65,0.59,0.6,0.59,0.6,0.58,0.58),
                           C_se=c(0.0075,0.025,0.01,
                                  0.01,0.015,0.01,
                                  0.035,0.05,0.015,0.04,0.02,0.055,0.025,
                                  0.01,0.01,0.005,0.005,0.01,0.01,0.01),
                           intercept=c(0.16,-0.93,0.1,
                                       0.20,1.12,0.13,
                                       0.12,0.27,0.15,0.22,0.16,0.24,0.19,
                                       -1.18,-2.34,-1.19,-1.28,-1.1,-1.83,-1.37
                                       
                           ),
                           slope=c(0.92,1.62,0.98,
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
  geom_text(aes(x=0.66, y=0.95,  label="0.88"), stat="unique", size=8, color="#599ec4", family="Roboto Condensed") +
  geom_text(aes(x=1.00, y=0.95,  label="0.62"), stat="unique", size=8, color="#ecc363", family="Roboto Condensed") +
  geom_text(aes(x=1.33, y=0.95,  label="0.89"), stat="unique", size=8, color="#7bae72", family="Roboto Condensed") +
  geom_text(aes(x=1.66, y=0.95,  label="0.86"), stat="unique", size=8, color="#599ec4", family="Roboto Condensed") +
  geom_text(aes(x=2.00, y=0.95,  label="0.60"), stat="unique", size=8, color="#ecc363", family="Roboto Condensed") +
  geom_text(aes(x=2.33, y=0.95,  label="0.87"), stat="unique", size=8, color="#7bae72", family="Roboto Condensed") +
  
  geom_text(aes(x=2.61, y=0.94,  label="0.61"), stat="unique", size=8, color="#c8526a", family="Roboto Condensed") +
  geom_text(aes(x=2.74, y=0.96,  label="0.54"), stat="unique", size=8, color="#599ec4", family="Roboto Condensed") +
  geom_text(aes(x=2.87, y=0.94,  label="0.58"), stat="unique", size=8, color="#ecc363", family="Roboto Condensed") +
  geom_text(aes(x=3.00, y=0.96,  label="0.55"), stat="unique", size=8, color="#8d828b", family="Roboto Condensed") +
  geom_text(aes(x=3.13, y=0.94,  label="0.57"), stat="unique", size=8, color="#d65f3d", family="Roboto Condensed") +
  geom_text(aes(x=3.26, y=0.96,  label="0.55"), stat="unique", size=8, color="#7bae72", family="Roboto Condensed") +
  geom_text(aes(x=3.39, y=0.94,  label="0.55"), stat="unique", size=8, color="#ab6758", family="Roboto Condensed") +
  
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
  scale_y_continuous(limits=c(0.4,1), breaks = seq(0.4, 1, by = 0.10)) +
  guides(color = guide_legend(title = "Predictors")) +
  theme(text = element_text(family="Roboto", face="bold", size=21),legend.title = element_text(size=23),legend.text = element_text(size = 23))
ggsave("~/Dropbox/Work/PPS/EU-GEI/Summary_181224.png",width = 42, height = 32, units = "cm")

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
dcaGen_PPS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/lr_PPS_dcacgen_171224.csv")
dcaGen_PRS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/lr_PRS_dcacgen_171224.csv")
dcaGen_PPSPRS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/lr_all_dcacgen_171224.csv")

dcaGen <- dcaGen_PPS %>% filter(label!="pred")
dcaGen_PPS <- dcaGen_PPS %>% filter(label=="pred")
dcaGen_PRS <- dcaGen_PRS %>% filter(label=="pred")
dcaGen_PPSPRS <- dcaGen_PPSPRS %>% filter(label=="pred")

dcaGen_PPS$variable <- "PPS"
dcaGen_PRS$variable <- "PRS"
dcaGen_PPSPRS$variable <- "PPS+PRS"
dcaGen$variable <- NA

dcaGen_PPS$label <- "PPS"
dcaGen_PRS$label <- "PRS"
dcaGen_PPSPRS$label <- "PPS+PRS"

dcaGen <- rbind(dcaGen, dcaGen_PPS)
dcaGen <- rbind(dcaGen, dcaGen_PRS)
dcaGen <- rbind(dcaGen, dcaGen_PPSPRS)

dcaGen$label <- factor(dcaGen$label,levels = c("Treat All", "Treat None", "PPS", "PRS", "PPS+PRS"))
dcaGen2 <- dcaGen %>% aggregate(net_benefit ~ label + threshold, FUN="mean")
dcaGen_wide <- pivot_wider(dcaGen2, 
                            names_from = label, 
                            values_from = net_benefit)
dcaGenSummary <- dcaGen_wide %>% group_by(threshold) %>% mutate(PPS_nb = case_when(`Treat All`>0 ~ PPS -`Treat All`,
                                                                                     TRUE ~ PPS),
                                                                  PRS_nb = case_when(`Treat All`>0 ~ PRS -`Treat All`,
                                                                                     TRUE ~ PRS),
                                                                  all_nb = case_when(`Treat All`>0 ~ `PPS+PRS` -`Treat All`,
                                                                                     TRUE ~ `PPS+PRS`))

dcaGenSummary2 <- dcaGenSummary %>% subset(select=c(threshold, PPS_nb, PRS_nb, all_nb)) %>% pivot_longer(cols = -threshold, 
                                                                                                        names_to = "label", 
                                                                                                        values_to = "net_benefit")
dcaGenSummary2 <- dcaGenSummary2 %>% mutate(standardised_net_benefit = net_benefit/0.017) %>% 
                                     filter(threshold==0 | threshold==0.01 |threshold==0.02 |threshold==0.03 |threshold==0.04 |threshold==0.05 |threshold==0.06 |threshold==0.07 |threshold==0.08 |threshold==0.09 |threshold==0.1)
write_csv(dcaGenSummary2,"~/Dropbox/Work/PPS/EU-GEI/dca_gen_results.csv")                               

ggplot(data=dcaGen2, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2) +
  coord_cartesian(ylim = c(-0.005, 0.02)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,0.03)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  scale_color_manual(values=c("gray80","#000000","#599ec4","#ecc363","#7bae72")) +
  theme(text = element_text(family="Roboto", face="bold", size=30),legend.title = element_text(size=23),legend.text = element_text(size = 23)) +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/dca_det_gen_summary_171224.png", width=20, height=15, scale=0.5)

dcaClin_PPS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/lr_PPS_dcaclin_171224.csv")
dcaClin_PRS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/lr_PRS_dcaclin_171224.csv")
dcaClin_PPSPRS <- read_csv("~/Dropbox/Work/PPS/EU-GEI/lr_all_dcaclin_171224.csv")

dcaClin <- dcaClin_PPS %>% filter(label!="pred")
dcaClin_PPS <- dcaClin_PPS %>% filter(label=="pred")
dcaClin_PRS <- dcaClin_PRS %>% filter(label=="pred")
dcaClin_PPSPRS <- dcaClin_PPSPRS %>% filter(label=="pred")

dcaClin_PPS$variable <- "PPS"
dcaClin_PRS$variable <- "PRS"
dcaClin_PPSPRS$variable <- "PPS+PRS"
dcaClin$variable <- NA

dcaClin_PPS$label <- "PPS"
dcaClin_PRS$label <- "PRS"
dcaClin_PPSPRS$label <- "PPS+PRS"

dcaClin <- rbind(dcaClin, dcaClin_PPS)
dcaClin <- rbind(dcaClin, dcaClin_PRS)
dcaClin <- rbind(dcaClin, dcaClin_PPSPRS)

dcaClin$label <- factor(dcaClin$label,levels = c("Treat All", "Treat None", "PPS", "PRS", "PPS+PRS"))
dcaClin2 <- dcaClin %>% aggregate(net_benefit ~ label + threshold, FUN="mean")
dcaClin2 <- dcaClin2 %>% mutate(standardised_net_benefit = net_benefit/0.192)
dcaClin_wide <- pivot_wider(dcaClin2, 
                            names_from = label, 
                            values_from = c(net_benefit,standardised_net_benefit))
dcaClinSummary <- dcaClin_wide %>% group_by(threshold) %>% mutate(PPS_nb = case_when(`net_benefit_Treat All`>0 ~ net_benefit_PPS -`net_benefit_Treat All`,
                                                                                     TRUE ~ net_benefit_PPS),
                                                                  PRS_nb = case_when(`net_benefit_Treat All`>0 ~ net_benefit_PRS -`net_benefit_Treat All`,
                                                                                     TRUE ~ net_benefit_PRS),
                                                                  all_nb = case_when(`net_benefit_Treat All`>0 ~ `net_benefit_PPS+PRS` -`net_benefit_Treat All`,
                                                                                     TRUE ~ `net_benefit_PPS+PRS`))


dcaClinSummary2 <- dcaClinSummary %>% subset(select=c(threshold, PPS_nb, PRS_nb, all_nb)) %>% pivot_longer(cols = -threshold, 
                                                                                                         names_to = "label", 
                                                                                                         values_to = "net_benefit")
dcaClinSummary2 <- dcaClinSummary2 %>% mutate(standardised_net_benefit = net_benefit/0.192) %>% 
  filter(threshold==0 | threshold==0.1 |threshold==0.2 |threshold==0.3 |threshold==0.4 |threshold==0.5)
write_csv(dcaClinSummary2,"~/Dropbox/Work/PPS/EU-GEI/dca_clin_results.csv")                               

ggplot(data=dcaClin2, aes(x = threshold, y = net_benefit, color = label)) +
  stat_smooth(method = "loess", se = FALSE, formula = "y ~ x", span = 0.2) +
  coord_cartesian(ylim = c(-0.005, 0.2)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits=c(0,0.5)) +
  labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
  scale_color_manual(values=c("gray80","#000000","#599ec4","#ecc363","#7bae72")) +
  theme(text = element_text(family="Roboto", face="bold", size=30),legend.title = element_text(size=23),legend.text = element_text(size = 23)) +
  theme_classic()
ggsave("~/Dropbox/Work/PPS/EU-GEI/Plots/dca_det_clin_summary_171224.png", width=20, height=15, scale=0.5)

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

##### Harmonised - Detection LR #####

# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

# Wrangle data
PPS_scored_all[,c(3:8,10:14,16:17)] <- lapply(PPS_scored_all[,c(3:8,10:14,16:17)], factor)

###### Run model fitting and internal validation ######

# Perform nested cross-validation with repeats
PPS_PRS <- PPS_scored_all %>% subset(select=c(Gender:Tobacco, PRS_resid, site, chr))

results_all <- harmonised_lr_repeated_nested_cv(
  combined_df = PPS_PRS, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

PPS_genes <- PPS_scored_all %>% subset(select=c(PRS_resid, site, chr))

results_genes <- harmonised_lr_repeated_nested_cv(
  combined_df = PPS_genes, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

PPS_PPS <- PPS_scored_all %>% subset(select=c(Gender:Tobacco, site, chr))

results_PPS <- harmonised_lr_repeated_nested_cv(
  combined_df = PPS_PPS, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)


###### Generate performance metrics summary table######

lr_results <- data.frame(model=c("PPS", "PRS", "All"),
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
                         ppv=c(paste0(round(mean(results_PPS$c_stat_nested$ppv, na.rm=TRUE)*100,1), "% (",
                                      round((mean(results_PPS$c_stat_nested$ppv, na.rm=TRUE) - 1.96*(sd(results_PPS$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS$c_stat_nested$ppv, na.rm=TRUE) + 1.96*(sd(results_PPS$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes$c_stat_nested$ppv, na.rm=TRUE)*100,1), "% (",
                                      round((mean(results_genes$c_stat_nested$ppv, na.rm=TRUE) - 1.96*(sd(results_genes$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes$c_stat_nested$ppv, na.rm=TRUE) + 1.96*(sd(results_genes$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all$c_stat_nested$ppv, na.rm=TRUE)*100,1), "% (",
                                      round((mean(results_all$c_stat_nested$ppv, na.rm=TRUE) - 1.96*(sd(results_all$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all$c_stat_nested$ppv, na.rm=TRUE) + 1.96*(sd(results_all$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%)")),
                         npv=c(paste0(round(mean(results_PPS$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_PPS$c_stat_nested$npv) - 1.96*(sd(results_PPS$c_stat_nested$npv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS$c_stat_nested$npv) + 1.96*(sd(results_PPS$c_stat_nested$npv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_genes$c_stat_nested$npv) - 1.96*(sd(results_genes$c_stat_nested$npv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes$c_stat_nested$npv) + 1.96*(sd(results_genes$c_stat_nested$npv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_all$c_stat_nested$npv) - 1.96*(sd(results_all$c_stat_nested$npv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all$c_stat_nested$npv) + 1.96*(sd(results_all$c_stat_nested$npv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%)")),
                         intercept=c(paste0(round(mean(results_PPS$calibration_slopes$intercept),2), " (",
                                            round((mean(results_PPS$calibration_slopes$intercept) - 1.96*(sd(results_PPS$calibration_slopes$intercept) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_PPS$calibration_slopes$intercept) + 1.96*(sd(results_PPS$calibration_slopes$intercept) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_genes$calibration_slopes$intercept),2), " (",
                                            round((mean(results_genes$calibration_slopes$intercept) - 1.96*(sd(results_genes$calibration_slopes$intercept) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_genes$calibration_slopes$intercept) + 1.96*(sd(results_genes$calibration_slopes$intercept) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_all$calibration_slopes$intercept),2), " (",
                                            round((mean(results_all$calibration_slopes$intercept) - 1.96*(sd(results_all$calibration_slopes$intercept) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_all$calibration_slopes$intercept) + 1.96*(sd(results_all$calibration_slopes$intercept) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                            ")")),
                         slope=c(paste0(round(mean(results_PPS$calibration_slopes$slope),2), " (",
                                        round((mean(results_PPS$calibration_slopes$slope) - 1.96*(sd(results_PPS$calibration_slopes$slope) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS$calibration_slopes$slope) + 1.96*(sd(results_PPS$calibration_slopes$slope) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes$calibration_slopes$slope),2), " (",
                                        round((mean(results_genes$calibration_slopes$slope) - 1.96*(sd(results_genes$calibration_slopes$slope) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes$calibration_slopes$slope) + 1.96*(sd(results_genes$calibration_slopes$slope) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all$calibration_slopes$slope),2), " (",
                                        round((mean(results_all$calibration_slopes$slope) - 1.96*(sd(results_all$calibration_slopes$slope) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all$calibration_slopes$slope) + 1.96*(sd(results_all$calibration_slopes$slope) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        ")")),
                         Brier=c(paste0(round(mean(results_PPS$calibration_slopes$brier),2), " (",
                                        round((mean(results_PPS$calibration_slopes$brier) - 1.96*(sd(results_PPS$calibration_slopes$brier) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS$calibration_slopes$brier) + 1.96*(sd(results_PPS$calibration_slopes$brier) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes$calibration_slopes$brier),2), " (",
                                        round((mean(results_genes$calibration_slopes$brier) - 1.96*(sd(results_genes$calibration_slopes$brier) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes$calibration_slopes$brier) + 1.96*(sd(results_genes$calibration_slopes$brier) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all$calibration_slopes$brier),2), " (",
                                        round((mean(results_all$calibration_slopes$brier) - 1.96*(sd(results_all$calibration_slopes$brier) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all$calibration_slopes$brier) + 1.96*(sd(results_all$calibration_slopes$brier) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        ")")),
                         ICI=c(paste0(round(mean(results_PPS$calibration_slopes$ICI),2), " (",
                                      round((mean(results_PPS$calibration_slopes$ICI) - 1.96*(sd(results_PPS$calibration_slopes$ICI) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_PPS$calibration_slopes$ICI) + 1.96*(sd(results_PPS$calibration_slopes$ICI) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_genes$calibration_slopes$ICI),2), " (",
                                      round((mean(results_genes$calibration_slopes$ICI) - 1.96*(sd(results_genes$calibration_slopes$ICI) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_genes$calibration_slopes$ICI) + 1.96*(sd(results_genes$calibration_slopes$ICI) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_all$calibration_slopes$ICI),2), " (",
                                      round((mean(results_all$calibration_slopes$ICI) - 1.96*(sd(results_all$calibration_slopes$ICI) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_all$calibration_slopes$ICI) + 1.96*(sd(results_all$calibration_slopes$ICI) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                      ")"))
                         
)
write_csv(lr_results, "~/Dropbox/Work/PPS/EU-GEI/lr_results_181224_harmonised.csv")

###### Generate coefficients #######

temp <- PPS_scored_all %>% select(all_of((predictors[[1]])), site)
temp[,c(1:13)] <- lapply(temp[,c(1:13)], factor)
temp_mat <- model.matrix(~.-1,temp)
temp <- as.data.frame(temp_mat)

### Compute Global Means ###
global_mean <- colMeans(temp[, 1:(ncol(temp)-2), drop = FALSE], na.rm = TRUE)

### Mean Offset Correction ###
batch_temp <- as.factor(temp$site)
temp_corrected <- temp # Start with the original data

for (b in levels(batch_temp)) {
  batch_indices <- which(batch_temp == b) # Indices for samples in batch `b`
  
  if (length(batch_indices) == 0) {
    warning(paste("Batch", b, "is empty. Skipping."))
    next
  }
  
  # Extract the subset of test for the current batch
  batch_data <- temp[batch_indices, 1:(ncol(temp)-2), drop = FALSE] # Exclude site and outcome
  
  # Compute means for each predictor in this batch
  batch_mean <- colMeans(batch_data, na.rm = TRUE)
  
  if (length(batch_mean) == 0) {
    warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
    next
  }
  
  # Compute the offset: batch mean - global mean
  offset <- batch_mean - global_mean
  
  # Subtract the offset to align the batch with the global mean
  test_corrected[batch_indices, 1:(ncol(temp)-2)] <- sweep(
    batch_data,
    1,
    offset,
    "-"
  )
}
temp_mat <- temp_corrected %>% subset(select=c(-site))

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

##### Harmonised - Prognosis Cox #####
# Run clinical model
clinical <- PPS_scored %>% subset(select=c(CAARMS, gafex01, gafex02, Transition, day_exit, site))

results_clin_surv <- harmonised_Cox_repeated_nested_cv(
  combined_df = clinical, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS model
PPS_PPS <- PPS_scored %>% subset(select=c(Gender:Tobacco, Transition, day_exit, site))

results_PPS_surv <- harmonised_Cox_repeated_nested_cv(
  combined_df = PPS_PPS, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PRS model
PPS_genes <- PPS_scored %>% subset(select=c(PRS_resid, Transition, day_exit, site))

results_PRS_surv <- harmonised_Cox_repeated_nested_cv(
  combined_df = PPS_genes, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS+clinical model
PPS_clin <- PPS_scored %>% subset(select=c(Gender:Tobacco, CAARMS, gafex01, gafex02, Transition, day_exit, site))

results_PPS_clin_surv <- harmonised_Cox_repeated_nested_cv(
  combined_df = PPS_clin, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PRS+clinical model
PRS_clin <- PPS_scored %>% subset(select=c(CAARMS, gafex01, gafex02, PRS_resid, Transition, day_exit, site))

results_PRS_clin_surv <- harmonised_Cox_repeated_nested_cv(
  combined_df = PRS_clin, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS+PRS model
PPS_PRS <- PPS_scored %>% subset(select=c(Gender:Tobacco, PRS_resid, Transition, day_exit, site))

results_PPS_PRS_surv <- harmonised_Cox_repeated_nested_cv(
  combined_df = PPS_PRS, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS+PRS+clinical model
PPS_PRS_clin <- PPS_scored %>% subset(select=c(Gender:Tobacco, CAARMS, gafex01, gafex02, PRS_resid, Transition, day_exit, site))

results_PPS_PRS_clin_surv <- harmonised_Cox_repeated_nested_cv(
  combined_df = PPS_PRS_clin, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

###### Generate performance metrics summary table######

cox_results <- data.frame(model=c("Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All"),
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
                          calibration_in_the_large=c(paste0(round(mean(results_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PRS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PRS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PRS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                                            ")")),
                          slope=c(paste0(round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                         ")")),
                          Brier=c(paste0(round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                         ")"))
                          
)
write_csv(cox_results, "~/Dropbox/Work/PPS/EU-GEI/harmonised_cox_results_131224.csv")

##### Harmonised Detection Random Forest ####
# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

###### Run model fitting and internal validation ######

# Perform nested cross-validation with repeats
PPS_PRS <- PPS_scored_all %>% subset(select=c(Gender:Tobacco, PRS_resid, chr))

# Define tuning grid
tune_grid <- expand.grid(
  mtry=c(3:(ncol(PPS_PRS)-2))
)

results_all <- harmonised_RF_repeated_nested_cv(
  combined_df = PPS_PRS, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

PPS_genes <- PPS_scored_all %>% subset(select=c(PRS_resid, chr))
tune_grid_genes <- expand.grid(
  mtry=1
)
results_genes <- harmonised_RF_repeated_nested_cv(
  combined_df = PPS_genes, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid_genes,
  seed = 231
)

PPS_PPS <- PPS_scored_all %>% subset(select=c(Gender:Tobacco, chr))
tune_grid_PPS <- expand.grid(
  mtry=c(3:(ncol(PPS_PPS)-2))
)
results_PPS <- harmonised_RF_repeated_nested_cv(
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
                         ppv=c(paste0(round(mean(results_PPS$c_stat_nested$ppv)*100,1), "% (",
                                      round((mean(results_PPS$c_stat_nested$ppv) - 1.96*(sd(results_PPS$c_stat_nested$ppv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS$c_stat_nested$ppv) + 1.96*(sd(results_PPS$c_stat_nested$ppv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes$c_stat_nested$ppv)*100,1), "% (",
                                      round((mean(results_genes$c_stat_nested$ppv) - 1.96*(sd(results_genes$c_stat_nested$ppv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes$c_stat_nested$ppv) + 1.96*(sd(results_genes$c_stat_nested$ppv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all$c_stat_nested$ppv)*100,1), "% (",
                                      round((mean(results_all$c_stat_nested$ppv) - 1.96*(sd(results_all$c_stat_nested$ppv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all$c_stat_nested$ppv) + 1.96*(sd(results_all$c_stat_nested$ppv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%)")),
                         npv=c(paste0(round(mean(results_PPS$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_PPS$c_stat_nested$npv) - 1.96*(sd(results_PPS$c_stat_nested$npv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS$c_stat_nested$npv) + 1.96*(sd(results_PPS$c_stat_nested$npv) / sqrt(nrow(results_PPS$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_genes$c_stat_nested$npv) - 1.96*(sd(results_genes$c_stat_nested$npv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes$c_stat_nested$npv) + 1.96*(sd(results_genes$c_stat_nested$npv) / sqrt(nrow(results_genes$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_all$c_stat_nested$npv) - 1.96*(sd(results_all$c_stat_nested$npv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all$c_stat_nested$npv) + 1.96*(sd(results_all$c_stat_nested$npv) / sqrt(nrow(results_all$c_stat_nested))))*100,1),
                                      "%)")),
                         intercept=c(paste0(round(mean(results_PPS$calibration_slopes$intercept),2), " (",
                                            round((mean(results_PPS$calibration_slopes$intercept) - 1.96*(sd(results_PPS$calibration_slopes$intercept) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_PPS$calibration_slopes$intercept) + 1.96*(sd(results_PPS$calibration_slopes$intercept) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_genes$calibration_slopes$intercept),2), " (",
                                            round((mean(results_genes$calibration_slopes$intercept) - 1.96*(sd(results_genes$calibration_slopes$intercept) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_genes$calibration_slopes$intercept) + 1.96*(sd(results_genes$calibration_slopes$intercept) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_all$calibration_slopes$intercept),2), " (",
                                            round((mean(results_all$calibration_slopes$intercept) - 1.96*(sd(results_all$calibration_slopes$intercept) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_all$calibration_slopes$intercept) + 1.96*(sd(results_all$calibration_slopes$intercept) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                            ")")),
                         slope=c(paste0(round(mean(results_PPS$calibration_slopes$slope),2), " (",
                                        round((mean(results_PPS$calibration_slopes$slope) - 1.96*(sd(results_PPS$calibration_slopes$slope) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS$calibration_slopes$slope) + 1.96*(sd(results_PPS$calibration_slopes$slope) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes$calibration_slopes$slope),2), " (",
                                        round((mean(results_genes$calibration_slopes$slope) - 1.96*(sd(results_genes$calibration_slopes$slope) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes$calibration_slopes$slope) + 1.96*(sd(results_genes$calibration_slopes$slope) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all$calibration_slopes$slope),2), " (",
                                        round((mean(results_all$calibration_slopes$slope) - 1.96*(sd(results_all$calibration_slopes$slope) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all$calibration_slopes$slope) + 1.96*(sd(results_all$calibration_slopes$slope) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        ")")),
                         Brier=c(paste0(round(mean(results_PPS$calibration_slopes$brier),2), " (",
                                        round((mean(results_PPS$calibration_slopes$brier) - 1.96*(sd(results_PPS$calibration_slopes$brier) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS$calibration_slopes$brier) + 1.96*(sd(results_PPS$calibration_slopes$brier) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes$calibration_slopes$brier),2), " (",
                                        round((mean(results_genes$calibration_slopes$brier) - 1.96*(sd(results_genes$calibration_slopes$brier) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes$calibration_slopes$brier) + 1.96*(sd(results_genes$calibration_slopes$brier) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all$calibration_slopes$brier),2), " (",
                                        round((mean(results_all$calibration_slopes$brier) - 1.96*(sd(results_all$calibration_slopes$brier) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all$calibration_slopes$brier) + 1.96*(sd(results_all$calibration_slopes$brier) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                        ")")),
                         ICI=c(paste0(round(mean(results_PPS$calibration_slopes$ICI),2), " (",
                                      round((mean(results_PPS$calibration_slopes$ICI) - 1.96*(sd(results_PPS$calibration_slopes$ICI) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_PPS$calibration_slopes$ICI) + 1.96*(sd(results_PPS$calibration_slopes$ICI) / sqrt(nrow(results_PPS$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_genes$calibration_slopes$ICI),2), " (",
                                      round((mean(results_genes$calibration_slopes$ICI) - 1.96*(sd(results_genes$calibration_slopes$ICI) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_genes$calibration_slopes$ICI) + 1.96*(sd(results_genes$calibration_slopes$ICI) / sqrt(nrow(results_genes$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_all$calibration_slopes$ICI),2), " (",
                                      round((mean(results_all$calibration_slopes$ICI) - 1.96*(sd(results_all$calibration_slopes$ICI) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_all$calibration_slopes$ICI) + 1.96*(sd(results_all$calibration_slopes$ICI) / sqrt(nrow(results_all$calibration_slopes)))),2),
                                      ")"))
)
write_csv(rf_results, "~/Dropbox/Work/PPS/EU-GEI/harmonised_rf_results_081224.csv")

##### Harmonised Prognosis Random Survival Forest ####
# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

###### Run model fitting and internal validation ######
# Clinical
clin_rf_surv <- PPS_scored %>% subset(select=c(CAARMS, gafex01, gafex02, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3,3,3)
)

results_clin_surv <- harmonised_RSF_repeated_nested_cv(
  combined_df = clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PPS

PPS_rf_surv <- PPS_scored %>% subset(select=c(Gender:Tobacco, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_rf_surv)-2)
)

results_PPS_surv <- harmonised_RSF_repeated_nested_cv(
  combined_df = PPS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PRS
PRS_rf_surv <- PPS_scored %>% subset(select=c(PRS_resid, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(1)
)

results_PRS_surv <- harmonised_RSF_repeated_nested_cv(
  combined_df = PRS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PPS

PPS_clin_rf_surv <- PPS_scored %>% subset(select=c(Gender:Tobacco, CAARMS, gafex01, gafex02, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_clin_rf_surv)-2)
)

results_PPS_clin_surv <- harmonised_RSF_repeated_nested_cv(
  combined_df = PPS_clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PRS

PRS_clin_rf_surv <- PPS_scored %>% subset(select=c(PRS_resid, CAARMS, gafex01, gafex02, Transition, day_exit))
PRS_clin_rf_surv <- as.data.frame(model.matrix(~.-1,PRS_clin_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_rf_surv)-2)
)

results_PRS_clin_surv <- harmonised_RSF_repeated_nested_cv(
  combined_df = PRS_clin_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PPS + PRS

PPS_PRS_rf_surv <- PPS_scored %>% subset(select=c(Gender:Tobacco, PRS_resid, Transition, day_exit))
PPS_PRS_rf_surv <- as.data.frame(model.matrix(~.-1,PPS_PRS_rf_surv))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_PRS_rf_surv)-2)
)

results_PPS_PRS_surv <- harmonised_RSF_repeated_nested_cv(
  combined_df = PPS_PRS_rf_surv, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PPS + PRS 

PPS_PRS_clin_rf_surv <- PPS_scored %>% subset(select=c(Gender:Tobacco, PRS_resid, CAARMS, gafex01, gafex02, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_PRS_clin_rf_surv)-2)
)

results_PPS_PRS_clin_surv <- harmonised_RSF_repeated_nested_cv(
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
                              calibration_in_the_large=c(paste0(round(mean(results_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PPS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PPS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PPS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PRS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PRS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PRS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PPS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PRS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                                                ")"),
                                                         paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large),2), " (",
                                                                round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                                                "-",
                                                                round(mean(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                                                ")")),
                              slope=c(paste0(round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                             ")")),
                              Brier=c(paste0(round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                             ")"),
                                      paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                             round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                             "-",
                                             round(mean(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                             ")")),
                              ICI=c(paste0(round(mean(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                           round(mean(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                           "-",
                                           round(mean(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_clin_surv$calibration_slopes))),2),
                                           ")"),
                                    paste0(round(mean(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                           round(mean(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                           "-",
                                           round(mean(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_surv$calibration_slopes))),2),
                                           ")"),
                                    paste0(round(mean(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                           round(mean(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                           "-",
                                           round(mean(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_surv$calibration_slopes))),2),
                                           ")"),
                                    paste0(round(mean(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                           round(mean(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                           "-",
                                           round(mean(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv$calibration_slopes))),2),
                                           ")"),
                                    paste0(round(mean(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                           round(mean(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                           "-",
                                           round(mean(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv$calibration_slopes))),2),
                                           ")"),
                                    paste0(round(mean(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                           round(mean(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                           "-",
                                           round(mean(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv$calibration_slopes))),2),
                                           ")"),
                                    paste0(round(mean(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                           round(mean(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                           "-",
                                           round(mean(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv$calibration_slopes))),2),
                                           ")"))
                              
)
write_csv(rf_surv_results, "~/Dropbox/Work/PPS/EU-GEI/harmonised_rf_surv_results_081224.csv")

##### White Ancestry only #####
PPS_scored_WA <- read_csv("~/Dropbox/Work/PPS/EU-GEI/Databases/PPS_SA_PRS.csv")
PPS_scored_all_WA <- PPS_scored_WA
PPS_scored_WA <- PPS_scored_WA %>% filter(chr==1)
PPS_scored_WA <- PPS_scored_WA %>% mutate(Transition = case_when(day_exit>1095 ~ 0,
                                                           TRUE ~ Transition),
                                    day_exit = case_when(day_exit>1095 ~ 1095,
                                                         TRUE ~ day_exit))

##### Detection LR ####

# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

###### Run model fitting and internal validation ######

# Run PPS model
PPS_PPS_WA <- PPS_scored_all_WA %>% subset(select=c(Gender:Tobacco, chr))

results_PPS_WA <- LR_repeated_nested_cv(
  combined_df = PPS_PPS_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PRS model
PPS_genes_WA <- PPS_scored_all_WA %>% subset(select=c(PRS_resid, chr))

results_genes_WA <- LR_repeated_nested_cv(
  combined_df = PPS_genes_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run combined PPS and PRS model
PPS_PRS_WA <- PPS_scored_all_WA %>% subset(select=c(Gender:Tobacco, PRS_resid, chr))

results_all_WA <- LR_repeated_nested_cv(
  combined_df = PPS_PRS_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)
###### Generate performance metrics summary table######

lr_results_WA <- data.frame(model=c("PPS", "PRS", "All"),
                         C=c(paste0(round(mean(results_PPS_WA$c_stat_nested$C),2), " (",
                                    round(mean(results_PPS_WA$c_stat_nested$C) - 1.96*(sd(results_PPS_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_WA$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_PPS_WA$c_stat_nested$C) + 1.96*(sd(results_PPS_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_WA$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_genes_WA$c_stat_nested$C),2), " (",
                                    round(mean(results_genes_WA$c_stat_nested$C) - 1.96*(sd(results_genes_WA$c_stat_nested$C) / sqrt(nrow(results_genes_WA$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_genes_WA$c_stat_nested$C) + 1.96*(sd(results_genes_WA$c_stat_nested$C) / sqrt(nrow(results_genes_WA$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_all_WA$c_stat_nested$C),2), " (",
                                    round(mean(results_all_WA$c_stat_nested$C) - 1.96*(sd(results_all_WA$c_stat_nested$C) / sqrt(nrow(results_all_WA$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_all_WA$c_stat_nested$C) + 1.96*(sd(results_all_WA$c_stat_nested$C) / sqrt(nrow(results_all_WA$c_stat_nested))),2),
                                    ")")),
                         BAC=c(paste0(round(mean(results_PPS_WA$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                      round((mean(results_PPS_WA$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_PPS_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS_WA$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_PPS_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes_WA$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                      round((mean(results_genes_WA$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_genes_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes_WA$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_genes_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all_WA$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                      round((mean(results_all_WA$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_all_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all_WA$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_all_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%)")),
                         sensitivity=c(paste0(round(mean(results_PPS_WA$c_stat_nested$sensitivity)*100,1), "% (",
                                              round((mean(results_PPS_WA$c_stat_nested$sensitivity) - 1.96*(sd(results_PPS_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_PPS_WA$c_stat_nested$sensitivity) + 1.96*(sd(results_PPS_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                              "%)"),
                                       paste0(round(mean(results_genes_WA$c_stat_nested$sensitivity)*100,1), "% (",
                                              round((mean(results_genes_WA$c_stat_nested$sensitivity) - 1.96*(sd(results_genes_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_genes_WA$c_stat_nested$sensitivity) + 1.96*(sd(results_genes_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                              "%)"),
                                       paste0(round(mean(results_all_WA$c_stat_nested$sensitivity)*100,1), "% (",
                                              round((mean(results_all_WA$c_stat_nested$sensitivity) - 1.96*(sd(results_all_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_all_WA$c_stat_nested$sensitivity) + 1.96*(sd(results_all_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                              "%)")),
                         specificity=c(paste0(round(mean(results_PPS_WA$c_stat_nested$specificity)*100,1), "% (",
                                              round((mean(results_PPS_WA$c_stat_nested$specificity) - 1.96*(sd(results_PPS_WA$c_stat_nested$specificity) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_PPS_WA$c_stat_nested$specificity) + 1.96*(sd(results_PPS_WA$c_stat_nested$specificity) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                              "%)"),
                                       paste0(round(mean(results_genes_WA$c_stat_nested$specificity)*100,1), "% (",
                                              round((mean(results_genes_WA$c_stat_nested$specificity) - 1.96*(sd(results_genes_WA$c_stat_nested$specificity) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_genes_WA$c_stat_nested$specificity) + 1.96*(sd(results_genes_WA$c_stat_nested$specificity) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                              "%)"),
                                       paste0(round(mean(results_all_WA$c_stat_nested$specificity)*100,1), "% (",
                                              round((mean(results_all_WA$c_stat_nested$specificity) - 1.96*(sd(results_all_WA$c_stat_nested$specificity) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_all_WA$c_stat_nested$specificity) + 1.96*(sd(results_all_WA$c_stat_nested$specificity) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                              "%)")),
                         ppv=c(paste0(round(mean(results_PPS_WA$c_stat_nested$ppv, na.rm=TRUE)*100,1), "% (",
                                      round((mean(results_PPS_WA$c_stat_nested$ppv, na.rm=TRUE) - 1.96*(sd(results_PPS_WA$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS_WA$c_stat_nested$ppv, na.rm=TRUE) + 1.96*(sd(results_PPS_WA$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes_WA$c_stat_nested$ppv, na.rm=TRUE)*100,1), "% (",
                                      round((mean(results_genes_WA$c_stat_nested$ppv, na.rm=TRUE) - 1.96*(sd(results_genes_WA$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes_WA$c_stat_nested$ppv, na.rm=TRUE) + 1.96*(sd(results_genes_WA$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all_WA$c_stat_nested$ppv, na.rm=TRUE)*100,1), "% (",
                                      round((mean(results_all_WA$c_stat_nested$ppv, na.rm=TRUE) - 1.96*(sd(results_all_WA$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all_WA$c_stat_nested$ppv, na.rm=TRUE) + 1.96*(sd(results_all_WA$c_stat_nested$ppv, na.rm=TRUE) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%)")),
                         npv=c(paste0(round(mean(results_PPS_WA$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_PPS_WA$c_stat_nested$npv) - 1.96*(sd(results_PPS_WA$c_stat_nested$npv) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS_WA$c_stat_nested$npv) + 1.96*(sd(results_PPS_WA$c_stat_nested$npv) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes_WA$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_genes_WA$c_stat_nested$npv) - 1.96*(sd(results_genes_WA$c_stat_nested$npv) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes_WA$c_stat_nested$npv) + 1.96*(sd(results_genes_WA$c_stat_nested$npv) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all_WA$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_all_WA$c_stat_nested$npv) - 1.96*(sd(results_all_WA$c_stat_nested$npv) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all_WA$c_stat_nested$npv) + 1.96*(sd(results_all_WA$c_stat_nested$npv) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%)")),
                         intercept=c(paste0(round(mean(results_PPS_WA$calibration_slopes$intercept),2), " (",
                                            round((mean(results_PPS_WA$calibration_slopes$intercept) - 1.96*(sd(results_PPS_WA$calibration_slopes$intercept) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_PPS_WA$calibration_slopes$intercept) + 1.96*(sd(results_PPS_WA$calibration_slopes$intercept) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_genes_WA$calibration_slopes$intercept),2), " (",
                                            round((mean(results_genes_WA$calibration_slopes$intercept) - 1.96*(sd(results_genes_WA$calibration_slopes$intercept) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_genes_WA$calibration_slopes$intercept) + 1.96*(sd(results_genes_WA$calibration_slopes$intercept) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_all_WA$calibration_slopes$intercept),2), " (",
                                            round((mean(results_all_WA$calibration_slopes$intercept) - 1.96*(sd(results_all_WA$calibration_slopes$intercept) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_all_WA$calibration_slopes$intercept) + 1.96*(sd(results_all_WA$calibration_slopes$intercept) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                            ")")),
                         slope=c(paste0(round(mean(results_PPS_WA$calibration_slopes$slope),2), " (",
                                        round((mean(results_PPS_WA$calibration_slopes$slope) - 1.96*(sd(results_PPS_WA$calibration_slopes$slope) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS_WA$calibration_slopes$slope) + 1.96*(sd(results_PPS_WA$calibration_slopes$slope) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes_WA$calibration_slopes$slope),2), " (",
                                        round((mean(results_genes_WA$calibration_slopes$slope) - 1.96*(sd(results_genes_WA$calibration_slopes$slope) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes_WA$calibration_slopes$slope) + 1.96*(sd(results_genes_WA$calibration_slopes$slope) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all_WA$calibration_slopes$slope),2), " (",
                                        round((mean(results_all_WA$calibration_slopes$slope) - 1.96*(sd(results_all_WA$calibration_slopes$slope) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all_WA$calibration_slopes$slope) + 1.96*(sd(results_all_WA$calibration_slopes$slope) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                        ")")),
                         Brier=c(paste0(round(mean(results_PPS_WA$calibration_slopes$brier),2), " (",
                                        round((mean(results_PPS_WA$calibration_slopes$brier) - 1.96*(sd(results_PPS_WA$calibration_slopes$brier) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS_WA$calibration_slopes$brier) + 1.96*(sd(results_PPS_WA$calibration_slopes$brier) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes_WA$calibration_slopes$brier),2), " (",
                                        round((mean(results_genes_WA$calibration_slopes$brier) - 1.96*(sd(results_genes_WA$calibration_slopes$brier) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes_WA$calibration_slopes$brier) + 1.96*(sd(results_genes_WA$calibration_slopes$brier) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all_WA$calibration_slopes$brier),2), " (",
                                        round((mean(results_all_WA$calibration_slopes$brier) - 1.96*(sd(results_all_WA$calibration_slopes$brier) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all_WA$calibration_slopes$brier) + 1.96*(sd(results_all_WA$calibration_slopes$brier) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                        ")")),
                         ICI=c(paste0(round(mean(results_PPS_WA$calibration_slopes$ICI),2), " (",
                                      round((mean(results_PPS_WA$calibration_slopes$ICI) - 1.96*(sd(results_PPS_WA$calibration_slopes$ICI) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_PPS_WA$calibration_slopes$ICI) + 1.96*(sd(results_PPS_WA$calibration_slopes$ICI) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_genes_WA$calibration_slopes$ICI),2), " (",
                                      round((mean(results_genes_WA$calibration_slopes$ICI) - 1.96*(sd(results_genes_WA$calibration_slopes$ICI) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_genes_WA$calibration_slopes$ICI) + 1.96*(sd(results_genes_WA$calibration_slopes$ICI) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_all_WA$calibration_slopes$ICI),2), " (",
                                      round((mean(results_all_WA$calibration_slopes$ICI) - 1.96*(sd(results_all_WA$calibration_slopes$ICI) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_all_WA$calibration_slopes$ICI) + 1.96*(sd(results_all_WA$calibration_slopes$ICI) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                      ")"))
                         
)
write_csv(lr_results_WA, "~/Dropbox/Work/PPS/EU-GEI/lr_results_WA_181224.csv")

#### Prognosis Cox ####

# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

###### Run model fitting and internal validation ######

# Run clinical_WA model
clinical_WA <- PPS_scored_WA %>% subset(select=c(CAARMS, gafex01, gafex02, Transition, day_exit))

results_clin_surv_WA <- Cox_repeated_nested_cv(
  combined_df = clinical_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS model
PPS_PPS_WA <- PPS_scored_WA %>% subset(select=c(Gender:Tobacco, Transition, day_exit))

results_PPS_surv_WA <- Cox_repeated_nested_cv(
  combined_df = PPS_PPS_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PRS model
PPS_genes_WA <- PPS_scored_WA %>% subset(select=c(PRS_resid, Transition, day_exit))

results_PRS_surv_WA <- Cox_repeated_nested_cv(
  combined_df = PPS_genes_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS+clinical_WA model
PPS_clin_WA <- PPS_scored_WA %>% subset(select=c(Gender:Tobacco, CAARMS, gafex01, gafex02, Transition, day_exit))

results_PPS_clin_surv_WA <- Cox_repeated_nested_cv(
  combined_df = PPS_clin_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PRS+clinical_WA model
PRS_clin_WA <- PPS_scored_WA %>% subset(select=c(CAARMS, gafex01, gafex02, PRS_resid, Transition, day_exit))

results_PRS_clin_surv_WA <- Cox_repeated_nested_cv(
  combined_df = PRS_clin_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS+PRS model
PPS_PRS_WA <- PPS_scored_WA %>% subset(select=c(Gender:Tobacco, PRS_resid, Transition, day_exit))

results_PPS_PRS_surv_WA <- Cox_repeated_nested_cv(
  combined_df = PPS_PRS_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

# Run PPS+PRS+clinical_WA model
PPS_PRS_clin_WA <- PPS_scored_WA %>% subset(select=c(Gender:Tobacco, CAARMS, gafex01, gafex02, PRS_resid, Transition, day_exit))

results_PPS_PRS_clin_WA_surv <- Cox_repeated_nested_cv(
  combined_df = PPS_PRS_clin_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  seed = 231
)

###### Generate performance metrics summary table######

cox_results <- data.frame(model=c("clinical_WA", "PPS", "PRS", "clinical_WA+PPS", "clinical_WA+PRS", "PPS+PRS", "All"),
                          C=c(paste0(round(mean(results_clin_surv_WA$c_stat_nested$C),2), " (",
                                     round(mean(results_clin_surv_WA$c_stat_nested$C) - 1.96*(sd(results_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_clin_surv_WA$c_stat_nested))),2),
                                     "-",
                                     round(mean(results_clin_surv_WA$c_stat_nested$C) + 1.96*(sd(results_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_clin_surv_WA$c_stat_nested))),2),
                                     ")"),
                              paste0(round(mean(results_PPS_surv_WA$c_stat_nested$C),2), " (",
                                     round(mean(results_PPS_surv_WA$c_stat_nested$C) - 1.96*(sd(results_PPS_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_surv_WA$c_stat_nested))),2),
                                     "-",
                                     round(mean(results_PPS_surv_WA$c_stat_nested$C) + 1.96*(sd(results_PPS_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_surv_WA$c_stat_nested))),2),
                                     ")"),
                              paste0(round(mean(results_PRS_surv_WA$c_stat_nested$C),2), " (",
                                     round(mean(results_PRS_surv_WA$c_stat_nested$C) - 1.96*(sd(results_PRS_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PRS_surv_WA$c_stat_nested))),2),
                                     "-",
                                     round(mean(results_PRS_surv_WA$c_stat_nested$C) + 1.96*(sd(results_PRS_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PRS_surv_WA$c_stat_nested))),2),
                                     ")"),
                              paste0(round(mean(results_PPS_clin_surv_WA$c_stat_nested$C),2), " (",
                                     round(mean(results_PPS_clin_surv_WA$c_stat_nested$C) - 1.96*(sd(results_PPS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_clin_surv_WA$c_stat_nested))),2),
                                     "-",
                                     round(mean(results_PPS_clin_surv_WA$c_stat_nested$C) + 1.96*(sd(results_PPS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_clin_surv_WA$c_stat_nested))),2),
                                     ")"),
                              paste0(round(mean(results_PRS_clin_surv_WA$c_stat_nested$C),2), " (",
                                     round(mean(results_PRS_clin_surv_WA$c_stat_nested$C) - 1.96*(sd(results_PRS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PRS_clin_surv_WA$c_stat_nested))),2),
                                     "-",
                                     round(mean(results_PRS_clin_surv_WA$c_stat_nested$C) + 1.96*(sd(results_PRS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PRS_clin_surv_WA$c_stat_nested))),2),
                                     ")"),
                              paste0(round(mean(results_PPS_PRS_clin_WA_surv$c_stat_nested$C),2), " (",
                                     round(mean(results_PPS_PRS_clin_WA_surv$c_stat_nested$C) - 1.96*(sd(results_PPS_PRS_clin_WA_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$c_stat_nested))),2),
                                     "-",
                                     round(mean(results_PPS_PRS_clin_WA_surv$c_stat_nested$C) + 1.96*(sd(results_PPS_PRS_clin_WA_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$c_stat_nested))),2),
                                     ")"),
                              paste0(round(mean(results_PPS_PRS_clin_WA_surv$c_stat_nested$C),2), " (",
                                     round(mean(results_PPS_PRS_clin_WA_surv$c_stat_nested$C) - 1.96*(sd(results_PPS_PRS_clin_WA_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$c_stat_nested))),2),
                                     "-",
                                     round(mean(results_PPS_PRS_clin_WA_surv$c_stat_nested$C) + 1.96*(sd(results_PPS_PRS_clin_WA_surv$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$c_stat_nested))),2),
                                     ")")),
                          calibration_in_the_large=c(paste0(round(mean(results_clin_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_clin_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_clin_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PRS_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PRS_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PRS_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_clin_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_clin_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_clin_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PRS_clin_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_PRS_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_PRS_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_PRS_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                                            ")"),
                                                     paste0(round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$calibration_in_large),2), " (",
                                                            round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_clin_WA_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$calibration_slopes))),2),
                                                            "-",
                                                            round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_clin_WA_surv$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$calibration_slopes))),2),
                                                            ")")),
                          slope=c(paste0(round(mean(results_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_WA_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_WA_surv$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$calibration_slopes))),2),
                                         ")")),
                          Brier=c(paste0(round(mean(results_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                         ")"),
                                  paste0(round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                         round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_WA_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$calibration_slopes))),2),
                                         "-",
                                         round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_WA_surv$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$calibration_slopes))),2),
                                         ")")),
                          ICI=c(paste0(round(mean(results_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PPS_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PPS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PPS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PPS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PPS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PPS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PPS_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PPS_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PPS_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                       ")"),
                                paste0(round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                       round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_WA_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$calibration_slopes))),2),
                                       "-",
                                       round(mean(results_PPS_PRS_clin_WA_surv$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_WA_surv$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_WA_surv$calibration_slopes))),2),
                                       ")"))
                          
)
write_csv(cox_results, "~/Dropbox/Work/PPS/EU-GEI/cox_results_WA_031224.csv")

##### Detection Random Forest ####
# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

###### Run model fitting and internal validation ######

# Perform nested cross-validation with repeats
PPS_PRS_WA <- PPS_scored_all_WA %>% subset(select=c(Gender:Tobacco, PRS_resid, chr))

# Define tuning grid
tune_grid <- expand.grid(
  mtry=c(3:(ncol(PPS_PRS_WA)-2))
)

results_all_WA <- RF_repeated_nested_cv(
  combined_df = PPS_PRS_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

PPS_genes_WA <- PPS_scored_all_WA %>% subset(select=c(PRS_resid, chr))
tune_grid_genes <- expand.grid(
  mtry=1
)
results_genes_WA <- RF_repeated_nested_cv(
  combined_df = PPS_genes_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid_genes,
  seed = 231
)

PPS_PPS_WA <- PPS_scored_all_WA %>% subset(select=c(Gender:Tobacco, chr))
tune_grid_PPS <- expand.grid(
  mtry=c(3:(ncol(PPS_PPS_WA)-2))
)
results_PPS_WA <- RF_repeated_nested_cv(
  combined_df = PPS_PPS_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid_PPS,
  seed = 231
)

###### Generate performance metrics summary table######

rf_results_WA <- data.frame(model=c("PPS", "PRS", "All"),
                         C=c(paste0(round(mean(results_PPS_WA$c_stat_nested$C),2), " (",
                                    round(mean(results_PPS_WA$c_stat_nested$C) - 1.96*(sd(results_PPS_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_WA$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_PPS_WA$c_stat_nested$C) + 1.96*(sd(results_PPS_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_WA$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_genes_WA$c_stat_nested$C),2), " (",
                                    round(mean(results_genes_WA$c_stat_nested$C) - 1.96*(sd(results_genes_WA$c_stat_nested$C) / sqrt(nrow(results_genes_WA$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_genes_WA$c_stat_nested$C) + 1.96*(sd(results_genes_WA$c_stat_nested$C) / sqrt(nrow(results_genes_WA$c_stat_nested))),2),
                                    ")"),
                             paste0(round(mean(results_all_WA$c_stat_nested$C),2), " (",
                                    round(mean(results_all_WA$c_stat_nested$C) - 1.96*(sd(results_all_WA$c_stat_nested$C) / sqrt(nrow(results_all_WA$c_stat_nested))),2),
                                    "-",
                                    round(mean(results_all_WA$c_stat_nested$C) + 1.96*(sd(results_all_WA$c_stat_nested$C) / sqrt(nrow(results_all_WA$c_stat_nested))),2),
                                    ")")),
                         BAC=c(paste0(round(mean(results_PPS_WA$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                      round((mean(results_PPS_WA$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_PPS_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS_WA$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_PPS_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes_WA$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                      round((mean(results_genes_WA$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_genes_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes_WA$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_genes_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all_WA$c_stat_nested$balanced_accuracy)*100,1), "% (",
                                      round((mean(results_all_WA$c_stat_nested$balanced_accuracy) - 1.96*(sd(results_all_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all_WA$c_stat_nested$balanced_accuracy) + 1.96*(sd(results_all_WA$c_stat_nested$balanced_accuracy) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%)")),
                         sensitivity=c(paste0(round(mean(results_PPS_WA$c_stat_nested$sensitivity)*100,1), "% (",
                                              round((mean(results_PPS_WA$c_stat_nested$sensitivity) - 1.96*(sd(results_PPS_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_PPS_WA$c_stat_nested$sensitivity) + 1.96*(sd(results_PPS_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                              "%)"),
                                       paste0(round(mean(results_genes_WA$c_stat_nested$sensitivity)*100,1), "% (",
                                              round((mean(results_genes_WA$c_stat_nested$sensitivity) - 1.96*(sd(results_genes_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_genes_WA$c_stat_nested$sensitivity) + 1.96*(sd(results_genes_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                              "%)"),
                                       paste0(round(mean(results_all_WA$c_stat_nested$sensitivity)*100,1), "% (",
                                              round((mean(results_all_WA$c_stat_nested$sensitivity) - 1.96*(sd(results_all_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_all_WA$c_stat_nested$sensitivity) + 1.96*(sd(results_all_WA$c_stat_nested$sensitivity) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                              "%)")),
                         specificity=c(paste0(round(mean(results_PPS_WA$c_stat_nested$specificity)*100,1), "% (",
                                              round((mean(results_PPS_WA$c_stat_nested$specificity) - 1.96*(sd(results_PPS_WA$c_stat_nested$specificity) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_PPS_WA$c_stat_nested$specificity) + 1.96*(sd(results_PPS_WA$c_stat_nested$specificity) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                              "%)"),
                                       paste0(round(mean(results_genes_WA$c_stat_nested$specificity)*100,1), "% (",
                                              round((mean(results_genes_WA$c_stat_nested$specificity) - 1.96*(sd(results_genes_WA$c_stat_nested$specificity) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_genes_WA$c_stat_nested$specificity) + 1.96*(sd(results_genes_WA$c_stat_nested$specificity) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                              "%)"),
                                       paste0(round(mean(results_all_WA$c_stat_nested$specificity)*100,1), "% (",
                                              round((mean(results_all_WA$c_stat_nested$specificity) - 1.96*(sd(results_all_WA$c_stat_nested$specificity) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                              "%-",
                                              round((mean(results_all_WA$c_stat_nested$specificity) + 1.96*(sd(results_all_WA$c_stat_nested$specificity) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                              "%)")),
                         ppv=c(paste0(round(mean(results_PPS_WA$c_stat_nested$ppv)*100,1), "% (",
                                      round((mean(results_PPS_WA$c_stat_nested$ppv) - 1.96*(sd(results_PPS_WA$c_stat_nested$ppv) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS_WA$c_stat_nested$ppv) + 1.96*(sd(results_PPS_WA$c_stat_nested$ppv) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes_WA$c_stat_nested$ppv)*100,1), "% (",
                                      round((mean(results_genes_WA$c_stat_nested$ppv) - 1.96*(sd(results_genes_WA$c_stat_nested$ppv) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes_WA$c_stat_nested$ppv) + 1.96*(sd(results_genes_WA$c_stat_nested$ppv) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all_WA$c_stat_nested$ppv)*100,1), "% (",
                                      round((mean(results_all_WA$c_stat_nested$ppv) - 1.96*(sd(results_all_WA$c_stat_nested$ppv) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all_WA$c_stat_nested$ppv) + 1.96*(sd(results_all_WA$c_stat_nested$ppv) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%)")),
                         npv=c(paste0(round(mean(results_PPS_WA$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_PPS_WA$c_stat_nested$npv) - 1.96*(sd(results_PPS_WA$c_stat_nested$npv) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_PPS_WA$c_stat_nested$npv) + 1.96*(sd(results_PPS_WA$c_stat_nested$npv) / sqrt(nrow(results_PPS_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_genes_WA$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_genes_WA$c_stat_nested$npv) - 1.96*(sd(results_genes_WA$c_stat_nested$npv) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_genes_WA$c_stat_nested$npv) + 1.96*(sd(results_genes_WA$c_stat_nested$npv) / sqrt(nrow(results_genes_WA$c_stat_nested))))*100,1),
                                      "%)"),
                               paste0(round(mean(results_all_WA$c_stat_nested$npv)*100,1), "% (",
                                      round((mean(results_all_WA$c_stat_nested$npv) - 1.96*(sd(results_all_WA$c_stat_nested$npv) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%-",
                                      round((mean(results_all_WA$c_stat_nested$npv) + 1.96*(sd(results_all_WA$c_stat_nested$npv) / sqrt(nrow(results_all_WA$c_stat_nested))))*100,1),
                                      "%)")),
                         intercept=c(paste0(round(mean(results_PPS_WA$calibration_slopes$intercept),2), " (",
                                            round((mean(results_PPS_WA$calibration_slopes$intercept) - 1.96*(sd(results_PPS_WA$calibration_slopes$intercept) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_PPS_WA$calibration_slopes$intercept) + 1.96*(sd(results_PPS_WA$calibration_slopes$intercept) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_genes_WA$calibration_slopes$intercept),2), " (",
                                            round((mean(results_genes_WA$calibration_slopes$intercept) - 1.96*(sd(results_genes_WA$calibration_slopes$intercept) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_genes_WA$calibration_slopes$intercept) + 1.96*(sd(results_genes_WA$calibration_slopes$intercept) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                            ")"),
                                     paste0(round(mean(results_all_WA$calibration_slopes$intercept),2), " (",
                                            round((mean(results_all_WA$calibration_slopes$intercept) - 1.96*(sd(results_all_WA$calibration_slopes$intercept) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                            "-",
                                            round((mean(results_all_WA$calibration_slopes$intercept) + 1.96*(sd(results_all_WA$calibration_slopes$intercept) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                            ")")),
                         slope=c(paste0(round(mean(results_PPS_WA$calibration_slopes$slope),2), " (",
                                        round((mean(results_PPS_WA$calibration_slopes$slope) - 1.96*(sd(results_PPS_WA$calibration_slopes$slope) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS_WA$calibration_slopes$slope) + 1.96*(sd(results_PPS_WA$calibration_slopes$slope) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes_WA$calibration_slopes$slope),2), " (",
                                        round((mean(results_genes_WA$calibration_slopes$slope) - 1.96*(sd(results_genes_WA$calibration_slopes$slope) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes_WA$calibration_slopes$slope) + 1.96*(sd(results_genes_WA$calibration_slopes$slope) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all_WA$calibration_slopes$slope),2), " (",
                                        round((mean(results_all_WA$calibration_slopes$slope) - 1.96*(sd(results_all_WA$calibration_slopes$slope) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all_WA$calibration_slopes$slope) + 1.96*(sd(results_all_WA$calibration_slopes$slope) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                        ")")),
                         Brier=c(paste0(round(mean(results_PPS_WA$calibration_slopes$brier),2), " (",
                                        round((mean(results_PPS_WA$calibration_slopes$brier) - 1.96*(sd(results_PPS_WA$calibration_slopes$brier) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_PPS_WA$calibration_slopes$brier) + 1.96*(sd(results_PPS_WA$calibration_slopes$brier) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_genes_WA$calibration_slopes$brier),2), " (",
                                        round((mean(results_genes_WA$calibration_slopes$brier) - 1.96*(sd(results_genes_WA$calibration_slopes$brier) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_genes_WA$calibration_slopes$brier) + 1.96*(sd(results_genes_WA$calibration_slopes$brier) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                        ")"),
                                 paste0(round(mean(results_all_WA$calibration_slopes$brier),2), " (",
                                        round((mean(results_all_WA$calibration_slopes$brier) - 1.96*(sd(results_all_WA$calibration_slopes$brier) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                        "-",
                                        round((mean(results_all_WA$calibration_slopes$brier) + 1.96*(sd(results_all_WA$calibration_slopes$brier) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                        ")")),
                         ICI=c(paste0(round(mean(results_PPS_WA$calibration_slopes$ICI),2), " (",
                                      round((mean(results_PPS_WA$calibration_slopes$ICI) - 1.96*(sd(results_PPS_WA$calibration_slopes$ICI) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_PPS_WA$calibration_slopes$ICI) + 1.96*(sd(results_PPS_WA$calibration_slopes$ICI) / sqrt(nrow(results_PPS_WA$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_genes_WA$calibration_slopes$ICI),2), " (",
                                      round((mean(results_genes_WA$calibration_slopes$ICI) - 1.96*(sd(results_genes_WA$calibration_slopes$ICI) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_genes_WA$calibration_slopes$ICI) + 1.96*(sd(results_genes_WA$calibration_slopes$ICI) / sqrt(nrow(results_genes_WA$calibration_slopes)))),2),
                                      ")"),
                               paste0(round(mean(results_all_WA$calibration_slopes$ICI),2), " (",
                                      round((mean(results_all_WA$calibration_slopes$ICI) - 1.96*(sd(results_all_WA$calibration_slopes$ICI) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                      "-",
                                      round((mean(results_all_WA$calibration_slopes$ICI) + 1.96*(sd(results_all_WA$calibration_slopes$ICI) / sqrt(nrow(results_all_WA$calibration_slopes)))),2),
                                      ")"))
)
write_csv(rf_results_WA, "~/Dropbox/Work/PPS/EU-GEI/rf_results_WA_081224.csv")

##### Prognosis Random Survival Forest ####
# Parameters for nested cross-validation
outerFolds <- 5
outerRepeats <- 10
innerFolds <- 5
innerRepeats <- 10

###### Run model fitting and internal validation ######
# Clinical
clin_rf_surv_WA <- PPS_scored_WA %>% subset(select=c(CAARMS, gafex01, gafex02, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3,3,3)
)

results_clin_surv_WA <- RSF_repeated_nested_cv(
  combined_df = clin_rf_surv_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PPS

PPS_rf_surv_WA <- PPS_scored_WA %>% subset(select=c(Gender:Tobacco, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_rf_surv_WA)-2)
)

results_PPS_surv_WA <- RSF_repeated_nested_cv(
  combined_df = PPS_rf_surv_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PRS
PRS_rf_surv_WA <- PPS_scored_WA %>% subset(select=c(PRS_resid, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PRS_rf_surv_WA)-2)
)

results_PRS_surv_WA <- RSF_repeated_nested_cv(
  combined_df = PRS_rf_surv_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PPS

PPS_clin_rf_surv_WA <- PPS_scored_WA %>% subset(select=c(Gender:Tobacco, CAARMS, gafex01, gafex02, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_clin_rf_surv_WA)-2)
)

results_PPS_clin_surv_WA <- RSF_repeated_nested_cv(
  combined_df = PPS_clin_rf_surv_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PRS

PRS_clin_rf_surv_WA <- PPS_scored_WA %>% subset(select=c(PRS_resid, CAARMS, gafex01, gafex02, Transition, day_exit))
PRS_clin_rf_surv_WA <- as.data.frame(model.matrix(~.-1,PRS_clin_rf_surv_WA))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_rf_surv_WA)-2)
)

results_PRS_clin_surv_WA <- RSF_repeated_nested_cv(
  combined_df = PRS_clin_rf_surv_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# PPS + PRS

PPS_PRS_rf_surv_WA <- PPS_scored_WA %>% subset(select=c(Gender:Tobacco, PRS_resid, Transition, day_exit))
PPS_PRS_rf_surv_WA <- as.data.frame(model.matrix(~.-1,PPS_PRS_rf_surv_WA))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_PRS_rf_surv_WA)-2)
)

results_PPS_PRS_surv_WA <- RSF_repeated_nested_cv(
  combined_df = PPS_PRS_rf_surv_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)

# Clinical + PPS + PRS 

PPS_PRS_clin_rf_surv_WA <- PPS_scored_WA %>% subset(select=c(Gender:Tobacco, PRS_resid, CAARMS, gafex01, gafex02, Transition, day_exit))

tune_grid <- expand.grid(
  mtry=c(3:ncol(PPS_PRS_clin_rf_surv_WA)-2)
)

results_PPS_PRS_clin_surv_WA <- RSF_repeated_nested_cv(
  combined_df = PPS_PRS_clin_rf_surv_WA, 
  outerFolds = outerFolds, 
  outerRepeats = outerRepeats, 
  innerFolds = innerFolds, 
  innerRepeats = innerRepeats, 
  tuneGrid = tune_grid,
  seed = 231
)
###### Generate performance metrics summary table######

rf_surv_WA_results <- data.frame(model=c("Clinical", "PPS", "PRS", "Clinical+PPS", "Clinical+PRS", "PPS+PRS", "All"),
                                 C=c(paste0(round(mean(results_clin_surv_WA$c_stat_nested$C),2), " (",
                                            round(mean(results_clin_surv_WA$c_stat_nested$C) - 1.96*(sd(results_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_clin_surv_WA$c_stat_nested))),2),
                                            "-",
                                            round(mean(results_clin_surv_WA$c_stat_nested$C) + 1.96*(sd(results_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_clin_surv_WA$c_stat_nested))),2),
                                            ")"),
                                     paste0(round(mean(results_PPS_surv_WA$c_stat_nested$C),2), " (",
                                            round(mean(results_PPS_surv_WA$c_stat_nested$C) - 1.96*(sd(results_PPS_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_surv_WA$c_stat_nested))),2),
                                            "-",
                                            round(mean(results_PPS_surv_WA$c_stat_nested$C) + 1.96*(sd(results_PPS_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_surv_WA$c_stat_nested))),2),
                                            ")"),
                                     paste0(round(mean(results_PRS_surv_WA$c_stat_nested$C),2), " (",
                                            round(mean(results_PRS_surv_WA$c_stat_nested$C) - 1.96*(sd(results_PRS_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PRS_surv_WA$c_stat_nested))),2),
                                            "-",
                                            round(mean(results_PRS_surv_WA$c_stat_nested$C) + 1.96*(sd(results_PRS_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PRS_surv_WA$c_stat_nested))),2),
                                            ")"),
                                     paste0(round(mean(results_PPS_clin_surv_WA$c_stat_nested$C),2), " (",
                                            round(mean(results_PPS_clin_surv_WA$c_stat_nested$C) - 1.96*(sd(results_PPS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_clin_surv_WA$c_stat_nested))),2),
                                            "-",
                                            round(mean(results_PPS_clin_surv_WA$c_stat_nested$C) + 1.96*(sd(results_PPS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_clin_surv_WA$c_stat_nested))),2),
                                            ")"),
                                     paste0(round(mean(results_PRS_clin_surv_WA$c_stat_nested$C),2), " (",
                                            round(mean(results_PRS_clin_surv_WA$c_stat_nested$C) - 1.96*(sd(results_PRS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PRS_clin_surv_WA$c_stat_nested))),2),
                                            "-",
                                            round(mean(results_PRS_clin_surv_WA$c_stat_nested$C) + 1.96*(sd(results_PRS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PRS_clin_surv_WA$c_stat_nested))),2),
                                            ")"),
                                     paste0(round(mean(results_PPS_PRS_clin_surv_WA$c_stat_nested$C),2), " (",
                                            round(mean(results_PPS_PRS_clin_surv_WA$c_stat_nested$C) - 1.96*(sd(results_PPS_PRS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$c_stat_nested))),2),
                                            "-",
                                            round(mean(results_PPS_PRS_clin_surv_WA$c_stat_nested$C) + 1.96*(sd(results_PPS_PRS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$c_stat_nested))),2),
                                            ")"),
                                     paste0(round(mean(results_PPS_PRS_clin_surv_WA$c_stat_nested$C),2), " (",
                                            round(mean(results_PPS_PRS_clin_surv_WA$c_stat_nested$C) - 1.96*(sd(results_PPS_PRS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$c_stat_nested))),2),
                                            "-",
                                            round(mean(results_PPS_PRS_clin_surv_WA$c_stat_nested$C) + 1.96*(sd(results_PPS_PRS_clin_surv_WA$c_stat_nested$C) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$c_stat_nested))),2),
                                            ")")),
                                 calibration_in_the_large=c(paste0(round(mean(results_clin_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                                   round(mean(results_clin_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                                                   "-",
                                                                   round(mean(results_clin_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                                                   ")"),
                                                            paste0(round(mean(results_PPS_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                                   round(mean(results_PPS_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                                                   "-",
                                                                   round(mean(results_PPS_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                                                   ")"),
                                                            paste0(round(mean(results_PRS_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                                   round(mean(results_PRS_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                                                   "-",
                                                                   round(mean(results_PRS_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                                                   ")"),
                                                            paste0(round(mean(results_PPS_clin_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                                   round(mean(results_PPS_clin_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                                                   "-",
                                                                   round(mean(results_PPS_clin_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                                                   ")"),
                                                            paste0(round(mean(results_PRS_clin_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                                   round(mean(results_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                                                   "-",
                                                                   round(mean(results_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                                                   ")"),
                                                            paste0(round(mean(results_PPS_PRS_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                                   round(mean(results_PPS_PRS_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                                                   "-",
                                                                   round(mean(results_PPS_PRS_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                                                   ")"),
                                                            paste0(round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$calibration_in_large),2), " (",
                                                                   round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) - 1.96*(sd(results_PPS_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$calibration_slopes))),2),
                                                                   "-",
                                                                   round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) + 1.96*(sd(results_PPS_PRS_clin_surv_WA$calibration_slopes$calibration_in_large) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$calibration_slopes))),2),
                                                                   ")")),
                                 slope=c(paste0(round(mean(results_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                                round(mean(results_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PPS_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                                round(mean(results_PPS_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PPS_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                                round(mean(results_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PPS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                                round(mean(results_PPS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PPS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                                round(mean(results_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PPS_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                                round(mean(results_PPS_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PPS_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE),2), " (",
                                                round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv_WA$calibration_slopes$slope, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$calibration_slopes))),2),
                                                ")")),
                                 Brier=c(paste0(round(mean(results_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                                round(mean(results_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PPS_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                                round(mean(results_PPS_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PPS_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                                round(mean(results_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PPS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                                round(mean(results_PPS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PPS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                                round(mean(results_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PPS_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                                round(mean(results_PPS_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PPS_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                                ")"),
                                         paste0(round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE),2), " (",
                                                round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$calibration_slopes))),2),
                                                "-",
                                                round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv_WA$calibration_slopes$brier, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$calibration_slopes))),2),
                                                ")")),
                                 ICI=c(paste0(round(mean(results_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                              round(mean(results_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                              "-",
                                              round(mean(results_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_clin_surv_WA$calibration_slopes))),2),
                                              ")"),
                                       paste0(round(mean(results_PPS_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                              round(mean(results_PPS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                              "-",
                                              round(mean(results_PPS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_surv_WA$calibration_slopes))),2),
                                              ")"),
                                       paste0(round(mean(results_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                              round(mean(results_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                              "-",
                                              round(mean(results_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_surv_WA$calibration_slopes))),2),
                                              ")"),
                                       paste0(round(mean(results_PPS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                              round(mean(results_PPS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                              "-",
                                              round(mean(results_PPS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_clin_surv_WA$calibration_slopes))),2),
                                              ")"),
                                       paste0(round(mean(results_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                              round(mean(results_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                              "-",
                                              round(mean(results_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PRS_clin_surv_WA$calibration_slopes))),2),
                                              ")"),
                                       paste0(round(mean(results_PPS_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                              round(mean(results_PPS_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                              "-",
                                              round(mean(results_PPS_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_surv_WA$calibration_slopes))),2),
                                              ")"),
                                       paste0(round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE),2), " (",
                                              round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) - 1.96*(sd(results_PPS_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$calibration_slopes))),2),
                                              "-",
                                              round(mean(results_PPS_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) + 1.96*(sd(results_PPS_PRS_clin_surv_WA$calibration_slopes$ICI, na.rm=TRUE) / sqrt(nrow(results_PPS_PRS_clin_surv_WA$calibration_slopes))),2),
                                              ")"))
                                 
)
write_csv(rf_surv_WA_results, "~/Dropbox/Work/PPS/EU-GEI/rf_surv_WA_results_181224.csv")
