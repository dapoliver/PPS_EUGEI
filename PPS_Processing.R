library(tidyverse)
library(readxl)
library(naniar)
library(gtsummary)

# Read in data
PPS_raw <- read_excel("~/Dropbox/Work/PPS/EU-GEI/Databases/EU-GEI_Dominic_Oliver_single_record_3.8.2020.xlsx")
PPS_ethnic_density <- read_excel("~/Dropbox/Work/PPS/EU-GEI/Databases/Ethnic Density.xlsx")
PPS_raw <- merge(PPS_raw,PPS_ethnic_density, by="st_subjid", all.x=TRUE)

# Add Parental mental health data
PPS_figs <- read_excel("~/Dropbox/Work/PPS/EU-GEI/Databases/EU_FIGS_REL-multi_record_28.02.2019.xlsx")
PPS_figs <- PPS_figs %>% group_by(st_subjid) %>% mutate(PSM = case_when(sum(total)>0 ~ 1,
                                                                        TRUE ~ 0))
PPS_figs <- PPS_figs %>% subset(select=c(st_subjid, PSM)) %>% group_by(st_subjid) %>% slice(1)
PPS_raw_figs <- merge(PPS_raw, PPS_figs, by="st_subjid", all.x=TRUE)

PPS_clean <- data.frame(ID=PPS_raw$st_subjid)

PPS_clean <- PPS_clean %>% mutate(site=PPS_raw$site,
                                  Age = PPS_raw$Age,
                                  Male = case_when(PPS_raw$Gender==1 ~ 1,
                                                   is.na(PPS_raw$Gender) ~ NA,
                                                   TRUE ~ 0),
                                  Gender = case_when(PPS_raw$Age<36 & PPS_raw$Age>24 & PPS_raw$Gender==1 ~ 1,
                                                     is.na(PPS_raw$Age) | is.na(PPS_raw$Gender) ~ NA,
                                                     TRUE ~ 0),
                                  Handedness = case_when(PPS_raw$mri_handedness==2  ~ 1,
                                                         is.na(PPS_raw$mri_handedness) ~ NA,
                                                         TRUE ~ 0),
                                  Ethnicity = case_when(PPS_raw$Ethnicity==1  ~ "White",
                                                        PPS_raw$Ethnicity==2  ~ "Black",
                                                        is.na(PPS_raw$Ethnicity) ~ NA,
                                                        TRUE  ~ "Other"),
                                  High_ethnic_density = PPS_raw$`High ethnic density`,
                                  Low_ethnic_density = PPS_raw$`Low ethnic density`,
                                  Ethnicity_ED = case_when(Ethnicity == "Black" & High_ethnic_density==1 ~ "BlackHighED",
                                                           Ethnicity == "Black" & Low_ethnic_density==1 ~ "BlackLowED",
                                                           (Ethnicity == "Black" & High_ethnic_density==0 & Low_ethnic_density==0) | (Ethnicity == "Black" & is.na(High_ethnic_density) & is.na(Low_ethnic_density)) ~ "BlackMediumED",
                                                           Ethnicity == "Other" & High_ethnic_density==1 ~ "OtherHighED",
                                                           Ethnicity == "Other" & Low_ethnic_density==1 ~ "OtherLowED",
                                                           (Ethnicity == "Other" & High_ethnic_density==0 & Low_ethnic_density==0) | (Ethnicity == "Other" & is.na(High_ethnic_density) & is.na(Low_ethnic_density)) ~ "OtherMediumED",
                                                           Ethnicity == "White" ~ "White"),
                                  FirstGenImmigrantNA = case_when(grepl("Alg", PPS_raw$mrc1_socde06)| grepl("Tun", PPS_raw$mrc1_socde06)| grepl("Egy", PPS_raw$mrc1_socde06)| grepl("Libya", PPS_raw$mrc1_socde06)| grepl("Western Sahara", PPS_raw$mrc1_socde06)| grepl("Sudan", PPS_raw$mrc1_socde06)| grepl("Morocco", PPS_raw$mrc1_socde06) | grepl("Maroc", PPS_raw$mrc1_socde06) ~ 1,
                                                                   TRUE ~ 0),
                                  FirstGenImmigrant = case_when(PPS_raw$site==2 & PPS_raw$mrc1_socde05==11 ~ 0,
                                                                PPS_raw$site==6 & PPS_raw$mrc1_socde05==9 ~ 0,
                                                                PPS_raw$site==20 & PPS_raw$mrc1_socde05==1 ~ 0,
                                                                PPS_raw$site==21 & PPS_raw$mrc1_socde05==8 ~ 0,
                                                                PPS_raw$site==22 & PPS_raw$mrc1_socde05==4 ~ 0,
                                                                PPS_raw$site==24 & PPS_raw$mrc1_socde05==13 ~ 0,
                                                                PPS_raw$site==26 & PPS_raw$mrc1_socde05==14 & grepl("Denmark", PPS_raw$mrc1_socde06) ~ 0, # Denmark
                                                                PPS_raw$site==27 & PPS_raw$mrc1_socde05==2 ~ 0,
                                                                PPS_raw$site==31 & PPS_raw$mrc1_socde05==3 ~ 0,
                                                                PPS_raw$site==32 & PPS_raw$mrc1_socde05==7 ~ 0,
                                                                PPS_raw$site==33 & PPS_raw$mrc1_socde05==12 ~ 0,
                                                                FirstGenImmigrantNA==1 ~ 0,
                                                                TRUE ~ 1),
                                  SecondGenImmigrantNA = case_when(grepl("Alg", PPS_raw$mrc1_socde09)| grepl("Tun", PPS_raw$mrc1_socde09)| grepl("Egy", PPS_raw$mrc1_socde09)| grepl("Libya", PPS_raw$mrc1_socde09)| grepl("Western Sahara", PPS_raw$mrc1_socde09)| grepl("Sudan", PPS_raw$mrc1_socde09)| grepl("Morocco", PPS_raw$mrc1_socde09) | grepl("Maroc", PPS_raw$mrc1_socde09)  |
                                                                   grepl("Alg", PPS_raw$mrc1_socde11)| grepl("Tun", PPS_raw$mrc1_socde11)| grepl("Egy", PPS_raw$mrc1_socde11)| grepl("Libya", PPS_raw$mrc1_socde11)| grepl("Western Sahara", PPS_raw$mrc1_socde11)| grepl("Sudan", PPS_raw$mrc1_socde11)| grepl("Morocco", PPS_raw$mrc1_socde11) | grepl("Maroc", PPS_raw$mrc1_socde11) ~ 1,
                                                                   TRUE ~ 0),
                                  SecondGenImmigrant = case_when(PPS_raw$site==2 &  PPS_raw$mrc1_socde08==11 & PPS_raw$mrc1_socde10==11 ~ 0,
                                                                 PPS_raw$site==6 &  PPS_raw$mrc1_socde08==9  & PPS_raw$mrc1_socde10==9 ~ 0,
                                                                 PPS_raw$site==20 & PPS_raw$mrc1_socde08==1  & PPS_raw$mrc1_socde10==1 ~ 0,
                                                                 PPS_raw$site==21 & PPS_raw$mrc1_socde08==8  & PPS_raw$mrc1_socde10==8 ~ 0,
                                                                 PPS_raw$site==22 & PPS_raw$mrc1_socde08==4  & PPS_raw$mrc1_socde10==4 ~ 0,
                                                                 PPS_raw$site==24 & PPS_raw$mrc1_socde08==13 & PPS_raw$mrc1_socde10==13 ~ 0,
                                                                 PPS_raw$site==26 & PPS_raw$mrc1_socde08==14 & PPS_raw$mrc1_socde10==14 & grepl("Denmark", PPS_raw$mrc1_socde09) & grepl("Denmark", PPS_raw$mrc1_socde11) ~ 0, # Denmark
                                                                 PPS_raw$site==27 & PPS_raw$mrc1_socde08==2  & PPS_raw$mrc1_socde10==2 ~ 0,
                                                                 PPS_raw$site==31 & PPS_raw$mrc1_socde08==3  & PPS_raw$mrc1_socde10==3 ~ 0,
                                                                 PPS_raw$site==32 & PPS_raw$mrc1_socde08==7  & PPS_raw$mrc1_socde10==7 ~ 0,
                                                                 PPS_raw$site==33 & PPS_raw$mrc1_socde08==12 & PPS_raw$mrc1_socde10==12 ~ 0,
                                                                 FirstGenImmigrant == 1 | SecondGenImmigrantNA==1 ~ 0,
                                                                 TRUE ~ 1),
                                  PaternalSES = case_when(PPS_raw$mrc1_socde12>6 ~ 1,
                                                          is.na(PPS_raw$mrc1_socde12) ~ NA,
                                                          TRUE ~ 0),
                                  Parental_SMI = case_when(PPS_raw_figs$PSM>0 ~ 1,
                                                          is.na(PPS_raw_figs$PSM) ~ NA,
                                                          TRUE ~ 0),
                                  AdultLifeEvents = case_when(PPS_raw$LTE_total_score>0 ~ 1,
                                                              is.na(PPS_raw$LTE_total_score) ~ NA,
                                                              TRUE ~ 0),
                                  Cannabis = case_when(PPS_raw$ceq15_9<3 ~ 1,
                                                       is.na(PPS_raw$ceq15_1) ~ NA,
                                                       TRUE ~ 0),
                                  ChildhoodTrauma = case_when(PPS_raw$imp_total_ctq_score>55 ~ 1,
                                                              is.na(PPS_raw$imp_total_ctq_score) ~ NA,
                                                              TRUE ~ 0),
                                  Anhedonia = case_when(PPS_raw$caarms_anh1>2 ~ 1,
                                                        is.na(PPS_raw$caarms_anh1) & is.na(PPS_raw$caarms_anh4) ~ NA,
                                                        TRUE ~ 0))


PPS_scored <- read.csv("~/Dropbox/Work/PPS/EU-GEI/Databases/PPS_scored_EUGEI_150221.csv")
PPS_scored <- PPS_scored %>% mutate(Urbanicity = case_when(Urbanicity==1 ~ 1,
                                                           is.na(Urbanicity) ~ NA,
                                                           TRUE ~ 0),
                                    Pollution = case_when(Pollution==2 ~ 1,
                                                          is.na(Pollution) ~ NA,
                                                          TRUE ~ 0),
                                    Tobacco = case_when(Tobacco==3 ~ 1,
                                                          is.na(Tobacco) ~ NA,
                                                          TRUE ~ 0))

PPS_scored <- PPS_scored %>% subset(select=c(ID, Urbanicity, Pollution, Tobacco, day_exit, ARMS_2yr))

PPS_all <- merge(PPS_clean, PPS_scored, by="ID", all.x=TRUE)

PPS_all <- PPS_all %>% mutate(Transition = PPS_raw$Transition_status)

# Integrate PRS
PRS <- read.table("~/Dropbox/Work/PPS/EU-GEI/Databases/PRS/eugei_sczprs22_pcs.txt", header=TRUE)
PRS2 <- PRS[c(1,3:15)] # subset to only include ID, PCs1:10, PRS_SCZ in all subjects
PRS2 <- PRS2 %>% dplyr::rename(ID=st_subjid) # Rename ID to match PPS_scored
resid <- lm(prs_scz ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=PRS2) # Regress out PCs
PRS2$PRS_resid <- resid$residuals # Save out residuals
PRS2 <- PRS2[c(1,2:4,15)] # subset to only include ID and residuals 
PPS_all <- join(PPS_all, PRS2, by="ID") # Add PRS columns to PPS_scored

# Add in clinical variables
clinical <- read_excel("~/Dropbox/Work/PPS/EU-GEI/Databases/eugei_data_all.xlsx")
clinical <- clinical %>% subset(select=c(st_subjid, caarms_unu4, caarms_unu5, caarms_non4, caarms_non5, caarms_per4, caarms_per5, caarms_dis4,caarms_dis5, gafex01,gafex02))
clinical[,c(2:11)] <- lapply(clinical[,c(2:11)], FUN=as.numeric)
clinical <- as.data.frame(clinical)

# Create overall CAARMS-P score with severity*frequency for each subscale
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
clinical <- clinical %>% dplyr::rename(ID=st_subjid) # Rename ID to match PPS_all
PPS_all <- join(PPS_all, clinical, by="ID") # Add clinical variables to PPS_all

PPS_inc <- read_excel("~/Dropbox/Work/PPS/EU-GEI/Databases/EU-GEI_Dominic_Oliver_single_record_3.8.2020.xlsx")
PPS_inc <- PPS_inc %>% subset(select=c(st_subjid, Subject_status))
PPS_all <- merge(PPS_all, PPS_inc, by.x="ID", by.y="st_subjid")
PPS_all <- PPS_all %>% mutate(chr = case_when(Subject_status==1 ~ 1,
                                              Subject_status==3 ~ 0))

table_HC <- tbl_summary(PPS_all,include=c(Age, Male, Handedness, Urbanicity, Pollution, Ethnicity, 
                                                 High_ethnic_density,  Low_ethnic_density, FirstGenImmigrantNA, 
                                                 FirstGenImmigrant, SecondGenImmigrantNA, SecondGenImmigrant,  
                                                 PaternalSES, Parental_SMI, AdultLifeEvents,
                                                 Tobacco, Cannabis,  
                                                 ChildhoodTrauma, Anhedonia, prs_scz, CAARMS, gafex01, gafex02),
                        by=chr,
                        statistic = all_continuous() ~ c("{mean}", "{sd}"),
                        type = list(Age ~ 'continuous2',
                                    prs_scz ~ 'continuous2',
                                    CAARMS ~ 'continuous2',
                                    gafex01 ~ 'continuous2',
                                    gafex02 ~ 'continuous2'),
                        digits = list(all_continuous() ~ c(1, 1),
                                      all_categorical() ~ c(0, 1)))  
table_HC %>% as_gt() %>% gt::gtsave("~/Dropbox/Work/PPS/EU-GEI/Table_HC.docx")

PPS_SA_PRS <- PPS_all %>% filter(Ethnicity=="White" & singleton==TRUE & gen_hom==TRUE)
PPS_SA_PRS <- PPS_SA_PRS %>% subset(select = c(-Age, -Male, -Ethnicity, -Low_ethnic_density, -High_ethnic_density, -Subject_status))
write_csv(PPS_all,"~/Dropbox/Work/PPS/EU-GEI/Databases/PPS_SA_PRS.csv")

PPS_all <- PPS_all %>% subset(select = c(-Age, -Male, -Ethnicity, -Low_ethnic_density, -High_ethnic_density, -Subject_status))
write_csv(PPS_all,"~/Dropbox/Work/PPS/EU-GEI/Databases/PPS_processed.csv")

# Visualise data missingness
PPS_missing <- PPS_all %>% subset(select=c(Gender, Handedness, Urbanicity, Pollution, Ethnicity_ED, 
                                          FirstGenImmigrantNA, 
                                           FirstGenImmigrant, SecondGenImmigrantNA, SecondGenImmigrant,  
                                           PaternalSES, Parental_SMI, AdultLifeEvents,
                                           Tobacco, Cannabis,  
                                           ChildhoodTrauma, Anhedonia, prs_scz, CAARMS, gafex01, gafex02))

png("~/Dropbox/Work/PPS/EU-GEI/Plots/gg_miss_var.png", width = 600, height = 600)
gg_miss_var(PPS_missing)
dev.off()

png("~/Dropbox/Work/PPS/EU-GEI/Plots/gg_miss_upset.png", width = 600, height = 600)
gg_miss_upset(PPS_missing) # Looks like there is no clear pattern
dev.off()

# Check harmonisation again