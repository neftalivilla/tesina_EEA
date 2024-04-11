##Prediction of CVD risk by VAT markers
#Analisis: Neftali Eduardo Antonio Villa
#Junio 2023

setwd("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INGER/Mexico City Cohort Study/Proyectos/VAT Validation")

#####Library#####

library(dplyr)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(readr)
library(survival)
library(jtools)
library(gtsummary)
library(data.table)
library(epiR)
library(riskRegression)
library(prodlim)
library(survival)
library(cmprsk)
library(lava)
library(ggplot2)
library(MASS)
library(treemapify)
library("maxstat")
library("survival")
library(survminer)
library("survival")
library(rms)


#####Import Data####

base <- read_csv("~/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INGER/Mexico City Cohort Study/Data/2021-004 MCPS BASELINE.csv")
base.mort <- read_csv("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INGER/Mexico City Cohort Study/Data/2022-012 MCPS MORTALITY.csv")
base.mrn <- read_csv("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INGER/Mexico City Cohort Study/Data/2022-012 MCPS BASELINE_NMR.csv")
base.rsv <- read_csv("/Users/nefoantonio/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INGER/Mexico City Cohort Study/Data/2022-012 MCPS RESURVEY.csv")
base.edu <- read_csv("~/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INGER/Mexico City Cohort Study/Data/2022-012_Neftali_Omar/2022-012 MCPS BASELINE_EDUGP.csv")
base.ids <- read_csv("~/Library/CloudStorage/OneDrive-UNIVERSIDADNACIONALAUTÓNOMADEMÉXICO/PROYECTOS/INGER/Mexico City Cohort Study/Data/2022-012_Neftali_Omar/2022-012 MCPS BASELINE_IDS.csv")

base <- base %>% 
  left_join(base.mort, by="PATID") %>% 
  left_join(base.mrn, by="PATID") %>% 
  left_join(base.rsv, by="PATID") %>% 
  left_join(base.edu, by="PATID") %>%
  left_join(base.ids, by="PATID")

# set.seed(123)
# base.equipo.3<-base%>%
#   dplyr::select(PATID,MALE,
#                 AGE,
#                 HEIGHT,
#                 WEIGHT,
#                 COYOACAN,
#                 EDU_LEVEL,
#                 INCOME,
#                 EVER_SMOK,
#                 CURR_SMOK,
#                 EVER_PHYS,
#                 DUR_PHYS,
#                 BASE_DIABETES,
#                 BASE_HYPERTENSION,
#                 R_HBA1C,
#                 PAS_PROM,
#                 PAD_PROM,
#                 D003)%>% distinct %>% sample_n(4000)
# 
# write.csv(base.equipo.3,"BASE_EQUIPO_3_MCPS.csv")

#####Estimation of Formulas######

#Mean SBP
base$PAS_PROM<-rowMeans(as.matrix(base%>%dplyr::select(SBP1,SBP2,SBP3)))

#Mean DBP
base$PAD_PROM<-rowMeans(as.matrix(base%>%dplyr::select(DBP1,DBP2,DBP3)))

#Diabetes
base$DIABETES_FINAL<-NULL
base$DIABETES_FINAL[base$BASE_DIABETES==1 | 
                      base$DRUG_D1 == 1 |
                      base$DRUG_D2 == 1 |
                      base$DRUG_D3 == 1 |
                      base$DRUG_D4 == 1 |
                      base$BASE_HBA1C >= 6.5]<-1
base$DIABETES_FINAL<-na.tools::na.replace(base$DIABETES_FINAL,0)

#Cardiovascular Disease
base$CVD_BASAL<-NULL
base$CVD_BASAL[base$BASE_HEARTATTACK==1 |
                 base$BASE_ANGINA==1 |
                 base$BASE_STROKE==1]<-1
base$CVD_BASAL<-na.tools::na.replace(base$CVD_BASAL,0)

#Hypertension Definition
base$HAS_FINAL<-NULL
base$HAS_FINAL[base$BASE_HYPERTENSION==1 | 
                 base$DRUG_A1 == 1 |
                 base$DRUG_A2 == 1 |
                 base$DRUG_A3 == 1 |
                 base$DRUG_A4 == 1 |
                 base$DRUG_A5 == 1 |
                 base$DRUG_A6 == 1 |
                 base$DRUG_A7 == 1 |
                 base$DRUG_A8 == 1 |
                 base$DRUG_A9 == 1 |
                 base$DRUG_A10 == 1 |
                 base$DRUG_A11 == 1 |
                 base$PAS_PROM >= 140 | base$PAD_PROM >= 90]<-1
base$HAS_FINAL<-na.tools::na.replace(base$HAS_FINAL,0)

#CVD DEATHS
base$VASC_NON_CARD<-NULL
base$VASC_NON_CARD[base$D011==1]<-1
base$VASC_NON_CARD[base$D012==1]<-1
base$VASC_NON_CARD[is.na(base$VASC_NON_CARD)]<-0

## Numeric Variables
base$IMC<-base$WEIGHT/((base$HEIGHT/100)^2)
base$ICE<-base$WAISTC/base$HEIGHT

#Transform Laboratories to mg/dl
base$Glc_mgdl<-base$Glc*18;base$Glc_mgdl[base$Glc_mgdl<=0]<-NA
base$Serum_TG_mgdl<-base$Serum_TG*88.57;base$Serum_TG_mgdl[base$Serum_TG_mgdl<=0]<-NA
base$HDL_C_mgdl<-base$HDL_C*38.67;base$HDL_C_mgdl[base$HDL_C_mgdl<=0]<-NA
base$LDL_C_mgdl<-base$LDL_C*38.67;base$LDL_C_mgdl[base$LDL_C_mgdl<=0]<-NA

## Log Transform Income

base$INCOME_LN<-log(base$INCOME+1)

#Education Recode
base$EDUGP_2<-NULL
base$EDUGP_2[base$EDUGP==1]<-1
base$EDUGP_2[base$EDUGP==2]<-1
base$EDUGP_2[base$EDUGP==3]<-1

base$EDUGP_2[base$EDUGP==4]<-2
base$EDUGP_2[base$EDUGP==5]<-2
base$EDUGP_2[base$EDUGP==6]<-2

base$EDUGP_2[base$EDUGP==7]<-3
base$EDUGP_2[base$EDUGP==8]<-3
base$EDUGP_2[base$EDUGP==9]<-3

base$EDUGP_2[base$EDUGP==10]<-4
base$EDUGP_2[base$EDUGP==11]<-4
base$EDUGP_2[base$EDUGP==12]<-4
base$EDUGP_2[base$EDUGP==13]<-4

base$EDAD_CAT<-NULL
base$EDAD_CAT[base$AGE<45]<-1
base$EDAD_CAT[base$AGE>=45 & base$AGE<65]<-2
base$EDAD_CAT[base$AGE>=65]<-3

## Transform selected numeric to factors
base$EDAD_CAT<-factor(base$EDAD_CAT,labels = c("<45","45-65",">65"))
base$COYOACAN<- factor(base$COYOACAN,levels = c(0,1),labels = c("Iztapalapa","Coyoacan"))
base$EDU_LEVEL<- factor(base$EDU_LEVEL,levels = c(1:4),labels = c("University/College","High-School","Elementary Only","Other"))
base$EDUGP<- factor(base$EDUGP,levels = c(1:13),labels = c("Illiterate",
                                                           "Knows how to read",
                                                           "Knows how to read and write",
                                                           "Incomplete elementary",
                                                           "Complete elementary",
                                                           "Tecnical Studies with complete elementary",
                                                           "Incomplete high School",
                                                           "Complete high School",
                                                           "Tecnical Studies with complete high school",
                                                           "Collegue",
                                                           "Tecnical Studies with complete collegue school",
                                                           "Incomplete univesity",
                                                           "Complete univesity"))


base$EDUGP_2<- factor(base$EDUGP_2,labels = c("Illiterate or Non-proper education",
                                              "Elementary",
                                              "High School",
                                              "Collegue"))

## Ocupation Categories
base$OCCUPATION_REC<-NULL

base$OCCUPATION_REC[base$OCCUPATION==10]<-1
base$OCCUPATION_REC[base$OCCUPATION==11]<-1
base$OCCUPATION_REC[base$OCCUPATION==1]<-1
base$OCCUPATION_REC[base$OCCUPATION==13]<-1
base$OCCUPATION_REC[base$OCCUPATION==14]<-1

base$OCCUPATION_REC[base$OCCUPATION==2]<-2
base$OCCUPATION_REC[base$OCCUPATION==3]<-2
base$OCCUPATION_REC[base$OCCUPATION==4]<-2
base$OCCUPATION_REC[base$OCCUPATION==6]<-2
base$OCCUPATION_REC[base$OCCUPATION==16]<-2
base$OCCUPATION_REC[base$OCCUPATION==19]<-2
base$OCCUPATION_REC[base$OCCUPATION==17]<-2
base$OCCUPATION_REC[base$OCCUPATION==18]<-2
base$OCCUPATION_REC[base$OCCUPATION==12]<-2

base$OCCUPATION_REC[base$OCCUPATION==7]<-3
base$OCCUPATION_REC[base$OCCUPATION==8]<-3
base$OCCUPATION_REC[base$OCCUPATION==5]<-3
base$OCCUPATION_REC[base$OCCUPATION==20]<-3

base$OCCUPATION_REC[base$OCCUPATION==15]<-4
base$OCCUPATION_REC[base$OCCUPATION==21]<-4

base$OCCUPATION_REC<-factor(base$OCCUPATION_REC,levels = c(1:4),labels = c("Private Employers and Professionals",
                                                                           "Blue-Collar Workers",
                                                                           "Public Sector Workers",
                                                                           "Retired or Unemployed"))

#Health Provider
base$HEALTH_PROVIDER_2<-NULL
base$HEALTH_PROVIDER_2[base$HEALTH_PROVIDER==1]<-1
base$HEALTH_PROVIDER_2[base$HEALTH_PROVIDER==2]<-1
base$HEALTH_PROVIDER_2[base$HEALTH_PROVIDER==3]<-1
base$HEALTH_PROVIDER_2[base$HEALTH_PROVIDER==4]<-1
base$HEALTH_PROVIDER_2[base$HEALTH_PROVIDER==5]<-1
base$HEALTH_PROVIDER_2[base$HEALTH_PROVIDER==6]<-1
base$HEALTH_PROVIDER_2[base$HEALTH_PROVIDER==7]<-1
base$HEALTH_PROVIDER_2[base$HEALTH_PROVIDER==8]<-2
base$HEALTH_PROVIDER_2[base$HEALTH_PROVIDER==9]<-2
base$HEALTH_PROVIDER_2[base$HEALTH_PROVIDER==10]<-3
base$HEALTH_PROVIDER_2[is.na(base$HEALTH_PROVIDER)]<-4
base$HEALTH_PROVIDER_2<-factor(base$HEALTH_PROVIDER_2,labels = c("Public HC","Private HC","Non-Specified","Missing"))


#Alcohol Intake
base$ALCGP<-factor(base$ALCGP,labels = c("Never","Former",">3 times a month", ">2 times a week", ">3 times a week"))

#Physical Activity
base$PHYSGP<-factor(base$PHYSGP,labels = c("None", ">2 times a week",">3 times a week"))

#Visceral Adiposiy Metrics

#METS_IR
base$METS_IR<-((log((2*base$Glc_mgdl)+base$Serum_TG_mgdl)*base$IMC))/(log(base$HDL_C_mgdl))

#METS_VF
base$METSVF<-(4.466+0.011*(log(base$METS_IR)^3)+3.239*(log(base$ICE)^3)+0.319*(base$MALE)
              +0.594*(log(base$AGE)))

#Grasa Visceral
base$VAT_METS<-exp(base$METSVF)

#Indices 
base<-base %>% 
  mutate(DAAT_index = if_else(MALE >= 1,
                              (-382.9+(1.09*base$WEIGHT)+(6.04*base$WAISTC)+(-2.29*base$IMC)),
                              (-278+(-0.86*base$WEIGHT)+(5.19*base$WAISTC))),
         Depres_index = if_else(MALE >= 1,
                                (-225.39 +(2.125*base$AGE)+(2.843*base$WAISTC)),
                                NA),
         VAI_index = if_else(MALE >= 1,
                             ((base$WAISTC/(39.68+(1.88*base$IMC)))*(base$Serum_TG/1.03)*(1.31/base$HDL_C)),
                             ((base$WAISTC/(36.58+(1.89*base$IMC)))*(base$Serum_TG/0.81)*(1.52/base$HDL_C))),
         VAI_GEA_index = if_else(MALE >= 1,
                                 ((base$WAISTC/(22.79+(2.68*base$IMC)))*(base$Serum_TG/1.37)*(1.19/base$HDL_C)),
                                 ((base$WAISTC/(24.02+(2.37*base$IMC)))*(base$Serum_TG/1.32)*(1.43/base$HDL_C))),
         LAAP_index = if_else(MALE >= 1,
                              ((base$WAISTC-65)*base$Serum_TG),
                              ((base$WAISTC-58)*base$Serum_TG)),
         EVA_index = if_else(MALE >= 1,
                             ((1.28*base$AGE)+(4.12*base$WAISTC)-(0.53*base$HDL_C_mgdl)+(0.14*base$Glc_mgdl)-319.50),
                             ((1.26*base$AGE)+(1.89*base$IMC)+(2.16*base$WAISTC)-(0.43*base$HDL_C_mgdl)+(0.11*base$LDL_C_mgdl)+(0.18*base$Glc_mgdl)-207.22)))

####Recoding labels####

setattr(base$AGE, "label", "Edad, (Años)")
setattr(base$MALE, "label", "Sexo, (%)")
setattr(base$COYOACAN, "label", "Municipio, (%)")
setattr(base$EDU_LEVEL, "label", "Educación, (%)")
setattr(base$INCOME, "label", "Ingreso, (pesos/month)")
setattr(base$HEALTH_PROVIDER_2, "label", "Provedor de Servicios de Salud, (%)")
setattr(base$OCCUPATION_REC, "label", "Ocupación, (%)")
setattr(base$EVER_SMOK, "label", "Tabaquismo, (%)")
setattr(base$ALCGP, "label", "Habitos de Consumo de Alcohol, (%)")
setattr(base$PHYSGP, "label", "Actividad Física, (%)")
setattr(base$IMC, "label", "Indice de Masa Corporal, (kg/m2)")
setattr(base$ICE, "label", "Waist-to-Height Ratio, (%)")
setattr(base$HIPC, "label", "Hip Circunference, (cm)")
setattr(base$PAS_PROM, "label", "PAS, (mmHg)")
setattr(base$PAD_PROM, "label", "TAD, (mmHg)")
setattr(base$BASE_HBA1C, "label", "HbA1c, (%)")
setattr(base$Glc_mgdl, "label", "Glucose, (mg/dl)")
setattr(base$Serum_TG_mgdl, "label", "Trigliceridos, (mg/dl)")
setattr(base$HDL_C_mgdl, "label", "HDL-C, (mg/dl)")
setattr(base$LDL_C_mgdl, "label", "LDL-C, (mg/dl)")

#####Dataset Depury####

base.2<-base%>%
  dplyr::filter(DIABETES_FINAL!=1)%>%
  dplyr::filter(CVD_BASAL!=1)%>%
  dplyr::filter(CVD_BASAL!=1)%>%
  dplyr::filter(AGE<80)%>%
  dplyr::filter(IMC<40)%>%
  dplyr::filter(IMC>=18.5)%>%
  dplyr::mutate(Depres.meanval = mean(Depres_index,na.rm = T), Depres.stdev = sd(Depres_index,na.rm = T),
                DAAT.meanval = mean(DAAT_index,na.rm = T), DAAT.stdev = sd(DAAT_index,na.rm = T),
                LAAP.meanval = mean(LAAP_index,na.rm = T), LAAP.stdev = sd(LAAP_index,na.rm = T),
                EVA.meanval = mean(EVA_index,na.rm = T), EVA.stdev = sd(EVA_index,na.rm = T),
                METS.meanval = mean(METSVF,na.rm = T), METS.stdev = sd(METSVF,na.rm = T),
                DAI.meanval = mean(VAI_GEA_index,na.rm = T), DAI.stdev = sd(VAI_GEA_index,na.rm = T))

nrow(base.2)

base.2<-base.2%>%
  dplyr::filter(!is.na(VAI_index))%>% 
  dplyr::filter(!is.na(EVA_index))%>% 
  dplyr::filter(!is.na(PERSON_YEARS))%>% 
  dplyr::filter(STATUS!="U")%>% 
  dplyr::filter(abs((DAAT_index-DAAT.meanval)/DAAT.stdev)<3)%>%
  dplyr::filter(abs((LAAP_index-LAAP.meanval)/LAAP.stdev)<3)%>%
  dplyr::filter(abs((EVA_index-EVA.meanval)/EVA.stdev)<3)%>%
  dplyr::filter(abs((METSVF-METS.meanval)/METS.stdev)<3)%>%
  dplyr::filter(abs((VAI_GEA_index-DAI.meanval)/DAI.stdev)<3)

nrow(base.2)

sum(is.na(base.2$Depres_index))
sum(is.na(base.2$DAAT_index))
sum(is.na(base.2$LAAP_index))
sum(is.na(base.2$VAI_index))
sum(is.na(base.2$EVA_index))
sum(is.na(base.2$METSVF))
sum(is.na(base.2$VAI_GEA_index))
sum(is.na(base.2$PERSON_YEARS))




#####Descriptive Characteristics of the Included Population (Table 1)#####

t.0<-base.2 %>% 
  dplyr::select(EDAD_CAT,AGE,COYOACAN,EDUGP_2,
                INCOME,OCCUPATION_REC,HEALTH_PROVIDER_2,
                EVER_SMOK,ALCGP,PHYSGP,
                IMC,ICE,HIPC,PAS_PROM,PAD_PROM,BASE_HBA1C,Glc_mgdl,Serum_TG_mgdl,HDL_C_mgdl,LDL_C_mgdl)%>%
  tbl_summary(missing = "no")%>%
  bold_labels()%>%
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))


t1.1<-base.2 %>% 
  dplyr::filter(MALE==1)%>%
  dplyr::select(EDAD_CAT,AGE,COYOACAN,EDUGP_2,
                INCOME,OCCUPATION_REC,HEALTH_PROVIDER_2,
                EVER_SMOK,ALCGP,PHYSGP,
                IMC,ICE,HIPC,PAS_PROM,PAD_PROM,BASE_HBA1C,Glc_mgdl,Serum_TG_mgdl,HDL_C_mgdl,LDL_C_mgdl)%>%
  tbl_summary(by = EDAD_CAT,missing = "no")%>%
  bold_labels()%>%
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))

t1.2<-base.2 %>% 
  dplyr::filter(MALE!=1)%>%
  dplyr::select(EDAD_CAT,AGE,COYOACAN,EDUGP_2,
                INCOME,OCCUPATION_REC,HEALTH_PROVIDER_2,
                EVER_SMOK,ALCGP,PHYSGP,
                IMC,ICE,HIPC,PAS_PROM,PAD_PROM,BASE_HBA1C,Glc_mgdl,Serum_TG_mgdl,HDL_C_mgdl,LDL_C_mgdl)%>%
  tbl_summary(by = EDAD_CAT,missing = "no")%>%
  bold_labels()%>%
  modify_spanning_header(all_stat_cols() ~ "**Overall Sample**")%>%
  modify_table_body(
    dplyr::mutate,
    label = ifelse(label == "N missing (% missing)",
                   "Unknown",
                   label))

tbl_merge(
  tbls = list(t.0,t1.1, t1.2),
  tab_spanner = c("Muestra Total (n=29,034)","**Hombres (n=10,593)**", "**Mujeres (n=18,441)**"))%>%
  as_flex_table()%>%
  flextable::save_as_docx(path="Table_1.docx")


#####Cause-Specific Mortality Models (Table 2)#####

#Muerte All Vascular
base.2$MUERTE_VASCULAR
base.2$MUERTE_VASCULAR[base.2$STATUS=="A" & base.2$D015==0]<-0
base.2$MUERTE_VASCULAR[base.2$STATUS=="D" & base.2$D015==1]<-1
base.2$MUERTE_VASCULAR[base.2$STATUS=="D" & base.2$D015==0]<-2

event<- base.2$MUERTE_VASCULAR
event<- factor(event, 0:2, labels=c("censor","Vascular","death_other"))

#METS-VF
fgfit1<-FGR(Hist(PERSON_YEARS,MUERTE_VASCULAR)~scale(METSVF),data=base.2,cause = 1)
summary(fgfit1)

#Deep-Abdominal Adipose Tissue index
fgfit2<-FGR(Hist(PERSON_YEARS,MUERTE_VASCULAR)~scale(DAAT_index),data=base.2,cause = 1)
summary(fgfit2)

#Visceral Adiposity Index
fgfit3<-FGR(Hist(PERSON_YEARS,MUERTE_VASCULAR)~scale(VAI_index),data=base.2,cause = 1)
summary(fgfit3)

#GEA-VAI
fgfit4<-FGR(Hist(PERSON_YEARS,MUERTE_VASCULAR)~scale(VAI_GEA_index),data=base.2,cause = 1)
summary(fgfit4)

#Lipid Accumulation Product
fgfit5<-FGR(Hist(PERSON_YEARS,MUERTE_VASCULAR)~scale(LAAP_index),data=base.2,cause = 1)
summary(fgfit5)

#Eva Index
fgfit6<-FGR(Hist(PERSON_YEARS,MUERTE_VASCULAR)~scale(EVA_index),data=base.2,cause = 1)
summary(fgfit6)

#Depress
fgfit7<-FGR(Hist(PERSON_YEARS,MUERTE_VASCULAR)~scale(Depres_index),data=base.2,cause = 1)
summary(fgfit7)

#Any Cardiac Death
base.2$MUERTE_CARDIO
base.2$MUERTE_CARDIO[base.2$STATUS=="A" & base.2$D003==0]<-0
base.2$MUERTE_CARDIO[base.2$STATUS=="D" & base.2$D003==1]<-1
base.2$MUERTE_CARDIO[base.2$STATUS=="D" & base.2$D003==0]<-2

#METS-VF
fg.fit.1.cardio<-FGR(Hist(PERSON_YEARS,MUERTE_CARDIO)~scale(METSVF),data=base.2,cause = 1)
summary(fg.fit.1.cardio)

#Deep-Abdominal Adipose Tissue index
fg.fit.2.cardio<-FGR(Hist(PERSON_YEARS,MUERTE_CARDIO)~scale(DAAT_index),data=base.2,cause = 1)
summary(fg.fit.2.cardio)

#Visceral Adiposity Index
fg.fit.3.cardio<-FGR(Hist(PERSON_YEARS,MUERTE_CARDIO)~scale(VAI_index),data=base.2,cause = 1)
summary(fg.fit.3.cardio)

#GEA-VAI
fg.fit.4.cardio<-FGR(Hist(PERSON_YEARS,MUERTE_CARDIO)~scale(VAI_GEA_index),data=base.2,cause = 1)
summary(fg.fit.4.cardio)

#Lipid Accumulation Product
fg.fit.5.cardio<-FGR(Hist(PERSON_YEARS,MUERTE_CARDIO)~scale(LAAP_index),data=base.2,cause = 1)
summary(fg.fit.5.cardio)

#Eva Index
fg.fit.6.cardio<-FGR(Hist(PERSON_YEARS,MUERTE_CARDIO)~scale(EVA_index),data=base.2,cause = 1)
summary(fg.fit.6.cardio)

#Depress
fg.fit.7.cardio<-FGR(Hist(PERSON_YEARS,MUERTE_CARDIO)~scale(Depres_index),data=base.2,cause = 1)
summary(fg.fit.7.cardio)

#Any Cerebrovascular Death
base.2$MUERTE_STROKE
base.2$MUERTE_STROKE[base.2$STATUS=="A" & base.2$D008==0]<-0
base.2$MUERTE_STROKE[base.2$STATUS=="D" & base.2$D008==1]<-1
base.2$MUERTE_STROKE[base.2$STATUS=="D" & base.2$D008==0]<-2

#METS-VF
fg.fit.1.stroke<-FGR(Hist(PERSON_YEARS,MUERTE_STROKE)~scale(METSVF),data=base.2,cause = 1)
summary(fg.fit.1.stroke)

#Deep-Abdominal Adipose Tissue index
fg.fit.2.stroke<-FGR(Hist(PERSON_YEARS,MUERTE_STROKE)~scale(DAAT_index),data=base.2,cause = 1)
summary(fg.fit.2.stroke)

#Visceral Adiposity Index
fg.fit.3.stroke<-FGR(Hist(PERSON_YEARS,MUERTE_STROKE)~scale(VAI_index),data=base.2,cause = 1)
summary(fg.fit.3.stroke)

#GEA-VAI
fg.fit.4.stroke<-FGR(Hist(PERSON_YEARS,MUERTE_STROKE)~scale(VAI_GEA_index),data=base.2,cause = 1)
summary(fg.fit.4.stroke)

#Lipid Accumulation Product
fg.fit.5.stroke<-FGR(Hist(PERSON_YEARS,MUERTE_STROKE)~scale(LAAP_index),data=base.2,cause = 1)
summary(fg.fit.5.stroke)

#Eva Index
fg.fit.6.stroke<-FGR(Hist(PERSON_YEARS,MUERTE_STROKE)~scale(EVA_index),data=base.2,cause = 1)
summary(fg.fit.6.stroke)

#Depress
fg.fit.7.stroke<-FGR(Hist(PERSON_YEARS,MUERTE_STROKE)~scale(Depres_index),data=base.2,cause = 1)
summary(fg.fit.7.stroke)

#Other CVD 
base.2$MUERTE_OTRAS
base.2$MUERTE_OTRAS[base.2$D009==1]<-1
base.2$MUERTE_OTRAS[base.2$D010==1]<-1
base.2$MUERTE_OTRAS[base.2$D011==1]<-1
base.2$MUERTE_OTRAS[base.2$D012==1]<-1
base.2$MUERTE_OTRAS<-na.tools::na.replace(base.2$MUERTE_OTRAS,0)

base.2$MUERTE_OTRAS_2[base.2$STATUS=="A" & base.2$MUERTE_OTRAS==0]<-0
base.2$MUERTE_OTRAS_2[base.2$STATUS=="D" & base.2$MUERTE_OTRAS==1]<-1
base.2$MUERTE_OTRAS_2[base.2$STATUS=="D" & base.2$MUERTE_OTRAS==0]<-2

#METS-VF
fg.fit.1.others<-FGR(Hist(PERSON_YEARS,MUERTE_OTRAS_2)~scale(METSVF),data=base.2,cause = 1)
summary(fg.fit.1.others)

  
#Deep-Abdominal Adipose Tissue index
fg.fit.2.others<-FGR(Hist(PERSON_YEARS,MUERTE_OTRAS_2)~scale(DAAT_index),data=base.2,cause = 1)
summary(fg.fit.2.others)

#Visceral Adiposity Index
fg.fit.3.others<-FGR(Hist(PERSON_YEARS,MUERTE_OTRAS_2)~scale(VAI_index),data=base.2,cause = 1)
summary(fg.fit.3.others)

#GEA-VAI
fg.fit.4.others<-FGR(Hist(PERSON_YEARS,MUERTE_OTRAS_2)~scale(VAI_GEA_index),data=base.2,cause = 1)
summary(fg.fit.4.others)

#Lipid Accumulation Product
fg.fit.5.others<-FGR(Hist(PERSON_YEARS,MUERTE_OTRAS_2)~scale(LAAP_index),data=base.2,cause = 1)
summary(fg.fit.5.others)

#Eva Index
fg.fit.6.others<-FGR(Hist(PERSON_YEARS,MUERTE_OTRAS_2)~scale(EVA_index),data=base.2,cause = 1)
summary(fg.fit.6.others)

#Depress
fg.fit.7.others<-FGR(Hist(PERSON_YEARS,MUERTE_OTRAS_2)~scale(Depres_index),data=base.2,cause = 1)
summary(fg.fit.7.others)


#####Overall AUROC Curve (Figure 2)#####

#ROC Metrics (Overall CVD Mortality)
ROC.df.1<-Score(list("METS-VF"=fgfit1,
                "DAAT"=fgfit2,
                "VAI"=fgfit3,
                "DAI"=fgfit4,
                "LAP"=fgfit5,
                "EVA"=fgfit6),
           formula=Hist(PERSON_YEARS,MUERTE_VASCULAR)~1,
           data=base.2,
           conf.int=TRUE,
           summary="risks",
           metrics="auc",
           plots="roc")

#ROC Metrics (Ischemic Heart Disease)
ROC.df.2<-Score(list("METS-VF"=fg.fit.1.cardio,
                     "DAAT"=fg.fit.2.cardio,
                     "VAI"=fg.fit.3.cardio,
                     "DAI"=fg.fit.4.cardio,
                     "LAP"=fg.fit.5.cardio,
                     "EVA"=fg.fit.6.cardio),
                formula=Hist(PERSON_YEARS,MUERTE_CARDIO)~1,
                data=base.2,
                conf.int=TRUE,
                summary="risks",
                metrics="auc",
                plots="roc")

#ROC Metrics (Stroke)
ROC.df.3<-Score(list("METS-VF"=fg.fit.1.stroke,
                     "DAAT"=fg.fit.2.stroke,
                     "VAI"=fg.fit.3.stroke,
                     "DAI"=fg.fit.4.stroke,
                     "LAP"=fg.fit.5.stroke,
                     "EVA"=fg.fit.6.stroke),
                formula=Hist(PERSON_YEARS,MUERTE_STROKE)~1,
                data=base.2,
                conf.int=TRUE,
                summary="risks",
                metrics="auc",
                plots="roc")

#ROC Metrics (Others)
ROC.df.4<-Score(list("METS-VF"=fg.fit.1.others,
                     "DAAT"=fg.fit.2.others,
                     "VAI"=fg.fit.3.others,
                     "DAI"=fg.fit.4.others,
                     "LAP"=fg.fit.5.others,
                     "EVA"=fg.fit.6.others),
                formula=Hist(PERSON_YEARS,MUERTE_OTRAS_2)~1,
                data=base.2,
                conf.int=TRUE,
                summary="risks",
                metrics="auc",
                plots="roc")

#ROC CURVE

Figure2A<-ggplotify::as.ggplot(~plotROC(ROC.df.1,col = c(ggsci::pal_jama("default")(6)),
                              xlab = "Tasa de Falsos Positivos (1-Especificidad)",
                              ylab = "Tasa de Verdaderos Positivos (Sensibilidad)"))+
  ggtitle("Todas las Causas Cardiovasculares")



Figure2B<-ggplotify::as.ggplot(~plotROC(ROC.df.2,col = c(ggsci::pal_jama("default")(6)),
                                        xlab = "Tasa de Falsos Positivos (1-Especificidad)",
                                        ylab = "Tasa de Verdaderos Positivos (Sensibilidad)"))+
  ggtitle("Afecciones Cardíacas")

Figure2C<-ggplotify::as.ggplot(~plotROC(ROC.df.3,col = c(ggsci::pal_jama("default")(6)),
                                        xlab = "Tasa de Falsos Positivos (1-Especificidad)",
                                        ylab = "Tasa de Verdaderos Positivos (Sensibilidad)"))+
  ggtitle("Enfermedades Cerebrovasculares")

Figure2D<-ggplotify::as.ggplot(~plotROC(ROC.df.4,col = c(ggsci::pal_jama("default")(6)),
                                        xlab = "Tasa de Falsos Positivos (1-Especificidad)",
                                        ylab = "Tasa de Verdaderos Positivos (Sensibilidad)"))+
  ggtitle("Otras Causas Cardiovasculares")

Figure2<-ggarrange(Figure2A,Figure2B,Figure2C,Figure2D,ncol = 2,nrow = 2)

ggsave(Figure2,
       filename = "Figura_2.png", 
       width = 35, 
       height = 30,
       units=c("cm"),
       dpi = 300,
       bg = "white",
       limitsize = FALSE)

#####AUROC Curve Time Dependent (Figure 3)#####

#ROC Metrics (Overall CVD Mortality)
ROC.time.df.1<-Score(list("METS-VF"=fgfit1,
                   "DAAT"=fgfit2,
                   "VAI"=fgfit3,
                   "DAI"=fgfit4,
                   "LAP"=fgfit5,
                   "EVA"=fgfit6),
              formula=Hist(PERSON_YEARS,MUERTE_VASCULAR)~1,
              data=base.2,
              times = seq(1,20,1),
              conf.int=TRUE,
              summary="risks",
              metrics="auc",
              plots="roc")

#ROC Metrics (Cardiac Causes)
ROC.time.df.2<-Score(list("METS-VF"=fg.fit.1.cardio,
                        "DAAT"=fg.fit.2.cardio,
                        "VAI"=fg.fit.3.cardio,
                        "DAI"=fg.fit.4.cardio,
                        "LAP"=fg.fit.5.cardio,
                        "EVA"=fg.fit.6.cardio),
                   formula=Hist(PERSON_YEARS,MUERTE_CARDIO)~1,
                   data=base.2,
                   times = seq(1,20,1),
                   conf.int=TRUE,
                   summary="risks",
                   metrics="auc",
                   plots="roc")

#ROC Metrics (Cerebrovascular Causes)
ROC.time.df.3<-Score(list("METS-VF"=fg.fit.1.stroke,
                          "DAAT"=fg.fit.2.stroke,
                          "VAI"=fg.fit.3.stroke,
                          "DAI"=fg.fit.4.stroke,
                          "LAP"=fg.fit.5.stroke,
                          "EVA"=fg.fit.6.stroke),
                     formula=Hist(PERSON_YEARS,MUERTE_STROKE)~1,
                     data=base.2,
                     times = seq(1,20,1),
                     conf.int=TRUE,
                     summary="risks",
                     metrics="auc",
                     plots="roc")

#ROC Metrics (Other Causes)
ROC.time.df.4<-Score(list("METS-VF"=fg.fit.1.others,
                          "DAAT"=fg.fit.2.others,
                          "VAI"=fg.fit.3.others,
                          "DAI"=fg.fit.4.others,
                          "LAP"=fg.fit.5.others,
                          "EVA"=fg.fit.6.others),
                     formula=Hist(PERSON_YEARS,MUERTE_OTRAS_2)~1,
                     data=base.2,
                     times = seq(1,20,1),
                     conf.int=TRUE,
                     summary="risks",
                     metrics="auc",
                     plots="roc")


Figure3A<-ggplot(ROC.time.df.1$AUC$score,
       aes(times,AUC,colour=model))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=lower,
                  ymax=upper,
                  colour=model),
              alpha=0.1,linetype=2)+
  ggsci::scale_colour_jama()+
  labs(colour="Indice")+
  theme_classic()+
  ylab("Area Bajo la Curva")+
  xlab("Tiempo de Seguimiento, (Años)")+
  ggtitle("Todas las Causas Cardiovasculares")+
  geom_hline(yintercept = 0.5, linetype=1, 
               color = "red", size=1)

Figure3B<-ggplot(ROC.time.df.2$AUC$score,
                 aes(times,AUC,colour=model))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=lower,
                  ymax=upper,
                  colour=model),
              alpha=0.1,linetype=2)+
  ggsci::scale_colour_jama()+
  labs(colour="Indice")+
  theme_classic()+
  ylab("Area Bajo la Curva")+
  xlab("Tiempo de Seguimiento, (Años)")+
  ggtitle("Afecciones Cardíacas")+
  geom_hline(yintercept = 0.5, linetype=1, 
             color = "red", size=1)

Figure3C<-ggplot(ROC.time.df.3$AUC$score,
                 aes(times,AUC,colour=model))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=lower,
                  ymax=upper,
                  colour=model),
              alpha=0.1,linetype=2)+
  ggsci::scale_colour_jama()+
  labs(colour="Indice")+
  theme_classic()+
  ylab("Area Bajo la Curva")+
  xlab("Tiempo de Seguimiento, (Años)")+
  ggtitle("Enfermedades Cerebrovasculares")+
  geom_hline(yintercept = 0.5, linetype=1, 
             color = "red", size=1)

Figure3D<-ggplot(ROC.time.df.4$AUC$score,
                 aes(times,AUC,colour=model))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=lower,
                  ymax=upper,
                  colour=model),
              alpha=0.1,linetype=2)+
  ggsci::scale_colour_jama()+
  labs(colour="Indice")+
  theme_classic()+
  ylab("Area Bajo la Curva")+
  xlab("Tiempo de Seguimiento, (Años)")+
  ggtitle("Otras Causas Cardiovasculares")+
  geom_hline(yintercept = 0.5, linetype=1, 
             color = "red", size=1)

Figure3<-ggarrange(Figure3A,Figure3B,Figure3C,Figure3D,ncol = 2,nrow = 2,labels = LETTERS[1:4])

ggsave(Figure3,
       filename = "Figura_3.png", 
       width = 30, 
       height = 25,
       units=c("cm"),
       dpi = 300,
       bg = "white",
       limitsize = FALSE)

#####Calibration Curve (Figure 4)#####

Figure4A<-ggplotify::as.ggplot(~calPlot(list("METS-VF"=fgfit1,
                                   "DAAT"=fgfit2,
                                   "LAP"=fgfit5,
                                   "EVA"=fgfit6),
                              time=c(20),cores = 3,col = c(ggsci::pal_jama("default")(6))[c(1:2,5:6)],
                              showPseudo=FALSE,
                              type="risk",
                              data=base.2,
                              xlim=c(0,0.20),
                              ylim=c(0,0.20),
                              xlab = c("Probabilidad del Evento Predicho"),
                              ylab = c("Frecuencia de Eventos Observados")))+
  ggtitle("Todas las Causas Cardiovasculares")
                       

Figure4B<-ggplotify::as.ggplot(~calPlot(list("METS-VF"=fg.fit.1.cardio,
             "DAAT"=fg.fit.2.cardio,
             "LAP"=fg.fit.5.cardio,
             "EVA"=fg.fit.6.cardio),
        time=c(20),cores = 3,col = c(ggsci::pal_jama("default")(6))[c(1:2,5:6)],
        showPseudo=FALSE,
        type="risk",
        data=base.2,
        xlim=c(0,0.20),
        ylim=c(0,0.20),
        xlab = c("Probabilidad del Evento Predicho"),
        ylab = c("Frecuencia de Eventos Observados")))+
  ggtitle("Afecciones Cardíacas")


Figure4C<-ggplotify::as.ggplot(~calPlot(list("METS-VF"=fg.fit.1.stroke,
             "DAAT"=fg.fit.2.stroke,
             "LAP"=fg.fit.5.stroke,
             "EVA"=fg.fit.6.stroke),
        time=c(15),cores = 3,col = c(ggsci::pal_jama("default")(6))[c(1:2,5:6)],
        showPseudo=FALSE,
        type="risk",
        data=base.2,
        xlim=c(0,0.05),
        ylim=c(0,0.05),
        xlab = c("Probabilidad del Evento Predicho"),
        ylab = c("Frecuencia de Eventos Observados")))+
  ggtitle("Enfermedades Cerebrovasculares")

Figure4<-ggarrange(Figure4A,Figure4B,Figure4C,ncol = 3,nrow = 1)

ggsave(Figure4,
       filename = "Figura_4.png", 
       width = 40, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       bg = "white",
       limitsize = FALSE)

#####Tressholds Evaluation (Table 3)######

#Muerte All Vascular
#METS-VF
res.cut.all.1 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D015",variables = c("METSVF"))
summary(res.cut.all.1)

#Deep-Abdominal Adipose Tissue index
res.cut.all.2 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D015",variables = c("DAAT_index"))
summary(res.cut.all.2)

#Visceral Adiposity Index
res.cut.all.3 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D015",variables = c("VAI_index"))
summary(res.cut.all.3)

#GEA-VAI
res.cut.all.4 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D015",variables = c("VAI_GEA_index"))
summary(res.cut.all.4)

#Lipid Accumulation Product
res.cut.all.5 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D015",variables = c("LAAP_index"))
summary(res.cut.all.5)

#Eva Index
res.cut.all.6 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D015",variables = c("EVA_index"))
summary(res.cut.all.6)

#Depress
res.cut.all.7 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D015",variables = c("Depres_index"))
summary(res.cut.all.7)

###Cardiac Deaths
#METS-VF
res.cut.cardiac.1 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D003",variables = c("METSVF"))
summary(res.cut.cardiac.1)

#Deep-Abdominal Adipose Tissue index
res.cut.cardiac.2 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D003",variables = c("DAAT_index"))
summary(res.cut.cardiac.2)

#Visceral Adiposity Index
res.cut.cardiac.3 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D003",variables = c("VAI_index"))
summary(res.cut.cardiac.3)

#GEA-VAI
res.cut.cardiac.4 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D003",variables = c("VAI_GEA_index"))
summary(res.cut.cardiac.4)

#Lipid Accumulation Product
res.cut.cardiac.5 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D003",variables = c("LAAP_index"))
summary(res.cut.cardiac.5)

#Eva Index
res.cut.cardiac.6 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D003",variables = c("EVA_index"))
summary(res.cut.cardiac.6)

#Depress
res.cut.cardiac.7 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D003",variables = c("Depres_index"))
summary(res.cut.cardiac.7)

##Any Cerebrovascular Deaths
#METS-VF
res.cut.cerebrovasc.1 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D008",variables = c("METSVF"))
summary(res.cut.cerebrovasc.1)

#Deep-Abdominal Adipose Tissue index
res.cut.cerebrovasc.2 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D008",variables = c("DAAT_index"))
summary(res.cut.cerebrovasc.2)

#Visceral Adiposity Index
res.cut.cerebrovasc.3 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D008",variables = c("VAI_index"))
summary(res.cut.cerebrovasc.3)

#GEA-VAI
res.cut.cerebrovasc.4 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D008",variables = c("VAI_GEA_index"))
summary(res.cut.cerebrovasc.4)

#Lipid Accumulation Product
res.cut.cerebrovasc.5 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D008",variables = c("LAAP_index"))
summary(res.cut.cerebrovasc.5)

#Eva Index
res.cut.cerebrovasc.6 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D008",variables = c("EVA_index"))
summary(res.cut.cerebrovasc.6)

#Depress
res.cut.cerebrovasc.7 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D008",variables = c("Depres_index"))
summary(res.cut.cerebrovasc.7)

#Other CVD
#METS-VF
res.cut.other.1 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D009",variables = c("METSVF"))
summary(res.cut.other.1)

#Deep-Abdominal Adipose Tissue index
res.cut.other.2 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D009",variables = c("DAAT_index"))
summary(res.cut.other.2)

#Visceral Adiposity Index
res.cut.other.3 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D009",variables = c("VAI_index"))
summary(res.cut.other.3)

#GEA-VAI
res.cut.other.4 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D009",variables = c("VAI_GEA_index"))
summary(res.cut.other.4)

#Lipid Accumulation Product
res.cut.other.5 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D009",variables = c("LAAP_index"))
summary(res.cut.other.5)

#Eva Index
res.cut.other.6 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D009",variables = c("EVA_index"))
summary(res.cut.other.6)

#Depress
res.cut.other.7 <- survminer::surv_cutpoint(base.2, time = "PERSON_YEARS", event = "D009",variables = c("Depres_index"))
summary(res.cut.other.7)




#####Tressholds Evaluation - Figure (Figure 6)#####

##Muerte All Vascular
res.cat.all.1 <- surv_categorize(res.cut.all.1)
res.cat.all.2 <- surv_categorize(res.cut.all.2)
res.cat.all.3 <- surv_categorize(res.cut.all.3)
res.cat.all.4 <- surv_categorize(res.cut.all.4)
res.cat.all.5 <- surv_categorize(res.cut.all.5)
res.cat.all.6 <- surv_categorize(res.cut.all.6)

cox.cat.all.1<-coxph(Surv(PERSON_YEARS, D015)~METSVF, data=res.cat.all.1)
cox.cat.all.2<-coxph(Surv(PERSON_YEARS, D015)~DAAT_index, data=res.cat.all.2)
cox.cat.all.3<-coxph(Surv(PERSON_YEARS, D015)~VAI_index, data=res.cat.all.3)
cox.cat.all.4<-coxph(Surv(PERSON_YEARS, D015)~VAI_GEA_index, data=res.cat.all.4)
cox.cat.all.5<-coxph(Surv(PERSON_YEARS, D015)~LAAP_index, data=res.cat.all.5)
cox.cat.all.6<-coxph(Surv(PERSON_YEARS, D015)~EVA_index, data=res.cat.all.6)

c.stat.all.1<-paste0(round(cox.cat.all.1$concordance[6],2)," (",round(cox.cat.all.1$concordance[6]-(cox.cat.all.1$concordance[7]*1.96),2),"-",round(cox.cat.all.1$concordance[6]+(cox.cat.all.1$concordance[7]*1.96),2),")")
c.stat.all.2<-paste0(round(cox.cat.all.2$concordance[6],2)," (",round(cox.cat.all.2$concordance[6]-(cox.cat.all.2$concordance[7]*1.96),2),"-",round(cox.cat.all.2$concordance[6]+(cox.cat.all.2$concordance[7]*1.96),2),")")
c.stat.all.3<-paste0(round(cox.cat.all.3$concordance[6],2)," (",round(cox.cat.all.3$concordance[6]-(cox.cat.all.3$concordance[7]*1.96),2),"-",round(cox.cat.all.3$concordance[6]+(cox.cat.all.3$concordance[7]*1.96),2),")")
c.stat.all.4<-paste0(round(cox.cat.all.4$concordance[6],2)," (",round(cox.cat.all.4$concordance[6]-(cox.cat.all.4$concordance[7]*1.96),2),"-",round(cox.cat.all.4$concordance[6]+(cox.cat.all.4$concordance[7]*1.96),2),")")
c.stat.all.5<-paste0(round(cox.cat.all.5$concordance[6],2)," (",round(cox.cat.all.5$concordance[6]-(cox.cat.all.5$concordance[7]*1.96),2),"-",round(cox.cat.all.5$concordance[6]+(cox.cat.all.5$concordance[7]*1.96),2),")")
c.stat.all.6<-paste0(round(cox.cat.all.6$concordance[6],2)," (",round(cox.cat.all.6$concordance[6]-(cox.cat.all.6$concordance[7]*1.96),2),"-",round(cox.cat.all.6$concordance[6]+(cox.cat.all.6$concordance[7]*1.96),2),")")

##Any Cardiac Death
res.cat.cardiac.1 <- surv_categorize(res.cut.cardiac.1)
res.cat.cardiac.2 <- surv_categorize(res.cut.cardiac.2)
res.cat.cardiac.3 <- surv_categorize(res.cut.cardiac.3)
res.cat.cardiac.4 <- surv_categorize(res.cut.cardiac.4)
res.cat.cardiac.5 <- surv_categorize(res.cut.cardiac.5)
res.cat.cardiac.6 <- surv_categorize(res.cut.cardiac.6)

cox.cat.cardiac.1<-coxph(Surv(PERSON_YEARS, D003)~METSVF, data=res.cat.cardiac.1)
cox.cat.cardiac.2<-coxph(Surv(PERSON_YEARS, D003)~DAAT_index, data=res.cat.cardiac.2)
cox.cat.cardiac.3<-coxph(Surv(PERSON_YEARS, D003)~VAI_index, data=res.cat.cardiac.3)
cox.cat.cardiac.4<-coxph(Surv(PERSON_YEARS, D003)~VAI_GEA_index, data=res.cat.cardiac.4)
cox.cat.cardiac.5<-coxph(Surv(PERSON_YEARS, D003)~LAAP_index, data=res.cat.cardiac.5)
cox.cat.cardiac.6<-coxph(Surv(PERSON_YEARS, D003)~EVA_index, data=res.cat.cardiac.6)

c.stat.cardiac.1<-paste0(round(cox.cat.cardiac.1$concordance[6],2)," (",round(cox.cat.cardiac.1$concordance[6]-(cox.cat.cardiac.1$concordance[7]*1.96),2),"-",round(cox.cat.cardiac.1$concordance[6]+(cox.cat.cardiac.1$concordance[7]*1.96),2),")")
c.stat.cardiac.2<-paste0(round(cox.cat.cardiac.2$concordance[6],2)," (",round(cox.cat.cardiac.2$concordance[6]-(cox.cat.cardiac.2$concordance[7]*1.96),2),"-",round(cox.cat.cardiac.2$concordance[6]+(cox.cat.cardiac.2$concordance[7]*1.96),2),")")
c.stat.cardiac.3<-paste0(round(cox.cat.cardiac.3$concordance[6],2)," (",round(cox.cat.cardiac.3$concordance[6]-(cox.cat.cardiac.3$concordance[7]*1.96),2),"-",round(cox.cat.cardiac.3$concordance[6]+(cox.cat.cardiac.3$concordance[7]*1.96),2),")")
c.stat.cardiac.4<-paste0(round(cox.cat.cardiac.4$concordance[6],2)," (",round(cox.cat.cardiac.4$concordance[6]-(cox.cat.cardiac.4$concordance[7]*1.96),2),"-",round(cox.cat.cardiac.4$concordance[6]+(cox.cat.cardiac.4$concordance[7]*1.96),2),")")
c.stat.cardiac.5<-paste0(round(cox.cat.cardiac.5$concordance[6],2)," (",round(cox.cat.cardiac.5$concordance[6]-(cox.cat.cardiac.5$concordance[7]*1.96),2),"-",round(cox.cat.cardiac.5$concordance[6]+(cox.cat.cardiac.5$concordance[7]*1.96),2),")")
c.stat.cardiac.6<-paste0(round(cox.cat.cardiac.6$concordance[6],2)," (",round(cox.cat.cardiac.6$concordance[6]-(cox.cat.cardiac.6$concordance[7]*1.96),2),"-",round(cox.cat.cardiac.6$concordance[6]+(cox.cat.cardiac.6$concordance[7]*1.96),2),")")

#Cerebrovascular
res.cat.cerebrovasc.1 <- surv_categorize(res.cut.cerebrovasc.1)
res.cat.cerebrovasc.2 <- surv_categorize(res.cut.cerebrovasc.2)
res.cat.cerebrovasc.3 <- surv_categorize(res.cut.cerebrovasc.3)
res.cat.cerebrovasc.4 <- surv_categorize(res.cut.cerebrovasc.4)
res.cat.cerebrovasc.5 <- surv_categorize(res.cut.cerebrovasc.5)
res.cat.cerebrovasc.6 <- surv_categorize(res.cut.cerebrovasc.6)

cox.cat.cerebrovasc.1<-coxph(Surv(PERSON_YEARS, D008)~METSVF, data=res.cat.cerebrovasc.1)
cox.cat.cerebrovasc.2<-coxph(Surv(PERSON_YEARS, D008)~DAAT_index, data=res.cat.cerebrovasc.2)
cox.cat.cerebrovasc.3<-coxph(Surv(PERSON_YEARS, D008)~VAI_index, data=res.cat.cerebrovasc.3)
cox.cat.cerebrovasc.4<-coxph(Surv(PERSON_YEARS, D008)~VAI_GEA_index, data=res.cat.cerebrovasc.4)
cox.cat.cerebrovasc.5<-coxph(Surv(PERSON_YEARS, D008)~LAAP_index, data=res.cat.cerebrovasc.5)
cox.cat.cerebrovasc.6<-coxph(Surv(PERSON_YEARS, D008)~EVA_index, data=res.cat.cerebrovasc.6)

c.stat.cerebrovasc.1<-paste0(round(cox.cat.cerebrovasc.1$concordance[6],2)," (",round(cox.cat.cerebrovasc.1$concordance[6]-(cox.cat.cerebrovasc.1$concordance[7]*1.96),2),"-",round(cox.cat.cerebrovasc.1$concordance[6]+(cox.cat.cerebrovasc.1$concordance[7]*1.96),2),")")
c.stat.cerebrovasc.2<-paste0(round(cox.cat.cerebrovasc.2$concordance[6],2)," (",round(cox.cat.cerebrovasc.2$concordance[6]-(cox.cat.cerebrovasc.2$concordance[7]*1.96),2),"-",round(cox.cat.cerebrovasc.2$concordance[6]+(cox.cat.cerebrovasc.2$concordance[7]*1.96),2),")")
c.stat.cerebrovasc.3<-paste0(round(cox.cat.cerebrovasc.3$concordance[6],2)," (",round(cox.cat.cerebrovasc.3$concordance[6]-(cox.cat.cerebrovasc.3$concordance[7]*1.96),2),"-",round(cox.cat.cerebrovasc.3$concordance[6]+(cox.cat.cerebrovasc.3$concordance[7]*1.96),2),")")
c.stat.cerebrovasc.4<-paste0(round(cox.cat.cerebrovasc.4$concordance[6],2)," (",round(cox.cat.cerebrovasc.4$concordance[6]-(cox.cat.cerebrovasc.4$concordance[7]*1.96),2),"-",round(cox.cat.cerebrovasc.4$concordance[6]+(cox.cat.cerebrovasc.4$concordance[7]*1.96),2),")")
c.stat.cerebrovasc.5<-paste0(round(cox.cat.cerebrovasc.5$concordance[6],2)," (",round(cox.cat.cerebrovasc.5$concordance[6]-(cox.cat.cerebrovasc.5$concordance[7]*1.96),2),"-",round(cox.cat.cerebrovasc.5$concordance[6]+(cox.cat.cerebrovasc.5$concordance[7]*1.96),2),")")
c.stat.cerebrovasc.6<-paste0(round(cox.cat.cerebrovasc.6$concordance[6],2)," (",round(cox.cat.cerebrovasc.6$concordance[6]-(cox.cat.cerebrovasc.6$concordance[7]*1.96),2),"-",round(cox.cat.cerebrovasc.6$concordance[6]+(cox.cat.cerebrovasc.6$concordance[7]*1.96),2),")")

#Other
res.cat.other.1 <- surv_categorize(res.cut.other.1)
res.cat.other.2 <- surv_categorize(res.cut.other.2)
res.cat.other.3 <- surv_categorize(res.cut.other.3)
res.cat.other.4 <- surv_categorize(res.cut.other.4)
res.cat.other.5 <- surv_categorize(res.cut.other.5)
res.cat.other.6 <- surv_categorize(res.cut.other.6)

cox.cat.other.1<-coxph(Surv(PERSON_YEARS, D009)~METSVF, data=res.cat.other.1)
cox.cat.other.2<-coxph(Surv(PERSON_YEARS, D009)~DAAT_index, data=res.cat.other.2)
cox.cat.other.3<-coxph(Surv(PERSON_YEARS, D009)~VAI_index, data=res.cat.other.3)
cox.cat.other.4<-coxph(Surv(PERSON_YEARS, D009)~VAI_GEA_index, data=res.cat.other.4)
cox.cat.other.5<-coxph(Surv(PERSON_YEARS, D009)~LAAP_index, data=res.cat.other.5)
cox.cat.other.6<-coxph(Surv(PERSON_YEARS, D009)~EVA_index, data=res.cat.other.6)

c.stat.other.1<-paste0(round(cox.cat.other.1$concordance[6],2)," (",round(cox.cat.other.1$concordance[6]-(cox.cat.other.1$concordance[7]*1.96),2),"-",round(cox.cat.other.1$concordance[6]+(cox.cat.other.1$concordance[7]*1.96),2),")")
c.stat.other.2<-paste0(round(cox.cat.other.2$concordance[6],2)," (",round(cox.cat.other.2$concordance[6]-(cox.cat.other.2$concordance[7]*1.96),2),"-",round(cox.cat.other.2$concordance[6]+(cox.cat.other.2$concordance[7]*1.96),2),")")
c.stat.other.3<-paste0(round(cox.cat.other.3$concordance[6],2)," (",round(cox.cat.other.3$concordance[6]-(cox.cat.other.3$concordance[7]*1.96),2),"-",round(cox.cat.other.3$concordance[6]+(cox.cat.other.3$concordance[7]*1.96),2),")")
c.stat.other.4<-paste0(round(cox.cat.other.4$concordance[6],2)," (",round(cox.cat.other.4$concordance[6]-(cox.cat.other.4$concordance[7]*1.96),2),"-",round(cox.cat.other.4$concordance[6]+(cox.cat.other.4$concordance[7]*1.96),2),")")
c.stat.other.5<-paste0(round(cox.cat.other.5$concordance[6],2)," (",round(cox.cat.other.5$concordance[6]-(cox.cat.other.5$concordance[7]*1.96),2),"-",round(cox.cat.other.5$concordance[6]+(cox.cat.other.5$concordance[7]*1.96),2),")")
c.stat.other.6<-paste0(round(cox.cat.other.6$concordance[6],2)," (",round(cox.cat.other.6$concordance[6]-(cox.cat.other.6$concordance[7]*1.96),2),"-",round(cox.cat.other.6$concordance[6]+(cox.cat.other.6$concordance[7]*1.96),2),")")

##Table of CutOffs

df.cut.1<-data.frame(index_1=c("METS-VF","DAAT","VAI","DAI","LAAP","EVA"),
                     cutoff_1=c(round(as.numeric(res.cut.all.1$cutpoint[1]),2),
                                round(as.numeric(res.cut.all.2$cutpoint[1]),2),
                                round(as.numeric(res.cut.all.3$cutpoint[1]),2),
                                round(as.numeric(res.cut.all.4$cutpoint[1]),2),
                                round(as.numeric(res.cut.all.5$cutpoint[1]),2),
                                round(as.numeric(res.cut.all.6$cutpoint[1]),2)),
                     statistic_1=c(round(as.numeric(res.cut.all.1$cutpoint[2]),2),
                                   round(as.numeric(res.cut.all.2$cutpoint[2]),2),
                                   round(as.numeric(res.cut.all.3$cutpoint[2]),2),
                                   round(as.numeric(res.cut.all.4$cutpoint[2]),2),
                                   round(as.numeric(res.cut.all.5$cutpoint[2]),2),
                                   round(as.numeric(res.cut.all.6$cutpoint[2]),2)),
                     CI_95_1=c(rbind(c.stat.all.1,
                                   c.stat.all.2,
                                   c.stat.all.3,
                                   c.stat.all.4,
                                   c.stat.all.5,
                                   c.stat.all.6)))

df.cut.2<-data.frame(index_2=c("METS-VF","DAAT","VAI","DAI","LAAP","EVA"),
                     cutoff_2=c(round(as.numeric(res.cut.cardiac.1$cutpoint[1]),2),
                                round(as.numeric(res.cut.cardiac.2$cutpoint[1]),2),
                                round(as.numeric(res.cut.cardiac.3$cutpoint[1]),2),
                                round(as.numeric(res.cut.cardiac.4$cutpoint[1]),2),
                                round(as.numeric(res.cut.cardiac.5$cutpoint[1]),2),
                                round(as.numeric(res.cut.cardiac.6$cutpoint[1]),2)),
                     statistic_2=c(round(as.numeric(res.cut.cardiac.1$cutpoint[2]),2),
                                   round(as.numeric(res.cut.cardiac.2$cutpoint[2]),2),
                                   round(as.numeric(res.cut.cardiac.3$cutpoint[2]),2),
                                   round(as.numeric(res.cut.cardiac.4$cutpoint[2]),2),
                                   round(as.numeric(res.cut.cardiac.5$cutpoint[2]),2),
                                   round(as.numeric(res.cut.cardiac.6$cutpoint[2]),2)),
                     CI_95_2=c(rbind(c.stat.cardiac.1,
                                   c.stat.cardiac.2,
                                   c.stat.cardiac.3,
                                   c.stat.cardiac.4,
                                   c.stat.cardiac.5,
                                   c.stat.cardiac.6)))

df.cut.3<-data.frame(index_3=c("METS-VF","DAAT","VAI","DAI","LAAP","EVA"),
                     cutoff_3=c(round(as.numeric(res.cut.cerebrovasc.1$cutpoint[1]),2),
                                round(as.numeric(res.cut.cerebrovasc.2$cutpoint[1]),2),
                                round(as.numeric(res.cut.cerebrovasc.3$cutpoint[1]),2),
                                round(as.numeric(res.cut.cerebrovasc.4$cutpoint[1]),2),
                                round(as.numeric(res.cut.cerebrovasc.5$cutpoint[1]),2),
                                round(as.numeric(res.cut.cerebrovasc.6$cutpoint[1]),2)),
                     statistic_3=c(round(as.numeric(res.cut.cerebrovasc.1$cutpoint[2]),2),
                                   round(as.numeric(res.cut.cerebrovasc.2$cutpoint[2]),2),
                                   round(as.numeric(res.cut.cerebrovasc.3$cutpoint[2]),2),
                                   round(as.numeric(res.cut.cerebrovasc.4$cutpoint[2]),2),
                                   round(as.numeric(res.cut.cerebrovasc.5$cutpoint[2]),2),
                                   round(as.numeric(res.cut.cerebrovasc.6$cutpoint[2]),2)),
                     CI_95_3=c(rbind(c.stat.cerebrovasc.1,
                                   c.stat.cerebrovasc.2,
                                   c.stat.cerebrovasc.3,
                                   c.stat.cerebrovasc.4,
                                   c.stat.cerebrovasc.5,
                                   c.stat.cerebrovasc.6)))

df.cut.4<-data.frame(index_4=c("METS-VF","DAAT","VAI","DAI","LAAP","EVA"),
                     cutoff_4=c(round(as.numeric(res.cut.other.1$cutpoint[1]),2),
                                round(as.numeric(res.cut.other.2$cutpoint[1]),2),
                                round(as.numeric(res.cut.other.3$cutpoint[1]),2),
                                round(as.numeric(res.cut.other.4$cutpoint[1]),2),
                                round(as.numeric(res.cut.other.5$cutpoint[1]),2),
                                round(as.numeric(res.cut.other.6$cutpoint[1]),2)),
                     statistic_4=c(round(as.numeric(res.cut.other.1$cutpoint[2]),2),
                                   round(as.numeric(res.cut.other.2$cutpoint[2]),2),
                                   round(as.numeric(res.cut.other.3$cutpoint[2]),2),
                                   round(as.numeric(res.cut.other.4$cutpoint[2]),2),
                                   round(as.numeric(res.cut.other.5$cutpoint[2]),2),
                                   round(as.numeric(res.cut.other.6$cutpoint[2]),2)),
                     CI_95_4=c(rbind(c.stat.other.1,
                                   c.stat.other.2,
                                   c.stat.other.3,
                                   c.stat.other.4,
                                   c.stat.other.5,
                                   c.stat.other.6)))


df.cut<-cbind(df.cut.1,df.cut.2[,-1],df.cut.3[,-1],df.cut.4[,-1])
df.cut<-df.cut %>% 
  rename("Indice de Grasa Visceral" = "index_1",
         "Punto de Corte\n(Cualquier Causa Cardiovasculares)" = "cutoff_1",
         "Estadístico-Log-Rank\n(Cualquier Causa Cardiovasculares)" = "statistic_1",
         "Estadístico-C\n(Cualquier Causa Cardiovasculares)" = "CI_95_1",
         "Punto de Corte\n(Afecciones Cardiacas)" = "cutoff_2",
         "Estadístico-Log-Rank\n(Afecciones Cardiacas)" = "statistic_2",
         "Estadístico-C\n(Afecciones Cardiacas)" = "CI_95_2",
         "Punto de Corte\n(Enfermedades Cerebrovasculares)" = "cutoff_3",
         "Estadístico-Log-Rank\n(Enfermedades Cerebrovasculares)" = "statistic_3",
         "Estadístico-C\n(Enfermedades Cerebrovasculares)" = "CI_95_3",
         "Punto de Corte\n(Otras Causas Cardiovasculares)" = "cutoff_4",
         "Estadístico-Log-Rank\n(Otras Causas Cardiovasculares)" = "statistic_4",
         "Estadístico-C\n(Otras Causas Cardiovasculares)" = "CI_95_4")

Table1_Flex<-flextable::align(flextable::flextable(df.cut,cwidth=7),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table1_Flex,path="Table_3.docx")


#####Distribution of Visceral Fat Indexes (Supplementary Figure 1)#####

#METS-VF
Sup.Fig.1A<-base.2 %>% 
  ggplot(aes(x = METSVF)) + 
  geom_histogram(aes(y=..density..),fill="#e79b33",col="white", binwidth = 0.05) +
  geom_density(col="#29498d")+
  labs(title = "Metabolic Score for Visceral Fat")+
  xlab ("METS-VF") +
  ylab ("Densidad de área") +
  theme_classic()

Sup.Fig.1A<-egg::tag_facet(Sup.Fig.1A, x = Inf, y = Inf, 
               hjust = 1, open = "", close = "",
               tag_pool = c(paste0("Media (D.E.)= ", paste0(round(mean(base.2$METSVF),2)," (",round(sd(base.2$METSVF),2),")"),
                                   "\nMediana (R.I.Q.)= ",  paste0(round(median(base.2$METSVF),2)," (IQR: ",round(quantile(base.2$METSVF,probs = c(0.25)),2),"-",round(quantile(base.2$METSVF,probs = c(0.75)),2),")"),
                                   "\nn= ",nrow(base.2))))

#DAAT
Sup.Fig.1B<-base.2 %>% 
  ggplot(aes(x = DAAT_index)) + 
  geom_histogram(aes(y=..density..),fill="#e79b33",col="white", binwidth = 10) +
  geom_density(col="#29498d")+
  labs(title = "Deep-abdominal-adipose-tissue")+
  xlab ("DAAT") +
  ylab ("Densidad de área") +
  theme_classic()


Sup.Fig.1B<-egg::tag_facet(Sup.Fig.1B, x = Inf, y = Inf, 
                           hjust = 1, open = "", close = "",
                           tag_pool = c(paste0("Media (D.E.)= ", paste0(round(mean(base.2$DAAT_index),2)," (",round(sd(base.2$DAAT_index),2),")"),
                                               "\nMediana (R.I.Q.)= ",  paste0(round(median(base.2$DAAT_index),2)," (IQR: ",round(quantile(base.2$DAAT_index,probs = c(0.25)),2),"-",round(quantile(base.2$DAAT_index,probs = c(0.75)),2),")"),
                                               "\nn= ",nrow(base.2))))

#VAI
Sup.Fig.1C<-base.2 %>% 
  ggplot(aes(x = VAI_index)) + 
  geom_histogram(aes(y=..density..),fill="#e79b33",col="white", binwidth = 0.1) +
  geom_density(col="#29498d")+
  labs(title = "Visceral Adiposity Index ")+
  xlab ("VAI") +
  ylab ("Densidad de área") +
  theme_classic()


Sup.Fig.1C<-egg::tag_facet(Sup.Fig.1C, x = Inf, y = Inf, 
                           hjust = 1, open = "", close = "",
                           tag_pool = c(paste0("Media (D.E.)= ", paste0(round(mean(base.2$VAI_index),2)," (",round(sd(base.2$VAI_index),2),")"),
                                               "\nMediana (R.I.Q.)= ",  paste0(round(median(base.2$VAI_index),2)," (IQR: ",round(quantile(base.2$VAI_index,probs = c(0.25)),2),"-",round(quantile(base.2$VAI_index,probs = c(0.75)),2),")"),
                                               "\nn= ",nrow(base.2))))

#VAI-GEA
Sup.Fig.1D<-base.2 %>% 
  ggplot(aes(x = VAI_GEA_index)) + 
  geom_histogram(aes(y=..density..),fill="#e79b33",col="white", binwidth = 0.1) +
  geom_density(col="#29498d")+
  labs(title = "Dysfunctional adiposity index ")+
  xlab ("DAI") +
  ylab ("Densidad de área") +
  theme_classic()


Sup.Fig.1D<-egg::tag_facet(Sup.Fig.1D, x = Inf, y = Inf, 
                           hjust = 1, open = "", close = "",
                           tag_pool = c(paste0("Media (D.E.)= ", paste0(round(mean(base.2$VAI_GEA_index),2)," (",round(sd(base.2$VAI_GEA_index),2),")"),
                                               "\nMediana (R.I.Q.)= ",  paste0(round(median(base.2$VAI_GEA_index),2)," (IQR: ",round(quantile(base.2$VAI_GEA_index,probs = c(0.25)),2),"-",round(quantile(base.2$VAI_GEA_index,probs = c(0.75)),2),")"),
                                               "\nn= ",nrow(base.2))))


#LAAP
Sup.Fig.1E<-base.2 %>% 
  ggplot(aes(x = LAAP_index)) + 
  geom_histogram(aes(y=..density..),fill="#e79b33",col="white", binwidth = 2.5) +
  geom_density(col="#29498d")+
  labs(title = "Deep-abdominal-adipose-tissue ")+
  xlab ("DAAT") +
  ylab ("Densidad de área") +
  theme_classic()


Sup.Fig.1E<-egg::tag_facet(Sup.Fig.1E, x = Inf, y = Inf, 
                           hjust = 1, open = "", close = "",
                           tag_pool = c(paste0("Media (D.E.)= ", paste0(round(mean(base.2$LAAP_index),2)," (",round(sd(base.2$LAAP_index),2),")"),
                                               "\nMediana (R.I.Q.)= ",  paste0(round(median(base.2$LAAP_index),2)," (IQR: ",round(quantile(base.2$LAAP_index,probs = c(0.25)),2),"-",round(quantile(base.2$LAAP_index,probs = c(0.75)),2),")"),
                                               "\nn= ",nrow(base.2))))

#EVA
Sup.Fig.1F<-base.2 %>% 
  ggplot(aes(x = EVA_index)) + 
  geom_histogram(aes(y=..density..),fill="#e79b33",col="white", binwidth = 5) +
  geom_density(col="#29498d")+
  labs(title = "Estimate of Visceral Adipose Tissue Area ")+
  xlab ("EVA") +
  ylab ("Densidad de área") +
  theme_classic()


Sup.Fig.1F<-egg::tag_facet(Sup.Fig.1F, x = Inf, y = Inf, 
                           hjust = 1, open = "", close = "",
                           tag_pool = c(paste0("Media (D.E.)= ", paste0(round(mean(base.2$EVA_index),2)," (",round(sd(base.2$EVA_index),2),")"),
                                               "\nMediana (R.I.Q.)= ",  paste0(round(median(base.2$EVA_index),2)," (IQR: ",round(quantile(base.2$EVA_index,probs = c(0.25)),2),"-",round(quantile(base.2$EVA_index,probs = c(0.75)),2),")"),
                                               "\nn= ",nrow(base.2))))

Sup.Fig.1<-ggarrange(Sup.Fig.1A,Sup.Fig.1B,Sup.Fig.1C,Sup.Fig.1D,Sup.Fig.1E,Sup.Fig.1F,ncol = 3,nrow = 2,labels = LETTERS[1:6])

ggsave(Sup.Fig.1,
       filename = "Figura_Suplemenetaria_1.png", 
       width = 35, 
       height = 20,
       units=c("cm"),
       dpi = 300,
       bg = "white",
       limitsize = FALSE)

#####Proportion of Deaths in MCPS (Supplementary Figure 2)#####

df<- data.frame(CLUSTERS=c("Afección Cardiaca",
                           "Afección Cardiaca",
                           "Enfermedad Cerebrovascular",
                           "Enfermedad Cerebrovascular",
                           "Enfermedad Cerebrovascular",
                           "Enfermedad Cerebrovascular",
                           "Otras Causas Cardiovasculares",
                           "Otras Causas Cardiovasculares",
                           "Otras Causas Cardiovasculares"),
                prevalence=c(table(base.2$D001)[2],
                             table(base.2$D002)[2],
                             table(base.2$D004)[2],
                             table(base.2$D005)[2],
                             table(base.2$D006)[2],
                             table(base.2$D007)[2],
                             table(base.2$D009)[2],
                             table(base.2$D010)[2],
                             table(base.2$D011)[2]),
                hemisphere=c("Afección Cardiaca",
                             "Afección Cardiaca",
                             "Enfermedad Cerebrovascular",
                             "Enfermedad Cerebrovascular",
                             "Enfermedad Cerebrovascular",
                             "Enfermedad Cerebrovascular",
                             "Otras Causas Cardiovasculares",
                             "Otras Causas Cardiovasculares",
                             "Otras Causas Cardiovasculares"),
                labels=c(paste0("Cardiopatía Isquémica\n n = ",table(base.2$D001)[2]),
                         paste0("Cardiaca No especificada\n n = ",table(base.2$D002)[2]),
                         paste0("EVC Isquémico\n n = ",table(base.2$D004)[2]),
                         paste0("EVC Hemorragico\n n = ",table(base.2$D005)[2]),
                         paste0("EVC Desconocido\n n = ",table(base.2$D006)[2]),
                         paste0("Otro Tipo de EVC\n n = ",table(base.2$D007)[2]),
                         paste0("Tromboembolismo\n n = ",table(base.2$D009)[2]),
                         paste0("Enf. Arterial Periférica\n n = ",table(base.2$D010)[2]),
                         paste0("Otras Muertes No Especificadas\n n = ",table(base.2$D011)[2])))


Sup.Figure.2<-ggplot2::ggplot(df,aes(area=prevalence,
                       fill=CLUSTERS,
                       label=labels,
                       subgroup = hemisphere))+
  treemapify::geom_treemap(layout="squarified")+
  geom_treemap_subgroup_border(colour = "white")+
  geom_treemap_text(place = "centre",size = 14,colour = "white")+
  labs(fill="Grupo de ECV",title = "Mortalidad por Causas Cardiovasculares",subtitle = "Cohorte Prospectiva de la Ciudad de México")+
  ggsci::scale_fill_jama()+
  theme_classic()+
  xlab("")+
  ylab("")+
  theme(legend.position = "top")

ggsave(Sup.Figure.2,
       filename = "Figura_Suplemenetaria_2.png", 
       width = 35, 
       height = 20,
       units=c("cm"),
       dpi = 300,
       bg = "white",
       limitsize = FALSE)

#Cardiac deaths

mcps_lexis <- Lexis(entry = list("period" = date_recruited, 
                                 "age" = AGE, "time_at_entry" = 0),
                    exit = list("period" = date_recruited + PERSON_YEARS), 
                    data = mcps_fin, exit.status = D003) ####

mcps_fin2<-splitLexis(mcps_lexis, breaks=c(35,40,45,50,55,60,65,70,75,80,85), time.scale = "age")
mcps_fin2$age.cat <- timeBand(mcps_fin2, "age", type = "factor")
mcps_fin2$premature<-ifelse(mcps_fin2$lex.Xst==T,1,0); mcps_fin2$premature[mcps_fin2$premature==1 & mcps_fin2$ageout>=75]<-0 

cor.test(base.2$AGE,base.2$PAS_PROM)
base.ex<-base.2%>%
  sample_n(1000)

ggscatter(base.ex, x = "AGE", y = "PAS_PROM", 
          cor.coef = FALSE, cor.method = "spearman",
          add = "reg.line", conf.int = TRUE,
          xlab = "Edad, (Años)", ylab = "Presión Arterial Sistólica, (mmHg)",
          title = "Asociación entre Presión Arterial Sistólica y Edad",
          palette = "lancet",color = "red",size = 1)

ggplot(base.ex, aes(AGE, PAS_PROM))+
  geom_point(col="black")+
  theme_classic()+
  xlab("Edad, (Años)")+
  ylab("Presión Arterial Sistólica, (mmHg)")+
  geom_smooth(method='lm', formula= y~x,col="red")+
  labs(title = "Asociación entre Presión Arterial Sistólica y Edad",
       subtitle = "Cohorte de la Ciudad de México",
       caption  = "n=1000")+  # Plot regression slope
  geom_segment(aes(xend = AGE, yend = predicted), alpha = 2.0,col="green") +  # alpha to fade lines
  geom_point() +
  geom_point(aes(y = predicted), shape = 1)
  
base.ex$AGE_10<-(base.ex$AGE/10)
lm1<-lm(PAS_PROM~AGE,data = base.ex)
summary(lm1)
confint(lm1)
base.ex$predicted <- predict(lm1)   # Save the predicted values
base.ex$residuals <- residuals(lm1) # Save the residual values


lm2<-lm(PAS_PROM~AGE_10,data = base.ex)
summary(lm2)


pred_lines = data.frame(x=pred_x,
                        y=predict(my_lm, data.frame(x=pred_x)),
                        obs_Or_Pred=rep(c("Obs","Pred"), each=2))

ggplot(pred_lines, aes(x, y, colour=obs_Or_Pred, shape=obs_Or_Pred, linetype=obs_Or_Pred)) +
  geom_point(data=data.frame(x,y, obs_Or_Pred="Obs"), size=3) +
  geom_line(size=1) +
  scale_shape_manual(values=c(16,NA)) +
  theme_bw()

base$HAS_FINAL<-factor(base$HAS_FINAL,labels = c("NO-HTA","HTA"))
base$D015<-factor(base$D015,labels = c("NO-MUERTE-CVD","MUERTE-CVD"))
table(base$HAS_FINAL,base$D015)
base$MALE
base$IMC
m1<-glm(D015~HAS_FINAL,data = base,family = "binomial")
summary(m1)
jtools::summ(m1,exp=T)

