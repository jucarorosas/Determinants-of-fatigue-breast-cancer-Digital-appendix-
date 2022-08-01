
##############################################################################
##########                                                      ##############
########## Determinants of fatigue and longitudinal changes up  ##############
############      to 2 years post-radiotherapy among       ###################
###############     breast cancer patients in REQUITE    #####################
#################                                       ######################
##################            Main R script           ########################
####################        Date: 01.08.2022       ###########################     
##############################################################################


Sys.setenv(LANG = "en")

library(tidyverse)
library(dplyr)
library(mice)
library(foreign)
library(pastecs)
library(car)
library(Amelia)
library(Zelig)
library(table1)
library(arsenal)
library(ggplot2)
library(PROscorer)
library(writexl)
library(readr)
library(readxl)
library(knitr)
library(janitor)
library(kableExtra)
library(MASS)
library(jtools)
library(Rcpp)
library(huxtable)
library(lme4)
library(gtsummary)
library(performance)
library(car)
library(broom)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(gridExtra)
library(forestmodel)
library(ggstance)
library(broom.mixed)


# Load the seven databases we are going to use 
MFI_data <- read_csv("MFI_2021-04-09.csv")
Baseline_data <- read_excel("B2a_Breast_patient_factors_baseline_2021-04-09.xlsx")
Breastclinical_data <- read_excel ("B3_Breast_clinical_and_treatment_data_form_2021-04-09.xlsx")
GPAQ_data <- read_excel("GPAQ_2021-04-09.xlsx")
C30_data <- read_excel("C30_2021-04-09.xlsx")
BR23_data <- read_excel("BR23_2021-04-09.xlsx")
SD_data <- read_excel("Subject_details_2021-04-09.xlsx")
Rad_variables_ntd <- read_excel("Radiation_Variables_ntd_2021-10-11.xlsx")
 

# create duplicate work-databases 
MFIdat <- data.frame(MFI_data)
BDdat <- data.frame(Baseline_data)
B3BCdat <- data.frame(Breastclinical_data)
GPAQdat <- data.frame(GPAQ_data)
C30dat <- data.frame(C30_data)
BR23dat <- data.frame(BR23_data)
SDdat <- data.frame(SD_data)
Radntd <- data.frame(Rad_variables_ntd)

# Reverse scoring for MFI_data

MFI_data_rev <- MFIdat
MFI_data_rev[, c("MFI_q02_do_little", "MFI_q05_tired", "MFI_q09_dread_doing_things", 
                 "MFI_q10_do_very_little", "MFI_q13_effort_concentration", 
                 "MFI_q14_bad_condition", "MFI_q16_tire_easily", 
                 "MFI_q17_get_little_done", "MFI_q18_dont_feel_do_anything", 
                 "MFI_q19_thoughts_wander")] <-  6 - MFIdat [, c("MFI_q02_do_little", "MFI_q05_tired", "MFI_q09_dread_doing_things", 
                                                                 "MFI_q10_do_very_little", "MFI_q13_effort_concentration", 
                                                                 "MFI_q14_bad_condition", "MFI_q16_tire_easily", 
                                                                 "MFI_q17_get_little_done", "MFI_q18_dont_feel_do_anything", 
                                                                 "MFI_q19_thoughts_wander")]

# Checking if the scores where correctly reversed in the new data.frame "MFI_data_rev"
head(MFI_data$MFI_q02_do_little)
head(MFI_data_rev$MFI_q02_do_little)


# Rename C30dat column names (it is necessary to rename the columns to use the 
# package "PROscorer")

C30dat <- setNames(C30dat, c("Subject.Id","Visit","event_date", "submitdate", "q1", "q2", 
                             "q3", "q4", "q5", "q6", "q7", "q8", "q9", "q10", "q11", 
                             "q12", "q13", "q14", "q15", "q16", "q17", "q18", "q19", 
                             "q20", "q21", "q22", "q23", "q24", "q25", "q26", 
                             "q27", "q28", "q29", "q30"))


## SELECT VARIABLES AND CREATE NEW DATAFRAME

# First, we are going to create a new MFI data.frame with only the 24-month 
# follow-up. For this we use the data frame with the reverse scores "MFI_data_rev"

MFI_data_rev %>% 
  filter(MFIdat$Visit == "24m follow-up") -> MFI_24month

length(MFI_24month$Subject.Id) # 1159 patients at 24mths with MFI scores

# Also, create new data.frames with only baseline for GPAQ, C30, and BR-23
GPAQdat %>% 
  filter(GPAQdat$Visit == "Baseline") -> GPAQ_baseline

length(GPAQ_baseline$Subject.Id) # 1261 patients with GPAQ at baseline

C30dat %>% 
  filter(C30dat$Visit == "Baseline") -> C30_baseline

length(C30_baseline$Subject.Id) # 2025 patients with C30 at baseline 

BR23dat %>% 
  filter(BR23dat$Visit == "Baseline") -> BR23_baseline

length(BR23_baseline$Subject.Id) # 2025 patients with BR23 at baseline

# Select the relevant variables from each data frame and create
# smaller data frames. 

GPAQ_baseline %>%
  dplyr::select(Subject.Id, GPAQ_q16_sitting_recline_min) -> GPAQ_small

C30_baseline %>%
  dplyr::select(Subject.Id, q8, q9, q11, q19, q29, q30) -> C30_small

BR23_baseline %>%
  dplyr::select(Subject.Id, BR23_q47_pain_arm_shoulder, 
                BR23_q48_swollen_arm_hand, 
         BR23_q49_difficult_raise_arm, BR23_q51_swollen_breast)-> BR23_small

SDdat %>%
  dplyr::select(Subject.Id, Mainsite)-> SDdat_small

Radntd %>%
  dplyr::select(Subject.Id, b3radio_breast_fractions_dose_pe, ntd_tc_fibrose,  
                b3radio_elec_boost_dose_precise_, b3radio_photon_boostdose_precise, 
                b3radio_breast_dose_Gy, b3radio_boost) -> Rad_small



## MODIFY VARIABLES 


# Cleaning Rad_small
# Change "." for "0" in "b3radio_elec_boost_dose_precise_" and "b3radio_photon_boostdose_precise"

Rad_small$b3radio_elec_boost_dose_precise_[Rad_small$b3radio_elec_boost_dose_precise_ == "."] <- "0"

Rad_small$b3radio_photon_boostdose_precise[Rad_small$b3radio_photon_boostdose_precise == "."] <- "0"

## Convert those two variables from character to numeric 

Rad_small$b3radio_elec_boost_dose_precise_ <- as.numeric(Rad_small$b3radio_elec_boost_dose_precise_)

Rad_small$b3radio_photon_boostdose_precise <- as.numeric(Rad_small$b3radio_photon_boostdose_precise)

str(Rad_small) ## Check if that worked 
  
# Convert b3radio_boost to factor

Rad_small$b3radio_boost <- as.factor(Rad_small$b3radio_boost)
levels(Rad_small$b3radio_boost)


# New dichotomous "fractionation scheme" variable

Rad_small$fractionation_bin <- rep(NA, nrow(Rad_small))

Rad_small$fractionation_bin <- ifelse((Rad_small$b3radio_breast_fractions_dose_pe <= 2), "Standard fractionation", 
                                      "Hypofractionation")

Rad_small$fractionation_bin <- as.factor(Rad_small$fractionation_bin)
table(Rad_small$fractionation_bin)

# Create new variable "ntd_tc_cat" from ntd_tc_fibrose (numeric). 
# 1= (Gy <50), 2= (50<=Gy<60), 3= (60<=Gy<70), 4= Gy >=70 

Rad_small$ntd_tc_cat <- rep(NA, nrow(Rad_small))
Rad_small$ntd_tc_cat[Rad_small$ntd_tc_fibrose <50] <- "< 50 Gy" 
Rad_small$ntd_tc_cat[Rad_small$ntd_tc_fibrose>=50 & 
                       Rad_small$ntd_tc_fibrose<60 ] <- "50<=Gy<60"
Rad_small$ntd_tc_cat[Rad_small$ntd_tc_fibrose>=60 & 
                       Rad_small$ntd_tc_fibrose<70 ] <- "60<=Gy<70"
Rad_small$ntd_tc_cat[Rad_small$ntd_tc_fibrose>=70] <- "Gy>=70"

Rad_small$ntd_tc_cat <- factor(Rad_small$ntd_tc_cat)
table(Rad_small$ntd_tc_cat)


# Create new variable "boost_frac" with "b3radio_boost" and "fractionation_bin" 

Rad_small$boost_frac <- rep(NA, nrow(Rad_small))

Rad_small$boost_frac <- ifelse((Rad_small$fractionation_bin  == "Standard fractionation") & (Rad_small$b3radio_boost == "0"), 0,
                        + ifelse((Rad_small$fractionation_bin  == "Hypofractionation") & (Rad_small$b3radio_boost == "0"), 1,
                        + ifelse((Rad_small$fractionation_bin  == "Standard fractionation") & (Rad_small$b3radio_boost == "1"), 2,
                        + ifelse((Rad_small$fractionation_bin  == "Hypofractionation") & (Rad_small$b3radio_boost == "1"), 3, NA))))

Rad_small$boost_frac <- as.factor(Rad_small$boost_frac) # convert to factor
table(Rad_small$boost_frac) # verify levels

levels(Rad_small$boost_frac) # verify levels�order

levels(Rad_small$boost_frac) <- c("No boost & standard", "No boost & hypo", 
                                  "Boost & standard", "Boost & hypo") # re-label levels 

table(Rad_small$boost_frac) # verify new labels 


# New variable "radio_total_dose_sum" adding up the total dose, the electron 
# boost dose and the proton boost dose

Rad_small$radio_total_dose_sum <- Rad_small$b3radio_breast_dose_Gy +
  Rad_small$b3radio_elec_boost_dose_precise_ + 
  Rad_small$b3radio_photon_boostdose_precise

summary(Rad_small$radio_total_dose_sum)

# New dichotomous "histology_di" variable
B3BCdat$histology_di <- car::recode(B3BCdat$b3tumour_histological_type, 
                                    "3 = 'in situ'; c(1, 2, 4, 5) = 'invasive'")
table(B3BCdat$b3tumour_histological_type, B3BCdat$histology_di)

fakeTable <- B3BCdat %>% 
  tabyl( b3tumour_histological_type, histology_di ) %>% 
  adorn_totals( where = c("row", "col") ) %>%
  adorn_percentages("row") %>%
  adorn_pct_formatting() %>%
  adorn_ns( position = "front" ) %>%
  adorn_title("combined")

fakeTable %>%
  kable() %>%
  kable_styling(bootstrap_options = c("condensed", "striped", "bordered")) 


# Create new variable "Ethnicity_bin" from Ethnicity (categorical). 
# "0" for white and "1" for non-white
class(BDdat$ethnicity)
table(BDdat$ethnicity)

BDdat$ethnicity_bin <- rep(NA, nrow(BDdat))
BDdat$ethnicity_bin[BDdat$ethnicity==1] <- "white"
BDdat$ethnicity_bin[BDdat$ethnicity %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                           13, 14, 15, 16, 17, 19)] <- "other" 
table(BDdat$ethnicity_bin)

BDdat$ethnicity_bin <- factor(BDdat$ethnicity_bin, levels = c("white", "other")) 
table(BDdat$ethnicity_bin)



# Create new variable "Age_cat" from Age (numeric). 
# 1= age <50 , 2= (50<=age<60), 3= (60<=age<70), 4= age>=70
BDdat$age_cat <- rep(NA, nrow(BDdat))
BDdat$age_cat[BDdat$age_at_radiotherapy_start_yrs<50] <- "<50" 
BDdat$age_cat[BDdat$age_at_radiotherapy_start_yrs>=50 & 
                BDdat$age_at_radiotherapy_start_yrs<60 ] <- "50-59"
BDdat$age_cat[BDdat$age_at_radiotherapy_start_yrs>=60 & 
                BDdat$age_at_radiotherapy_start_yrs<70 ] <- "60-69"
BDdat$age_cat[BDdat$age_at_radiotherapy_start_yrs>=70] <- ">=70"

BDdat$age_cat <- factor(BDdat$age_cat, levels = c("<50", "50-59", "60-69", 
                                                  ">=70"))
table(BDdat$age_cat)


# Create new variable "household_bin" from household_members (numeric). 
# "0" for alone and "1" for the 2 or more
class(BDdat$household_members)
table(BDdat$household_members)

BDdat$household_bin <- rep(NA, nrow(BDdat))
BDdat$household_bin[BDdat$household_members==1] <- "living alone"
BDdat$household_bin[BDdat$household_members %in% c(2, 3, 4, 5, 6, 7, 8)] <- "not living alone" 

table(BDdat$household_bin)

# Create new variable "educational_lmh" from education_profession (numeric). 
# "1" for low, "2" for medium,  "3" for high, "4" for others.  
class(BDdat$education_profession)
table(BDdat$education_profession)

table(BDdat$education_profession, BDdat$secondary_school_selection)

BDdat$education_lmh <- rep(NA, nrow(BDdat))
BDdat$education_lmh[BDdat$education_profession %in% 1 | 
                      BDdat$secondary_school_selection %in% c(5, 13)] <- "low" 

BDdat$education_lmh[BDdat$education_profession==3 |
                      BDdat$education_profession==5 |
                      BDdat$secondary_school_selection 
                    %in% c(1, 3, 4, 6, 8, 9, 11, 14, 16, 17)] <- "medium"

BDdat$education_lmh[BDdat$education_profession==4 | 
                      BDdat$secondary_school_selection %in% 
                      c(2, 7, 10, 12)] <- "high"

BDdat$education_lmh <- factor(BDdat$education_lmh, levels = c("low", "medium", "high"))

table(SDdat$Mainsite, BDdat$secondary_school_selection)
table(SDdat$Mainsite, BDdat$education_profession)
table(BDdat$education_lmh, BDdat$secondary_school_selection)
table(BDdat$education_profession)
table(BDdat$education_profession, BDdat$education_lmh)

# Create bmi variable from weight and height. 
BDdat$bmi <- BDdat$weight_at_cancer_diagnosis_kg/(BDdat$height_cm/100)^2 

class(BDdat$bmi)
head(BDdat$bmi)
summary(BDdat$bmi)



# Create new variable "bmi_cat" from bmi (numeric). 
# "1" for BMI<25 kg/m², "2" for (25,0<=BMI<=29,9 kg/m²),  
# "3" for (BMI >=30,0 kg/m²)
BDdat$bmi_cat <- rep(NA, nrow(BDdat))
BDdat$bmi_cat[BDdat$bmi<25] <- "<25" 
BDdat$bmi_cat[BDdat$bmi>=25 & BDdat$bmi<30 ] <- "25-29"
BDdat$bmi_cat[BDdat$bmi>=30] <- ">=30"

BDdat$bmi_cat <- factor(BDdat$bmi_cat, levels = c("<25", "25-29", ">=30"))

table(BDdat$bmi_cat)


# Create new variable "smoking_bin" from smoker (numeric). 
# "1" for former/never, "2" for current. 

class(BDdat$smoker)
table(BDdat$smoker)

BDdat$smoker_bin <- rep(NA, nrow(BDdat))
BDdat$smoker_bin[BDdat$smoker %in% c(0, 1)] <- "never/former"
BDdat$smoker_bin[BDdat$smoker %in% c(2,3)] <- "current"

table(BDdat$smoker_bin)

# For alcohol, first, we are going to replace all NAs in 
# "alcohol_current_consumption" to 0 if they answered "never" 
# for the variable "alcohol_intake".

BDdat$alcohol_current_consumption <- coalesce(BDdat$alcohol_current_consumption, 
                                              BDdat$alcohol_previous_consumption)

BDdat$alcohol_current_consumption <- ifelse(BDdat$alcohol_intake == 0 & 
                                              is.na(BDdat$alcohol_current_consumption), 
                                            0, BDdat$alcohol_current_consumption)   


# Create new variable "alcohol_bin" from alcohol_current_consumption

class(BDdat$alcohol_current_consumption)
summary(BDdat$alcohol_current_consumption)

BDdat$alcohol_bin <- rep(NA, nrow(BDdat))
BDdat$alcohol_bin[BDdat$alcohol_current_consumption<8] <- "low/moderate" 
BDdat$alcohol_bin[BDdat$alcohol_current_consumption>=8] <- "high"

BDdat$alcohol_bin <- factor(BDdat$alcohol_bin, levels = c("low/moderate", "high")) 

table(BDdat$alcohol_bin)
table(BDdat$alcohol_intake, BDdat$alcohol_current_consumption) #  Only 212 (from 878) people 
# that said they "never" drink alcohol, answered "0" drinks for the other question. 
summary(BDdat$alcohol_bin)



# Create new variable "systemic_cat" from sys_treatment, tamoxifen and aromatase

B3BCdat$systemic_cat <- rep(NA, nrow(B3BCdat))

B3BCdat$systemic_cat <- ifelse((B3BCdat$b3sys_treatment == "0"), 0,
                               + ifelse((B3BCdat$b3sys_tamoxifen == "0") & 
                                          (B3BCdat$b3sys_aromatase == "0"), 0,
                                        + ifelse((B3BCdat$b3sys_treatment == "1") & 
                                                   (B3BCdat$b3sys_tamoxifen == "1") & 
                                                   (B3BCdat$b3sys_aromatase == "0"), 1,
                                                 + ifelse((B3BCdat$b3sys_treatment == "1") & 
                                                            (B3BCdat$b3sys_tamoxifen == "0") & 
                                                            (B3BCdat$b3sys_aromatase == "1"), 2, NA))))

table(B3BCdat$systemic_cat)

B3BCdat$systemic_cat <- as.factor(B3BCdat$systemic_cat)
levels(B3BCdat$systemic_cat) <- c("none", "tamoxifen", "aromatase inhibitors")

table(B3BCdat$systemic_cat)


# Create new variable "pato_t_cat" from "b3pathologic_uicc_t_stage"

B3BCdat$b3pathologic_uicc_t_stage <- coalesce(B3BCdat$b3pathologic_uicc_t_stage, 
                                              B3BCdat$b3neo_adj_uicc_t_stage)

B3BCdat$pato_t_cat <- rep(NA, nrow(B3BCdat))

B3BCdat$pato_t_cat <- ifelse((B3BCdat$b3pathologic_uicc_t_stage == "Tis") |
                               (B3BCdat$b3pathologic_uicc_t_stage == "T1a") | 
                                          (B3BCdat$b3pathologic_uicc_t_stage == "T1b") | 
                                          (B3BCdat$b3pathologic_uicc_t_stage == "T1c"), 0,
                                                 + ifelse((B3BCdat$b3pathologic_uicc_t_stage == "T2") |
                                                   (B3BCdat$b3pathologic_uicc_t_stage == "T3") | 
                                                            (B3BCdat$b3pathologic_uicc_t_stage == "T4a") | 
                                                            (B3BCdat$b3pathologic_uicc_t_stage == "T4b"), 1, NA))

table(B3BCdat$pato_t_cat, B3BCdat$b3pathologic_uicc_t_stage)

# Create new variable "pato_n_cat" from "b3pathologic_uicc_n_stage"

B3BCdat$pato_n_cat <- rep(NA, nrow(B3BCdat))

B3BCdat$pato_n_cat <- ifelse((B3BCdat$b3pathologic_uicc_n_stage == "N0"), 0,
                             + ifelse((B3BCdat$b3pathologic_uicc_n_stage == "N1") |
                                        (B3BCdat$b3pathologic_uicc_n_stage == "N2") |
                                        (B3BCdat$b3pathologic_uicc_n_stage == "N3"), 1, NA))

table(B3BCdat$pato_n_cat, B3BCdat$b3pathologic_uicc_n_stage)


table(is.na(B3BCdat$b3pathologic_uicc_n_stage))
table(is.na(B3BCdat$b3neo_adj_uicc_n_stage))

# Converting to factors 

## BDdat variables 
BDdat$diabetes <-  as.factor(BDdat$diabetes)

BDdat$hypertension <-  as.factor(BDdat$hypertension)

BDdat$history_of_heart_disease <-  as.factor(BDdat$history_of_heart_disease)

BDdat$ra <-  as.factor(BDdat$ra)

BDdat$depression <-  as.factor(BDdat$depression)

BDdat$breast_cancer_family_history_1st_degree <-  as.factor(BDdat$breast_cancer_family_history_1st_degree)

## B3BCdat variables 
B3BCdat$b3tumour_side_of_primary <-  as.factor(B3BCdat$b3tumour_side_of_primary)
levels(B3BCdat$b3tumour_side_of_primary) <- c("left", "right")

B3BCdat$b3surgery_type <-  as.factor(B3BCdat$b3surgery_type)
levels(B3BCdat$b3surgery_type) <- c("Segmentectomy/Quadrantectomy", "Wide local excision")

B3BCdat$b3sys_treatment <-  as.factor(B3BCdat$b3sys_treatment)

B3BCdat$pato_t_cat <-  as.factor(B3BCdat$pato_t_cat)
levels(B3BCdat$pato_t_cat) <- c("Tis-T1", "T2-T4")

B3BCdat$pato_n_cat <-  as.factor(B3BCdat$pato_n_cat)
levels(B3BCdat$pato_n_cat) <- c("N0", "N+")

# Change NAs to 0 for those patients who answered 0 for the variable systemic treatment. 
B3BCdat$b3sys_antiher2 <- ifelse(B3BCdat$b3sys_treatment == 0 & 
                                   is.na(B3BCdat$b3sys_antiher2), 0, B3BCdat$b3sys_antiher2)
B3BCdat$b3sys_antiher2 <-  as.factor(B3BCdat$b3sys_antiher2)
  


B3BCdat$b3chemo_neo_adjuvant <-  as.factor(B3BCdat$b3chemo_neo_adjuvant)
B3BCdat$b3chemo_neoadjuvant_anthracycline <-  as.factor(B3BCdat$b3chemo_neoadjuvant_anthracycline)
B3BCdat$b3chemo_neoadjuvant_nonanthracycline <-  as.factor(B3BCdat$b3chemo_neoadjuvant_nonanthracycline)

B3BCdat$b3chemo_adjuvant <-  as.factor(B3BCdat$b3chemo_adjuvant)
B3BCdat$b3chemo_adjuvant_anthracycline <-  as.factor(B3BCdat$b3chemo_adjuvant_anthracycline)
B3BCdat$b3chemo_adjuvant_nonanthracycline <-  as.factor(B3BCdat$b3chemo_adjuvant_nonanthracycline)

B3BCdat$b3radio_imrt <-  as.factor(B3BCdat$b3radio_imrt)

B3BCdat$b3radio_boost <-  as.factor(B3BCdat$b3radio_boost) 

## create a variable chemo_mix and then create a variable chemo_bin 
B3BCdat <-  transform(B3BCdat, chemo_mix=paste(b3chemo_adjuvant, b3chemo_neo_adjuvant))

B3BCdat$chemo_bin <- rep(NA, nrow(B3BCdat))
B3BCdat$chemo_bin[B3BCdat$chemo_mix=="0 0"] <- "No"
B3BCdat$chemo_bin[B3BCdat$chemo_mix %in% c("1 0", "0 1")] <- "Yes"  

table(B3BCdat$chemo_bin)

## C30_small variables 

###FIRST OPTION###

### We are going to create new binary variables 1-2 and 3-4, and converting to
### factors at the same time, but only for the dyspnea, pain and insomnia questions.
### The QoL questions (29 and 30) stay the same, at least for now. 

table(C30_small$q8)
table(C30_small$q9)
table(C30_small$q11)
table(C30_small$q19)
table(C30_small$q29)
table(C30_small$q30)

C30_small$q8_bin <- rep(NA, nrow(C30_small))
C30_small$q8_bin[C30_small$q8 %in% c(1, 2)] <- "0"
C30_small$q8_bin[C30_small$q8 %in% c(3, 4)] <- "1"
table(C30_small$q8_bin)

C30_small$q9_bin <- rep(NA, nrow(C30_small))
C30_small$q9_bin[C30_small$q9 %in% c(1, 2)] <- "0"
C30_small$q9_bin[C30_small$q9 %in% c(3, 4)] <- "1"
table(C30_small$q9_bin)

C30_small$q11_bin <- rep(NA, nrow(C30_small))
C30_small$q11_bin[C30_small$q11 %in% c(1, 2)] <- "0"
C30_small$q11_bin[C30_small$q11 %in% c(3, 4)] <- "1"
table(C30_small$q11_bin)

C30_small$q19_bin <- rep(NA, nrow(C30_small))
C30_small$q19_bin[C30_small$q19 %in% c(1, 2)] <- "0"
C30_small$q19_bin[C30_small$q19 %in% c(3, 4)] <- "1"
table(C30_small$q19_bin)

###SECOND OPTION###

### We are going to use the sub-scale transformed scores from the  QLQ-C30 (QL, PA, DY, and SL).
### First we use the qlq_c30 function from PROscore to calculate the transformed scores for every
### subscale of the QLQ_30. 

c30scores <- qlq_c30(C30_baseline, iprefix = "q" )

### Then, we merge the scores with the initial data frame so that we have the Subject IDs with the
### new scores in the same data frame. 
C30dat_scores <- merge(C30_baseline, c30scores, by = 0)

### Then, we select only the sub-scales of interest (QL2, PA, DY, and SL). 
C30dat_scores %>%
  dplyr::select(Subject.Id, QL, PA, DY, SL) -> C30dat_scores_small


#######Scoring the BR23 dimensions)################

# Create new data.frame with only the BR23 questions that we need

BR23_baseline %>%
  dplyr::select(Subject.Id, BR23_q47_pain_arm_shoulder, BR23_q48_swollen_arm_hand, 
         BR23_q49_difficult_raise_arm, BR23_q50_pain_breast,
         BR23_q51_swollen_breast, BR23_q52_oversensitive_breast, 
         BR23_q53_skin_problems)-> BR23_7Q

# Add the BR23 questions to C30_baseline

C30_baseline %>%
  left_join(BR23_7Q, by="Subject.Id") -> c30_br23

# Remove the C30 questions that are going to be replaced by the BR23 questions
# (i.e. remove 10, 12, 18, 21, 22, 23, 24)

c30_br23$q10 <- NULL
c30_br23$q12 <- NULL
c30_br23$q18 <- NULL
c30_br23$q21 <- NULL
c30_br23$q22 <- NULL
c30_br23$q23 <- NULL
c30_br23$q24 <- NULL

# Rename the BR23 questions with the names of the previously
# removed C30 questions


c30_br23 <- setNames(c30_br23, c("Subject.Id","Visit","event_date", "submitdate", "q1", "q2", 
                                 "q3", "q4", "q5", "q6", "q7", "q8", "q9", "q11", 
                                 "q13", "q14", "q15", "q16", "q17", "q19", 
                                 "q20", "q25", "q26", "q27", "q28", "q29", "q30", 
                                 "q10", "q12", "q18", "q21", "q22", "q23", "q24"))
str(c30_br23)

# Order the questions

c30_br23_ordered <- c30_br23[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 28, 14, 29, 
                                 15, 16, 17, 18, 19, 30, 20, 21, 31, 32, 33, 34, 22, 
                                 23, 24, 25, 26, 27)]
str(c30_br23_ordered)

# Now we use the function qlq_c30 from PROscorer to calculate the "C30" scores
# (we do this to calculate the BR23 Breast and Arm symptoms dimensions)

c30scores_cheat <- qlq_c30(c30_br23_ordered, iprefix = "q" )

# Merge the scores with the initial data frame so that we have the 
# Subject IDs with the new scores in the same data frame. 
c30scores_cheat_id <- merge(c30scores_cheat, c30_br23_ordered, by = 0)

# Select only the sub-scales of interest. In this case, "EF" is BRBS and 
# "FA" is BRAS. 

c30scores_cheat_id %>%
  dplyr::select(Subject.Id, FA, EF ) -> BR23_small_wrong_labels

# Change the names according to the BR23 dimensions "Breast symptoms" (BRBS)
# and "Arm symptoms" (BRAS)

str(BR23_small_wrong_labels)

BR23_final <- setNames(BR23_small_wrong_labels, c("Subject.Id","BRAS", "BRBS"))

BR23_final$BRBS <- 100 - BR23_final$BRBS 

## Now we can create BDdat_small 

BDdat %>%
  dplyr::select(Subject.Id, age_at_radiotherapy_start_yrs, age_cat, ethnicity_bin, 
                household_bin, education_lmh, 
        bmi, bmi_cat, diabetes, hypertension, history_of_heart_disease, 
         ra, depression, breast_cancer_family_history_1st_degree, 
         smoker_bin, alcohol_bin) -> BDdat_small

BDdat_small$bmi <- round(BDdat_small$bmi,digits=1)

# We can also create B3BC_small 

B3BCdat %>%
  dplyr::select(Subject.Id, b3tumour_side_of_primary, b3surgery_type, 
                histology_di, b3sys_treatment, 
         systemic_cat, b3sys_antiher2, chemo_bin, b3radio_breast_dose_Gy, 
         b3radio_imrt, b3radio_boost, pato_t_cat, 
         pato_n_cat) -> B3BC_small

# And we also create Rad_small2

Rad_small %>%
  dplyr::select(Subject.Id, fractionation_bin, boost_frac, 
                radio_total_dose_sum) -> Rad_small2


# Then, we are going to left join BDdat to MFI_24month. This way, we only keep the  
# subjects that have MFI-20 measurements at 24 months and add the baseline variables. 
MFI_24month %>%
  left_join(BDdat_small, by="Subject.Id") -> df1

## left join SDdat_small to df1. 

df1 %>%
  left_join(SDdat_small, by="Subject.Id") -> df2

## left join GPAQ_small to df2

df2 %>%
  left_join(GPAQ_small, by="Subject.Id") -> df3

## left join c30dat_scores_small to df3

df3 %>%
  left_join(C30dat_scores_small, by="Subject.Id") -> df4

## left join BR23_small to df4

df4 %>%
  left_join(BR23_final, by="Subject.Id") -> df5

## left join B3BC_small to df5

df5 %>%
  left_join(B3BC_small, by="Subject.Id") -> df6

## left join Rad_small2 to df6 

df6 %>%
  left_join(Rad_small2, by="Subject.Id") -> df6


## Create descriptive summary statistics table with table1

table1::label(df6$age_cat) <- "Age"
table1::label(df6$ethnicity_bin) <- "Ethnicity"
table1::label(df6$household_bin) <- "Living status"
table1::label(df6$education_lmh) <- "Educational level"
table1::label(df6$bmi_cat) <- "BMI"
table1::label(df6$diabetes) <- "Diabetes"
table1::label(df6$hypertension) <- "Hypertension"
table1::label(df6$history_of_heart_disease) <- "Heart disease"
table1::label(df6$ra) <- "Rheumatoid arthritis"
table1::label(df6$depression) <- "Depression"
table1::label(df6$breast_cancer_family_history_1st_degree) <- "Breast cancer 1st degree"
table1::label(df6$smoker_bin) <- "Smoking status"
table1::label(df6$alcohol_bin) <- "Current alcohol consumption"
table1::label(df6$GPAQ_q16_sitting_recline_min) <- "Time spent sitting (minutes per day)"
table1::label(df6$QL) <- "Global Health Status / QoL (QLQ-C30)"
table1::label(df6$PA) <- "Pain (QLQ-C30)"
table1::label(df6$SL) <- "Insomnia (QLQ-C30)"
table1::label(df6$DY) <- "Dyspnea (QLQ-C30)"
table1::label(df6$BRAS) <- "Arm symptoms-BR23"
table1::label(df6$BRBS) <- "Breast symptoms-BR23"
table1::label(df6$histology_di) <- "Histologic type"
table1::label(df6$b3tumour_side_of_primary) <- "Side of primary tumor"
table1::label(df6$b3surgery_type) <- "Surgery type"
table1::label(df6$b3sys_treatment) <- "Other systemic treatment"
table1::label(df6$systemic_cat) <- "Systemic hormonal treatment (CA related)"
table1::label(df6$b3sys_antiher2) <- "Anti-HER2 therapie"
table1::label(df6$chemo_bin) <- "Chemotherapie"
table1::label(df6$b3radio_breast_dose_Gy) <- "Radiotherapie dose"
table1::label(df6$radio_total_dose_sum) <- "Radiotherapie total dose (including boosts)"
table1::label(df6$b3radio_boost) <- "Radiotherapie boost"
table1::label(df6$fractionation_bin) <- "Fractionation scheme"
table1::label(df6$boost_frac) <- "Fractionation and boost scheme"
table1::label(df6$b3radio_imrt) <- "IMRT"
table1::label(df6$pato_t_cat) <- "Pathologic T-stage"
table1::label(df6$pato_n_cat) <- "Pathologic N-status"

table1::table1(~age_cat + ethnicity_bin + household_bin 
               + education_lmh + bmi_cat + diabetes + hypertension +
                 history_of_heart_disease + ra + depression + 
                 breast_cancer_family_history_1st_degree + 
                 smoker_bin + alcohol_bin + 
                 GPAQ_q16_sitting_recline_min + QL + PA + SL + DY + 
                 BRAS + BRBS +
                 b3tumour_side_of_primary + b3surgery_type + histology_di + 
                 systemic_cat + b3sys_treatment + b3sys_antiher2 + chemo_bin + 
                 b3radio_breast_dose_Gy + radio_total_dose_sum +
                 b3radio_boost + fractionation_bin + boost_frac +
                 b3radio_imrt + pato_t_cat + pato_n_cat, data = df6)


# combine the individual scores for each of the MFI dimensions 
# and create the dimensions' variables
df6$MFI_general <- MFI_24month$MFI_q01_fit + MFI_24month$MFI_q05_tired + 
  MFI_24month$MFI_q12_rested + MFI_24month$MFI_q16_tire_easily

df6$MFI_physic <- MFI_24month$MFI_q02_do_little + 
  MFI_24month$MFI_q08_take_a_lot_physical + MFI_24month$MFI_q14_bad_condition + 
  MFI_24month$MFI_q20_excellent_condition

df6$MFI_mental <- MFI_24month$MFI_q07_keep_thoughts + 
  MFI_24month$MFI_q11_can_concentrate + MFI_24month$MFI_q13_effort_concentration +
  MFI_24month$MFI_q19_thoughts_wander

df6$MFI_redact <- MFI_24month$MFI_q03_very_active + MFI_24month$MFI_q06_do_lot +
  MFI_24month$MFI_q10_do_very_little + MFI_24month$MFI_q17_get_little_done

df6$MFI_redmot <- MFI_24month$MFI_q04_doing_all_sorts + 
  MFI_24month$MFI_q09_dread_doing_things + MFI_24month$MFI_q15_plans + 
  MFI_24month$MFI_q18_dont_feel_do_anything

df6 %>%
  dplyr::select(age_cat, MFI_general, MFI_physic, MFI_mental,
                MFI_redact, MFI_redmot) -> df_MFI



######## MFI descriptives 

MFI_datx <- data.frame(MFI_data_rev)

MFI_datx$MFI_general <- MFI_data_rev$MFI_q01_fit + MFI_data_rev$MFI_q05_tired +
  MFI_data_rev$MFI_q12_rested + MFI_data_rev$MFI_q16_tire_easily

MFI_datx$MFI_physic <- MFI_data_rev$MFI_q02_do_little + 
  MFI_data_rev$MFI_q08_take_a_lot_physical + MFI_data_rev$MFI_q14_bad_condition +
  MFI_data_rev$MFI_q20_excellent_condition

MFI_datx$MFI_mental <- MFI_data_rev$MFI_q07_keep_thoughts + 
  MFI_data_rev$MFI_q11_can_concentrate + MFI_data_rev$MFI_q13_effort_concentration +
  MFI_data_rev$MFI_q19_thoughts_wander

MFI_datx$MFI_redact <- MFI_data_rev$MFI_q03_very_active + MFI_data_rev$MFI_q06_do_lot +
  MFI_data_rev$MFI_q10_do_very_little + MFI_data_rev$MFI_q17_get_little_done

MFI_datx$MFI_redmot <- MFI_data_rev$MFI_q04_doing_all_sorts + 
  MFI_data_rev$MFI_q09_dread_doing_things + MFI_data_rev$MFI_q15_plans +
  MFI_data_rev$MFI_q18_dont_feel_do_anything

MFI_datx$MFI_total <- MFI_datx$MFI_general + MFI_datx$MFI_physic + 
  MFI_datx$MFI_mental +
  MFI_datx$MFI_redact + MFI_datx$MFI_redmot

# Filter by visit 
MFI_datx %>% 
  filter(MFI_datx$Visit == "Baseline") -> MFI_baseline

MFI_datx %>% 
  filter(MFI_datx$Visit == "Post-RT") -> MFI_postRT

MFI_datx %>% 
  filter(MFI_datx$Visit == "12m follow-up") -> MFI_12m

MFI_datx %>% 
  filter(MFI_datx$Visit == "24m follow-up") -> MFI_24m

MFI_datx %>% 
  filter(MFI_datx$Visit == "3m follow-up") -> MFI_3m

MFI_datx$Visit <- factor(MFI_datx$Visit, levels = 
                           c("Baseline", "Post-RT", "3m follow-up", 
                                                    "12m follow-up", "24m follow-up"))


# Mean fatigue scores and SD by follow-up 
tab1 <- tableby(Visit ~ MFI_general + MFI_physic + MFI_mental + MFI_redmot + MFI_redact, data=MFI_datx)
summary(tab1)

############# Missing values (outcome) with AMELIA ############# 

MFIdat %>% 
  filter(MFIdat$Visit == "24m follow-up") -> M24

MFIdat %>% 
  filter(MFIdat$Visit == "12m follow-up") -> M12

MFIdat %>% 
  filter(MFIdat$Visit == "Baseline") -> Base

MFIdat %>% 
  filter(MFIdat$Visit == "Post-RT") -> Post

Base$Subject.Id <- NULL
Base$Visit <- NULL
Base$event_date <- NULL
Base$submitdate <- NULL
Base$MFI_q21_tired_past_week <- NULL
Base$MFI_q22_tired_year_pre_diag <- NULL

Post$MFI_q22_tired_year_pre_diag <- NULL
Post$Subject.Id <- NULL
Post$Visit <- NULL
Post$event_date <- NULL
Post$submitdate <- NULL
Post$MFI_q21_tired_past_week <- NULL

M12$MFI_q22_tired_year_pre_diag <- NULL
M12$Subject.Id <- NULL
M12$Visit <- NULL
M12$event_date <- NULL
M12$submitdate <- NULL
M12$MFI_q21_tired_past_week <- NULL

M24$MFI_q22_tired_year_pre_diag <- NULL
M24$Subject.Id <- NULL
M24$Visit <- NULL
M24$event_date <- NULL
M24$submitdate <- NULL
M24$MFI_q21_tired_past_week <- NULL

## Removing all rows with only NAs (for MFI-20 questions only) and creating
## new data.frames with rows that at least have one question without NAs

Base_full <- Base[rowSums(is.na(Base)) != ncol(Base), ]
Post_full <- Post[rowSums(is.na(Post)) != ncol(Post), ]
M12_full <- M12[rowSums(is.na(M12)) != ncol(M12), ]
M24_full <- M24[rowSums(is.na(M24)) != ncol(M24), ]

missmap(Base_full, main = "Missingness map MFI-20 Baseline")
missmap(Post_full, main = "Missingness map MFI-20 Post-RT")
missmap(M12_full, main = "Missingness map MFI-20 12 months")
missmap(M24_full, main = "Missingness map MFI-20 24 months")


## Want to know # of patients at 24 months have complete data for 
## the general fatigue dimension of the MFI-20. We create M24_full2.
## Then we obtain the total number of rows with NAs in the 
## general fatigue column, and also for the other dimensions.  

#General fatigue
M24_full2$MFI_general <- M24_full$MFI_q01_fit + 
  M24_full$MFI_q05_tired + M24_full$MFI_q12_rested + 
  M24_full$MFI_q16_tire_easily

sum(is.na(M24_full2$MFI_general)) 
      # The total number of patients at 24 months with complete General fatigue
      # data is 1098 - 51 = 1047

#Physical fatigue
M24_full2$MFI_physic <- M24_full$MFI_q02_do_little + 
  M24_full$MFI_q08_take_a_lot_physical + 
  M24_full$MFI_q14_bad_condition + 
  M24_full$MFI_q20_excellent_condition


sum(is.na(M24_full2$MFI_physic)) 
      # The total number of patients at 24 months with complete General fatigue
      # data is 1098 - 42 = 1056

#Mental fatigue
M24_full2$MFI_mental <- M24_full$MFI_q07_keep_thoughts + 
  M24_full$MFI_q11_can_concentrate + 
  M24_full$MFI_q13_effort_concentration + 
  M24_full$MFI_q19_thoughts_wander


sum(is.na(M24_full2$MFI_mental)) 
      # The total number of patients at 24 months with complete General fatigue
      # data is 1098 - 52 = 1046

#Reduced activity
M24_full2$MFI_redact <- M24_full$MFI_q03_very_active + 
  M24_full$MFI_q06_do_lot + 
  M24_full$MFI_q10_do_very_little + 
  M24_full$MFI_q17_get_little_done


sum(is.na(M24_full2$MFI_redact)) 
      # The total number of patients at 24 months with complete General fatigue
      # data is 1098 - 41 = 1057

#Reduced motivation
M24_full2$MFI_redmot <- M24_full$MFI_q04_doing_all_sorts + 
  M24_full$MFI_q09_dread_doing_things + 
  M24_full$MFI_q15_plans + 
  M24_full$MFI_q18_dont_feel_do_anything


sum(is.na(M24_full2$MFI_redmot)) 
      # The total number of patients at 24 months with complete General fatigue
      # data is 1098 - 48 = 1050



MFI_baseline_2 <- MFI_baseline

MFI_baseline_2$MFI_q01_fit[is.na(MFI_baseline_2$MFI_q01_fit)] <- 0
MFI_baseline_2$MFI_q02_do_little[is.na(MFI_baseline_2$MFI_q02_do_little)] <- 0
MFI_baseline_2$MFI_q03_very_active[is.na(MFI_baseline_2$MFI_q03_very_active)] <- 0
MFI_baseline_2$MFI_q04_doing_all_sorts[is.na(MFI_baseline_2$MFI_q04_doing_all_sorts)] <- 0
MFI_baseline_2$MFI_q05_tired[is.na(MFI_baseline_2$MFI_q05_tired)] <- 0
MFI_baseline_2$MFI_q01_fit[is.na(MFI_baseline_2$MFI_q01_fit)] <- 0
MFI_baseline_2$MFI_q06_do_lot[is.na(MFI_baseline_2$MFI_q06_do_lot)] <- 0
MFI_baseline_2$MFI_q07_keep_thoughts[is.na(MFI_baseline_2$MFI_q07_keep_thoughts)] <- 0
MFI_baseline_2$MFI_q08_take_a_lot_physical[is.na(MFI_baseline_2$MFI_q08_take_a_lot_physical)] <- 0
MFI_baseline_2$MFI_q09_dread_doing_things[is.na(MFI_baseline_2$MFI_q09_dread_doing_things)] <- 0
MFI_baseline_2$MFI_q10_do_very_little[is.na(MFI_baseline_2$MFI_q10_do_very_little)] <- 0
MFI_baseline_2$MFI_q11_can_concentrate[is.na(MFI_baseline_2$MFI_q11_can_concentrate)] <- 0
MFI_baseline_2$MFI_q12_rested[is.na(MFI_baseline_2$MFI_q12_rested)] <- 0
MFI_baseline_2$MFI_q13_effort_concentration[is.na(MFI_baseline_2$MFI_q13_effort_concentration)] <- 0
MFI_baseline_2$MFI_q14_bad_condition[is.na(MFI_baseline_2$MFI_q14_bad_condition)] <- 0
MFI_baseline_2$MFI_q15_plans[is.na(MFI_baseline_2$MFI_q15_plans)] <- 0
MFI_baseline_2$MFI_q16_tire_easily[is.na(MFI_baseline_2$MFI_q16_tire_easily)] <- 0
MFI_baseline_2$MFI_q17_get_little_done[is.na(MFI_baseline_2$MFI_q17_get_little_done)] <- 0
MFI_baseline_2$MFI_q18_dont_feel_do_anything[is.na(MFI_baseline_2$MFI_q18_dont_feel_do_anything)] <- 0
MFI_baseline_2$MFI_q19_thoughts_wander[is.na(MFI_baseline_2$MFI_q19_thoughts_wander)] <- 0
MFI_baseline_2$MFI_q20_excellent_condition[is.na(MFI_baseline_2$MFI_q20_excellent_condition)] <- 0

## Sum all questions and obtain a column with the total. Those rows with 0 
## are the rows where no questions were answered. Those are the ones
## we want to delete. 
MFI_baseline_2$sumall <- MFI_baseline_2$MFI_q01_fit + MFI_baseline_2$MFI_q02_do_little+
                               MFI_baseline_2$MFI_q03_very_active + MFI_baseline_2$MFI_q04_doing_all_sorts+
                               MFI_baseline_2$MFI_q05_tired+MFI_baseline_2$MFI_q06_do_lot+
                               MFI_baseline_2$MFI_q07_keep_thoughts+MFI_baseline_2$MFI_q08_take_a_lot_physical+
                               MFI_baseline_2$MFI_q09_dread_doing_things+MFI_baseline_2$MFI_q10_do_very_little+
                               MFI_baseline_2$MFI_q11_can_concentrate+MFI_baseline_2$MFI_q12_rested+
                               MFI_baseline_2$MFI_q13_effort_concentration+MFI_baseline_2$MFI_q14_bad_condition+
                               MFI_baseline_2$MFI_q15_plans+MFI_baseline_2$MFI_q16_tire_easily+
                               MFI_baseline_2$MFI_q17_get_little_done+MFI_baseline_2$MFI_q18_dont_feel_do_anything+
                               MFI_baseline_2$MFI_q19_thoughts_wander+MFI_baseline_2$MFI_q20_excellent_condition


MFI_baseline_3 <- MFI_baseline_2[MFI_baseline_2$sumall!= 0, ]

## Now we want the same but for MFI_24m
MFI_24m_2 <- MFI_24m

MFI_24m_2$MFI_q01_fit[is.na(MFI_24m_2$MFI_q01_fit)] <- 0
MFI_24m_2$MFI_q02_do_little[is.na(MFI_24m_2$MFI_q02_do_little)] <- 0
MFI_24m_2$MFI_q03_very_active[is.na(MFI_24m_2$MFI_q03_very_active)] <- 0
MFI_24m_2$MFI_q04_doing_all_sorts[is.na(MFI_24m_2$MFI_q04_doing_all_sorts)] <- 0
MFI_24m_2$MFI_q05_tired[is.na(MFI_24m_2$MFI_q05_tired)] <- 0
MFI_24m_2$MFI_q01_fit[is.na(MFI_24m_2$MFI_q01_fit)] <- 0
MFI_24m_2$MFI_q06_do_lot[is.na(MFI_24m_2$MFI_q06_do_lot)] <- 0
MFI_24m_2$MFI_q07_keep_thoughts[is.na(MFI_24m_2$MFI_q07_keep_thoughts)] <- 0
MFI_24m_2$MFI_q08_take_a_lot_physical[is.na(MFI_24m_2$MFI_q08_take_a_lot_physical)] <- 0
MFI_24m_2$MFI_q09_dread_doing_things[is.na(MFI_24m_2$MFI_q09_dread_doing_things)] <- 0
MFI_24m_2$MFI_q10_do_very_little[is.na(MFI_24m_2$MFI_q10_do_very_little)] <- 0
MFI_24m_2$MFI_q11_can_concentrate[is.na(MFI_24m_2$MFI_q11_can_concentrate)] <- 0
MFI_24m_2$MFI_q12_rested[is.na(MFI_24m_2$MFI_q12_rested)] <- 0
MFI_24m_2$MFI_q13_effort_concentration[is.na(MFI_24m_2$MFI_q13_effort_concentration)] <- 0
MFI_24m_2$MFI_q14_bad_condition[is.na(MFI_24m_2$MFI_q14_bad_condition)] <- 0
MFI_24m_2$MFI_q15_plans[is.na(MFI_24m_2$MFI_q15_plans)] <- 0
MFI_24m_2$MFI_q16_tire_easily[is.na(MFI_24m_2$MFI_q16_tire_easily)] <- 0
MFI_24m_2$MFI_q17_get_little_done[is.na(MFI_24m_2$MFI_q17_get_little_done)] <- 0
MFI_24m_2$MFI_q18_dont_feel_do_anything[is.na(MFI_24m_2$MFI_q18_dont_feel_do_anything)] <- 0
MFI_24m_2$MFI_q19_thoughts_wander[is.na(MFI_24m_2$MFI_q19_thoughts_wander)] <- 0
MFI_24m_2$MFI_q20_excellent_condition[is.na(MFI_24m_2$MFI_q20_excellent_condition)] <- 0

## Sum all questions and obtain a column with the total. Those rows with 0 
## are the rows where no questions were answered. Those are the ones
## we want to delete. 
MFI_24m_2$sumall <- MFI_24m_2$MFI_q01_fit + MFI_24m_2$MFI_q02_do_little+
  MFI_24m_2$MFI_q03_very_active + MFI_24m_2$MFI_q04_doing_all_sorts+
  MFI_24m_2$MFI_q05_tired+MFI_24m_2$MFI_q06_do_lot+
  MFI_24m_2$MFI_q07_keep_thoughts+MFI_24m_2$MFI_q08_take_a_lot_physical+
  MFI_24m_2$MFI_q09_dread_doing_things+MFI_24m_2$MFI_q10_do_very_little+
  MFI_24m_2$MFI_q11_can_concentrate+MFI_24m_2$MFI_q12_rested+
  MFI_24m_2$MFI_q13_effort_concentration+MFI_24m_2$MFI_q14_bad_condition+
  MFI_24m_2$MFI_q15_plans+MFI_24m_2$MFI_q16_tire_easily+
  MFI_24m_2$MFI_q17_get_little_done+MFI_24m_2$MFI_q18_dont_feel_do_anything+
  MFI_24m_2$MFI_q19_thoughts_wander+MFI_24m_2$MFI_q20_excellent_condition

MFI_24m_3 <- MFI_24m_2[MFI_24m_2$sumall!= 0, ]


## Now let's create a new column with all the same values. This step is required
## before we can have in one data.frame the subjects divided by those 
## present at baseline and those lost to follow-up
MFI_24m_3$follow_up <- rep(NA, nrow(MFI_24m_3))

MFI_24m_3$follow_up <- "yes"

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

################## LOGISTIC REGRESSION ANALYSIS ###########################

####### Preparing the data #######

## 1. Create data.frame with only the IDs and the 5 dimension's 
## scores. 

MFI_baseline_3 %>%
  dplyr::select(Subject.Id, MFI_general, MFI_physic, MFI_mental, 
         MFI_redact, MFI_redmot) -> MFI_baseline_dimensions

## 2. Change the names of these dimensions (to differentiate them 
## from the 24mth measurements) and then left joint them to df7

str(MFI_baseline_dimensions)
MFI_baseline_dimensions <- setNames(MFI_baseline_dimensions, 
                                    c("Subject.Id", "MFI_general_b", "MFI_physic_b", 
                                      "MFI_mental_b", "MFI_redact_b", "MFI_redmot_b"))

## 3. Create new binary variables for every MFI-20 dimension (Baseline)

##General Fatigue
MFI_baseline_dimensions$MFI_general_bin_base <- rep(NA, nrow(MFI_baseline_dimensions))
MFI_baseline_dimensions$MFI_general_bin_base[MFI_baseline_dimensions$MFI_general_b<=12] <- "0" 
MFI_baseline_dimensions$MFI_general_bin_base[MFI_baseline_dimensions$MFI_general_b>12] <- "1"

MFI_baseline_dimensions$MFI_general_bin_base <- factor(MFI_baseline_dimensions$MFI_general_bin_base,
                                                       levels = c("0", "1"))

table(MFI_baseline_dimensions$MFI_general_bin_base)

##Physical Fatigue
MFI_baseline_dimensions$MFI_phy_bin_base <- rep(NA, nrow(MFI_baseline_dimensions))
MFI_baseline_dimensions$MFI_phy_bin_base[MFI_baseline_dimensions$MFI_physic_b<=12] <- "0" 
MFI_baseline_dimensions$MFI_phy_bin_base[MFI_baseline_dimensions$MFI_physic_b>12] <- "1"

MFI_baseline_dimensions$MFI_phy_bin_base <- factor(MFI_baseline_dimensions$MFI_phy_bin_base,
                                                       levels = c("0", "1"))

table(MFI_baseline_dimensions$MFI_phy_bin_base)

##Mental Fatigue
MFI_baseline_dimensions$MFI_mental_bin_base <- rep(NA, nrow(MFI_baseline_dimensions))
MFI_baseline_dimensions$MFI_mental_bin_base[MFI_baseline_dimensions$MFI_mental_b<=12] <- "0" 
MFI_baseline_dimensions$MFI_mental_bin_base[MFI_baseline_dimensions$MFI_mental_b>12] <- "1"

MFI_baseline_dimensions$MFI_mental_bin_base <- factor(MFI_baseline_dimensions$MFI_mental_bin_base,
                                                   levels = c("0", "1"))

table(MFI_baseline_dimensions$MFI_mental_bin_base)

##Reduced activity 
MFI_baseline_dimensions$MFI_redact_bin_base <- rep(NA, nrow(MFI_baseline_dimensions))
MFI_baseline_dimensions$MFI_redact_bin_base[MFI_baseline_dimensions$MFI_redact_b<=12] <- "0" 
MFI_baseline_dimensions$MFI_redact_bin_base[MFI_baseline_dimensions$MFI_redact_b>12] <- "1"

MFI_baseline_dimensions$MFI_redact_bin_base <- factor(MFI_baseline_dimensions$MFI_redact_bin_base,
                                                      levels = c("0", "1"))

table(MFI_baseline_dimensions$MFI_redact_bin_base)

##Reduced motivation 
MFI_baseline_dimensions$MFI_redmot_bin_base <- rep(NA, nrow(MFI_baseline_dimensions))
MFI_baseline_dimensions$MFI_redmot_bin_base[MFI_baseline_dimensions$MFI_redmot_b<=12] <- "0" 
MFI_baseline_dimensions$MFI_redmot_bin_base[MFI_baseline_dimensions$MFI_redmot_b>12] <- "1"

MFI_baseline_dimensions$MFI_redmot_bin_base <- factor(MFI_baseline_dimensions$MFI_redmot_bin_base,
                                                      levels = c("0", "1"))

table(MFI_baseline_dimensions$MFI_redmot_bin_base)


## 4. Create data.frame with only binary baseline dimensions

str(MFI_baseline_dimensions)

MFI_baseline_dimensions %>%
  dplyr::select(Subject.Id, MFI_general_bin_base, MFI_phy_bin_base, 
         MFI_mental_bin_base, MFI_redact_bin_base, MFI_redmot_bin_base, MFI_general_b) -> MFI_baseline_bin

## 5. Add MFI binary baseline variables to df7

df6 %>%
  left_join(MFI_baseline_bin, by="Subject.Id") -> df7

## 6. Create new binary variables for every MFI-20 dimension (24 months)
## Moderate-severe fatigue (>12)

##General Fatigue
df7$MFI_general_bin <- rep(NA, nrow(df7))
df7$MFI_general_bin[df7$MFI_general<=12] <- "0" 
df7$MFI_general_bin[df7$MFI_general>12] <- "1"

df7$MFI_general_bin <- factor(df7$MFI_general_bin, levels = c("0", "1"))

table(df7$MFI_general_bin)

##Physical Fatigue
df7$MFI_phy_bin <- rep(NA, nrow(df7))
df7$MFI_phy_bin[df7$MFI_physic<=12] <- "0" 
df7$MFI_phy_bin[df7$MFI_physic>12] <- "1"

df7$MFI_phy_bin <- factor(df7$MFI_phy_bin, levels = c("0", "1"))

table(df7$MFI_phy_bin)

##Mental Fatigue
df7$MFI_mental_bin <- rep(NA, nrow(df7))
df7$MFI_mental_bin[df7$MFI_mental<=12] <- "0" 
df7$MFI_mental_bin[df7$MFI_mental>12] <- "1"

df7$MFI_mental_bin <- factor(df7$MFI_mental_bin, levels = c("0", "1"))

table(df7$MFI_mental_bin)

##Reduced activity 
df7$MFI_redact_bin <- rep(NA, nrow(df7))
df7$MFI_redact_bin[df7$MFI_redact<=12] <- "0" 
df7$MFI_redact_bin[df7$MFI_redact>12] <- "1"

df7$MFI_redact_bin <- factor(df7$MFI_redact_bin, levels = c("0", "1"))

table(df7$MFI_redact_bin)

##Reduced motivation 
df7$MFI_redmot_bin <- rep(NA, nrow(df7))
df7$MFI_redmot_bin[df7$MFI_redmot<=12] <- "0" 
df7$MFI_redmot_bin[df7$MFI_redmot>12] <- "1"

df7$MFI_redmot_bin <- factor(df7$MFI_redmot_bin, levels = c("0", "1"))

table(df7$MFI_redmot_bin)

## 7. Create new binary variables for every MFI-20 dimension (Post RT)
## Moderate-severe fatigue (>12)

##General Fatigue
MFI_postRT$MFI_general_bin <- rep(NA, nrow(MFI_postRT))
MFI_postRT$MFI_general_bin[MFI_postRT$MFI_general<=12] <- "0" 
MFI_postRT$MFI_general_bin[MFI_postRT$MFI_general>12] <- "1"

MFI_postRT$MFI_general_bin <- factor(MFI_postRT$MFI_general_bin, levels = c("0", "1"))

table(MFI_postRT$MFI_general_bin)

### 2. Create new data.frames for every dimension by filtering df7

## First, create new continuous variable with 10-unit intervals
df7$QL_10 <- df7$QL / 10
df7$PA_10 <- df7$PA / 10
df7$SL_10 <- df7$SL / 10
df7$DY_10 <- df7$DY / 10
df7$BRAS_10 <- df7$BRAS / 10
df7$BRBS_10 <- df7$BRBS / 10

## And same for age
df7$age_5 <- df7$age_at_radiotherapy_start_yrs / 5


### Just converting "Mainsite" to factor 

df7$Mainsite <- as.factor(df7$Mainsite)


##df7_general
df7 %>%
  dplyr::select(Subject.Id, age_cat, age_at_radiotherapy_start_yrs, age_5, ethnicity_bin, household_bin, education_lmh, 
         bmi, bmi_cat, diabetes, hypertension, history_of_heart_disease, ra, 
         depression, breast_cancer_family_history_1st_degree, smoker_bin, 
         alcohol_bin, GPAQ_q16_sitting_recline_min, QL, PA, SL, DY, BRAS,
         BRBS, histology_di, b3tumour_side_of_primary, 
         b3surgery_type, b3sys_treatment, systemic_cat, b3sys_antiher2, chemo_bin,
         b3radio_breast_dose_Gy, radio_total_dose_sum, fractionation_bin, 
         boost_frac, b3radio_boost, b3radio_imrt, pato_t_cat, 
         pato_n_cat, MFI_general_bin_base, MFI_general_bin, MFI_general, MFI_general_b) -> df7_general

df7_general_complete <- na.omit(df7_general)

## Complete cases when removing some
## selected variables 

df7 %>%
  dplyr::select(Subject.Id, age_at_radiotherapy_start_yrs, age_5, ethnicity_bin, Mainsite, 
         bmi, diabetes, hypertension, history_of_heart_disease, ra, 
         depression, breast_cancer_family_history_1st_degree, smoker_bin, 
         alcohol_bin, QL_10, PA_10, SL_10, DY_10,  BRAS_10, BRBS_10, histology_di, 
         b3tumour_side_of_primary, b3surgery_type, b3sys_treatment, 
         b3sys_antiher2, chemo_bin, radio_total_dose_sum, fractionation_bin, b3radio_boost, 
         boost_frac, b3radio_imrt, pato_t_cat, MFI_general_bin_base, 
         MFI_general_bin) -> df7_general_sin_edu_recline_hous_n

df7_general_complete4 <- na.omit(df7_general_sin_edu_recline_hous_n)


##df7_physical

df7 %>%
  dplyr::select(Subject.Id, age_at_radiotherapy_start_yrs, age_5, ethnicity_bin, Mainsite, 
                bmi, diabetes, hypertension, history_of_heart_disease, ra, 
                depression, breast_cancer_family_history_1st_degree, smoker_bin, 
                alcohol_bin, QL_10, PA_10, SL_10, DY_10,  BRAS_10, BRBS_10, histology_di, 
                b3tumour_side_of_primary, b3surgery_type, b3sys_treatment, 
                b3sys_antiher2, chemo_bin, radio_total_dose_sum, fractionation_bin, b3radio_boost, 
                boost_frac, b3radio_imrt, pato_t_cat, MFI_phy_bin_base,
               MFI_phy_bin) -> df7_phy_sin_edu_recline_hous_n

df7_phy_complete4 <- na.omit(df7_phy_sin_edu_recline_hous_n)



##df7_mental 

df7 %>%
  dplyr::select(Subject.Id, age_at_radiotherapy_start_yrs, age_5, ethnicity_bin, Mainsite, 
                bmi, diabetes, hypertension, history_of_heart_disease, ra, 
                depression, breast_cancer_family_history_1st_degree, smoker_bin, 
                alcohol_bin, QL_10, PA_10, SL_10, DY_10,  BRAS_10, BRBS_10, histology_di, 
                b3tumour_side_of_primary, b3surgery_type, b3sys_treatment, 
                b3sys_antiher2, chemo_bin, radio_total_dose_sum, fractionation_bin, b3radio_boost, 
                boost_frac, b3radio_imrt, pato_t_cat, MFI_mental_bin_base,
         MFI_mental_bin) -> df7_mental_sin_edu_recline_hous_n

df7_mental_complete4 <- na.omit(df7_mental_sin_edu_recline_hous_n)

##df7_redmot

df7 %>%
  dplyr::select(Subject.Id, age_at_radiotherapy_start_yrs, age_5, ethnicity_bin, Mainsite, 
                bmi, diabetes, hypertension, history_of_heart_disease, ra, 
                depression, breast_cancer_family_history_1st_degree, smoker_bin, 
                alcohol_bin, QL_10, PA_10, SL_10, DY_10,  BRAS_10, BRBS_10, histology_di, 
                b3tumour_side_of_primary, b3surgery_type, b3sys_treatment, 
                b3sys_antiher2, chemo_bin, radio_total_dose_sum, fractionation_bin, b3radio_boost, 
                boost_frac, b3radio_imrt, pato_t_cat, MFI_redmot_bin_base,
         MFI_redmot_bin) -> df7_redmot_sin_edu_recline_hous_n

df7_redmot_complete4 <- na.omit(df7_redmot_sin_edu_recline_hous_n)

##df7_redact

df7 %>%
  dplyr::select(Subject.Id, age_at_radiotherapy_start_yrs, age_5, ethnicity_bin, Mainsite, 
                bmi, diabetes, hypertension, history_of_heart_disease, ra, 
                depression, breast_cancer_family_history_1st_degree, smoker_bin, 
                alcohol_bin, QL_10, PA_10, SL_10, DY_10,  BRAS_10, BRBS_10, histology_di, 
                b3tumour_side_of_primary, b3surgery_type, b3sys_treatment, 
                b3sys_antiher2, chemo_bin, radio_total_dose_sum, fractionation_bin, b3radio_boost, 
                boost_frac, b3radio_imrt, pato_t_cat, MFI_redact_bin_base,
         MFI_redact_bin) -> df7_redact_sin_edu_recline_hous_n

df7_redact_complete4 <- na.omit(df7_redact_sin_edu_recline_hous_n)



################  Frequency tables by fatigue status at 24months   ############
################          Bi-variate associations                  ############


# Function to compute the p-values (t-test for continuous and chi-square for categorical variables)
  # Taken from: https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}


# Change level names for MFI binary variables (With and Without fatigue)

levels(df7_general$MFI_general_bin) <- c("Without general fatigue", "With general fatigue")
levels(df7$MFI_general_bin) <- c("Without general fatigue", "With general fatigue")
levels(df7$MFI_phy_bin) <- c("Without physical fatigue", "With physical fatigue")
levels(df7$MFI_mental_bin) <- c("Without mental fatigue", "With mental fatigue")
levels(df7$MFI_redact_bin) <- c("Without reduced activity", "With reduced activity")
levels(df7$MFI_redmot_bin) <- c("Without reduced motivation", "With reduced motivation")


#Table 1.1.0 General fatigue Baseline

table1::label(MFI_baseline_bin$MFI_mental_bin_base) <- "Mental"

table1::table1(~MFI_mental_bin_base | MFI_general_bin_base, data = MFI_baseline_bin, 
               overall=F, extra.col=list("P-value"= pvalue))


## Creating new categorical variable with the longitudinal fatigue groups for 
## general fatigue 

df7_general$fatigue_longi <- rep(NA, nrow(df7_general))
df7_general$fatigue_longi[df7_general$MFI_general_b<=12 &
                            df7_general$MFI_general<=12] <- "no fatigue"
df7_general$fatigue_longi[df7_general$MFI_general_b<=12 &
                            df7_general$MFI_general>12] <- "worsening of fatigue" 
df7_general$fatigue_longi[df7_general$MFI_general_b>12 &
                            df7_general$MFI_general<=12] <- "improvement of fatigue"
df7_general$fatigue_longi[df7_general$MFI_general_b>12 &
                            df7_general$MFI_general>12] <- "persistent high fatigue"


df7_general$fatigue_longi <- factor(df7_general$fatigue_longi, levels = c("no fatigue", 
                                                                          "worsening of fatigue", 
                                                                          "improvement of fatigue", 
                                                                          "persistent high fatigue"))
table(df7_general$fatigue_longi)



df7 %>%
  dplyr::select(Subject.Id, Mainsite, age_cat, age_at_radiotherapy_start_yrs, 
                age_5, ethnicity_bin, household_bin, education_lmh, 
                bmi, bmi_cat, diabetes, hypertension, history_of_heart_disease, ra, 
                depression, breast_cancer_family_history_1st_degree, smoker_bin, 
                alcohol_bin, GPAQ_q16_sitting_recline_min, QL, PA, SL, DY, BRAS,
                BRBS, histology_di, b3tumour_side_of_primary, 
                b3surgery_type, b3sys_treatment, systemic_cat, 
                b3sys_antiher2, chemo_bin,
                b3radio_breast_dose_Gy, radio_total_dose_sum, fractionation_bin, 
                boost_frac, b3radio_boost, b3radio_imrt, pato_t_cat, 
                pato_n_cat, MFI_general_bin_base, MFI_general_bin, 
                MFI_general, MFI_general_b) -> df7_general2



#Table 1.1.1 General fatigue 24 months 

table1::label(df7_general2$Mainsite) <- "Site"
table1::label(df7_general2$MFI_general_bin_base) <- "General fatigue at baseline"
table1::label(df7_general2$age_at_radiotherapy_start_yrs) <- "Age (years)"
table1::label(df7_general2$age_cat) <- "Age categories (years)"
table1::label(df7_general2$household_bin) <- "Living status"
table1::label(df7_general2$education_lmh) <- "Educational level"
table1::label(df7_general2$bmi) <- "BMI (kg/m2)"
table1::label(df7_general2$bmi_cat) <- "BMI categories (kg/m2)"
table1::label(df7_general2$smoker_bin) <- "Smoking status"
table1::label(df7_general2$depression) <- "Depression"
table1::label(df7_general2$GPAQ_q16_sitting_recline_min) <- "Time spent sitting (minutes per day)"
table1::label(df7_general2$QL) <- "Global Health Status / QoL (QLQ-C30)"
table1::label(df7_general2$PA) <- "Pain (QLQ-C30)"
table1::label(df7_general2$SL) <- "Insomnia (QLQ-C30)"
table1::label(df7_general2$b3radio_boost)<- "Radiotherapie boost"
table1::label(df7_general2$fractionation_bin) <- "Fractionation scheme"
table1::label(df7_general2$b3radio_imrt) <- "IMRT"
table1::label(df7_general2$histology_di) <- "Histologic type"
table1::label(df7_general2$b3surgery_type) <- "Surgery type"
table1::label(df7_general2$b3sys_treatment) <- "Hormonal cancer treatment"
table1::label(df7_general2$chemo_bin) <- "Chemotherapy"
table1::label(df7_general2$pato_t_cat) <- "Tumor-stage"
table1::label(df7_general2$pato_n_cat) <- "Pathologic N-status"

table1::table1(~Mainsite + MFI_general_bin_base + age_at_radiotherapy_start_yrs + age_cat + household_bin +
                 education_lmh + bmi + bmi_cat + 
                 smoker_bin + depression +  
                 GPAQ_q16_sitting_recline_min + QL + PA + SL + b3radio_boost+
                 fractionation_bin + b3radio_imrt + histology_di +
                 b3surgery_type + b3sys_treatment + chemo_bin + pato_t_cat 
                  | MFI_general_bin , overall=T, data = df7_general2)



#Table 1.1.2 General fatigue 24 months (Overall without testing)

table1::label(df7_general$MFI_general_bin_base) <- "General fatigue at baseline"
table1::label(df7_general$age_cat) <- "Age"
table1::label(df7_general$ethnicity_bin) <- "Ethnicity"
table1::label(df7_general$household_bin) <- "Living status"
table1::label(df7_general$education_lmh) <- "Educational level"
table1::label(df7_general$bmi_cat) <- "BMI"
table1::label(df7_general$diabetes) <- "Diabetes"
table1::label(df7_general$hypertension) <- "Hypertension"
table1::label(df7_general$history_of_heart_disease) <- "Heart disease"
table1::label(df7_general$ra) <- "Rheumatoid arthritis"
table1::label(df7_general$depression) <- "Depression"
table1::label(df7_general$breast_cancer_family_history_1st_degree) <- "Breast cancer 1st degree"
table1::label(df7_general$smoker_bin) <- "Smoking status"
table1::label(df7_general$alcohol_bin) <- "Current alcohol consumption"
table1::label(df7_general$GPAQ_q16_sitting_recline_min) <- "Time spent sitting (minutes per day)"
table1::label(df7_general$QL) <- "Global Health Status / QoL (QLQ-C30)"
table1::label(df7_general$PA) <- "Pain (QLQ-C30)"
table1::label(df7_general$SL) <- "Insomnia (QLQ-C30)"
table1::label(df7_general$DY) <- "Dyspnea (QLQ-C30)"
table1::label(df7_general$BRAS) <- "Arm symptoms-BR23"
table1::label(df7_general$BRBS) <- "Breast symptoms-BR23"
table1::label(df7_general$histology_di) <- "Histologic type"
table1::label(df7_general$b3tumour_side_of_primary) <- "Side of primary tumor"
table1::label(df7_general$b3surgery_type) <- "Surgery type"
table1::label(df7_general$b3sys_treatment) <- "Other systemic treatment"
table1::label(df7_general$b3sys_antiher2) <- "Anti-HER2 therapie"
table1::label(df7_general$chemo_bin) <- "Chemotherapie"
table1::label(df7_general$b3radio_breast_dose_Gy) <- "Radiotherapie dose"
table1::label(df7_general$radio_total_dose_sum) <- "Radiotherapie total dose (including boosts)"
table1::label(df7_general$b3radio_boost) <- "Radiotherapie boost"
table1::label(df7_general$fractionation_bin) <- "Fractionation scheme"
table1::label(df7_general$boost_frac) <- "Fractionation and boost scheme"
table1::label(df7_general$b3radio_imrt) <- "IMRT"
table1::label(df7_general$pato_t_cat) <- "Pathologic T-stage"
table1::label(df7_general$pato_n_cat) <- "Pathologic N-status"

table1::table1(~MFI_general_bin_base + age_cat + ethnicity_bin + household_bin +
                 education_lmh + bmi_cat + diabetes + hypertension +
                 history_of_heart_disease + ra + depression + 
                 breast_cancer_family_history_1st_degree + 
                 smoker_bin + alcohol_bin + 
                 GPAQ_q16_sitting_recline_min + QL + PA + SL + DY + 
                 BRAS + BRBS + histology_di +
                 b3tumour_side_of_primary + b3surgery_type + systemic_cat + 
                 b3sys_treatment + b3sys_antiher2 + chemo_bin + 
                 b3radio_breast_dose_Gy + radio_total_dose_sum + b3radio_boost +
                 fractionation_bin + boost_frac + b3radio_imrt + pato_t_cat + 
                 pato_n_cat | MFI_general_bin, data = df7_general)


#Table 1.1.3 General fatigue by longi 

## Only keep those rows with enough data at baseline and 24 months for 
## general fatigue 

df7_general <- df7_general[!is.na(df7_general$fatigue_longi),]


rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x))
  what <- switch(name,
                 age_at_radiotherapy_start_yrs  = c(.="Mean (SD)"),
                 bmi= c(.="Mean (SD)"),
                 GPAQ_q16_sitting_recline_min = c(.="Mean (SD)", .="Median [Min, Max]"), 
                 QL = "Mean (SD)",
                 PA  = "Mean (SD)", 
                 SL  = "Mean (SD)")
  parse.abbrev.render.code(c("", what))(x)
}


table1::label(df7_general$age_at_radiotherapy_start_yrs) <- "Age (years)"
table1::label(df7_general$age_cat) <- "Age categories (years)"
table1::label(df7_general$household_bin) <- "Living status"
table1::label(df7_general$education_lmh) <- "Educational level"
table1::label(df7_general$bmi) <- "BMI (kg/m2)"
table1::label(df7_general$bmi_cat) <- "BMI categories (kg/m2)"
table1::label(df7_general$smoker_bin) <- "Smoking status"
table1::label(df7_general$GPAQ_q16_sitting_recline_min) <- "Time spent sitting (minutes per day)"
table1::label(df7_general$QL) <- "Global Health Status / QoL (QLQ-C30)"
table1::label(df7_general$PA) <- "Pain (QLQ-C30)"
table1::label(df7_general$SL) <- "Insomnia (QLQ-C30)"
table1::label(df7_general$chemo_bin) <- "Chemotherapy"
table1::label(df7_general$boost_frac) <- "Fractionation and boost scheme"
table1::label(df7_general$pato_t_cat) <- "Pathologic T-stage"
table1::label(df7_general$pato_n_cat) <- "Pathologic N-status"

table1::table1(~age_at_radiotherapy_start_yrs + age_cat + household_bin +
                 education_lmh + bmi + bmi_cat + 
                 smoker_bin +  
                 GPAQ_q16_sitting_recline_min + QL + PA + SL + 
                 chemo_bin +  boost_frac + pato_t_cat + 
                 pato_n_cat | fatigue_longi, data = df7_general, render=rndr)




#Table 1.2 Physical fatigue 24 months 

table1::label(df7$MFI_phy_bin_base) <- "Physical fatigue at baseline"
table1::label(df7$age_cat) <- "Age"
table1::label(df7) <- "Ethnicity"
table1::label(df7$household_bin) <- "Living status"
table1::label(df7$education_lmh) <- "Educational level"
table1::label(df7$bmi_cat) <- "BMI"
table1::label(df7$diabetes) <- "Diabetes"
table1::label(df7$hypertension) <- "Hypertension"
table1::label(df7$history_of_heart_disease) <- "Heart disease"
table1::label(df7$ra) <- "Rheumatoid arthritis"
table1::label(df7$depression) <- "Depression"
table1::label(df7$breast_cancer_family_history_1st_degree) <- "Breast cancer 1st degree"
table1::label(df7$smoker_bin) <- "Smoking status"
table1::label(df7$alcohol_bin) <- "Current alcohol consumption"
table1::label(df7$GPAQ_q16_sitting_recline_min) <- "Time spent sitting (minutes per day)"
table1::label(df7$QL) <- "Global Health Status / QoL (QLQ-C30)"
table1::label(df7$PA) <- "Pain (QLQ-C30)"
table1::label(df7$SL) <- "Insomnia (QLQ-C30)"
table1::label(df7$DY) <- "Dyspnea (QLQ-C30)"
table1::label(df7$BRAS) <- "Arm symptoms-BR23"
table1::label(df7$BRBS) <- "Breast symptoms-BR23"
table1::label(df7$histology_di) <- "Histologic type"
table1::label(df7$b3tumour_side_of_primary) <- "Side of primary tumor"
table1::label(df7$b3surgery_type) <- "Surgery type"
table1::label(df7$b3sys_treatment) <- "Other systemic treatment"
table1::label(df7$b3sys_antiher2) <- "Anti-HER2 therapie"
table1::label(df7$chemo_bin) <- "Chemotherapie"
table1::label(df7$b3radio_breast_dose_Gy) <- "Radiotherapie dose"
table1::label(df7$radio_total_dose_sum) <- "Radiotherapie total dose (including boosts)"
table1::label(df7$b3radio_boost) <- "Radiotherapie boost"
table1::label(df7$fractionation_bin) <- "Fractionation scheme"
table1::label(df7$boost_frac) <- "Fractionation and boost scheme"
table1::label(df7$b3radio_imrt) <- "IMRT"
table1::label(df7$pato_t_cat) <- "Pathologic T-stage"
table1::label(df7$pato_n_cat) <- "Pathologic N-status"

table1::table1(~MFI_phy_bin_base + age_cat + ethnicity_bin + household_bin +
                 education_lmh + bmi_cat + diabetes + hypertension +
                 history_of_heart_disease + ra + depression + 
                 breast_cancer_family_history_1st_degree + 
                 smoker_bin + alcohol_bin + 
                 GPAQ_q16_sitting_recline_min + QL + PA + SL + DY + 
                 BRAS + BRBS + histology_di +
                 b3tumour_side_of_primary + b3surgery_type + systemic_cat + 
                 b3sys_treatment + b3sys_antiher2 + chemo_bin + 
                 b3radio_breast_dose_Gy + radio_total_dose_sum + b3radio_boost +
                 fractionation_bin + boost_frac + b3radio_imrt + pato_t_cat + 
                 pato_n_cat | MFI_phy_bin, data = df7, overall=F, extra.col=list("P-value"= pvalue))


#Table 1.3 Mental fatigue 24 months 

table1::label(df7$MFI_mental_bin_base) <- "Mental fatigue at baseline"
table1::label(df7$age_cat) <- "Age"
table1::label(df7) <- "Ethnicity"
table1::label(df7$household_bin) <- "Living status"
table1::label(df7$education_lmh) <- "Educational level"
table1::label(df7$bmi_cat) <- "BMI"
table1::label(df7$diabetes) <- "Diabetes"
table1::label(df7$hypertension) <- "Hypertension"
table1::label(df7$history_of_heart_disease) <- "Heart disease"
table1::label(df7$ra) <- "Rheumatoid arthritis"
table1::label(df7$depression) <- "Depression"
table1::label(df7$breast_cancer_family_history_1st_degree) <- "Breast cancer 1st degree"
table1::label(df7$smoker_bin) <- "Smoking status"
table1::label(df7$alcohol_bin) <- "Current alcohol consumption"
table1::label(df7$GPAQ_q16_sitting_recline_min) <- "Time spent sitting (minutes per day)"
table1::label(df7$QL) <- "Global Health Status / QoL (QLQ-C30)"
table1::label(df7$PA) <- "Pain (QLQ-C30)"
table1::label(df7$SL) <- "Insomnia (QLQ-C30)"
table1::label(df7$DY) <- "Dyspnea (QLQ-C30)"
table1::label(df7$BRAS) <- "Arm symptoms-BR23"
table1::label(df7$BRBS) <- "Breast symptoms-BR23"
table1::label(df7$histology_di) <- "Histologic type"
table1::label(df7$b3tumour_side_of_primary) <- "Side of primary tumor"
table1::label(df7$b3surgery_type) <- "Surgery type"
table1::label(df7$b3sys_treatment) <- "Other systemic treatment"
table1::label(df7$b3sys_antiher2) <- "Anti-HER2 therapie"
table1::label(df7$chemo_bin) <- "Chemotherapie"
table1::label(df7$b3radio_breast_dose_Gy) <- "Radiotherapie dose"
table1::label(df7$radio_total_dose_sum) <- "Radiotherapie total dose (including boosts)"
table1::label(df7$b3radio_boost) <- "Radiotherapie boost"
table1::label(df7$fractionation_bin) <- "Fractionation scheme"
table1::label(df7$boost_frac) <- "Fractionation and boost scheme"
table1::label(df7$b3radio_imrt) <- "IMRT"
table1::label(df7$pato_t_cat) <- "Pathologic T-stage"
table1::label(df7$pato_n_cat) <- "Pathologic N-status"

table1::table1(~MFI_mental_bin_base + age_cat + ethnicity_bin + household_bin +
                 education_lmh + bmi_cat + diabetes + hypertension +
                 history_of_heart_disease + ra + depression + 
                 breast_cancer_family_history_1st_degree + 
                 smoker_bin + alcohol_bin + 
                 GPAQ_q16_sitting_recline_min + QL + PA + SL + DY + 
                 BRAS + BRBS + histology_di +
                 b3tumour_side_of_primary + b3surgery_type + systemic_cat + 
                 b3sys_treatment + b3sys_antiher2 + chemo_bin + 
                 b3radio_breast_dose_Gy + radio_total_dose_sum + b3radio_boost +
                 fractionation_bin + boost_frac + b3radio_imrt + pato_t_cat + 
                 pato_n_cat | MFI_mental_bin, data = df7, overall=F, extra.col=list("P-value"= pvalue))

#Table 1.4 Reduced activity 24 months 

table1::label(df7$MFI_redact_bin_base) <- "Reduced activity at baseline"
table1::label(df7$age_cat) <- "Age"
table1::label(df7) <- "Ethnicity"
table1::label(df7$household_bin) <- "Living status"
table1::label(df7$education_lmh) <- "Educational level"
table1::label(df7$bmi_cat) <- "BMI"
table1::label(df7$diabetes) <- "Diabetes"
table1::label(df7$hypertension) <- "Hypertension"
table1::label(df7$history_of_heart_disease) <- "Heart disease"
table1::label(df7$ra) <- "Rheumatoid arthritis"
table1::label(df7$depression) <- "Depression"
table1::label(df7$breast_cancer_family_history_1st_degree) <- "Breast cancer 1st degree"
table1::label(df7$smoker_bin) <- "Smoking status"
table1::label(df7$alcohol_bin) <- "Current alcohol consumption"
table1::label(df7$GPAQ_q16_sitting_recline_min) <- "Time spent sitting (minutes per day)"
table1::label(df7$QL) <- "Global Health Status / QoL (QLQ-C30)"
table1::label(df7$PA) <- "Pain (QLQ-C30)"
table1::label(df7$SL) <- "Insomnia (QLQ-C30)"
table1::label(df7$DY) <- "Dyspnea (QLQ-C30)"
table1::label(df7$BRAS) <- "Arm symptoms-BR23"
table1::label(df7$BRBS) <- "Breast symptoms-BR23"
table1::label(df7$histology_di) <- "Histologic type"
table1::label(df7$b3tumour_side_of_primary) <- "Side of primary tumor"
table1::label(df7$b3surgery_type) <- "Surgery type"
table1::label(df7$b3sys_treatment) <- "Other systemic treatment"
table1::label(df7$b3sys_antiher2) <- "Anti-HER2 therapie"
table1::label(df7$chemo_bin) <- "Chemotherapie"
table1::label(df7$b3radio_breast_dose_Gy) <- "Radiotherapie dose"
table1::label(df7$radio_total_dose_sum) <- "Radiotherapie total dose (including boosts)"
table1::label(df7$b3radio_boost) <- "Radiotherapie boost"
table1::label(df7$fractionation_bin) <- "Fractionation scheme"
table1::label(df7$boost_frac) <- "Fractionation and boost scheme"
table1::label(df7$b3radio_imrt) <- "IMRT"
table1::label(df7$pato_t_cat) <- "Pathologic T-stage"
table1::label(df7$pato_n_cat) <- "Pathologic N-status"

table1::table1(~MFI_redact_bin_base + age_cat + ethnicity_bin + household_bin +
                 education_lmh + bmi_cat + diabetes + hypertension +
                 history_of_heart_disease + ra + depression + 
                 breast_cancer_family_history_1st_degree + 
                 smoker_bin + alcohol_bin + 
                 GPAQ_q16_sitting_recline_min + QL + PA + SL + DY + 
                 BRAS + BRBS + histology_di +
                 b3tumour_side_of_primary + b3surgery_type + 
                 b3sys_treatment + b3sys_antiher2 + chemo_bin + 
                 b3radio_breast_dose_Gy + radio_total_dose_sum + b3radio_boost +
                 fractionation_bin + boost_frac + b3radio_imrt + pato_t_cat + 
                 pato_n_cat | MFI_redact_bin, data = df7, overall=F, extra.col=list("P-value"= pvalue))


#Table 1.4 Reduced motivation 24 months 

table1::label(df7$MFI_redmot_bin_base) <- "Reduced motivation at baseline"
table1::label(df7$age_cat) <- "Age"
table1::label(df7) <- "Ethnicity"
table1::label(df7$household_bin) <- "Living status"
table1::label(df7$education_lmh) <- "Educational level"
table1::label(df7$bmi_cat) <- "BMI"
table1::label(df7$diabetes) <- "Diabetes"
table1::label(df7$hypertension) <- "Hypertension"
table1::label(df7$history_of_heart_disease) <- "Heart disease"
table1::label(df7$ra) <- "Rheumatoid arthritis"
table1::label(df7$depression) <- "Depression"
table1::label(df7$breast_cancer_family_history_1st_degree) <- "Breast cancer 1st degree"
table1::label(df7$smoker_bin) <- "Smoking status"
table1::label(df7$alcohol_bin) <- "Current alcohol consumption"
table1::label(df7$GPAQ_q16_sitting_recline_min) <- "Time spent sitting (minutes per day)"
table1::label(df7$QL) <- "Global Health Status / QoL (QLQ-C30)"
table1::label(df7$PA) <- "Pain (QLQ-C30)"
table1::label(df7$SL) <- "Insomnia (QLQ-C30)"
table1::label(df7$DY) <- "Dyspnea (QLQ-C30)"
table1::label(df7$BRAS) <- "Arm symptoms-BR23"
table1::label(df7$BRBS) <- "Breast symptoms-BR23"
table1::label(df7$histology_di) <- "Histologic type"
table1::label(df7$b3tumour_side_of_primary) <- "Side of primary tumor"
table1::label(df7$b3surgery_type) <- "Surgery type"
table1::label(df7$b3sys_treatment) <- "Other systemic treatment"
table1::label(df7$b3sys_antiher2) <- "Anti-HER2 therapie"
table1::label(df7$chemo_bin) <- "Chemotherapie"
table1::label(df7$b3radio_breast_dose_Gy) <- "Radiotherapie dose"
table1::label(df7$radio_total_dose_sum) <- "Radiotherapie total dose (including boosts)"
table1::label(df7$b3radio_boost) <- "Radiotherapie boost"
table1::label(df7$fractionation_bin) <- "Fractionation scheme"
table1::label(df7$boost_frac) <- "Fractionation and boost scheme"
table1::label(df7$b3radio_imrt) <- "IMRT"
table1::label(df7$pato_t_cat) <- "Pathologic T-stage"
table1::label(df7$pato_n_cat) <- "Pathologic N-status"

table1::table1(~MFI_redmot_bin_base + age_cat + ethnicity_bin + household_bin +
                 education_lmh + bmi_cat + diabetes + hypertension +
                 history_of_heart_disease + ra + depression + 
                 breast_cancer_family_history_1st_degree + 
                 smoker_bin + alcohol_bin + 
                 GPAQ_q16_sitting_recline_min + QL + PA + SL + DY + 
                 BRAS + BRBS + histology_di +
                 b3tumour_side_of_primary + b3surgery_type + 
                 b3sys_treatment + b3sys_antiher2 + chemo_bin + 
                 b3radio_breast_dose_Gy + radio_total_dose_sum + b3radio_boost +
                 fractionation_bin + boost_frac + b3radio_imrt + pato_t_cat + 
                 pato_n_cat | MFI_redmot_bin, data = df7, overall=F, extra.col=list("P-value"= pvalue))



table1::label(df7_general_complete4$MFI_general_bin) <- "General Fatigue"
table1::label(df7_general_complete4$Mainsite) <- "Mainsite"

table1::table1(~Mainsite | MFI_general_bin, data = df7_general_complete4, 
               overall=F, extra.col=list("P-value"= pvalue))



table1::label(df7_general_complete4$MFI_general_bin) <- "General Fatigue"
table1::label(df7_general_complete4$Mainsite) <- "Mainsite"

table1::table1(~MFI_general_bin | Mainsite, data = df7_general_complete4, 
               overall=F, extra.col=list("P-value"= pvalue))

# Chemotherapy by baseline general fatigue 
table1::label(df7_general_complete4$MFI_general_bin_base) <- "Baseline General Fatigue"
table1::label(df7_general_complete4$MFI_general_bin) <- "General Fatigue"
table1::label(df7_general_complete4$chemo_bin) <- "Chemotherapy"
table1::label(df7_general_complete4$b3sys_treatment) <- "Other systemic tx"


# Cross tables for fatigue (base and 24 months), chemo and systemic
table1::table1(~MFI_general_bin_base | chemo_bin, 
               data = df7_general, 
               overall=F, extra.col=list("P-value"= pvalue))

table1::table1(~MFI_general_bin_base | b3sys_treatment, 
               data = df7_general, 
               overall=F, extra.col=list("P-value"= pvalue))

table1::table1(~MFI_general_bin | chemo_bin, 
               data = df7_general, 
               overall=F, extra.col=list("P-value"= pvalue))

table1::table1(~MFI_general_bin | b3sys_treatment, 
               data = df7_general, 
               overall=F, extra.col=list("P-value"= pvalue))

table1::table1(~chemo_bin + b3sys_treatment | MFI_general_bin_base, 
               data = df7_general, 
               overall=F, extra.col=list("P-value"= pvalue))

table1::table1(~chemo_bin + b3sys_treatment | MFI_general_bin, 
               data = df7_general, 
               overall=F, extra.col=list("P-value"= pvalue))


table1::table1(~MFI_general_b | chemo_bin, 
               data = df7_general, 
               overall=F, extra.col=list("P-value"= pvalue))

table1::table1(~MFI_general | b3sys_treatment, 
               data = df7_general, 
               overall=F, extra.col=list("P-value"= pvalue))

table1::table1(~MFI_general_b | b3sys_treatment, 
               data = df7_general, 
               overall=F, extra.col=list("P-value"= pvalue))

table1::table1(~MFI_general | b3sys_treatment, 
               data = df7_general, 
               overall=F, extra.col=list("P-value"= pvalue))


######
#####
###
##
###################  LOGISTIC REGRESSION (BY CHEMO STATUS) ############
## 
###
####
######


##1. GENERAL FATIGUE 

df7_general_complete4 %>% 
  dplyr::filter(df7_general_complete4$chemo_bin == "Yes") -> df7_chemo_yes

df7_general_complete4 %>% 
  dplyr::filter(df7_general_complete4$chemo_bin == "No") -> df7_chemo_no


df7_general_complete4


## 1.1.a General Fatigue (chemo yes): Logistic Regression 

mod_general_chemo_yes <- glm(MFI_general_bin~MFI_general_bin_base + 
                           age_5 + bmi + 
                           history_of_heart_disease + depression +
                           smoker_bin + alcohol_bin + 
                           QL_10 + PA_10 + SL_10 + DY_10 + 
                           BRAS_10 + BRBS_10 + 
                           b3surgery_type + b3sys_treatment + 
                           fractionation_bin+ b3radio_boost + b3radio_imrt +
                           pato_t_cat, family="binomial", 
                         data = df7_chemo_yes)

summary(mod_general_chemo_yes)

## 1.2.a General Fatigue (chemo yes): Backward selection with stepAIC

step_gen_chemo_yes <- stepAIC(mod_general_chemo_yes, direction = "both", 
                           trace = F)

summary(step_gen_chemo_yes)

## 1.3.a Include adjusting variables 
gen_chemo_yes_final <- glm(MFI_general_bin~MFI_general_bin_base + 
                             age_5 + bmi + 
                             QL_10 + SL_10 + 
                             BRBS_10 + b3radio_boost + b3radio_imrt, family="binomial", 
                           data = df7_chemo_yes)


summary(gen_chemo_yes_final)

## 1.4.a Nice table with OR and CI
gen_chemo_yes_final %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_general_bin_base", "age_5","bmi",
               "QL_10", "SL_10", "BRBS_10", "b3radio_boost", "b3radio_imrt"), 
    show_single_row = c("MFI_general_bin_base", "b3radio_boost"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 

## 1.1.b General Fatigue (chemo no): Logistic Regression 

mod_general_chemo_no <- glm(MFI_general_bin~MFI_general_bin_base + 
                              age_5 + bmi + 
                              history_of_heart_disease + depression +
                              smoker_bin + alcohol_bin + 
                              QL_10 + PA_10 + SL_10 + DY_10 + 
                              BRAS_10 + BRBS_10 + 
                              b3surgery_type + b3sys_treatment + 
                              fractionation_bin+ b3radio_boost + b3radio_imrt +
                              pato_t_cat, family="binomial", 
                             data = df7_chemo_no)

summary(mod_general_chemo_no)

## 1.2.b General Fatigue (chemo no: Backward selection with stepAIC

step_gen_chemo_no <- stepAIC(mod_general_chemo_no, direction = "both", 
                              trace = F)

summary(step_gen_chemo_no)

## 1.3.b Include adjusting variables 
gen_chemo_no_final <- glm(MFI_general_bin~MFI_general_bin_base + 
                            age_5 + depression + 
                            PA_10 + SL_10 + 
                            DY_10 + pato_t_cat,           
                            family="binomial", data = df7_chemo_no)

## 1.4.b Nice table with OR and CI
gen_chemo_no_final %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_general_bin_base", "age_5", "depression",
               "PA_10", "SL_10","DY_10", "pato_t_cat"), 
    show_single_row = c("MFI_general_bin_base","depression"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 


##2. PHYSICAL FATIGUE 

df7_phy_complete4 %>% 
  filter(df7_phy_complete4$chemo_bin == "Yes") -> df7_phy_chemo_yes

df7_phy_complete4 %>% 
  filter(df7_phy_complete4$chemo_bin == "No") -> df7_phy_chemo_no



## 2.1.a Physical Fatigue (chemo yes): Logistic Regression 

mod_phy_chemo_yes <- glm(MFI_phy_bin~MFI_phy_bin_base + 
                           age_5 + bmi + 
                           history_of_heart_disease + depression +
                           smoker_bin + alcohol_bin + 
                           QL_10 + PA_10 + SL_10 + DY_10 + 
                           BRAS_10 + BRBS_10 + 
                           b3surgery_type + b3sys_treatment + 
                           fractionation_bin+ b3radio_boost + b3radio_imrt +
                           pato_t_cat, family="binomial", 
                             data = df7_phy_chemo_yes)

summary(mod_phy_chemo_yes)

## 2.2.a Physical Fatigue (chemo yes): Backward selection with stepAIC

step_phy_chemo_yes <- stepAIC(mod_phy_chemo_yes, direction = "both", 
                              trace = F)

summary(step_phy_chemo_yes)

## 2.3.a Include adjusting variables 
phy_chemo_yes_final <- glm(MFI_phy_bin~MFI_phy_bin_base + age_5 + bmi + SL_10 + DY_10 +
                             BRBS_10 + b3sys_treatment + fractionation_bin, family="binomial", 
                           data = df7_phy_chemo_yes)


summary(phy_chemo_yes_final)

## 2.4.a Nice table with OR and CI
phy_chemo_yes_final %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_phy_bin_base", "age_5","bmi", 
               "SL_10", "DY_10", "BRBS_10", "b3sys_treatment", "fractionation_bin"), 
    show_single_row = c("MFI_phy_bin_base", "b3sys_treatment", "fractionation_bin"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 

## 2.1.b Physical Fatigue (chemo no): Logistic Regression 

mod_phy_chemo_no <- glm(MFI_phy_bin~MFI_phy_bin_base + age_5 + bmi + 
                          history_of_heart_disease + depression +
                          smoker_bin + alcohol_bin + 
                          QL_10 + PA_10 + SL_10 + DY_10 + 
                          BRAS_10 + BRBS_10 + 
                          b3surgery_type + b3sys_treatment + 
                          fractionation_bin+ b3radio_boost + b3radio_imrt +
                          pato_t_cat, family="binomial", 
                         data = df7_phy_chemo_no)

summary(mod_phy_chemo_no)

## 2.2.b Physical Fatigue (chemo no): Backward selection with stepAIC

step_phy_chemo_no <- stepAIC(mod_phy_chemo_no, direction = "both", 
                              trace = F)

summary(step_phy_chemo_no)

## 2.3.b Include adjusting variables 
phy_chemo_no_final <- glm(MFI_phy_bin~MFI_phy_bin_base + 
                             age_5 + bmi  
                          + depression + alcohol_bin + PA_10 + SL_10 + DY_10 + BRBS_10,
                          family="binomial", 
                           data = df7_phy_chemo_no)


summary(phy_chemo_no_final)

## 2.4.b Nice table with OR and CI
phy_chemo_no_final %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_phy_bin_base", "age_5","bmi", "depression","alcohol_bin", "PA_10",
               "SL_10", "DY_10", "BRBS_10"), 
    show_single_row = c("MFI_phy_bin_base", "depression", 
                        "alcohol_bin"), 
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 

##3. MENTAL FATIGUE 

df7_mental_complete4 %>% 
  filter(df7_mental_complete4$chemo_bin == "Yes") -> df7_mental_chemo_yes

df7_mental_complete4 %>% 
  filter(df7_mental_complete4$chemo_bin == "No") -> df7_mental_chemo_no



## 3.1.a Mental Fatigue (chemo yes): Logistic Regression 

mod_mental_chemo_yes <- glm(MFI_mental_bin~MFI_mental_bin_base + age_5 + bmi + 
                              history_of_heart_disease + depression +
                              smoker_bin + alcohol_bin + 
                              QL_10 + PA_10 + SL_10 + DY_10 + 
                              BRAS_10 + BRBS_10 + 
                              b3surgery_type + b3sys_treatment + 
                              fractionation_bin+ b3radio_boost + b3radio_imrt +
                              pato_t_cat, family="binomial", 
                         data = df7_mental_chemo_yes)

summary(mod_mental_chemo_yes)

## 3.2.a Mental Fatigue (chemo yes): Backward selection with stepAIC

step_mental_chemo_yes <- stepAIC(mod_mental_chemo_yes, direction = "both", 
                              trace = F)

summary(step_mental_chemo_yes)

## 3.3.a Include adjusting variables 
mental_chemo_yes_final <- glm(MFI_mental_bin~MFI_mental_bin_base + 
                             age_5 + PA_10, family="binomial", 
                           data = df7_mental_chemo_yes)


summary(mental_chemo_yes_final)

## 3.4.a Nice table with OR and CI
mental_chemo_yes_final %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_mental_bin_base", "age_5", "PA_10"), 
    show_single_row = c("MFI_mental_bin_base"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 

## 3.1.b Mental Fatigue (chemo no): Logistic Regression 

mod_mental_chemo_no <- glm(MFI_mental_bin~MFI_mental_bin_base + age_5 + bmi + 
                             history_of_heart_disease + depression +
                             smoker_bin + alcohol_bin + 
                             QL_10 + PA_10 + SL_10 + DY_10 + 
                             BRAS_10 + BRBS_10 + 
                             b3surgery_type + b3sys_treatment + 
                             fractionation_bin+ b3radio_boost + b3radio_imrt +
                             pato_t_cat, family="binomial", 
                            data = df7_mental_chemo_no)

summary(mod_mental_chemo_no)

## 3.2.b Mental Fatigue (chemo no): Backward selection with stepAIC

step_mental_chemo_no <- stepAIC(mod_mental_chemo_no, direction = "both", 
                                 trace = F)

summary(step_mental_chemo_no)

## 3.3.b Include adjusting variables 
mental_chemo_no_final <- glm(MFI_mental_bin~MFI_mental_bin_base + 
                                age_5 + depression + alcohol_bin +
                                QL_10 + SL_10 + b3radio_boost, family="binomial", 
                                data = df7_mental_chemo_no)


summary(mental_chemo_no_final)

## 3.4.b Nice table with OR and CI
mental_chemo_no_final %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_mental_bin_base", "age_5","depression", "alcohol_bin", 
               "QL_10", "SL_10", "b3radio_boost"), 
    show_single_row = c("MFI_mental_bin_base", "depression", "alcohol_bin", 
                        "b3radio_boost"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 



##4. REDUCED MOTIVATION

df7_redmot_complete4 %>% 
  filter(df7_redmot_complete4$chemo_bin == "Yes") -> df7_redmot_chemo_yes

df7_redmot_complete4 %>% 
  filter(df7_redmot_complete4$chemo_bin == "No") -> df7_redmot_chemo_no



## 4.1.a Reduced motivation (chemo yes): Logistic Regression 

mod_redmot_chemo_yes <- glm(MFI_redmot_bin~MFI_redmot_bin_base + 
                              age_5 + bmi + 
                              history_of_heart_disease + depression +
                              smoker_bin + alcohol_bin + 
                              QL_10 + PA_10 + SL_10 + DY_10 + 
                              BRAS_10 + BRBS_10 + 
                              b3surgery_type + b3sys_treatment + 
                              fractionation_bin+ b3radio_boost + b3radio_imrt +
                              pato_t_cat, family="binomial", 
                            data = df7_redmot_chemo_yes)

summary(mod_redmot_chemo_yes)

## 4.2.a Reduced motivation (chemo yes): Backward selection with stepAIC

step_redmot_chemo_yes <- stepAIC(mod_redmot_chemo_yes, direction = "both", 
                                 trace = F)

summary(step_redmot_chemo_yes)

## 4.3.a Include adjusting variables 
redmot_chemo_yes_final <- glm(MFI_redmot_bin~MFI_redmot_bin_base + 
                                age_5 + b3surgery_type + 
                                depression, family="binomial", 
                              data = df7_redmot_chemo_yes)


summary(redmot_chemo_yes_final)

## 4.4.a Nice table with OR and CIhypertension
redmot_chemo_yes_final %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_redmot_bin_base", "age_5", "b3surgery_type",
               "depression"), 
    show_single_row = c("MFI_redmot_bin_base", "b3surgery_type", 
                        "depression"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 

## 4.1.b Reduced motivation (chemo no): Logistic Regression 

mod_redmot_chemo_no <- glm(MFI_redmot_bin~MFI_redmot_bin_base + age_5 + bmi + 
                             history_of_heart_disease + depression +
                             smoker_bin + alcohol_bin + 
                             QL_10 + PA_10 + SL_10 + DY_10 + 
                             BRAS_10 + BRBS_10 + 
                             b3surgery_type + b3sys_treatment + 
                             fractionation_bin+ b3radio_boost + b3radio_imrt +
                             pato_t_cat, family="binomial", 
                            data = df7_redmot_chemo_no)

summary(mod_redmot_chemo_no)

## 4.2.b Reduced motivation (chemo no): Backward selection with stepAIC

step_redmot_chemo_no <- stepAIC(mod_redmot_chemo_no, direction = "both", 
                                 trace = F)

summary(step_redmot_chemo_no)

## 4.3.b Include adjusting variables 
redmot_chemo_no_final <- glm(MFI_redmot_bin~MFI_redmot_bin_base + 
                                age_5 + QL_10 + DY_10 +
                                b3sys_treatment, family="binomial", 
                              data = df7_redmot_chemo_no)


summary(redmot_chemo_no_final)

## 4.4.b Nice table with OR and CI
redmot_chemo_no_final %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_redmot_bin_base", "age_5", "QL_10", "DY_10",
               "b3sys_treatment"), 
    show_single_row = c("MFI_redmot_bin_base", "b3sys_treatment"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels()



##5. REDUCED ACTIVITY 

df7_redact_complete4 %>% 
  filter(df7_redact_complete4$chemo_bin == "Yes") -> df7_redact_chemo_yes

df7_redact_complete4 %>% 
  filter(df7_redact_complete4$chemo_bin == "No") -> df7_redact_chemo_no



## 5.1.a Reduced activity (chemo yes): Logistic Regression 

mod_redact_chemo_yes <- glm(MFI_redact_bin~MFI_redact_bin_base + age_5 + bmi + 
                              history_of_heart_disease + depression +
                              smoker_bin + alcohol_bin + 
                              QL_10 + PA_10 + SL_10 + DY_10 + 
                              BRAS_10 + BRBS_10 + 
                              b3surgery_type + b3sys_treatment + 
                              fractionation_bin+ b3radio_boost + b3radio_imrt +
                              pato_t_cat, family="binomial", 
                            data = df7_redact_chemo_yes)

summary(mod_redact_chemo_yes)

table(df7_redact_chemo_yes$MFI_redact_bin)
table(df7_redact_chemo_yes$MFI_redact_bin, df7_redact_chemo_yes$hypertension)
table(df7_redact_chemo_yes$MFI_redact_bin, df7_redact_chemo_yes$depression)
table(df7_redact_chemo_yes$MFI_redact_bin, df7_redact_chemo_yes$diabetes)
table(df7_redact_chemo_yes$MFI_redact_bin, df7_redact_chemo_yes$breast_cancer_family_history_1st_degree)
table(df7_redact_chemo_yes$MFI_redact_bin, df7_redact_chemo_yes$b3surgery_type)
table(df7_redact_chemo_yes$MFI_redact_bin, df7_redact_chemo_yes$b3radio_boost)


## 5.2.a Reduced activity (chemo yes): Backward selection with stepAIC

step_redact_chemo_yes <- stepAIC(mod_redact_chemo_yes, direction = "both", 
                                 trace = F)

summary(step_redact_chemo_yes)

## 5.3.a Include adjusting variables 
redact_chemo_yes_final <- glm(MFI_redact_bin~MFI_redact_bin_base + 
                                age_5 + bmi + QL_10 + PA_10 + SL_10 + DY_10
                              , family="binomial", 
                              data = df7_redact_chemo_yes)


summary(redact_chemo_yes_final)

## 5.4.a Nice table with OR and CI
redact_chemo_yes_final %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_redact_bin_base", "age_5","bmi", "QL_10", "PA_10", "SL_10", 
               "DY_10"), 
    show_single_row = c("MFI_redact_bin_base"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 

## 5.1.b Reduced activity (chemo no): Logistic Regression 

mod_redact_chemo_no <- glm(MFI_redact_bin~MFI_redact_bin_base + age_5 + bmi + 
                             history_of_heart_disease + depression +
                             smoker_bin + alcohol_bin + 
                             QL_10 + PA_10 + SL_10 + DY_10 + 
                             BRAS_10 + BRBS_10 + 
                             b3surgery_type + b3sys_treatment + 
                             fractionation_bin+ b3radio_boost + b3radio_imrt +
                             pato_t_cat, family="binomial", 
                            data = df7_redact_chemo_no)

summary(mod_redact_chemo_no)

## 5.2.b Reduced activity (chemo no): Backward selection with stepAIC

step_redact_chemo_no <- stepAIC(mod_redact_chemo_no, direction = "both", 
                                 trace = F)

summary(step_redact_chemo_no)

## 5.3.b Include adjusting variables 
redact_chemo_no_final <- glm(MFI_redact_bin~MFI_redact_bin_base + 
                                age_5 + depression + QL_10 + DY_10 +
                                b3radio_imrt, family="binomial", 
                              data = df7_redact_chemo_no)


summary(redact_chemo_no_final)

## 5.4.b Nice table with OR and CI
redact_chemo_no_final %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_redact_bin_base", "age_5","depression", "QL_10", "DY_10", 
               "b3radio_imrt"), 
    show_single_row = c("MFI_redact_bin_base", "depression", "b3radio_imrt"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 

table(df7_redact_chemo_no$MFI_redact_bin, df7_redact_chemo_no$pato_t_cat)


plot_yes <- plot_summs(gen_chemo_yes_final, phy_chemo_yes_final, 
           mental_chemo_yes_final, redmot_chemo_yes_final,
           redact_chemo_yes_final, model.names = c("General fatigue", "Physical fatigue", 
                                                   "Mental fatigue", "Reduced motivation", 
                                                    "Reduced activity"))


plot_no <- plot_summs(gen_chemo_no_final, phy_chemo_no_final, 
           mental_chemo_no_final, redmot_chemo_no_final,
           redact_chemo_no_final, model.names = c("General fatigue", "Physical fatigue", 
                                                   "Mental fatigue", "Reduced motivation", 
                                                   "Reduced activity"))
#
####
#########
############## LOGISTIC REGRESSION BY CHEMO STATUS (FULL COHORT) ############
#########
####
#

#
### GENERAL FATIGUE
#

# First create new dataframe with only the final variables from the previous
# logistic models

df7 %>%
  dplyr::select(Subject.Id, chemo_bin, age_at_radiotherapy_start_yrs, age_5, 
                bmi, QL_10, SL_10, BRBS_10, fractionation_bin, b3radio_boost, 
                boost_frac, MFI_general_bin_base, 
                MFI_general_bin) -> df7_gf_final_only1

df7_gf_complete_final_only1 <- na.omit(df7_gf_final_only1)


df7 %>%
  dplyr::select(Subject.Id, chemo_bin, age_at_radiotherapy_start_yrs, age_5, 
                depression, PA_10, SL_10, DY_10, MFI_general_bin_base, 
                MFI_general_bin) -> df7_gf_final_only2

df7_gf_complete_final_only2 <- na.omit(df7_gf_final_only2)

# Divide "df7_gf_complete_final_only1" into two chemo dataframes 

df7_gf_complete_final_only1 %>% 
  dplyr::filter(df7_gf_complete_final_only1$chemo_bin == "Yes") -> df7_chemo_only_yes1 #320

df7_gf_complete_final_only2 %>% 
  dplyr::filter(df7_gf_complete_final_only2$chemo_bin == "No") -> df7_chemo_only_no2 # 646



# Now test the SEVENTH logistic models with the previously identified 
# variables and with the full number of subjects for those variables


## Chemo= yes 
gen_chemo_yes_only_final1 <- glm(MFI_general_bin~MFI_general_bin_base + 
                             age_5 + bmi + 
                             QL_10 + SL_10 + 
                             BRBS_10 + b3radio_boost, family="binomial", 
                           data = df7_chemo_only_yes1)


summary(gen_chemo_yes_only_final1)

## 1.4.a Nice table with OR and CI
gen_chemo_yes_only_final1 %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_general_bin_base", "age_5","bmi",
               "QL_10", "SL_10", "BRBS_10", "b3radio_boost"), 
    show_single_row = c("MFI_general_bin_base", "b3radio_boost"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 


## Chemo= no
gen_chemo_no_only_final2 <- glm(MFI_general_bin~MFI_general_bin_base + 
                                   age_5  + depression +
                                   PA_10 + SL_10 + DY_10, family="binomial", 
                                 data = df7_chemo_only_no2)


summary(gen_chemo_no_only_final2)

## 1.4.a Nice table with OR and CI
gen_chemo_no_only_final2 %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_general_bin_base", "age_5",
               "PA_10", "SL_10", "DY_10",  "depression"), 
    show_single_row = c("MFI_general_bin_base", "depression"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 


#
### PHYSICAL FATIGUE
#

# First create new data frame with only the final variables from the previous
# logistic models

df7 %>%
  dplyr::select(Subject.Id, chemo_bin, age_at_radiotherapy_start_yrs, age_5, 
                bmi , SL_10 , DY_10 ,BRBS_10, b3sys_treatment, 
                fractionation_bin, MFI_phy_bin_base, 
                MFI_phy_bin) -> df7_pf_final_only1

df7_pf_complete_final_only1 <- na.omit(df7_pf_final_only1)


df7 %>%
  dplyr::select(Subject.Id, chemo_bin, age_at_radiotherapy_start_yrs, age_5, 
                bmi, depression, alcohol_bin, PA_10, SL_10, DY_10, BRBS_10, 
                MFI_phy_bin_base, 
                MFI_phy_bin) -> df7_pf_final_only2

df7_pf_complete_final_only2 <- na.omit(df7_pf_final_only2)

# Filter by chemo status for the data frames:
# "df7_pf_complete_final_only1" and 
# "df7_pf_complete_final_only2"

df7_pf_complete_final_only1 %>% 
  dplyr::filter(df7_pf_complete_final_only1$chemo_bin == "Yes") -> df7_chemo_pf_only_yes1 #

df7_pf_complete_final_only2 %>% 
  dplyr::filter(df7_pf_complete_final_only2$chemo_bin == "No") -> df7_chemo_pf_only_no2 #



# Now test the SEVENTH logistic models with the previously identified 
# variables and with the full number of subjects for those variables


## Chemo= yes 
pf_chemo_yes_only_final1 <- glm(MFI_phy_bin~MFI_phy_bin_base + age_5 + bmi + SL_10 + DY_10 +
                                  BRBS_10 + b3sys_treatment + fractionation_bin, family="binomial", 
                                data = df7_chemo_pf_only_yes1)


summary(pf_chemo_yes_only_final1)

## Nice table with OR and CI
pf_chemo_yes_only_final1 %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_phy_bin_base", "age_5","bmi", 
               "SL_10", "DY_10", "BRBS_10", "b3sys_treatment", "fractionation_bin"), 
    show_single_row = c("MFI_phy_bin_base", "b3sys_treatment", "fractionation_bin"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 


## Chemo= no
pf_chemo_no_only_final2 <- glm(MFI_phy_bin~MFI_phy_bin_base + 
                                 age_5 + bmi + depression + alcohol_bin + 
                                 PA_10 + SL_10 + DY_10 + BRBS_10,
                               family="binomial", 
                               data = df7_chemo_pf_only_no2)


summary(pf_chemo_no_only_final2)

## Nice table with OR and CI
pf_chemo_no_only_final2 %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_phy_bin_base", "age_5","bmi", "depression","alcohol_bin", "PA_10",
               "SL_10", "DY_10", "BRBS_10"), 
    show_single_row = c("MFI_phy_bin_base", "depression", 
                        "alcohol_bin"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 


#
### MENTAL FATIGUE
#

# First create new data frame with only the final variables from the previous
# logistic models

df7 %>%
  dplyr::select(Subject.Id, chemo_bin, age_at_radiotherapy_start_yrs, age_5, 
                PA_10, MFI_mental_bin_base, 
                MFI_mental_bin) -> df7_mf_final_only1

df7_mf_complete_final_only1 <- na.omit(df7_mf_final_only1)


df7 %>%
  dplyr::select(Subject.Id, chemo_bin, age_at_radiotherapy_start_yrs, age_5, 
                depression, alcohol_bin, QL_10, SL_10, b3radio_boost, 
                MFI_mental_bin_base, 
                MFI_mental_bin) -> df7_mf_final_only2

df7_mf_complete_final_only2 <- na.omit(df7_mf_final_only2)

# Filter by chemo status for the data frames:
# "df7_mf_complete_final_only1" and 
# "df7_mf_complete_final_only2"

df7_mf_complete_final_only1 %>% 
  dplyr::filter(df7_mf_complete_final_only1$chemo_bin == "Yes") -> df7_chemo_mf_only_yes1 #

df7_mf_complete_final_only2 %>% 
  dplyr::filter(df7_mf_complete_final_only2$chemo_bin == "No") -> df7_chemo_mf_only_no2 #



# Now test the SEVENTH logistic models with the previously identified 
# variables and with the full number of subjects for those variables


## Chemo= yes 
mf_chemo_yes_only_final1 <- glm(MFI_mental_bin~MFI_mental_bin_base + 
                                  age_5 + PA_10, family="binomial", 
                                data = df7_chemo_mf_only_yes1)


summary(mf_chemo_yes_only_final1)

## Nice table with OR and CI
mf_chemo_yes_only_final1 %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_mental_bin_base", "age_5", "PA_10"), 
    show_single_row = c("MFI_mental_bin_base"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 


## Chemo= no
mf_chemo_no_only_final2 <- glm(MFI_mental_bin~MFI_mental_bin_base + 
                                 age_5 + depression + alcohol_bin +
                                 QL_10 + SL_10 + b3radio_boost, family="binomial", 
                               data = df7_chemo_mf_only_no2)


summary(mf_chemo_no_only_final2)

## Nice table with OR and CI
mf_chemo_no_only_final2 %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_mental_bin_base", "age_5","depression", "alcohol_bin", 
               "QL_10", "SL_10", "b3radio_boost"), 
    show_single_row = c("MFI_mental_bin_base", "depression", "alcohol_bin", 
                        "b3radio_boost"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 



#
### REDUCED MOTIVATION
#

# First create new data frame with only the final variables from the previous
# logistic models

df7 %>%
  dplyr::select(Subject.Id, chemo_bin, age_at_radiotherapy_start_yrs, age_5, 
                b3surgery_type, depression, MFI_redmot_bin_base, 
                MFI_redmot_bin) -> df7_rm_final_only1

df7_rm_complete_final_only1 <- na.omit(df7_rm_final_only1)


df7 %>%
  dplyr::select(Subject.Id, chemo_bin, age_at_radiotherapy_start_yrs, age_5, 
                QL_10, DY_10, b3sys_treatment, MFI_redmot_bin_base, 
                MFI_redmot_bin) -> df7_rm_final_only2

df7_rm_complete_final_only2 <- na.omit(df7_rm_final_only2)

# Filter by chemo status for the data frames:
# "df7_mf_complete_final_only1" and 
# "df7_mf_complete_final_only2"

df7_rm_complete_final_only1 %>% 
  dplyr::filter(df7_rm_complete_final_only1$chemo_bin == "Yes") -> df7_chemo_rm_only_yes1 #

df7_rm_complete_final_only2 %>% 
  dplyr::filter(df7_rm_complete_final_only2$chemo_bin == "No") -> df7_chemo_rm_only_no2 #



# Now test the SEVENTH logistic models with the previously identified 
# variables and with the full number of subjects for those variables


## Chemo= yes 
rm_chemo_yes_only_final1 <- glm(MFI_redmot_bin~MFI_redmot_bin_base + 
                                  age_5 + b3surgery_type + 
                                  depression, family="binomial", 
                                data = df7_chemo_rm_only_yes1)


summary(mf_chemo_yes_only_final1)

## Nice table with OR and CI
rm_chemo_yes_only_final1 %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_redmot_bin_base", "age_5", "b3surgery_type",
               "depression"), 
    show_single_row = c("MFI_redmot_bin_base", "b3surgery_type", 
                        "depression"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 


## Chemo= no
rm_chemo_no_only_final2 <- glm(MFI_redmot_bin~MFI_redmot_bin_base + 
                                 age_5 + QL_10 + DY_10 +
                                 b3sys_treatment, family="binomial", 
                               data = df7_chemo_rm_only_no2)


summary(mf_chemo_no_only_final2)

## Nice table with OR and CI
rm_chemo_no_only_final2 %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_redmot_bin_base", "age_5", "QL_10", "DY_10",
               "b3sys_treatment"), 
    show_single_row = c("MFI_redmot_bin_base", "b3sys_treatment"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 


#
### REDUCED ACTIVITY
#

# First create new data frame with only the final variables from the previous
# logistic models

df7 %>%
  dplyr::select(Subject.Id, chemo_bin, age_at_radiotherapy_start_yrs, age_5, 
                bmi, QL_10, PA_10, SL_10, DY_10, MFI_redact_bin_base, 
                MFI_redact_bin) -> df7_ra_final_only1

df7_ra_complete_final_only1 <- na.omit(df7_ra_final_only1)


df7 %>%
  dplyr::select(Subject.Id, chemo_bin, age_at_radiotherapy_start_yrs, age_5, 
                depression, QL_10, DY_10, b3radio_imrt, MFI_redact_bin_base, 
                MFI_redact_bin) -> df7_ra_final_only2

df7_ra_complete_final_only2 <- na.omit(df7_ra_final_only2)

# Filter by chemo status for the data frames:
# "df7_mf_complete_final_only1" and 
# "df7_mf_complete_final_only2"

df7_ra_complete_final_only1 %>% 
  dplyr::filter(df7_ra_complete_final_only1$chemo_bin == "Yes") -> df7_chemo_ra_only_yes1 #

df7_ra_complete_final_only2 %>% 
  dplyr::filter(df7_ra_complete_final_only2$chemo_bin == "No") -> df7_chemo_ra_only_no2 #



# Now test the SEVENTH logistic models with the previously identified 
# variables and with the full number of subjects for those variables


## Chemo= yes 
ra_chemo_yes_only_final1 <- glm(MFI_redact_bin~MFI_redact_bin_base + 
                                  age_5 + bmi + QL_10 + PA_10 + SL_10 + DY_10
                                , family="binomial", 
                                data = df7_chemo_ra_only_yes1)


summary(ra_chemo_yes_only_final1)

## Nice table with OR and CI
ra_chemo_yes_only_final1 %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_redact_bin_base", "age_5","bmi", "QL_10", "PA_10", "SL_10", 
               "DY_10"), 
    show_single_row = c("MFI_redact_bin_base"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 


## Chemo= no
ra_chemo_no_only_final2 <- glm(MFI_redact_bin~MFI_redact_bin_base + 
                                 age_5 + depression + QL_10 + DY_10 +
                                 b3radio_imrt, family="binomial", 
                               data = df7_chemo_ra_only_no2)


summary(ra_chemo_no_only_final2)

## Nice table with OR and CI
ra_chemo_no_only_final2 %>%
  tbl_regression(
    exponentiate = TRUE,
    include= c("MFI_redact_bin_base", "age_5","depression", "QL_10", "DY_10", 
               "b3radio_imrt"), 
    show_single_row = c("MFI_redact_bin_base", "depression", "b3radio_imrt"),
    add_estimate_to_reference_rows = T,
    pvalue_fun = ~style_pvalue(.x, digits = 3), ) %>% 
  add_global_p() %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() 



###############################################################################


# FOREST PLOT FOR MULTIPLE REGRESSION MODELS 
# https://strengejacke.github.io/sjPlot/reference/plot_models.html

plot_models(gen_chemo_yes_only_final1, pf_chemo_yes_only_final1, 
            mf_chemo_yes_only_final1, rm_chemo_yes_only_final1, 
            ra_chemo_yes_only_final1, grid = T,   m.labels = c("General fatigue", 
                                                               "Physical fatigue", 
                                                               "Mental fatigue", 
                                                               "Reduced motivation", 
                                                               "Reduced activity"),
            
            show.values = FALSE, show.p = FALSE, p.shape = TRUE)



plot_chemo_all <- plot_models(gen_chemo_yes_only_final1, pf_chemo_yes_only_final1, 
                              mf_chemo_yes_only_final1, rm_chemo_yes_only_final1, 
                              ra_chemo_yes_only_final1, grid = T, axis.lim= c(0.1, 10),  
                              m.labels = c("General fatigue", 
                                           
                                           "Physical fatigue", 
                                           
                                           "Mental fatigue", 
                                           
                                           "Reduced motivation", 
                                           
                                           "Reduced activity"),
                              axis.labels = c("Reduced activity at baseline", "Depression", 
                                              "Wide local excision", "Reduced motivation at baseline", 
                                              "Pain (10 unit increment)", "Mental fatigue at baseline", 
                                              "Standard fractionation", "Hormone cancer therapy",
                                              "Dyspnea (10 unit increment)", "Physical fatigue at baseline", 
                                              "RT boost", 
                                              "Breast symptoms (10 unit increment)", 
                                              "Insomnia (10 unit increment)",
                                              "Global Health Status (10 unit increment)", 
                                              "BMI (continuous, 1 unit increment", 
                                              "Age (continuous, 5 unit increment) ",
                                              "General fatigue at baseline"),
                              vline.color = "snow3", show.values = T, value.size = 3, 
                              show.p = T, p.shape = F)


plot_chemo_all$data <- plot_chemo_all$data %>% arrange(desc(plot_chemo_all$data$estimate))


plot_models(gen_chemo_no_only_final2, pf_chemo_no_only_final2, 
            mf_chemo_no_only_final2, rm_chemo_no_only_final2, 
            ra_chemo_no_only_final2, grid = T, 
            axis.lim= c(0.1, 10),  m.labels = c("General fatigue", 
                                                "Physical fatigue", 
                                                "Mental fatigue", 
                                                "Reduced motivation", 
                                                "Reduced activity"), 
            show.values = FALSE, show.p = FALSE, p.shape = TRUE)





plot_chemo_no_all <- plot_models(gen_chemo_no_only_final2, pf_chemo_no_only_final2, 
                                 mf_chemo_no_only_final2, rm_chemo_no_only_final2, 
                                 ra_chemo_no_only_final2, grid = T, 
                                 axis.lim= c(0.1, 10),  m.labels = c("General fatigue", 
                                                                     "Physical fatigue", 
                                                                     "Mental fatigue", 
                                                                     "Reduced motivation", 
                                                                     "Reduced activity"),
                                 axis.labels = c( "IMRT", "Reduced activity at baseline", 
                                                  "Hormone cancer therapy",
                                                  "Reduced motivation at baseline", "RT boost", 
                                                  "Global Health Status (10 unit increment)",
                                                  "Mental fatigue at baseline", "Breast symptoms (10 unit increment)", 
                                                  "High alcohol consumption",  "BMI (continuous, 1 unit increment",
                                                  "Physical fatigue at baseline", 
                                                  "Dyspnea (10 unit increment)", 
                                                  "Insomnia (10 unit increment)", "Pain (10 unit increment)",
                                                  "Depression", "Age (continuous, 5 unit increment) ", 
                                                  "General fatigue at baseline"),
                                 vline.color = "snow3", show.values = T, value.size = 3, 
                                 show.p = T, p.shape = F)







