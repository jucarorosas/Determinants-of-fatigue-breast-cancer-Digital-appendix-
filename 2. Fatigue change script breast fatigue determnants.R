
##############################################################################
##########                                                      ##############
########## Determinants of fatigue and longitudinal changes up  ##############
############      to 2 years post-radiotherapy among       ###################
###############     breast cancer patients in REQUITE    #####################
#################                                       ######################
##################         Fatigue change R script    ########################
####################        Date: 01.08.2022       ###########################     
##############################################################################
library(lcsm)
library(lcmm)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(lme4)
library(afex)
library(lattice) #xyplot
library(data.table) #setnames
library(reshape2) #melt
library(ggplot2)
library(texreg)
library(plyr)


## Create a data.frame in long format with the IDs, visits, and scores for
## the five dimensions. 

MFI_datx3 <- MFI_datx

MFI_datx3 %>%
  dplyr::select(Subject.Id, Visit, MFI_general, 
                MFI_physic, MFI_mental, MFI_redact, MFI_redmot)-> MFI_long_growth


MFI_long_growth %>% 
  filter(MFI_long_growth$Visit == "Baseline" | MFI_long_growth$Visit == "Post-RT" |
           MFI_long_growth$Visit == "12m follow-up" |
           MFI_long_growth$Visit == "24m follow-up" ) -> MFI_long_growth_2

MFI_drop <- droplevels(MFI_long_growth_2)

## Select variables from SDdat (Mainsite)
SDdat %>%
  dplyr::select(Subject.Id, Mainsite) -> SDdat_for_drop

## Left join SDdat_for_drop to MFI_drop

MFI_drop %>%
  left_join(SDdat_for_drop, by="Subject.Id") -> MFI_dropy

str(MFI_dropy)

MFI_drop %>%
  left_join(SDdat_for_drop, by="Subject.Id") -> MFI_dropx

str(MFI_dropx)

## Left join BDdat_chemo to MFI_dropx (to include chemo variable)

MFI_dropx %>%
  left_join(BDdat_chemo, by="Subject.Id") -> MFI_dropx

MFI_dropy %>%
  left_join(BDdat_chemo, by="Subject.Id") -> MFI_dropy


## Convert Mainsite, chemo and Subject.Id to factor
MFI_dropx$Mainsite <- as.factor(MFI_dropx$Mainsite)
MFI_dropx$Subject.Id <- as.factor(MFI_dropx$Subject.Id)
MFI_dropx$chemo_bin <- as.factor(MFI_dropx$chemo_bin)

MFI_dropy$Mainsite <- as.factor(MFI_dropy$Mainsite)
MFI_dropy$Subject.Id <- as.factor(MFI_dropy$Subject.Id)
MFI_dropy$chemo_bin <- as.factor(MFI_dropy$chemo_bin)

## Null model 
m_null<-lmer(MFI_general~1+(1|Subject.Id),MFI_dropx,REML=F)


## a. Unconditional Growth Model: General fatigue

 # Re-order Visit so that Post-RT is the comparison dummy variable
   MFI_dropy$Visit<-relevel(MFI_dropy$Visit, "Post-RT")
   
   MFI_dropz <- MFI_dropy
   MFI_dropz$Visit<-relevel(MFI_dropz$Visit, "12m follow-up")
   
## a.1 Unconditional Growth Model
m_un_growthx<-lmer(MFI_general ~ Visit + (1 | Subject.Id), data=MFI_dropx)
summary(m_un_growthx)

# a.2 Unconditional Growth Model (Re ordered levels for Visit) 
m_un_growthy<-lmer(MFI_general ~ Visit + (1 | Subject.Id), data=MFI_dropy)
summary(m_un_growthy)

## b. Unconditional Growth Model: Physical fatigue
m_un_growthx_phy<-lmer(MFI_physic ~ Visit + (1 | Subject.Id), data=MFI_dropx)
summary(m_un_growthx_phy)

## c. Unconditional Growth Model: Mental fatigue
m_un_growthx_men<-lmer(MFI_mental ~ Visit + (1 | Subject.Id), data=MFI_dropx)
summary(m_un_growthx_men)

## d. Unconditional Growth Model: Reduced activity
m_un_growthx_redact<-lmer(MFI_redact ~ Visit + (1 | Subject.Id), data=MFI_dropx)
summary(m_un_growthx_redact)

## e. Unconditional Growth Model: Reduced motivation 
m_un_growthx_redmot<-lmer(MFI_redmot ~ Visit + (1 | Subject.Id), data=MFI_dropx)
summary(m_un_growthx_redmot)

                            
                    #############################
             #######                            ########################
################## Mixed models by chemo groups ##############################
             ######                             ########################
                   ##############################

## First create two separate data.frames by chemo status (do it twice 
## because we re-ordered so that RT is the comparison in the "dropy" dataframes)

MFI_dropx %>% 
  filter(MFI_dropx$chemo_bin == "Yes") -> MFI_dropx_chemo_yes

MFI_dropx %>% 
  filter(MFI_dropx$chemo_bin == "No") -> MFI_dropx_chemo_no

MFI_dropy %>% 
  filter(MFI_dropy$chemo_bin == "Yes") -> MFI_dropy_chemo_yes

MFI_dropy %>% 
  filter(MFI_dropy$chemo_bin == "No") -> MFI_dropy_chemo_no

MFI_dropz %>% 
  filter(MFI_dropy$chemo_bin == "Yes") -> MFI_dropz_chemo_yes

MFI_dropz %>% 
  filter(MFI_dropy$chemo_bin == "No") -> MFI_dropz_chemo_no

# GENERAL FATIGUE 

## a.1 Unconditional Growth Model: MFI general, Chemo yes, comparison=baseline
growthx_chemo_yes<-lmer(MFI_general ~ Visit + (1 | Subject.Id), data=MFI_dropx_chemo_yes)
summary(growthx_chemo_yes)

# a.2 Unconditional Growth Model: MFI general, Chemo yes, comparison=Post_RT
growthy_chemo_yes<-lmer(MFI_general ~ Visit + (1 | Subject.Id), data=MFI_dropy_chemo_yes)
summary(growthy_chemo_yes)

## a.3 Unconditional Growth Model: MFI general, Chemo no, comparison=baseline
growthx_chemo_no<-lmer(MFI_general ~ Visit + (1 | Subject.Id), data=MFI_dropx_chemo_no)
summary(growthx_chemo_no)

# a.4 Unconditional Growth Model: MFI general, Chemo no, comparison=Post_RT
growthy_chemo_no<-lmer(MFI_general ~ Visit + (1 | Subject.Id), data=MFI_dropy_chemo_no)
summary(growthy_chemo_no)


######## # ######### # 


# PHYSICAL FATIGUE 
## a.1 Unconditional Growth Model: MFI general, Chemo yes, comparison=baseline
growthx_chemo_yes_pf<-lmer(MFI_physic ~ Visit + (1 | Subject.Id), data=MFI_dropx_chemo_yes)
summary(growthx_chemo_yes_pf)

# a.2 Unconditional Growth Model: MFI general, Chemo yes, comparison=Post_RT
growthy_chemo_yes_pf<-lmer(MFI_physic ~ Visit + (1 | Subject.Id), data=MFI_dropy_chemo_yes)
summary(growthy_chemo_yes_pf)

## a.3 Unconditional Growth Model: MFI general, Chemo no, comparison=baseline
growthx_chemo_no_pf<-lmer(MFI_physic ~ Visit + (1 | Subject.Id), data=MFI_dropx_chemo_no)
summary(growthx_chemo_no_pf)

# a.4 Unconditional Growth Model: MFI general, Chemo no, comparison=Post_RT
growthy_chemo_no_pf<-lmer(MFI_physic ~ Visit + (1 | Subject.Id), data=MFI_dropy_chemo_no)
summary(growthy_chemo_no_pf)


# MENTAL FATIGUE 
## a.1 Unconditional Growth Model: MFI mental, Chemo yes, comparison=baseline
growthx_chemo_yes_mf<-lmer(MFI_mental ~ Visit + (1 | Subject.Id), data=MFI_dropx_chemo_yes)
summary(growthx_chemo_yes_mf)

# a.2 Unconditional Growth Model: MFI mental, Chemo yes, comparison=Post_RT
growthy_chemo_yes_mf<-lmer(MFI_mental ~ Visit + (1 | Subject.Id), data=MFI_dropy_chemo_yes)
summary(growthy_chemo_yes_mf)

## a.3 Unconditional Growth Model: MFI mental, Chemo no, comparison=baseline
growthx_chemo_no_mf<-lmer(MFI_mental ~ Visit + (1 | Subject.Id), data=MFI_dropx_chemo_no)
summary(growthx_chemo_no_mf)

# a.4 Unconditional Growth Model: MFI mental, Chemo no, comparison=Post_RT
growthy_chemo_no_mf<-lmer(MFI_mental ~ Visit + (1 | Subject.Id), data=MFI_dropy_chemo_no)
summary(growthy_chemo_no_mf)


# REDUCED MOTIVATION
## a.1 Unconditional Growth Model: MFI redmot, Chemo yes, comparison=baseline
growthx_chemo_yes_rm<-lmer(MFI_redmot ~ Visit + (1 | Subject.Id), data=MFI_dropx_chemo_yes)
summary(growthx_chemo_yes_rm)

# a.2 Unconditional Growth Model: MFI redmot, Chemo yes, comparison=Post_RT
growthy_chemo_yes_rm<-lmer(MFI_redmot ~ Visit + (1 | Subject.Id), data=MFI_dropy_chemo_yes)
summary(growthy_chemo_yes_rm)

## a.3 Unconditional Growth Model: MFI redmot, Chemo no, comparison=baseline
growthx_chemo_no_rm<-lmer(MFI_redmot ~ Visit + (1 | Subject.Id), data=MFI_dropx_chemo_no)
summary(growthx_chemo_no_rm)

# a.4 Unconditional Growth Model: MFI redmot, Chemo no, comparison=Post_RT
growthy_chemo_no_rm<-lmer(MFI_redmot ~ Visit + (1 | Subject.Id), data=MFI_dropy_chemo_no)
summary(growthy_chemo_no_rm)



# REDUCED ACTIVITY
## a.1 Unconditional Growth Model: MFI redact, Chemo yes, comparison=baseline
growthx_chemo_yes_ra<-lmer(MFI_redact ~ Visit + (1 | Subject.Id), data=MFI_dropx_chemo_yes)
summary(growthx_chemo_yes_ra)

# a.2 Unconditional Growth Model: MFI redact, Chemo yes, comparison=Post_RT
growthy_chemo_yes_ra<-lmer(MFI_redact ~ Visit + (1 | Subject.Id), data=MFI_dropy_chemo_yes)
summary(growthy_chemo_yes_ra)

## a.3 Unconditional Growth Model: MFI redact, Chemo no, comparison=baseline
growthx_chemo_no_ra<-lmer(MFI_redact ~ Visit + (1 | Subject.Id), data=MFI_dropx_chemo_no)
summary(growthx_chemo_no_ra)

# a.4 Unconditional Growth Model: MFI redact, Chemo no, comparison=Post_RT
growthy_chemo_no_ra<-lmer(MFI_redact ~ Visit + (1 | Subject.Id), data=MFI_dropy_chemo_no)
summary(growthy_chemo_no_ra)


## Diagnostics 

plot(m_un_growthx)

qqnorm(resid(m_un_growthx))

MFI_dropx %>%
  dplyr::select(Subject.Id, Visit, MFI_general) -> MFI_dropx_gen

MFI_dropx_gen_complete <-   MFI_dropx_gen[complete.cases(MFI_dropx_gen), ]

MFI_dropx_gen_complete$residuals_gen <- residuals(m_un_growthx)


residuals(m_un_growthx)

#plotting RESIDUAL intraindividual change
ggplot(data = MFI_dropx_gen_complete, aes(x = Visit, y = residuals_gen, group = Subject.Id)) +
  ggtitle("Fixed Linear, Random Intercept") +
  #  geom_point() + 
  geom_line() +
  xlab("Visit") + 
  ylab("RESIDUAL WISC MFI_gen score") + ylim(0,20) #Note the removal of limits on y-axis


Plot.Model.F.Linearity<-plot(resid(m_un_growthx),MFI_dropx_gen_complete$MFI_general)

# Testing homogeneity of variance 
#MFI_dropx_gen_complete$Model.F.Res<- residuals(Model.F) #extracts the residuals and places them in a new column in our original data table
MFI_dropx_gen_complete$abs.growtx.F.Res <-abs(MFI_dropx_gen_complete$residuals_gen) #creates a new column with the absolute value of the residuals
MFI_dropx_gen_complete$Model.growthx.Res2 <- MFI_dropx_gen_complete$abs.growtx.F.Res^2 #squares the absolute values of the residuals to provide the more robust estimate
Levene.Model.F <- lm(Model.growthx.Res2 ~ Subject.Id, data=MFI_dropx_gen_complete) #ANOVA of the squared residuals
anova(Levene.Model.F) #displays the results

Plot.Model.F <- plot(m_un_growthx) #creates a fitted vs residual plot
Plot.Model.F



###############  WORSENING OF FATIGUE LEVELS #########################

# Reshape to wide format
MFI_gen_wide <-  reshape(MFI_datx_gen, idvar = "Subject.Id", 
                         timevar = "Visit", direction = "wide")

# Select the relvant variables
MFI_gen_wide %>%
  dplyr::select(Subject.Id, MFI_general.Baseline, 
                `MFI_general.24m follow-up`) -> MFI_gen_wide2
# Omit NAs
MFI_gen_wide3 <- na.omit(MFI_gen_wide2)

# Create new "difference" fatigue variable
MFI_gen_wide3$diff_gen <- MFI_gen_wide3$MFI_general.Baseline -
  MFI_gen_wide3$`MFI_general.24m follow-up`

# Count total number of rows 
nrow(MFI_gen_wide3)

# Count the number of rows with < 0
nrow(MFI_gen_wide3[MFI_gen_wide3$diff_gen < 0, ]) # 390 (39%)

# Count the number of rows with <= 2
nrow(MFI_gen_wide3[MFI_gen_wide3$diff_gen <= -2, ]) # 306 (31%)


summary(MFI_gen_wide3)

hist(MFI_gen_wide3$MFI_general.Baseline, breaks = 17, col = "lightblue", border = "salmon")
hist(MFI_gen_wide$`MFI_general.Post-RT`, breaks = 17, col = "lightblue", border = "salmon")
hist(MFI_gen_wide$`MFI_general.12m follow-up`, breaks = 17, col = "lightblue", border = "salmon")
hist(MFI_gen_wide3$`MFI_general.24m follow-up`, breaks = 17, col = "lightblue", border = "salmon")




########################  GROWTH MIXTURE MODELING WITH LCMM ####################

# 5.1 Package lcmm

## 5.1.1 LCGA with lcmm

# Models with 1-3 classes

MFI_dropy$Subject.Id2 <- as.numeric(MFI_dropy$Subject.Id)

MFI_dropy %>%
  dplyr::select(Subject.Id, Subject.Id2, Visit, MFI_general)-> MFI_dropy_2

MFI_dropy_f <-  MFI_dropy_2[complete.cases(MFI_dropy_2), ]

MFI_dropy_f$Visit_n <- as.numeric(MFI_dropy_f$Visit)

# Re-order Visit so that Baseline is the comparison dummy variable
MFI_dropy_f$Visit<-relevel(MFI_dropy_f$Visit, "Baseline")

# Create again the numeric visit variable (so that baseline is 1)
MFI_dropy_f$Visit_n <- as.numeric(MFI_dropy_f$Visit)


# Set seed so that the random procedures are reproducible 
set.seed(153)

lcga1 <- hlme(MFI_general ~ Visit, subject = "Subject.Id2", ng = 1, data = MFI_dropy_f) 

lcga2 <- gridsearch(rep = 50, maxiter = 10, minit = lcga1, 
                    hlme(MFI_general ~ Visit, subject = "Subject.Id2", ng = 2, data = MFI_dropy_f, mixture = ~ Visit)) 


lcga3 <- gridsearch(rep = 50, maxiter = 10, minit = lcga1, 
                    hlme(MFI_general ~ Visit, subject = "Subject.Id2", ng = 3, data = MFI_dropy_f, mixture = ~ Visit))

lcga4 <- gridsearch(rep = 50, maxiter = 10, minit = lcga1, 
                    hlme(MFI_general ~ Visit, subject = "Subject.Id2", ng = 4, data = MFI_dropy_f, mixture = ~ Visit))

lcga5 <- gridsearch(rep = 50, maxiter = 10, minit = lcga1, 
                    hlme(MFI_general ~ Visit, subject = "Subject.Id2", ng = 5, data = MFI_dropy_f, mixture = ~ Visit))


# make table with results for the 3 models: 
summarytable(lcga1, lcga2, lcga3)

summary(lcga1)
summary(lcga2)
summary(lcga3)

# 5.1.2 GMM-1 with lcmm

set.seed(153)
gmm1 <- hlme(MFI_general ~ Visit, subject = "Subject.Id2", random=~1, ng = 1, data = MFI_dropy_f)

gmm2 <- gridsearch(rep = 50, maxiter = 10, minit = gmm1, 
                   hlme(MFI_general ~ Visit, subject = "Subject.Id2", random=~1,
                        ng = 2, data = MFI_dropy_f, mixture = ~ Visit, nwg=T))

gmm3 <- gridsearch(rep = 50, maxiter = 10, minit = gmm1,
                   hlme(MFI_general ~ Visit, subject = "Subject.Id2", random=~1,
                        ng = 3, data = MFI_dropy_f, mixture = ~ Visit, nwg=T))


gmm4 <- gridsearch(rep = 50, maxiter = 10, minit = gmm1,
                   hlme(MFI_general ~ Visit, subject = "Subject.Id2", random=~1,
                        ng = 4, data = MFI_dropy_f, mixture = ~ Visit, nwg=T))


gmm5 <- gridsearch(rep = 50, maxiter = 10, minit = gmm1,
                   hlme(MFI_general ~ Visit, subject = "Subject.Id2", random=~1,
                        ng = 5, data = MFI_dropy_f, mixture = ~ Visit, nwg=T))


# make table with results for the 3 models: 
summarytable(gmm1, gmm2, gmm3)
summary(gmm5)

# 5.1.3 GMM-2 with lcmm

set.seed(153)
gmmr1 <- hlme(MFI_general ~ Visit, subject = "Subject.Id2", random=~1 + Visit, ng = 1, data = MFI_dropy_f)


gmmr2  <- gridsearch(rep = 50, maxiter = 10, minit = gmmr1, 
                        hlme(MFI_general ~ Visit, subject = "Subject.Id2", random=~1 + Visit,
                             ng = 2, data = MFI_dropy_f, mixture = ~ Visit, nwg=T))
 

gmmr3  <- gridsearch(rep = 50, maxiter = 10, minit = gmmr1,
                   hlme(MFI_general ~ Visit, subject = "Subject.Id2", random=~1 + Visit,
                        ng = 3, data = MFI_dropy_f, mixture = ~ Visit, nwg=T))

gmmr4  <- gridsearch(rep = 50, maxiter = 10, minit = gmmr1,
                     hlme(MFI_general ~ Visit, subject = "Subject.Id2", random=~1 + Visit,
                          ng = 4, data = MFI_dropy_f, mixture = ~ Visit, nwg=T))


gmmr5  <- gridsearch(rep = 50, maxiter = 10, minit = gmmr1,
                     hlme(MFI_general ~ Visit, subject = "Subject.Id2", random=~1 + Visit,
                          ng = 5, data = MFI_dropy_f, mixture = ~ Visit, nwg=T))

# make table with results for the 4 models: 
summarytable(lcga1, lcga2, lcga3, lcga4, lcga5, 
             gmm1, gmm2, gmm3, gmm4, gmm5,
             gmmr1, gmmr2, gmmr3, gmmr4, gmmr5,  
             which = c("G", "loglik", "conv", "npm", 
                       "AIC", "BIC", "SABIC", "entropy", "%class"))

# Summary plots 
  #1. lcga 
summaryplot(lcga1, lcga2, lcga3, lcga4, lcga5, which = c("AIC", "SABIC", "entropy"))

  #2. gmm 
summaryplot(gmm1, gmm2, gmm3, gmm4, gmm5, which = c("AIC", "SABIC", "entropy"))

  #3. gmmr 
summaryplot(gmmr1, gmmr2, gmmr3, gmmr4, gmmr5, which = c("AIC", "SABIC", "entropy"))

summary(lcga4)

postprob(gmm5)

postprob(lcga4)

str(MFI_dropy_f)

plot(lcga4, which="fit",var.time = "Visit_n", marg=T, shades=T,
     legend.loc="right")

plot(lcga4, which="fit",var.time = "Visit_n", marg=F, break.times=10, shades=T,
     legend.loc="right", legend=NULL)

plot(gmm5, which="fit",var.time = "Visit_n", marg=F, break.times=10, shades=T,
     legend.loc="right", legend=NULL)

plot(lcga4, var.time = "Visit_n", which="postprob")

# Computing and ploting predictions

newdatos<-data.frame(Time=seq(1,4,length=1494))


plot(predictY(gmmr5,newdatos,var.time = "Time"),legend.loc="right",bty="l")


################################################################################
    #       FREQUENCY TABLES BY LATENT CLASSES (using lcga4)        #


#Within the large hlme object "lcga4", we have the data.frame "lcga4$pprob" with 
#the assigned classes for every subject id. Let's assign that object to a new 
#data.frame

table_4classes <- lcga4$pprob

#Select only the column with IDs and the one with classes 

table_4classes <- table_4classes %>%
  select(Subject.Id2, class)

#Left join table_4classes to MFI_dropy_f. We do this because this these two 
#tables share the Subject.Id2 variable (necessary intermediate step)

MFI_dropy_f %>%
  left_join(table_4classes, by="Subject.Id2") -> MFI_dropy_4classes

table(MFI_dropy_4classes$class)

#filter and leave only "baseline" measurements 

MFI_dropy_4classes %>% 
  filter(MFI_dropy_4classes$Visit == "Baseline") -> MFI_dropy_4classes

#select only the ID columns and the class column 

MFI_dropy_4classes %>%
  dplyr::select(Subject.Id, Subject.Id2, class) -> MFI_dropy_4classes


#We already a data.frame (MFI_baseline) with only the IDs for those with 
#baseline MFI-20 measurements. We are going to select from this table
#only the Subject.Id and Visit. 

MFI_baseline %>%
  dplyr::select(Subject.Id, Visit) -> Baseline_ids #table with 1475 subjects 

#left-join "MFI_dropy_4classes" to "Baseline_ids. We do this to have in one
#data.frame all the baseline subjects with their respective latent class. 

Baseline_ids %>%
  left_join(MFI_dropy_4classes, by="Subject.Id") -> Baseline_classes 

#leave complete cases only 

Baseline_classes_1 <- na.omit(Baseline_classes) #table with 1415 subjects 


###LEFT-JOIN VARIABLES TO TABLE WITH CLASSES###

#1. chemo_bin 

#Baseline_classes_1 %>%
#  left_join(BDdat_chemo, by="Subject.Id") -> Baseline_classes_2

#table(Baseline_classes_2$chemo_bin, Baseline_classes_2$class)

#2.  
Baseline_classes_1 %>%
  left_join(BDdat_small, by="Subject.Id") -> Baseline_classes_2


#3.  
Baseline_classes_2 %>%
  left_join(SDdat_small, by="Subject.Id") -> Baseline_classes_3

#4.  
Baseline_classes_3 %>%
  left_join(GPAQ_small, by="Subject.Id") -> Baseline_classes_4

#5.  
Baseline_classes_4 %>%
  left_join(C30dat_scores_small, by="Subject.Id") -> Baseline_classes_5

#6.  
Baseline_classes_5 %>%
  left_join(BR23_final, by="Subject.Id") -> Baseline_classes_6

#7.  
Baseline_classes_6 %>%
  left_join(B3BC_small, by="Subject.Id") -> Baseline_classes_7

#8.  
Baseline_classes_7 %>%
  left_join(Rad_small2, by="Subject.Id") -> Base_class8

#Covert variable "class" to factor 
Base_class8$class <- as.factor(Base_class8$class)

levels(Base_class8$class)
Base_class8$class <- revalue(Base_class8$class, c("1"="Constant high", "2"="Constant low", 
                                                  "3"="Constant moderate", "4"="Decreasing"))

table(Base_class8$class)

# Function to compute the p-values (t-test for continuous and chi-square for categorical variables)
# Taken from: https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html

pval <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- oneway.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}


#Table 1: Descriptives and p-values by latent class 

table1::label(Base_class8$age_at_radiotherapy_start_yrs) <- "Mean age"
table1::label(Base_class8$age_cat) <- "Age (categorical)"
table1::label(Base_class8$bmi) <- "Mean BMI"
table1::label(Base_class8$bmi_cat) <- "BMI (categorical)"
table1::label(Base_class8$household_bin) <- "Living status"
table1::label(Base_class8$education_lmh) <- "Educational level"
table1::label(Base_class8$diabetes) <- "Diabetes"
table1::label(Base_class8$hypertension) <- "Hypertension"
table1::label(Base_class8$history_of_heart_disease) <- "Heart disease"
table1::label(Base_class8$ra) <- "Rheumatoid arthritis"
table1::label(Base_class8$depression) <- "Depression"
table1::label(Base_class8$breast_cancer_family_history_1st_degree) <- "Breast cancer 1st degree"
table1::label(Base_class8$smoker_bin) <- "Smoking status"
table1::label(Base_class8$alcohol_bin) <- "Current alcohol consumption"
table1::label(Base_class8$GPAQ_q16_sitting_recline_min) <- "Time spent sitting (minutes per day)"
table1::label(Base_class8$QL) <- "Global Health Status / QoL (QLQ-C30)"
table1::label(Base_class8$PA) <- "Pain (QLQ-C30)"
table1::label(Base_class8$SL) <- "Insomnia (QLQ-C30)"
table1::label(Base_class8$DY) <- "Dyspnea (QLQ-C30)"
table1::label(Base_class8$BRAS) <- "Arm symptoms-BR23"
table1::label(Base_class8$BRBS) <- "Breast symptoms-BR23"
table1::label(Base_class8$histology_di) <- "Histologic type"
table1::label(Base_class8$b3tumour_side_of_primary) <- "Side of primary tumor"
table1::label(Base_class8$b3surgery_type) <- "Surgery type"
table1::label(Base_class8$b3sys_treatment) <- "Other systemic treatment"
table1::label(Base_class8$b3sys_antiher2) <- "Anti-HER2 therapie"
table1::label(Base_class8$chemo_bin) <- "Chemotherapie"
table1::label(Base_class8$b3radio_breast_dose_Gy) <- "Radiotherapie dose"
table1::label(Base_class8$radio_total_dose_sum) <- "Radiotherapie total dose (including boosts)"
table1::label(Base_class8$b3radio_boost) <- "Radiotherapie boost"
table1::label(Base_class8$fractionation_bin) <- "Fractionation scheme"
table1::label(Base_class8$boost_frac) <- "Fractionation and boost scheme"
table1::label(Base_class8$b3radio_imrt) <- "IMRT"
table1::label(Base_class8$pato_t_cat) <- "Pathologic T-stage"
table1::label(Base_class8$pato_n_cat) <- "Pathologic N-status"

table1::table1(~age_at_radiotherapy_start_yrs + age_cat + bmi + 
                 bmi_cat + household_bin +
                 education_lmh + diabetes + hypertension +
                 history_of_heart_disease + ra + depression + 
                 breast_cancer_family_history_1st_degree + 
                 smoker_bin + alcohol_bin + 
                 GPAQ_q16_sitting_recline_min + QL + PA + SL + DY + 
                 BRAS + BRBS + histology_di + b3tumour_side_of_primary + 
                 b3surgery_type + b3sys_treatment + b3sys_antiher2 + chemo_bin + 
                 b3radio_breast_dose_Gy + radio_total_dose_sum + b3radio_boost +
                 fractionation_bin + boost_frac + b3radio_imrt + 
                 pato_t_cat + pato_n_cat | class, data = Base_class8, 
                 overall=F, extra.col=list("P-value"= pval))



################################################################################
#       MULTINOMIAL LOGISTIC REGRESSION (using lcga4 classes)        #

# Tutorial by: https://stats.oarc.ucla.edu/r/dae/multinomial-logistic-regression/#:~:text=Multinomial%20logistic%20regression%20is%20used,combination%20of%20the%20predictor%20variables.

install.packages("nnet")
library(nnet)
library(gtsummary)

# Choose the reference level 
Base_class8$class <- relevel(Base_class8$class, ref = "Constant low")

# Create new continuous variable with 10-unit intervals
Base_class8$QL_10 <- Base_class8$QL / 10
Base_class8$PA_10 <- Base_class8$PA / 10
Base_class8$SL_10 <- Base_class8$SL / 10
Base_class8$DY_10 <- Base_class8$DY / 10
Base_class8$BRAS_10 <- Base_class8$BRAS / 10
Base_class8$BRBS_10 <- Base_class8$BRBS / 10

# Same for age
Base_class8$age_5 <- Base_class8$age_at_radiotherapy_start_yrs / 5

# Run the multinomial model with the previously selected variables 
# from the frequency tables
multinom_mod1 <- multinom(class ~ age_5 + bmi_cat + 
                            smoker_bin + QL_10 + PA_10 + SL_10 + DY_10 + 
                            BRAS_10 + BRBS_10 + chemo_bin + b3radio_boost + 
                            pato_t_cat, data = Base_class8)

multinom_mod1 <- multinom(class ~ age_5 + bmi_cat + 
                            history_of_heart_disease + depression +
                            smoker_bin + alcohol_bin + 
                            QL_10 + PA_10 + SL_10 + DY_10 + 
                            BRAS_10 + BRBS_10 + 
                            b3surgery_type + b3sys_treatment + chemo_bin +
                            fractionation_bin + b3radio_boost + b3radio_imrt + 
                            pato_t_cat, data = Base_class8)

summary(multinom_mod1)

multinom_mod1%>%
  tbl_regression(exponentiate = T)


glance(multinom_mod1)

z <- summary(multinom_mod1)$coefficients/summary(multinom_mod1)$standard.errors
z

# 2-tailed z test
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p

# extract coefficients and exponentiate
exp(coef(multinom_mod1))





## Code from  https://gist.github.com/ddsjoberg/a55afa74ac58e1f895862fcabab62406
## for nice results table for multinomial logistic regression 

set.seed(100)
library(gtsummary)
library(magrittr)

multinom_pivot_wider <- function(x) {
  # check inputs match expectatations
  if (!inherits(x, "tbl_regression") || !inherits(x$model_obj, "multinom")) {
    stop("`x=` must be class 'tbl_regression' summary of a `nnet::multinom()` model.")
  }
  
  # create tibble of results
  df <- tibble::tibble(outcome_level = unique(x$table_body$groupname_col))
  df$tbl <- 
    purrr::map(
      df$outcome_level,
      function(lvl) {
        gtsummary::modify_table_body(
          x, 
          ~dplyr::filter(.x, .data$groupname_col %in% lvl) %>%
            dplyr::ungroup() %>%
            dplyr::select(-.data$groupname_col)
        )
      }
    )
  
  tbl_merge(df$tbl, tab_spanner = paste0("**", df$outcome_level, "**"))
}


# multinom model tabulated with gtsummary
tbl <-
  nnet::multinom(class ~ age_5 + bmi_cat + 
                   depression + hypertension +
                   smoker_bin +
                   QL_10 + PA_10 + SL_10 + DY_10 + 
                   BRAS_10 + BRBS_10 + 
                   b3sys_treatment + histology_di + chemo_bin +
                   b3radio_boost + b3radio_imrt +
                   pato_t_cat, data = Base_class8) %>%
  tbl_regression(exponentiate = TRUE) %>%
  multinom_pivot_wider()

