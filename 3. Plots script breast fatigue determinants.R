########## PLOTS MFI-20 DATA ############

library(lcsm)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(gridExtra)
library(ggpubr)




## First, create a data.frame with only the general fatigue dimension scores and
## all follow-up points. 

MFI_datx2 <- MFI_datx

MFI_datx2 %>%
  dplyr::select(Subject.Id, Visit, MFI_general, )-> MFI_datx_gen


MFI_datx_gen %>% 
  filter(MFI_datx_gen$Visit == "Baseline" | MFI_datx_gen$Visit == "Post-RT" |
           MFI_datx_gen$Visit == "12m follow-up" |
           MFI_datx_gen$Visit == "24m follow-up" ) -> MFI_datx_gen2


MFI_datx_gen2$Visit <- droplevels(MFI_datx_gen2$Visit)




# Rename the levels for "Visit"
MFI_datx_gen2$Visit <- factor(MFI_datx_gen2$Visit)
levels(MFI_datx_gen2$Visit) <- c("T1", "T2", "T3", "T4")
levels(MFI_datx_gen2$Visit)

# Reshape "MFI_datx_gen2" to wide format 
MFI_gen_wide <- reshape(MFI_datx_gen2, idvar = "Subject.Id", timevar = "Visit", 
                        direction = "wide")

# Concatenate the different measurements in a list 
x_var_list <- c( "MFI_general.T1", "MFI_general.T2", "MFI_general.T3", "MFI_general.T4")

## 1.1 Plot all the data for General fatigue 
plot_trajectories(data = MFI_gen_wide,
                  id_var = "Subject.Id", 
                  var_list = x_var_list,
                  xlab = "Time", ylab = "General fatigue score",
                  connect_missing = FALSE, 
                  title_n = TRUE) 

## 1.2 Plot a small random sample for General fatigue only
plot_trajectories(data = MFI_gen_wide,
                  id_var = "Subject.Id", 
                  var_list = x_var_list,
                  xlab = "Time", ylab = "General fatigue score",
                  random_sample_frac = 0.018, 
                  connect_missing = FALSE, 
                  title_n = TRUE)


## Same plot but only with participants with all measurements for general 
## fatigue 

MFI_gen_wide_full <- na.omit(MFI_gen_wide)

plot_trajectories(data = MFI_gen_wide_full,
                  id_var = "Subject.Id", 
                  var_list = x_var_list,
                  xlab = "Time", ylab = "General fatigue score",
                  connect_missing = FALSE, 
                  title_n = TRUE)



#####################################################

## CHEMO: YES 

timex=c("T1. Baseline","T2. End of RT","T3. 12 months","T4. 24 months")

General_fatigue_ch <- c(12.6, 12.77, 11.18, 11.1)
Physical_fatigue_ch <- c(12.03, 11.64, 10.23, 9.89)
Mental_fatigue_ch <- c(9.27, 9.46, 9.04, 8.88)
Reduced_motivation_ch <- c(8.97, 9.06, 8.2, 7.96)
Reduced_activity_ch <-  c(10.85, 10.81, 9.17, 8.76)

dat_ch <- data.frame(timex, General_fatigue_ch, Physical_fatigue_ch, Mental_fatigue_ch, 
                    Reduced_motivation_ch, Reduced_activity_ch)


dat_ch_plot <- melt(dat_ch, id=c("timex")) 

dat_ch_plot$Dimensions <- dat_ch_plot$variable

levels(dat_ch_plot$Dimensions)

# Re-label dimension levels 

levels(dat_ch_plot$Dimensions) <- c("General fatigue", "Physical fatigue", 
                                  "Mental fatigue", "Reduced motivation", "Reduced activity")


# plot
plot_MFI_ch <- ggplot(dat_ch_plot, aes(x=timex, y=value, color=Dimensions, group=Dimensions))+
  geom_line(linetype="dotted", size=1)+
  labs(title = "Mean fatigue scores (with chemotherapy)", 
       x = "Time point", 
       y = "MFI-20 score") +
  theme_light() +
  scale_color_manual(values=wes_palette(n=5, name="Darjeeling1"))+
  theme(plot.background = element_rect(fill = colors()[248]))


# plot with whole range of values in the y-axis 
plot_MFI_ch2 <- plot_MFI_ch + scale_y_continuous(breaks = seq(7, 14, 1))


print(plot_MFI_ch2)


# Plot with y axis limitis from 7 to 13, only grid lines and labels for integers
plot_MFI_ch3 <- plot_MFI_ch + scale_y_continuous(breaks=c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
                                                          15, 16, 17, 18, 19, 20),
                                           limits =c(4, 20)) +
  theme(panel.grid.minor = element_blank())+
  geom_line(linetype="dotted")+
  geom_point()
  

print(plot_MFI_ch3)

## CHEMO: NO

timex=c("T1. Baseline","T2. End of RT","T3. 12 months","T4. 24 months")

General_fatigue <- c(10.54, 12.14, 10.8, 10.72)
Physical_fatigue <- c(9.97, 10.86, 9.69, 9.62)
Mental_fatigue <- c(8.28, 8.78, 8.39, 8.34)
Reduced_motivation <- c(8.74, 9.51, 8.49, 8.6)
Reduced_activity <-  c(9.45, 10.41, 8.97, 8.94)

dat_not_ch <- data.frame(timex, General_fatigue, Physical_fatigue, 
                         Mental_fatigue, Reduced_motivation, Reduced_activity)


dat_not_ch_plot <- melt(dat_not_ch, id=c("timex")) 

dat_not_ch_plot$Dimensions <- dat_not_ch_plot$variable

# Re-label dimension levels 

levels(dat_not_ch_plot$Dimensions) <- c("General fatigue", "Physical fatigue", 
                                    "Mental fatigue", "Reduced motivation", "Reduced activity")


# plot
plot_MFI_n <- ggplot(dat_not_ch_plot, aes(x=timex, y=value, color=Dimensions, group=Dimensions))+
  geom_point(size=2)+
  geom_path(size=1)+
  labs(title = "Mean fatigue scores (without chemotherapy)", 
       x = "Time point", 
       y = "MFI-20 score") +
  theme_light() +
  scale_color_manual(values=wes_palette(n=5, name="Darjeeling1"))+
  theme(plot.background = element_rect(fill = colors()[248]))


# plot with whole range of values in the y-axis 
plot_MFI_n2 <- plot_MFI_n + scale_y_continuous(breaks = seq(7, 14, 1))


print(plot_MFI_n2)


# Plot with y axis limitis from 7 to 13, only grid lines and labels for integers
plot_MFI_n3 <- plot_MFI_n + scale_y_continuous(breaks=c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
                                                        15, 16, 17, 18, 19, 20),
                                               limits =c(4, 20)) +
  theme(panel.grid.minor = element_blank())

print(plot_MFI_n3)

##### MUTIPLE PLOTS IN ONE ######

theme_set(theme_pubr())

figure1 <- ggarrange(plot_MFI_n3, plot_MFI_ch3,
                     labels = c("A", "B"), 
                     ncol = 2, nrow = 1, 
                     common.legend = T, legend="bottom")

figure1




############### Longitudinal fatigue by chemotherapy #######################

## GENERAL FATIGUE ##

# Left join chemo variable to MFI general
MFI_datx_gen %>%
  left_join(B3BC_small, by="Subject.Id") -> gen_chemo

# Filter the four relevant time-points
gen_chemo %>% 
  filter(gen_chemo$Visit == "Baseline" | gen_chemo$Visit == "Post-RT" |
           gen_chemo$Visit == "12m follow-up" |
           gen_chemo$Visit == "24m follow-up" ) -> gen_chemo_2

# Drop unused levels 
gen_chemo_2$Visit <- droplevels(gen_chemo_2$Visit)


# Select only the necessary variables
gen_chemo_2 %>%
  dplyr::select(Subject.Id, Visit, MFI_general, chemo_bin) -> gen_chemo_3

gen_chemo_4 <- na.omit(gen_chemo_3)

# Plot longitudinal means by chemo_bin

gg.base.gen <- ggplot(gen_chemo_4, aes(x=Visit, y=MFI_general))+
  stat_summary(aes(group = chemo_bin, color = chemo_bin),
               geom = "line", fun = mean, size = 2)+
  labs(title = "Mean general fatigue scores by chemotherapy status", 
       x = "Time point", 
       y = "General fatigue score") +
       scale_y_continuous(breaks=c(7, 8, 9, 10, 11, 12, 13, 14, 15)) +
       expand_limits(y = c(7, 15)) +
       geom_errorbar(aes(ymin=lower, ymax=upper), width=.3, position=pd)
       theme_light() +
       theme(panel.grid.minor = element_blank())

print(gg.base.gen)

table(gen_chemo_4$Visit, gen_chemo_4$chemo_bin)

## PHYSICAL FATIGUE ##

MFI_datx %>%
  dplyr::select(Subject.Id, Visit, MFI_physic) -> MFI_datx_physic

# Left join chemo variable to MFI redact
MFI_datx_physic %>%
  left_join(B3BC_small, by="Subject.Id") -> physic_chemo

# Filter the four relevant time-points
physic_chemo %>% 
  filter(physic_chemo$Visit == "Baseline" | physic_chemo$Visit == "Post-RT" |
           physic_chemo$Visit == "12m follow-up" |
           physic_chemo$Visit == "24m follow-up" ) -> physic_chemo_2

# Drop unused levels 
physic_chemo_2$Visit <- droplevels(physic_chemo_2$Visit)


# Select only the necessary variables
physic_chemo_2 %>%
  dplyr::select(Subject.Id, Visit, MFI_physic, chemo_bin) -> physic_chemo_3

physic_chemo_4 <- na.omit(physic_chemo_3)


# Plot longitudinal means by chemo_bin
gg.base.physic <- ggplot(physic_chemo_4, aes(x=Visit, y=MFI_physic))+
  stat_summary(aes(group = chemo_bin, color = chemo_bin),
               geom = "line", fun = mean, size = 2)+
  labs(title = "Mean physical fatigue scores by chemotherapy status", 
       x = "Time point", 
       y = "Physical fatigue score") +
  scale_y_continuous(breaks=c(7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  expand_limits(y = c(7, 15)) +
  theme_light() +
  theme(panel.grid.minor = element_blank())

print(gg.base.physic)

table(physic_chemo_4$Visit, physic_chemo_4$chemo_bin)

## MENTAL FATIGUE ##

MFI_datx %>%
  dplyr::select(Subject.Id, Visit, MFI_mental) -> MFI_datx_mental

# Left join chemo variable to MFI redact
MFI_datx_mental %>%
  left_join(B3BC_small, by="Subject.Id") -> mental_chemo

# Filter the four relevant time-points
mental_chemo %>% 
  filter(mental_chemo$Visit == "Baseline" | mental_chemo$Visit == "Post-RT" |
           mental_chemo$Visit == "12m follow-up" |
           mental_chemo$Visit == "24m follow-up" ) -> mental_chemo_2

# Drop unused levels 
mental_chemo_2$Visit <- droplevels(mental_chemo_2$Visit)


# Select only the necessary variables
mental_chemo_2 %>%
  dplyr::select(Subject.Id, Visit, MFI_mental, chemo_bin) -> mental_chemo_3

mental_chemo_4 <- na.omit(mental_chemo_3)


# Plot longitudinal means by chemo_bin
gg.base.mental <- ggplot(mental_chemo_4, aes(x=Visit, y=MFI_mental))+
  stat_summary(aes(group = chemo_bin, color = chemo_bin),
               geom = "line", fun = mean, size = 2)+
  labs(title = "Mean mental fatigue scores by chemotherapy status", 
       x = "Time point", 
       y = "Mental fatigue score") +
  scale_y_continuous(breaks=c(7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  expand_limits(y = c(7, 15)) +
  theme_light() +
  theme(panel.grid.minor = element_blank())

print(gg.base.mental)

table(mental_chemo_4$Visit, mental_chemo_4$chemo_bin)


## REDUCED ACTIVITY ##

MFI_datx %>%
  dplyr::select(Subject.Id, Visit, MFI_redact) -> MFI_datx_redact

# Left join chemo variable to MFI redact
MFI_datx_redact %>%
  left_join(B3BC_small, by="Subject.Id") -> redact_chemo

# Filter the four relevant time-points
redact_chemo %>% 
  filter(redact_chemo$Visit == "Baseline" | redact_chemo$Visit == "Post-RT" |
           redact_chemo$Visit == "12m follow-up" |
           redact_chemo$Visit == "24m follow-up" ) -> redact_chemo_2

# Drop unused levels 
redact_chemo_2$Visit <- droplevels(redact_chemo_2$Visit)


# Select only the necessary variables
redact_chemo_2 %>%
  dplyr::select(Subject.Id, Visit, MFI_redact, chemo_bin) -> redact_chemo_3

redact_chemo_4 <- na.omit(redact_chemo_3)


# Plot longitudinal means by chemo_bin
gg.base.redact <- ggplot(redact_chemo_4, aes(x=Visit, y=MFI_redact))+
  stat_summary(aes(group = chemo_bin, color = chemo_bin),
               geom = "line", fun = mean, size = 2)+
  labs(title = "Mean reduced activity scores by chemotherapy status", 
       x = "Time point", 
       y = "Reduced activity score") +
  scale_y_continuous(breaks=c(7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  expand_limits(y = c(7, 15)) +
  theme_light() +
  theme(panel.grid.minor = element_blank())

print(gg.base.redact)

table(redact_chemo_4$Visit, redact_chemo_4$chemo_bin)



## REDUCED MOTIVATION ##

MFI_datx %>%
  dplyr::select(Subject.Id, Visit, MFI_redmot) -> MFI_datx_redmot

# Left join chemo variable to MFI redact
MFI_datx_redmot %>%
  left_join(B3BC_small, by="Subject.Id") -> redmot_chemo

# Filter the four relevant time-points
redmot_chemo %>% 
  filter(redmot_chemo$Visit == "Baseline" | redmot_chemo$Visit == "Post-RT" |
           redmot_chemo$Visit == "12m follow-up" |
           redmot_chemo$Visit == "24m follow-up" ) -> redmot_chemo_2

# Drop unused levels 
redmot_chemo_2$Visit <- droplevels(redmot_chemo_2$Visit)


# Select only the necessary variables
redmot_chemo_2 %>%
  dplyr::select(Subject.Id, Visit, MFI_redmot, chemo_bin) -> redmot_chemo_3

redmot_chemo_4 <- na.omit(redmot_chemo_3)


# Plot longitudinal means by chemo_bin
gg.base.redmot <- ggplot(redmot_chemo_4, aes(x=Visit, y=MFI_redmot))+
  stat_summary(aes(group = chemo_bin, color = chemo_bin),
               geom = "line", fun = mean, size = 2)+
  labs(title = "Mean reduced motivation scores by chemotherapy status", 
       x = "Time point", 
       y = "Reduced motivation score") +
  scale_y_continuous(breaks=c(7, 8, 9, 10, 11, 12, 13, 14, 15)) +
  expand_limits(y = c(7, 15)) +
  theme_light() +
  theme(panel.grid.minor = element_blank())

print(gg.base.redmot)

table(redmot_chemo_4$Visit, redmot_chemo_4$chemo_bin)

################## CALCULATING MFI-20 OVERALL NUMBERS  ################

# Filter the 4 relevant time-points 
MFI_datx %>%
  dplyr::select(Subject.Id, Visit, MFI_total)-> MFI_all_total

MFI_all_total %>% 
  filter(MFI_all_total$Visit == "Baseline" | MFI_all_total$Visit == "Post-RT" |
           MFI_all_total$Visit == "12m follow-up" |
           MFI_all_total$Visit == "24m follow-up" ) -> MFI_all_total_4visits

# Drop unused "Visit" levels 
MFI_all_total_4visits$Visit <- droplevels(MFI_all_total_4visits$Visit)

levels(MFI_all_total_4visits$Visit)

## Remove all NAs (leave only complete cases)

MFI_all_total_4visits[complete.cases(MFI_all_total_4visits), ] ->MFI_4visits_complete

# See the frequency of completed MFI-20 questionares for every time-point 
table(MFI_4visits_complete$Visit)


##### MUTIPLE PLOTS IN ONE ######

theme_set(theme_pubr())

figure1 <- ggarrange(plot_MFI_n3, plot_MFI_ch3,
                     labels = c("A", "B"), 
                     ncol = 2, nrow = 1, 
                     common.legend = T, legend="bottom")

figure1

