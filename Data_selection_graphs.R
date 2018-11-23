# plots over time and space of raw and selected datasets in Brazil Microcepahly analysis
rm(list = ls())

require(ggplot2)
require(RColorBrewer)

setwd("/Users/eideobra/Desktop/Brazil_data/BEST")

## 01 load data

# load original data
load("INDIVIDUAL/Indiv_reg_df_NOV18_170718.RData")
orig = model_df
rm(model_df)
# add an SDI index to it
sociocomps <- orig[, c("mun_Income_pc",
                        "mun_Mean_years_edu",
                        "mun_Fert_rate")]
# standardize NAs
sociocomps[is.na(rowSums(sociocomps)), ] = NA
SDInaind <- is.na(sociocomps[, 1])
#Scale
sociocomps[!SDInaind, ] = apply(sociocomps[!SDInaind, ], 2, function(x) (x - min(x)) / (max(x) - min(x)))
SDI = rowSums(sociocomps) / 3
orig$SDI = SDI

# and processed data (for the three separate analyses)
# Summary exposure dataset
load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_SumEXP.RData")
sum_dat = model_df
rm(model_df)

# Full exposure dataset
load("INDIVIDUAL/NOV18_issue/FES_compiled_confirmed_microcephaly.RData")
full_dat = dat
rm(dat)




## 02 distill to only relevant columns 
# (space / time / outcome / mother age / sex of baby / SDI)
orig = orig[, c("confirmed_microcephaly", "week", "mother_age",
                "baby_sex", "case_sta", "SDI")]

sum_dat = sum_dat[, c("confirmed_microcephaly", "week", "mother_age",
                      "baby_sex", "case_sta", "SDI")]
full_dat = full_dat[, c("MC", "week", "mother_age",
                        "baby_sex", "case_sta", "SDI")]
colnames(full_dat)[1] = "confirmed_microcephaly"


### 03 plot overt time
# first aggregate data sources
pdat = data.frame(week = 0:124, 
                  orig = rep(NA, 125),
                  sum_dat = rep(NA, 125),
                  full_dat = rep(NA, 125))
for(i in 1:nrow(pdat)){
  pdat$orig[i] = sum(orig$confirmed_microcephaly[orig$week == pdat[i, 1]])
  pdat$sum_dat[i] = sum(sum_dat$confirmed_microcephaly[sum_dat$week == pdat[i, 1]])
  pdat$full_dat[i] = sum(full_dat$confirmed_microcephaly[full_dat$week == pdat[i, 1]])
}

pdat = data.frame(Week = rep(0:124, 3), 
                  Cases = c(pdat$orig,
                            pdat$sum_dat,
                            pdat$full_dat),
                  Dataset = c(rep("Live Births", 125),
                              rep("Summary exposure history", 125),
                              rep("Full exposure history", 125)))
pdat$Dataset <- factor(pdat$Dataset, levels = c("Live Births",
                                                "Summary exposure history",
                                                "Full exposure history"))

timeplot <- ggplot(pdat) + 
  geom_line(aes(x = Week, y = Cases, col = Dataset)) +
  scale_color_manual(values = brewer.pal(3, "Dark2"))
timeplot




### 04 plot by state

pdat = data.frame(state = unique(orig$case_sta), 
                  orig = rep(NA, 28),
                  sum_dat = rep(NA, 28),
                  full_dat = rep(NA, 28))
pdat = pdat[!is.na(pdat[, 1]), ]

for(i in 1:nrow(pdat)){
  pdat$orig[i] = sum(orig$confirmed_microcephaly[as.character(orig$case_sta) == as.character(pdat[i, 1])], na.rm = T)
  pdat$sum_dat[i] = sum(sum_dat$confirmed_microcephaly[sum_dat$case_sta == pdat[i, 1]])
  pdat$full_dat[i] = sum(full_dat$confirmed_microcephaly[full_dat$case_sta == pdat[i, 1]])
}

pdat = data.frame(State = rep(pdat[, 1], 3), 
                  Cases = c(pdat$orig,
                            pdat$sum_dat,
                            pdat$full_dat),
                  Dataset = c(rep("Live Births", 27),
                              rep("Summary exposure history", 27),
                              rep("Full exposure history", 27)))

case_state = aggregate(pdat$Cases[pdat$Dataset == "Live Births"], 
                       list(pdat$State[pdat$Dataset == "Live Births"]), FUN = sum)

pdat$Dataset <- factor(pdat$Dataset, levels = rev(c("Live Births",
                                                "Summary exposure history",
                                                "Full exposure history")))

stateplot <- ggplot(pdat, aes(x = State, y = Cases, col = Dataset)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(limits=case_state[rev(order(case_state[, 2])), 1]) +
  scale_color_manual(values = rev(brewer.pal(3, "Dark2")))
stateplot


### 05 mother age
pdat = data.frame(Age = c(orig$mother_age[((orig$mother_age <= 55) & (orig$mother_age >= 10))],
                          sum_dat$mother_age,
                          full_dat$mother_age),
                  Dataset = c(rep("Live Births", length(orig$mother_age[((orig$mother_age <= 55) & (orig$mother_age >= 10))])),
                              rep("Summary exposure history", nrow(sum_dat)),
                              rep("Full exposure history", nrow(full_dat))))
pdat = pdat[!is.na(pdat[, 1]), ]

pdat$Dataset <- factor(pdat$Dataset, levels = c("Live Births",
                                                "Summary exposure history",
                                                "Full exposure history"))

ageplot = ggplot(pdat, aes(x = Dataset, y = Age, col = Dataset))+
  geom_boxplot() +
  scale_color_manual(values = brewer.pal(3, "Dark2"))
ageplot



### 06 baby sex
pdat = data.frame(Male_proportion = c(sum(orig$baby_sex == "M", na.rm = T) / sum(orig$baby_sex %in% c("M", "F"), na.rm = T),
                                      sum(sum_dat$baby_sex == "M") / sum(sum_dat$baby_sex %in% c("M", "F")),
                                      sum(full_dat$baby_sex == "M") / sum(full_dat$baby_sex %in% c("M", "F"))),
                  Dataset = c("Live Births",
                              "Summary exposure history",
                              "Full exposure history"))

pdat$Dataset <- factor(pdat$Dataset, levels = c("Live Births",
                                                "Summary exposure history",
                                                "Full exposure history"))

sexplot = ggplot(pdat, aes(x = Dataset, y = Male_proportion, col = Dataset))+
  geom_boxplot() + 
  ylim(0.50, 0.52) + 
  ylab("Proportion of babies male") +
  scale_color_manual(values = brewer.pal(3, "Dark2"))
  
sexplot




### 07 SDI
pdat = data.frame(SDI = c(orig$SDI,
                          sum_dat$SDI,
                          full_dat$SDI),
                  Dataset = c(rep("Live Births", nrow(orig)),
                              rep("Summary exposure history", nrow(sum_dat)),
                              rep("Full exposure history", nrow(full_dat))))
pdat = pdat[!is.na(pdat[, 1]), ]

pdat$Dataset <- factor(pdat$Dataset, levels = c("Live Births",
                                                "Summary exposure history",
                                                "Full exposure history"))

SDIplot = ggplot(pdat, aes(x = Dataset, y = SDI, col = Dataset))+
  geom_boxplot() +
  scale_color_manual(values = brewer.pal(3, "Dark2"))
SDIplot


### 08 saving plots

ggsave("FIGURES/Data_diff_timeplot.png", timeplot, width = 15, height = 8, units = "cm")
ggsave("FIGURES/Data_diff_stateplot.png", stateplot, width = 20, height = 10, units = "cm")
ggsave("FIGURES/Data_diff_ageplot.png", ageplot, width = 20, height = 10, units = "cm")
ggsave("FIGURES/Data_diff_sexplot.png", sexplot, width = 20, height = 10, units = "cm")
ggsave("FIGURES/Data_diff_SDIplot.png", SDIplot, width = 20, height = 10, units = "cm")










