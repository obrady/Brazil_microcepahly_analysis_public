# variable examination and processing

rm(list = ls())

setwd("/Users/eideobra/Desktop/Brazil_data/BEST")
load("INDIVIDUAL/Indiv_reg_df_NOV18_170718.RData")
load("INDIVIDUAL/FULL_TS/Zika_sus_FES_NOV18.RData")
load("INDIVIDUAL/Detailed_birth_defects_NOV18.RData")

# add region
sta_rs <- read.csv("STA_COVS/STA_regions.csv")
model_df$Region <- sta_rs[match(model_df$case_sta, sta_rs$State), "Region"]

# processing functions
ZIT <- function(x){
  log(x + 0.5 * min(x[x > 0], na.rm = T))
}

data.process <- function(dat, FES, DBD , case_conf = "sus"){
  
  # 01 add a Socio Development index
  sociocomps <- dat[, c("mun_Income_pc",
                        "mun_Mean_years_edu",
                        "mun_Fert_rate")]
  # standardize NAs
  sociocomps[is.na(rowSums(sociocomps)), ] = NA
  SDInaind <- is.na(sociocomps[, 1])
  #Scale
  sociocomps[!SDInaind, ] = apply(sociocomps[!SDInaind, ], 2, function(x) (x - min(x)) / (max(x) - min(x)))
  SDI = rowSums(sociocomps) / 3
  dat$SDI = SDI
  
  # 02 adjust birth date covariates
  dat$weeks_gestation[dat$weeks_gestation > 45] = 40
  dat$weeks_gestation[dat$weeks_gestation < 25] = 40
  dat$weeks_gestation[is.na(dat$weeks_gestation)] = 40
  
  # 04 Zero inflation transformation
  
  dat$mun_sus_zik_incid_PREG = ZIT(dat$mun_sus_zik_incid_PREG)
  dat$mun_con_zik_incid_PREG = ZIT(dat$mun_con_zik_incid_PREG)
  
  dat$mun_sus_den_incid_PREG = ZIT(dat$mun_sus_den_incid_PREG)
  dat$mun_con_den_incid_PREG = ZIT(dat$mun_con_den_incid_PREG)
  dat$mun_sus_den_incid_PRE12 = ZIT(dat$mun_sus_den_incid_PRE12)
  dat$mun_con_den_incid_PRE12 = ZIT(dat$mun_con_den_incid_PRE12)
  
  dat$mun_sus_chik_incid_PREG = ZIT(dat$mun_sus_chik_incid_PREG)
  dat$mun_con_chik_incid_PREG = ZIT(dat$mun_con_chik_incid_PREG)
  dat$mun_sus_chik_incid_PRE12 = ZIT(dat$mun_sus_chik_incid_PRE12)
  dat$mun_con_chik_incid_PRE12 = ZIT(dat$mun_con_chik_incid_PRE12)
  
  # 05 mother race recategorization - 180768 NAs
  dat$mother_race = as.character(dat$mother_race)
  dat$mother_race[dat$mother_race == "IGNORADO"] = NA
  dat$mother_race = as.factor(dat$mother_race)
  
  # 06 mother age constraints - 263 NAs
  dat$mother_age[dat$mother_age > 55] = NA
  dat$mother_age[dat$mother_age < 10] = NA
  
  # 07 sex of baby - 1312 NAs
  dat$baby_sex = as.character(dat$baby_sex)
  dat$baby_sex[dat$baby_sex == "I"] = NA
  dat$baby_sex = as.factor(dat$baby_sex)
  
  # 08 cattle density
  #hist(dat$mun_cattle_den)
  dat$mun_cattle_den = ZIT(dat$mun_cattle_den)
  
  # 09 farm density
  dat$mun_farm_den = ZIT(dat$mun_farm_den)
  
  # 10 water vulnerability
  dat$mun_water_vul[dat$mun_water_vul > 1] = 1
  dat$mun_water_vul = dat$mun_water_vul ^ (1/3)
  
  # 11 second NA trim (now only NA trim)
  if(case_conf == "sus"){
    # now excluding previous dengue hypothesis
    exposure_covlist <- c("mun_sus_zik_incid_PREG", "mun_sus_den_incid_PREG",
                          "mun_sus_chik_incid_PREG", "YF_Vcov", "mun_farm_den", 
                          "mun_water_vul")
    individual_covlist <- c("mother_age", "baby_sex", "SDI")
  }else{
    exposure_covlist <- c("mun_con_zik_incid_PREG", "mun_con_den_incid_PREG",
                          "mun_con_chik_incid_PREG", "YF_Vcov", "mun_farm_den", 
                          "mun_water_vul")
    individual_covlist <- c("mother_age", "baby_sex", "SDI")
  }
  
  
  # 12 trim model df
  if(case_conf == "sus"){
    summary_dat = dat[, c("suspected_microcephaly",
                          "case_sta",
                          "case_mun",
                          "Region",
                          "notified_birth_defect",
                          "BD_BRAIN",
                          "BD_EYE",
                          "BD_MSK",
                          "BD_HAN",
                          "week",
                          exposure_covlist,
                          individual_covlist)]
    complete_dat = dat[, c("suspected_microcephaly",
                           "case_sta",
                           "case_mun",
                           "Region",
                           "week",
                           "weeks_gestation",
                           individual_covlist)]
  }else{
    summary_dat = dat[, c("confirmed_microcephaly",
                          "case_sta",
                          "case_mun",
                          "Region",
                          "notified_birth_defect",
                          "BD_BRAIN",
                          "BD_EYE",
                          "BD_MSK",
                          "BD_HAN",
                          "week",
                          exposure_covlist,
                          individual_covlist)]
    complete_dat = dat[, c("confirmed_microcephaly",
                           "case_sta",
                           "case_mun",
                           "Region",
                           "week",
                           "weeks_gestation",
                           individual_covlist)]
  }
  
  # 13 - FINAL trim
  NArows_summary <- apply(summary_dat[, c(exposure_covlist, individual_covlist)], 1, function(x) any(is.na(x)))
  NArows_complete <- apply(complete_dat[, individual_covlist], 1, function(x) any(is.na(x)))
  # for methdos figure calculation:
  #missdatSUM <- apply(dat_HT[, c(hypothesis_covlist, "Region")], 2, function(x) sum(is.na(x))
  
  summary_dat = summary_dat[!NArows_summary, ]
  DBD = DBD[!NArows_summary, ]
  
  FES = FES[!NArows_complete, ]
  complete_dat = complete_dat[!NArows_complete, ]

  
  # 14- save outputs
  if(case_conf == "sus"){
    model_df = summary_dat
    save(model_df, file = "INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_SumEXP.RData")
    save(FES, file = "INDIVIDUAL/NOV18_issue/Zika_FES_processed_suspected_FullEXP.RData")
    save(DBD, file = "INDIVIDUAL/NOV18_issue/Detailed_birth_defects_processed_suspected_SumEXP.RData")
    save(complete_dat, file = "INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_FullEXP.RData")
  }else{
    model_df = summary_dat
    save(model_df, file = "INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_SumEXP.RData")
    save(FES, file = "INDIVIDUAL/NOV18_issue/Zika_FES_processed_confirmed_FullEXP.RData")
    save(DBD, file = "INDIVIDUAL/NOV18_issue/Detailed_birth_defects_processed_confirmed_SumEXP.RData")
    save(complete_dat, file = "INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_FullEXP.RData")
  }
}

# suspected cases
data.process(dat = model_df,
             FES = Zika_sus_FES,
             DBD = BD_mat,
             case_conf = "sus")

# confirmed cases
rm(Zika_sus_FES)
load("INDIVIDUAL/FULL_TS/Zika_con_FES_NOV18.RData")
data.process(dat = model_df,
             FES = Zika_con_FES,
             DBD = BD_mat,
             case_conf = "conf")










###################
# special edition for race sensitivity analysis
###################
data.process <- function(dat, FES, DBD , case_conf = "conf"){
  
  # 01 add a Socio Development index
  sociocomps <- dat[, c("mun_Income_pc",
                        "mun_Mean_years_edu",
                        "mun_Fert_rate")]
  # standardize NAs
  sociocomps[is.na(rowSums(sociocomps)), ] = NA
  SDInaind <- is.na(sociocomps[, 1])
  #Scale
  sociocomps[!SDInaind, ] = apply(sociocomps[!SDInaind, ], 2, function(x) (x - min(x)) / (max(x) - min(x)))
  SDI = rowSums(sociocomps) / 3
  dat$SDI = SDI
  
  # 02 adjust birth date covariates
  dat$weeks_gestation[dat$weeks_gestation > 45] = 40
  dat$weeks_gestation[dat$weeks_gestation < 25] = 40
  dat$weeks_gestation[is.na(dat$weeks_gestation)] = 40
  
  
  # 03 check for implausible values for  microcephaly incidence
  # microcephaly max allowed 20%
  if(case_conf == "sus"){
    mcrate <- as.matrix(t(table(dat$suspected_microcephaly, dat$case_mun)))
  }else{
    mcrate <- as.matrix(t(table(dat$confirmed_microcephaly, dat$case_mun)))
  }
  mcrate = cbind(mcrate, mcrate[, 2] / (mcrate[, 1] + mcrate[, 2]))
  # trim
  dat = dat[mcrate[, 3] <= 0.2, ] # removes 3 municipalities with 3 cases total
  DBD = DBD[mcrate[, 3] <= 0.2, ] #  or none in confirmed dataset
  FES = FES[mcrate[, 3] <= 0.2, ]
  
  # 04 Zero inflation transformation
  
  dat$mun_sus_zik_incid_PREG = ZIT(dat$mun_sus_zik_incid_PREG)
  dat$mun_con_zik_incid_PREG = ZIT(dat$mun_con_zik_incid_PREG)
  
  dat$mun_sus_den_incid_PREG = ZIT(dat$mun_sus_den_incid_PREG)
  dat$mun_con_den_incid_PREG = ZIT(dat$mun_con_den_incid_PREG)
  dat$mun_sus_den_incid_PRE12 = ZIT(dat$mun_sus_den_incid_PRE12)
  dat$mun_con_den_incid_PRE12 = ZIT(dat$mun_con_den_incid_PRE12)
  
  dat$mun_sus_chik_incid_PREG = ZIT(dat$mun_sus_chik_incid_PREG)
  dat$mun_con_chik_incid_PREG = ZIT(dat$mun_con_chik_incid_PREG)
  dat$mun_sus_chik_incid_PRE12 = ZIT(dat$mun_sus_chik_incid_PRE12)
  dat$mun_con_chik_incid_PRE12 = ZIT(dat$mun_con_chik_incid_PRE12)
  
  # 05 mother race recategorization - 180768 NAs
  dat$mother_race = as.character(dat$mother_race)
  dat$mother_race[dat$mother_race == "IGNORADO"] = NA
  dat$mother_race = as.factor(dat$mother_race)
  
  # 06 mother age constraints - 263 NAs
  dat$mother_age[dat$mother_age > 55] = NA
  dat$mother_age[dat$mother_age < 10] = NA
  
  # 07 sex of baby - 1312 NAs
  dat$baby_sex = as.character(dat$baby_sex)
  dat$baby_sex[dat$baby_sex == "I"] = NA
  dat$baby_sex = as.factor(dat$baby_sex)
  
  # 08 cattle density
  #hist(dat$mun_cattle_den)
  dat$mun_cattle_den = ZIT(dat$mun_cattle_den)
  
  # 09 farm density
  dat$mun_farm_den = ZIT(dat$mun_farm_den)
  
  # 10 water vulnerability
  dat$mun_water_vul[dat$mun_water_vul > 1] = 1
  dat$mun_water_vul = dat$mun_water_vul ^ (1/3)
  
  # 11 second NA trim (now only NA trim)
  if(case_conf == "sus"){
    hypothesis_covlist <- c("mun_sus_zik_incid_PREG", "mun_sus_den_incid_PREG",
                            "mun_sus_chik_incid_PREG", "mun_sus_den_incid_PRE12",
                            "YF_Vcov", "mother_age", "baby_sex","mun_farm_den", 
                            "mun_water_vul", "SDI", "mother_race")
    BD_covlist <- c("mun_sus_zik_incid_PREG", "mother_age", "baby_sex","SDI", "weeks_gestation")
  }else{
    hypothesis_covlist <- c("mun_con_zik_incid_PREG", "mun_con_den_incid_PREG",
                            "mun_con_chik_incid_PREG", "mun_con_den_incid_PRE12",
                            "YF_Vcov", "mother_age", "mun_farm_den", "mun_water_vul", 
                            "baby_sex", "SDI", "mother_race")
    BD_covlist <- c("mun_con_zik_incid_PREG", "mother_age","baby_sex", "SDI", "weeks_gestation")
  }
  
  
  # 12 trim model df
  if(case_conf == "sus"){
    dat_HT = dat[, c("suspected_microcephaly",
                     "case_sta",
                     "case_mun",
                     "Region",
                     "BD_BRAIN",
                     "BD_EYE",
                     "BD_MSK",
                     "BD_HAN",
                     "week",
                     hypothesis_covlist)]
    dat_BD = dat[, c("suspected_microcephaly",
                     "notified_birth_defect",
                     "case_sta",
                     "case_mun",
                     "Region",
                     "BD_BRAIN",
                     "BD_EYE",
                     "BD_MSK",
                     "BD_HAN",
                     "week",
                     BD_covlist)]
  }else{
    dat_HT = dat[, c("confirmed_microcephaly",
                     "case_sta",
                     "case_mun",
                     "Region",
                     "BD_BRAIN",
                     "BD_EYE",
                     "BD_MSK",
                     "BD_HAN",
                     "week",
                     hypothesis_covlist)]
    dat_BD = dat[, c("confirmed_microcephaly",
                     "notified_birth_defect",
                     "case_sta",
                     "case_mun",
                     "Region",
                     "BD_BRAIN",
                     "BD_EYE",
                     "BD_MSK",
                     "BD_HAN",
                     "week",
                     BD_covlist)]
  }
  
  # 13 - FINAL trim
  NArows_HT <- apply(dat_HT[, hypothesis_covlist], 1, function(x) any(is.na(x)))
  # for methdos figure calculation:
  #missdatSUM <- apply(dat_HT[, c(hypothesis_covlist, "Region")], 2, function(x) sum(is.na(x)))
  NArows_BD <- apply(dat_BD[, BD_covlist], 1, function(x) any(is.na(x)))
  
  dat_HT = dat_HT[!NArows_HT, ]
  dat_BD = dat_BD[!NArows_BD, ]
  FES = FES[!NArows_BD, ]
  DBD = DBD[!NArows_BD, ]
  
  # 14- save outputs
  if(case_conf == "sus"){
    model_df = dat_HT
    save(model_df, file = "INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_RACE.RData")
    #model_df = dat_BD
    #save(model_df, file = "INDIVIDUAL/Indiv_reg_df_JUL18_170718_processed_suspected_ZIKONLY.RData")
    #save(FES, file = "INDIVIDUAL/FULL_TS/Zika_FES_processed_suspected_ZIKONLY.RData")
    #save(DBD, file = "INDIVIDUAL/Detailed_birth_defects_processed_suspected_ZIKONLY.RData")
  }else{
    model_df = dat_HT
    save(model_df, file = "INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_RACE.RData")
    #model_df = dat_BD
    #save(model_df, file = "INDIVIDUAL/Indiv_reg_df_JUL18_170718_processed_confirmed_ZIKONLY.RData")
    #save(FES, file = "INDIVIDUAL/FULL_TS/Zika_FES_processed_confirmed_ZIKONLY.RData")
    #save(DBD, file = "INDIVIDUAL/Detailed_birth_defects_processed_confirmed_ZIKONLY.RData")
  }
}


###################
# data summariser for data sharing
###################
load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_SumEXP.RData")
model_df_sus = model_df
load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_SumEXP.RData")
model_df_con = model_df


require(dplyr)
require(lubridate)

# make a month variable (week too identifiable)
model_df_con$week2 = week(as.Date("2015-01-01") + model_df_con$week * 7)
model_df_con$month = month(as.Date("2015-01-01") + model_df_con$week * 7)
model_df_con$year = year(as.Date("2015-01-01") + model_df_con$week * 7)
model_df_con$yearweek = model_df_con$year * 100 + model_df_con$month

sum_dat <- model_df_con[, c("case_mun",
                            "yearweek",
                            "confirmed_microcephaly",
                            "notified_birth_defect",
                            "BD_BRAIN",
                            "BD_EYE",
                            "BD_MSK",
                            "BD_HAN")]
sum_dat[is.na(sum_dat)] = 0
sum_dat = data.frame(sum_dat)
sum_dat$births = 1

sds <- sum_dat %>%
  group_by(case_mun, yearweek) %>% 
  summarise_all(funs(sum))

# now other environemntal variables

sum_dat2  <- model_df_con[, c("case_mun",
                              "yearweek",
                              "YF_Vcov",
                              "mun_farm_den",
                              "mun_water_vul",
                              "mother_age",
                              "baby_sex",
                              "SDI")]

sum_dat2 = data.frame(sum_dat2)
# un ZIT farm density
sum_dat2$mun_farm_den = exp(sum_dat2$mun_farm_den) - min(exp(sum_dat2$mun_farm_den))
# baby sex => proportion of babies male
sum_dat2$baby_sex = as.numeric(sum_dat2$baby_sex == "M")

sds2 <- sum_dat2 %>%
  group_by(case_mun, yearweek) %>% 
  summarise_all(funs(mean))

# bring back together
pub_dat = data.frame(sds, sds2)

# drop duplicate mun and week columns
pub_dat = pub_dat[, !(colnames(pub_dat) %in% c("case_mun.1", "yearweek.1"))]
# give columns sensible names and add state, region and municipality population size
pub_dat = data.frame(IBGE_municipality = pub_dat$case_mun,
                     IBGE_mun_state = mun_demo$IBGE_STATE_CODE[match(pub_dat$case_mun, mun_demo$IBGE_CODE)],
                     IBGE_mun_region = mun_demo$Region[match(pub_dat$case_mun, mun_demo$IBGE_CODE)],
                     Year_week = pub_dat$yearweek,
                     MWSD = pub_dat$confirmed_microcephaly,
                     Any_birth_defect = pub_dat$notified_birth_defect,
                     Brain_birth_defect = pub_dat$BD_BRAIN,
                     Eye_birth_defect = pub_dat$BD_EYE,
                     Muscoloskeletal_birth_defect = pub_dat$BD_MSK,
                     Head_and_neck_birth_defect = pub_dat$BD_HAN,
                     #total_births = pub_dat$births, # makes the data too identifiable
                     Yellow_fever_vaccine_coverage = pub_dat$YF_Vcov,
                     Farm_density = pub_dat$mun_farm_den,
                     Water_vulnerability = pub_dat$mun_water_vul,
                     Mean_mother_age = pub_dat$mother_age,
                     Proportion_babies_male = pub_dat$baby_sex,
                     SDI = pub_dat$SDI)

# save output
write.csv(pub_dat, file = "INDIVIDUAL/FOR_PUBLIC/MC_data_public_release_V1.csv")


