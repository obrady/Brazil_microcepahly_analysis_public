
# examine individual birth defects and their association with zika infection
rm(list = ls())
gc()

require(mgcv)
require(lmtest)
require(logistf)
require(parallel)
require(pryr)

# choose suspected or confirmed microcepahly?
case_type <- "confirmed_microcephaly"#"suspected_microcephaly"# "confirmed_microcephaly"

testrun = FALSE

# automatically fills in the type for arbovirus darta
if(case_type == "suspected_microcephaly"){
  zik = "mun_sus_zik_incid_PREG"
  micro = "suspected_microcephaly"
}else{
  zik = "mun_con_zik_incid_PREG"
  micro = "confirmed_microcephaly"
}


# load in dataframe
setwd("/home/admin/OBRADY/Brazil_data")
setwd("/Users/eideobra/Desktop/Brazil_data/BEST")
if(case_type == "suspected_microcephaly"){
  load("INDIVIDUAL/NOV18_issue/Detailed_birth_defects_processed_suspected_SumEXP.RData")
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_SumEXP.RData")
}else{
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_SumEXP.RData")
  load("INDIVIDUAL/NOV18_issue/Detailed_birth_defects_processed_confirmed_SumEXP.RData")
}


# data subsetter for testruns
if(testrun){
  # for testing: sliming down dataset to get accurate times
  prop_micro <- sum(model_df[, 1]) / nrow(model_df)
  # aiming for a df of 100k
  new_micro <- round(100000 * prop_micro, 0)
  micro_df = model_df[model_df[, 1] == 1, ]
  normal_df = model_df[model_df[, 1] == 0, ]
  micro_sam = sample(1:nrow(micro_df), new_micro)
  normal_sam = sample(1:nrow(normal_df), 100000 - new_micro)
  micro_df = micro_df[micro_sam, ]
  normal_df = normal_df[normal_sam, ]
  model_df = rbind(normal_df, micro_df)
  DBD = DBD[c(micro_sam, normal_sam), ]
}

######################
# Part 1: fitting the models
######################
HT.model.wrapper <- function(f_formula, mdf2, modrtn = FALSE){
  if(modrtn){
    return(logistf(f_formula, data = mdf2, dataout = FALSE))
  }else{
    return(extractAIC(logistf(f_formula, data = mdf2, dataout = FALSE))[2])
  }
}

# hypothesised baseline risk factors
base_risk <- c("mother_age", "baby_sex","week", "Region", "SDI")

# different outcomes to test
colnames(DBD) = sapply(colnames(DBD), function(x) paste0("Q", x))
outcome_opts <- colnames(DBD)

# add a merged birth defect as misclassified
outcome_opts = c(outcome_opts, "Q690Q699")
DBD = data.frame(DBD)
DBD$Q690Q699 = as.numeric((DBD$Q699 + DBD$Q690) > 0)

#pate together outcomes and model df then check nothing irrelevant in model_df
model_df = data.frame(model_df, DBD)
mdf = model_df[, c(zik, base_risk, outcome_opts)]
rm(model_df)
rm(DBD)

# build model formulae to test
b_modform_options = c(as.list(data.frame(combn(base_risk, 5), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 4), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 3), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 2), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 1), stringsAsFactors = F)))

# now formulae with different outcomes
cand_mod_form <- list()
for(i in 1:length(outcome_opts)){
  no_covs = as.formula(paste(outcome_opts[i], "~", zik))
  cand_mod_form = c(cand_mod_form, c(no_covs,
                                     lapply(b_modform_options, 
                                            function(x, outcome, exposures){as.formula(paste(outcome, "~", paste(c(exposures, x), collapse = "+")))},
                                            outcome = outcome_opts[i],
                                            exposures = zik)))
}

cl <- makeCluster(8, type = "FORK")
clusterEvalQ(cl, library(logistf))

# big different in processing with 1 risk factor and processing with all 5, so best to process 
# in 5 chunks with equal length
form_length <- rbind(c(1, 1),
                     cbind(unname(unlist(lapply(b_modform_options, length))) + 1,
                           2:(length(b_modform_options) + 1)))
varlength = length(b_modform_options) + 1

cand_mod_form_0R = list()
cand_mod_form_1R = list()
cand_mod_form_2R = list()
cand_mod_form_3R = list()
cand_mod_form_4R = list()
cand_mod_form_5R = list()


for(i in 1:length(outcome_opts)){
  cand_mod_form_0R = c(cand_mod_form_0R, cand_mod_form[((i - 1) * varlength) + form_length[form_length[, 1] == 1, 2]])
  cand_mod_form_1R = c(cand_mod_form_1R, cand_mod_form[((i - 1) * varlength) + form_length[form_length[, 1] == 2, 2]])
  cand_mod_form_2R = c(cand_mod_form_2R, cand_mod_form[((i - 1) * varlength) + form_length[form_length[, 1] == 3, 2]])
  cand_mod_form_3R = c(cand_mod_form_3R, cand_mod_form[((i - 1) * varlength) + form_length[form_length[, 1] == 4, 2]])
  cand_mod_form_4R = c(cand_mod_form_4R, cand_mod_form[((i - 1) * varlength) + form_length[form_length[, 1] == 5, 2]])
  cand_mod_form_5R = c(cand_mod_form_5R, cand_mod_form[((i - 1) * varlength) + form_length[form_length[, 1] == 6, 2]])
}


# now run the models
modAICs_0R = parLapply(cl, cand_mod_form_0R, HT.model.wrapper, mdf2 = mdf)
modAICs_1R = parLapply(cl, cand_mod_form_1R, HT.model.wrapper, mdf2 = mdf)
modAICs_2R = parLapply(cl, cand_mod_form_2R, HT.model.wrapper, mdf2 = mdf)
modAICs_3R = parLapply(cl, cand_mod_form_3R, HT.model.wrapper, mdf2 = mdf)
modAICs_4R = parLapply(cl, cand_mod_form_4R, HT.model.wrapper, mdf2 = mdf)
modAICs_5R = parLapply(cl, cand_mod_form_5R, HT.model.wrapper, mdf2 = mdf)

# runtime estimation
#system.time({mod <- HT.model.wrapper(cand_mod_form_0R[[1]], mdf2 = mdf)})
# with 100k dataset:
# no covs = 1.914, length = 8
# 1 cov = 7.714, length = 40
# 2 cov = 11.125, length = 80
# 3 cov = 15.049 , length = 80
# 4 cov = 4.000 , length = 40
# 5 cov = 6.416 , length = 8

# 4 cores probably optimal (with 8 core machine)
# for 100k rows would take:
# 2629.12 / 4 cores = 11 mins, actual runtime with 2 cores = 14 mins
# for full dataset on 4 cores:
# 6.18 hours

# recompiling AICs and selecting the best model
master_list = list() # is n(outcomes) long
for(i in 1:length(outcome_opts)){
  col_opts_0R = length(modAICs_0R) / length(outcome_opts)
  col_opts_1R = length(modAICs_1R) / length(outcome_opts)
  col_opts_2R = length(modAICs_2R) / length(outcome_opts)
  col_opts_3R = length(modAICs_3R) / length(outcome_opts)
  col_opts_4R = length(modAICs_4R) / length(outcome_opts)
  col_opts_5R = length(modAICs_5R) / length(outcome_opts)
  
  # temporary list for within outcome
  tmplist <- list()
  tmplist[[1]] = modAICs_0R[((i - 1) * col_opts_0R + 1):(i * col_opts_0R)]
  tmplist[[2]] = modAICs_1R[((i - 1) * col_opts_1R + 1):(i * col_opts_1R)]
  tmplist[[3]] = modAICs_2R[((i - 1) * col_opts_2R + 1):(i * col_opts_2R)]
  tmplist[[4]] = modAICs_3R[((i - 1) * col_opts_3R + 1):(i * col_opts_3R)]
  tmplist[[5]] = modAICs_4R[((i - 1) * col_opts_4R + 1):(i * col_opts_4R)]
  tmplist[[6]] = modAICs_5R[((i - 1) * col_opts_5R + 1):(i * col_opts_5R)]
  # assign to the masterlist
  master_list[[i]] = unname(unlist(tmplist))
}


# identify best model
exp_eval <- list()
iden_RFs <- list()
selected_RFs <- list()
for(i in 1:length(outcome_opts)){
  # AICs of all models
  exp_eval[[i]] = master_list[[i]]
  # selected optimal formula
  selected_RFs[[i]] = cand_mod_form[(32 * (i - 1) + 1): (32 * i)][which.min(exp_eval[[i]])][[1]]
  # identified final risk factors
  iden_RFs[[i]] = c("", b_modform_options)[which.min(exp_eval[[i]])][[1]]
}

# Refit final models now we know which risk factors shoudl be included
final_mods = parLapply(cl, selected_RFs, HT.model.wrapper, mdf2 = mdf, modrtn = TRUE)

stopCluster(cl)

# now save models

for(i in 1:length(final_mods)){
  mod <- list(iden_RFs[[i]],
              final_mods[[i]])
  
  save(mod, file = paste0("MODELS/", case_type, "/Model_", outcome_opts[i], "_",
                          case_type, ".Rdata"))
}



















######################
# Part 2: building the output table
######################



# FUNCTIONS

predict.logisf <- function(coefs, newdata){
  cnames <- rownames(coefs)
  # trim factors
  whichfac <- colnames(newdata)[sapply(newdata[1, ], is.factor)]
  if(length(whichfac) > 0){
    for(i in 1:length(whichfac)){cnames = gsub(whichfac[i], "", cnames)}
    rownames(coefs) = cnames
    # add refernce valeus for factor covariates
    misvars <- vector()
    for(i in 1:length(whichfac)){
      if(length(newdata[, whichfac[i]]) > 1){
        misvars[i] = as.character(newdata[, whichfac[i]][!unique(newdata[, whichfac[i]]) %in% rownames(coefs)])
      }else{
        misvars[i] = as.character(newdata[, whichfac[i]])
      }
    }
    misvars = misvars[!is.na(misvars)]
    tmat <- matrix(1, ncol = 3, nrow = length(misvars))
    rownames(tmat) = misvars
    colnames(tmat) = colnames(coefs)
    coefs = rbind(coefs, tmat)
  }
  # prediction
  pred_df <- data.frame(pred = rep(NA, nrow(newdata)),
                        pred_lower = rep(NA, nrow(newdata)),
                        pred_upper = rep(NA, nrow(newdata)))
  for(i in 1:nrow(newdata)){
    # start with the intercept
    val <- as.numeric(coefs[1, ])
    # then add or subtract other values
    for(j in 1:ncol(newdata)){
      if(is.numeric(newdata[i, j])){
        val = val + as.numeric(coefs[match(names(newdata)[j], rownames(coefs)), ]) * newdata[i, j]
      }else{
        val = val + as.numeric(coefs[match(newdata[i, j], rownames(coefs)), ])
      }
    }
    # exp to reverse log transform
    #pred_df[i, ] = exp(val)
    pred_df[i, ] = exp(val)/(1+exp(val))
  }
  return(pred_df)
}

RR.pred3 <- function(model, exposure, covariates){
  # 1- refit with new alpha to reflect multiple hypothesis testing
  model = logistf(model$formula, BD_mat, alpha = 0.0042) # 12 hypothesis tests
  coefs1 <- data.frame(coefficients(model), model$ci.lower, model$ci.upper)
  # if there are interactions, add an interation column to the original data
  if(any(grepl(":", exposure))){
    ints <- exposure[grep(":", exposure)]
    ints2 <- strsplit(ints, ":")[[1]]
    BD_mat$x = BD_mat[, ints2[1]] * BD_mat[, ints2[2]]
    colnames(BD_mat)[ncol(BD_mat)] = gsub(":", "X", ints)
    exposure[grep(":", exposure)] = gsub(":", "X", exposure[grep(":", exposure)])
    rownames(coefs1) = gsub(":", "X", rownames(coefs1))
  }
  
  # zerovalue of zero inflated log transform
  if(length(exposure) > 1){
    zeroval = apply(BD_mat[, exposure], 2, min)
  }else{zeroval <- min(BD_mat[, exposure])}
  
  
  # variable list
  if(all(is.na(covariates))){
    varlist <- exposure
  }else{varlist= c(exposure, covariates)}
  
  # editing for sex as a factor
  if(any(grepl("sexM", varlist))){
    varlist[grepl("sexM", varlist)] = "baby_sex"
  }
  
  # editing for Region as a factor
  if(any(grepl("Region", varlist))){
    varlist[grepl("Region", varlist)] = NA
    varlist = c(varlist, "Region")
  }
  varlist = varlist[!is.na(varlist)]
  
  newdat1 = data.frame(NA)
  for(i in 1:length(varlist)){
    if(varlist[i] %in% exposure){
      newdat1 = cbind(newdat1, zeroval[which.max(exposure %in% varlist[i])])
    }else{
      if(is.factor(BD_mat[, varlist[i]])){
        tt <- table(BD_mat[, varlist[i]])
        newdat1 = cbind(newdat1, as.character(names(tt[which.max(tt)])))
      }else{
        newdat1 = cbind(newdat1, median(BD_mat[, varlist[i]], na.rm = T))
      }
    }
  }
  # remove the NA
  
  newdat1 = newdat1[, 2:ncol(newdat1)]
  if(length(varlist) == 1){names(newdat1) = varlist}else{colnames(newdat1) = varlist}
  
  # data with exposure
  newdat2 = newdat1
  if(length(varlist) == 1){
    newdat2 = median(BD_mat[BD_mat[, exposure] > zeroval, exposure])
    names(newdat2) = exposure
  }else{
    if(length(exposure) == 1){
      newdat2[, exposure] = median(BD_mat[BD_mat[, exposure] > zeroval, exposure], na.rm = T)
    }else{
      for(i in 1:length(exposure)){
        newdat2[, exposure[i]] = median(BD_mat[BD_mat[, exposure[i]] > zeroval[i], exposure[i]], na.rm = T)
      }
    }
  }
  
  
  # model predictions
  p0zik <- predict.logisf(coefs1, data.frame(newdat1))
  pMzik <- predict.logisf(coefs1, data.frame(newdat2))
  RR = (pMzik) / p0zik
  return(unlist(c(RR, model$prob[2])))
}

# reload data
if(case_type == "suspected_microcephaly"){
  load("INDIVIDUAL/NOV18_issue/Detailed_birth_defects_processed_suspected_SumEXP.RData")
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_SumEXP.RData")
}else{
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_SumEXP.RData")
  load("INDIVIDUAL/NOV18_issue/Detailed_birth_defects_processed_confirmed_SumEXP.RData")
}
colnames(DBD) = sapply(colnames(DBD), function(x) paste0("Q", x))

BD_mat = data.frame(zik = model_df[, zik], 
                    mother_age = model_df$mother_age,
                    Region = model_df$Region,
                    week = model_df$week,
                    baby_sex = model_df$baby_sex, # just needed for no covs. RR pred
                    SDI = model_df$SDI,
                    DBD)
colnames(BD_mat)[1] = zik
rm(model_df)

# colnames change
colnames(BD_mat)[grepl("X", colnames(BD_mat))] = gsub("X", "Q", colnames(BD_mat)[grepl("X", colnames(BD_mat))])

# sample sizes for testing each BD - now all pre-processed
# BDs start in coulumn 8 and finish in coulmn 147 = 140 total
#BD_SS = colSums(BD_mat[, 8:147])
#sum(BD_SS > 5) # 96
#sum(BD_SS > 10) # 81
# incidence greater than 0.5 in 10,000
#sum(BD_SS > (5 * nrow(BD_mat) / 100000)) # 13 = pcrit = 0.003846154

# main list we will investigate
#BD_SS_focus = names(BD_SS[ BD_SS > (5 * nrow(BD_mat) / 100000)])
# remove Q000 - NA code
#BD_SS_focus = BD_SS_focus[BD_SS_focus != "Q000"]
confounders <- c("mother_age", "week", "Region", "SDI", "baby_sex")
BD_SS_focus = colnames(BD_mat)[7:14]
# will take about 2h each - 16h total

#for(i in 1:length(BD_SS_focus)){
#  mod <- step.custom(response = BD_SS_focus[i],
#                     covariates = c(zik, confounders),
#                     selectable = c(FALSE, rep(TRUE, length(confounders))),
#                     glmdata = BD_mat)
#  save(mod, file = paste("WORKSPACE/MODELS/Model_", 
#                         BD_SS_focus[i], 
#                         "_", 
#                         case_type,
#                         ".RData", 
#                         sep = ""))
#  rm(mod)
#  print(i)
#}

# special- combine extra fingers with polydactyly, unspecified
BD_mat$Q690Q699 = as.numeric((BD_mat$Q690 + BD_mat$Q699) > 0)

#mod <- step.custom(response = "Q690Q699",
#                   covariates = c("mun_con_zika_incid_TRI1", confounders),
#                   selectable = c(FALSE, rep(TRUE, length(confounders))),
#                   glmdata = BD_mat)
#save(mod, file = paste("WORKSPACE/MODELS/Model_", "Q690Q699_", case_type, ".RData", sep = ""))

BD_SS_focus = c(BD_SS_focus, "Q690Q699")



# assemble the data frame with the final results


resdf <- data.frame(ICD10 = BD_SS_focus,
                    Znum = rep(NA, length(BD_SS_focus)),
                    Zrate = rep(NA, length(BD_SS_focus)),
                    NZnum = rep(NA, length(BD_SS_focus)),
                    NZrate = rep(NA, length(BD_SS_focus)),
                    CrudeRR = rep(NA, length(BD_SS_focus)),
                    AdjRR_mid = rep(NA, length(BD_SS_focus)),
                    AdjRR_low = rep(NA, length(BD_SS_focus)),
                    AdjRR_high = rep(NA, length(BD_SS_focus)),
                    Pval = rep(NA, length(BD_SS_focus)))

for(i in 1:length(BD_SS_focus)){
  
  # split birth outcome data into Zika and no Zika
  nozikval = min(BD_mat[, zik])
  Nozik = BD_mat[BD_mat[, zik] == nozikval, BD_SS_focus[i]]
  zik2 = BD_mat[BD_mat[, zik] > nozikval, BD_SS_focus[i]]
  
  # asign crude values
  resdf$Znum[i] = sum(zik2)
  resdf$Zrate[i] = 100000 * sum(zik2) / length(zik2)
  
  resdf$NZnum[i] = sum(Nozik)
  resdf$NZrate[i] = 100000 * sum(Nozik) / length(Nozik)
  
  resdf$CrudeRR[i] = resdf$Zrate[i] / resdf$NZrate[i]
  
  # adusted values
  load(paste("WORKSPACE/MODELS/", case_type, "/Model_", BD_SS_focus[i], "_", case_type, ".RData", sep = ""))
  vals = RR.pred3(mod[[2]], zik, 
                  mod[[1]][!grepl("zik", mod[[1]])])
  
  resdf$AdjRR_mid[i] = as.numeric(vals[1])
  resdf$AdjRR_low[i] = as.numeric(vals[2])
  resdf$AdjRR_high[i] = as.numeric(vals[3])
  
  resdf$Pval[i] = as.numeric(vals[4])
  
  # update me
  print(i)
}

# significat ones: Q174 (1.58)- misplaced ear,  and Q059 (1.25) - Spina bifida
write.csv(resdf, file = paste0("FIGURES/Diff_outcomes_tab_INDIV_",
                              case_type,
                              ".csv"))










#### EYE investigation
source("REFERENCE/BD_ICD10s.R")

BD_mat = BD_mat[, c(rep(TRUE, 7), colnames(BD_mat)[8:147] %in% EYE)]

# split into Zika and no Zika
BD_mat_NZ = BD_mat[BD_mat$mun_con_zika_incid_TRI1 == min(BD_mat$mun_con_zika_incid_TRI1), ]
BD_mat_Z = BD_mat[BD_mat$mun_con_zika_incid_TRI1 > min(BD_mat$mun_con_zika_incid_TRI1), ]

resdf <- data.frame(ICD10 = names(apply(BD_mat_Z[, 8:ncol(BD_mat_Z)], 2, sum) / nrow(BD_mat_Z)),
                    Rate_Z = round(10000 * apply(BD_mat_Z[, 8:ncol(BD_mat_Z)], 2, sum) / nrow(BD_mat_Z), 3),
                    Rate_NZ = round(10000 * apply(BD_mat_NZ[, 8:ncol(BD_mat_NZ)], 2, sum) / nrow(BD_mat_NZ), 3))
resdf$diff = resdf$Rate_Z - resdf$Rate_NZ
resdf[order(resdf$diff), ]

# major ones (descending):
# Q112 microphthalmos, 
# Q103 Other congenital malformations of eyelid, 
# Q120 Congenital cataract







