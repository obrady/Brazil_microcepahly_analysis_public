
rm(list = ls())
gc()

require(mgcv)
require(lmtest)
require(logistf)
require(parallel)
require(pryr)

# choose suspected or confirmed microcepahly?
case_type <- "suspected_microcephaly"# "confirmed_microcephaly"
# test run for timings
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
setwd("/Users/eideobra/Desktop/Brazil_data/BEST")
setwd("/home/admin/OBRADY/Brazil_data")
if(case_type == "suspected_microcephaly"){
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_SumEXP.RData")
}else{
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_SumEXP.RData")
}

# data subsetter for testruns
if(testrun){
  # for testing: sliming down dataset to get accurate times
  prop_micro <- sum(model_df[, 1]) / nrow(model_df)
  # aiming for a df of 10k
  new_micro <- round(10000 * prop_micro, 0)
  micro_df = model_df[model_df[, 1] == 1, ]
  normal_df = model_df[model_df[, 1] == 0, ]
  micro_df = micro_df[sample(1:nrow(micro_df), new_micro), ]
  normal_df = normal_df[sample(1:nrow(normal_df), 10000 - new_micro), ]
  model_df = rbind(normal_df, micro_df)
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
outcome_opts <- c("notified_birth_defect", "BD_BRAIN", "BD_EYE", "BD_MSK", "BD_HAN")

#check nothing irrelevant in model_df
mdf = model_df[, c(zik, base_risk, outcome_opts)]
rm(model_df)

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

cl <- makeCluster(10, type = "FORK")
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
Sys.time()
modAICs_0R = parLapply(cl, cand_mod_form_0R, HT.model.wrapper, mdf2 = mdf)
modAICs_1R = parLapply(cl, cand_mod_form_1R, HT.model.wrapper, mdf2 = mdf)
modAICs_2R = parLapply(cl, cand_mod_form_2R, HT.model.wrapper, mdf2 = mdf)
modAICs_3R = parLapply(cl, cand_mod_form_3R, HT.model.wrapper, mdf2 = mdf)
modAICs_4R = parLapply(cl, cand_mod_form_4R, HT.model.wrapper, mdf2 = mdf)
modAICs_5R = parLapply(cl, cand_mod_form_5R, HT.model.wrapper, mdf2 = mdf)
Sys.time()

# runtime estimation
#system.time({mod <- HT.model.wrapper(cand_mod_form_1R[[1]], mdf2 = mdf)})
# with 100k dataset:
# no covs = 2.855, length = 5
# 1 cov = 5.134, length = 25
# 2 cov = 34.739, length = 50
# 3 cov = 66.351 , length = 50
# 4 cov = 75.019, length = 25
# 5 cov = 2.903, length = 5

# 5 cores probably optimal
# for 100k rows with 2 cores it takes:
# observed time = 1285.786 / 21 minutes
# for full dataset with 4 cores
# 9.5 hours

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
mod_names <- c(paste0("NBDMOD_", case_type, ".RData"),
               paste0("BDBRAINMOD_", case_type, ".RData"),
               paste0("BDEYEMOD_", case_type, ".RData"),
               paste0("BDMSKMOD_", case_type, ".RData"),
               paste0("BDHANMOD_", case_type, ".RData"))

for(i in 1:length(final_mods)){
  mod <- list(iden_RFs[[i]],
              final_mods[[i]])
  
  save(mod, file = paste0("MODELS/", case_type, "/", mod_names[i], ".Rdata"))
}













######################
# Part 2: building the output table
######################
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
    # exp to logit transform
    pred_df[i, ] = exp(val)/(1+exp(val))
    #pred_df[i, ] = exp(val)
  }
  return(pred_df)
}

RR.pred3 <- function(model, exposure, covariates){
  # 1- refit with new alpha to reflect multiple hypothesis testing
  model = logistf(model$formula, model_df, alpha = 0.0042) # 12 hypothesis tests
  coefs1 <- data.frame(coefficients(model), model$ci.lower, model$ci.upper)
  # if there are interactions, add an interation column to the original data
  if(any(grepl(":", exposure))){
    ints <- exposure[grep(":", exposure)]
    ints2 <- strsplit(ints, ":")[[1]]
    model_df$x = model_df[, ints2[1]] * model_df[, ints2[2]]
    colnames(model_df)[ncol(model_df)] = gsub(":", "X", ints)
    exposure[grep(":", exposure)] = gsub(":", "X", exposure[grep(":", exposure)])
    rownames(coefs1) = gsub(":", "X", rownames(coefs1))
  }
  
  # zerovalue of zero inflated log transform
  if(length(exposure) > 1){
    zeroval = apply(model_df[, exposure], 2, min)
  }else{zeroval <- min(model_df[, exposure])}
  
  
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
      if(is.factor(model_df[, varlist[i]])){
        tt <- table(model_df[, varlist[i]])
        newdat1 = cbind(newdat1, as.character(names(tt[which.max(tt)])))
      }else{
        newdat1 = cbind(newdat1, median(model_df[, varlist[i]], na.rm = T))
      }
    }
  }
  # remove the NA
  
  newdat1 = newdat1[, 2:ncol(newdat1)]
  if(length(varlist) == 1){names(newdat1) = varlist}else{colnames(newdat1) = varlist}
  
  # data with exposure
  newdat2 = newdat1
  if(length(varlist) == 1){
    newdat2 = median(model_df[model_df[, exposure] > zeroval, exposure])
    names(newdat2) = exposure
  }else{
    if(length(exposure) == 1){
      newdat2[, exposure] = median(model_df[model_df[, exposure] > zeroval, exposure], na.rm = T)
    }else{
      for(i in 1:length(exposure)){
        newdat2[, exposure[i]] = median(model_df[model_df[, exposure[i]] > zeroval[i], exposure[i]], na.rm = T)
      }
    }
  }
  
  
  # model predictions
  p0zik <- predict.logisf(coefs1, data.frame(newdat1))
  pMzik <- predict.logisf(coefs1, data.frame(newdat2))
  RR = (pMzik) / p0zik
  return(unlist(c(RR, model$prob[2])))
}

# reload model df
if(case_type == "suspected_microcephaly"){
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_SumEXP.RData")
}else{
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_SumEXP.RData")
}

# split dataset into some zika exposure and no zika exposure

zeroval2 <- min(round(model_df[, zik], 5))
#zerovalDEN <- min(round(model_df$mun_true_den_incid_PRE6, 5))
noZikcon <- model_df[round(model_df[, zik], 5) == zeroval2, ]
Zikcon <- model_df[round(model_df[, zik], 5) != zeroval2, ]

# exclude Mc from all birth defects and from NTDs
Zikcon$notified_birth_defect[Zikcon[, micro] == 1] = 0

noZikcon$notified_birth_defect[noZikcon[, micro] == 1] = 0


# choose exposure
expo <- Zikcon
noexpo <- noZikcon

# produce table
restab <- data.frame(outcome = c("NBD", "BRAIN", "EYE", "MSK", "HAN"),
                     num_expo = c(sum(expo$notified_birth_defect, na.rm = T),
                                  sum(expo$BD_BRAIN, na.rm = T),
                                  sum(expo$BD_EYE, na.rm = T),
                                  sum(expo$BD_MSK, na.rm = T),
                                  sum(expo$BD_HAN, na.rm = T)),
                     rate_expo = 100000 * c(sum(expo$notified_birth_defect, na.rm = T),
                                   sum(expo$BD_BRAIN, na.rm = T),
                                   sum(expo$BD_EYE, na.rm = T),
                                   sum(expo$BD_MSK, na.rm = T),
                                   sum(expo$BD_HAN, na.rm = T)) / nrow(expo),
                     num_noexpo = c(sum(noexpo$notified_birth_defect, na.rm = T),
                                  sum(noexpo$BD_BRAIN, na.rm = T),
                                  sum(noexpo$BD_EYE, na.rm = T),
                                  sum(noexpo$BD_MSK, na.rm = T),
                                  sum(noexpo$BD_HAN, na.rm = T)),
                     rate_noexpo = 100000 * c(sum(noexpo$notified_birth_defect, na.rm = T),
                                            sum(noexpo$BD_BRAIN, na.rm = T),
                                            sum(noexpo$BD_EYE, na.rm = T),
                                            sum(noexpo$BD_MSK, na.rm = T),
                                            sum(noexpo$BD_HAN, na.rm = T)) / nrow(noexpo))
restab$RR_unadjusted <- restab$rate_expo / restab$rate_noexpo

# 01 simple adjusted RR
confounders <- base_risk


# reload models (if necssary)
load(paste("WORKSPACE/MODELS/", case_type,"/NBDMOD_", case_type, ".RData", sep = ""))
minmod_NBD = mod
load(paste("WORKSPACE/MODELS/", case_type, "/BDBRAINMOD_", case_type, ".RData", sep = ""))
minmod_BDBRAIN = mod
load(paste("WORKSPACE/MODELS/", case_type, "/BDEYEMOD_", case_type, ".RData", sep = ""))
minmod_BDEYE = mod
load(paste("WORKSPACE/MODELS/", case_type, "/BDMSKMOD_", case_type, ".RData", sep = ""))
minmod_BDMSK = mod
load(paste("WORKSPACE/MODELS/", case_type, "/BDHANMOD_", case_type, ".RData", sep = ""))
minmod_BDHAN = mod

restab <- cbind(restab, rbind(RR.pred3(minmod_NBD[[2]], zik, minmod_NBD[[1]][!grepl("zik", minmod_NBD[[1]])]),
                              RR.pred3(minmod_BDBRAIN[[2]], zik, minmod_BDBRAIN[[1]][!grepl("zik", minmod_BDBRAIN[[1]])]),
                              RR.pred3(minmod_BDEYE[[2]], zik, minmod_BDEYE[[1]][!grepl("zik", minmod_BDEYE[[1]])]),
                              RR.pred3(minmod_BDMSK[[2]], zik, minmod_BDMSK[[1]][!grepl("zik", minmod_BDMSK[[1]])]),
                              RR.pred3(minmod_BDHAN[[2]], zik, minmod_BDHAN[[1]][!grepl("zik", minmod_BDHAN[[1]])])))
colnames(restab)[10] = "p_val"
# save file
write.csv(restab, file = paste("FIGURES/Diff_outcomes_tab", case_type, ".csv", sep = ""))


# specific tests
summary(minmod_NBD[[2]]) # p = 0.4352024749
summary(minmod_BDBRAIN[[2]]) # p = 0.0156525995
summary(minmod_BDEYE[[2]]) # p = 6.870039e-02
summary(minmod_BDMSK[[2]]) # p = 8.242330e-01
summary(minmod_BDHAN[[2]]) # p = 6.389602e-04






##### plots over time
plot(aggregate(model_df$BD_BRAIN, list(model_df$week), sum, na.rm = T))
plot(aggregate(model_df$BD_EYE, list(model_df$week), sum, na.rm = T))
plot(aggregate(model_df$BD_MSK, list(model_df$week), sum, na.rm = T))
plot(aggregate(model_df$BD_HAN, list(model_df$week), sum, na.rm = T))

# table weekly cases 
BRAIN = t(table(as.numeric(model_df$BD_BRAIN > 0), model_df$week))
BRAIN = 10000 * BRAIN[, 2] / (BRAIN[, 1] + BRAIN[, 2])

EYE = t(table(as.numeric(model_df$BD_EYE > 0), model_df$week))
EYE = 10000 * EYE[, 2] / (EYE[, 1] + EYE[, 2])

MSK = t(table(as.numeric(model_df$BD_MSK > 0), model_df$week))
MSK = 10000 * MSK[, 2] / (MSK[, 1] + MSK[, 2])

HAN = t(table(as.numeric(model_df$BD_HAN > 0), model_df$week))
HAN = 10000 * HAN[, 2] / (HAN[, 1] + HAN[, 2])


require(lubridate)

pdf("FIGURES/BDs_over_time.pdf", width = 10, height = 5)
par(mar = c(5, 4, 6, 2) + 0.1)
plot(as.Date("2015-01-01") + as.numeric(names(BRAIN[1:70])) * 7, 
     BRAIN[1:70], type = "l", xlab = "Date", ylab = "Cases per 10,000 births",
     ylim = c(0, 40), main = "", col = "red", cex.axis = 1.5, cex.lab = 1.5, lwd = 3)
axis(1, at = round_date(as.Date("2015-01-01") + as.numeric(names(EYE[1:70])) * 7, unit = "month"), 
     labels = F)

lines(as.Date("2015-01-01") + as.numeric(names(EYE[1:70])) * 7, EYE[1:70], col = "dark green", lwd = 3)
lines(as.Date("2015-01-01") + as.numeric(names(MSK[1:70])) * 7, MSK[1:70], col = "orange", lwd = 3)
lines(as.Date("2015-01-01") + as.numeric(names(HAN[1:70])) * 7, HAN[1:70], col = "dark blue", lwd = 3)

legend(16750, 50, legend = c("Brain", "Eye", "Musculoskeletal", "Head and neck"),
       horiz = T, xpd = T, col= c("red", "dark green", "orange", "dark blue"), lty = 1,
       cex = 1.2, lwd = 3)
dev.off()


































