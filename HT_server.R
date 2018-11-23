# server version of Hypothesis testing code for Brazil microcepahly paper
# parralellises model selection


rm(list = ls())
gc()

require(lmtest)
require(logistf)
require(lubridate)
#require(snowfall)
require(parallel)
require(pryr)

# choose suspected or confirmed microcepahly?
case_type <- "suspected_microcephaly"# "suspected_microcephaly"# "confirmed_microcephaly"
# test run for timings
testrun = TRUE

# automatically fills in the type for arbovirus darta
if(case_type == "suspected_microcephaly"){
  arbo_type <- list(zika_PREG = "mun_sus_zik_incid_PREG",
                    dengue_PREG = "mun_sus_den_incid_PREG",
                    chik_PREG = "mun_sus_chik_incid_PREG")
}else{
  arbo_type <- list(zika_PREG = "mun_con_zik_incid_PREG",
                    dengue_PREG = "mun_con_den_incid_PREG",
                    chik_PREG = "mun_con_chik_incid_PREG")
}



# load in dataframe
setwd("/home/admin/OBRADY/Brazil_data")
setwd("/Users/Admin/Desktop/Brazil_data/BEST")


if(case_type == "suspected_microcephaly"){
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_SumEXP.RData")
}else{
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_SumEXP.RData")
}




###################
# 01 fitting the models
###################

mdf.process<- function(mdf, testrun, interactions = FALSE, integer_vals = FALSE){
  
  # data subsampler (testrun = TRUE only for runtime estimation)
  # nrow(model_df) = 5360851
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
  
  # finer spatial resolution modification
  # states with > 10 MWSD per 10k live births =
  # PB, RN, SE, AM
  # states with > 9 MWSD per 10k live births =
  # BA, PE, RR, PI
  model_df$Region_plus = as.character(model_df$Region)
  model_df$Region_plus[model_df$case_sta == "PB"] = "Northeast_PB"
  model_df$Region_plus[model_df$case_sta == "RN"] = "Northeast_RN"
  model_df$Region_plus[model_df$case_sta == "SE"] = "Northeast_SE"
  model_df$Region_plus[model_df$case_sta == "AM"] = "North_AM"
  #model_df$Region_plus[model_df$case_sta == "BA"] = "Northeast_BA"
  #model_df$Region_plus[model_df$case_sta == "PE"] = "Northeast_PE"
  #model_df$Region_plus[model_df$case_sta == "RR"] = "North_RR"
  #model_df$Region_plus[model_df$case_sta == "PI"] = "Northeast_PI"
  model_df$Region_plus = as.factor(model_df$Region_plus)
  
  # time period modifier for finer temporal resolution
  # Key period is October 2-15 - July 2016 (when national incidence first exceeded and finally
  # descended below 4 MWSD per 10k live births)
  # use simple year effect for the rest of the points
  md_month <- as.Date("2015-01-01") + model_df$week * 7
  md_month2 <- year(md_month) * 100 + month(md_month)
  model_df$TF = year(md_month)
  model_df$TF[md_month2 %in% c(201510:201512, 201601:201607)] = md_month2[md_month2 %in% c(201510:201512, 201601:201607)]
  model_df$TF = as.factor(model_df$TF)
  
  # interactions
  if(interactions){
    # adding interaction columns
    test_ints <- c(paste(unlist(arbo_type)[1:2], collapse = "x"),
                   paste(unlist(arbo_type)[c(1, 3)], collapse = "x"),
                   paste(c(unlist(arbo_type)[1], "YF_Vcov"), collapse = "x"))
    model_df = data.frame(model_df,
                          model_df[, unlist(arbo_type)[1]] * model_df[, unlist(arbo_type)[2]],
                          model_df[, unlist(arbo_type)[1]] * model_df[, unlist(arbo_type)[3]],
                          model_df[, unlist(arbo_type)[1]] * model_df[, "YF_Vcov"])
    colnames(model_df)[(ncol(model_df) - 2):ncol(model_df)] = test_ints
  }
  
  # whether to convert to integer to save memory
  if(integer_vals){
    model_df$confirmed_microcephaly = as.integer(model_df$confirmed_microcephaly)
    model_df$mother_age = as.integer(model_df$mother_age)
    model_df$mun_con_zik_incid_PREG = as.integer(model_df$mun_con_zik_incid_PREG * 10^5)
    model_df$mun_con_den_incid_PREG = as.integer(model_df$mun_con_den_incid_PREG * 10^5)
    model_df$mun_con_chik_incid_PREG = as.integer(model_df$mun_con_chik_incid_PREG * 10^5)
    model_df$YF_Vcov = as.integer(model_df$YF_Vcov * 10^5)
    model_df$mun_farm_den = as.integer(model_df$mun_farm_den * 10^5)
    model_df$SDI = as.integer(model_df$SDI * 10^5)
  }
  
  return(model_df)
}

#model_df = mdf.process(model_df, testrun, integer_vals = TRUE)
model_df = mdf.process(model_df, testrun, integer_vals = FALSE, interactions = TRUE)


# define baseline riskfactors
#base_risk <- c("mother_age", "baby_sex","TF", "Region_plus", "SDI")
base_risk <- c("mother_age", "baby_sex","week", "Region", "SDI")
# tested exposure hypotheses
# if doing interactions:
test_ints <- c(paste(unlist(arbo_type)[1:2], collapse = "x"),
               paste(unlist(arbo_type)[c(1, 3)], collapse = "x"),
               paste(c(unlist(arbo_type)[1], "YF_Vcov"), collapse = "x"))
hyp_expos <- c(unlist(arbo_type), "mun_farm_den", "mun_water_vul", test_ints)
#hyp_expos <- c(unlist(arbo_type), "mun_farm_den", "mun_water_vul")

#check nothing irrelevant in model_df
model_df = model_df[, c(case_type, base_risk, hyp_expos)]


# variable selection options (31 total)
b_modform_options = c(as.list(data.frame(combn(base_risk, 5), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 4), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 3), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 2), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 1), stringsAsFactors = F)))



# wrapper funciton for model running and extracting AIC

#f_formula = cand_mod_form[[100]]
#mdf2 = mdf

HT.model.wrapper <- function(f_formula, mdf2, modrtn = FALSE){
  # 01 remove any unnecessary data
  #mdf2 = mdf2[, as.logical(sapply(colnames(mdf2), function(x) any(grepl(x, as.character(f_formula)))))]
  
  # 02 fit a logistf model
  #f_mod = logistf(f_formula, data = mdf2)
  
  # 03 extract and return model AIC
  if(modrtn){
    return(logistf(f_formula, data = mdf2, dataout = FALSE))
  }else{
    return(extractAIC(logistf(f_formula, data = mdf2, dataout = FALSE))[2])
  }
}



# Begin master workflow
outcome = case_type
exposures = hyp_expos
riskFactors = base_risk
mdf = model_df
rm(model_df)


# 01 compile lists of model formulae
# baseline model options first
no_covs = as.formula(paste(outcome, "~ week")) # - not used, but helps maintain consistent spacing and order in formula list
cand_mod_form <- c(no_covs,
                   lapply(b_modform_options, 
                          function(x, outcome){as.formula(paste(outcome, "~", paste(x, collapse = "+")))},
                          outcome = outcome))
# now formulae with candidate exposures
for(i in 1:length(exposures)){
  no_covs = as.formula(paste(outcome, "~", exposures[i]))
  cand_mod_form = c(cand_mod_form, c(no_covs,
                                     lapply(b_modform_options, 
                                            function(x, outcome, exposures){as.formula(paste(outcome, "~", paste(c(exposures, x), collapse = "+")))},
                                            outcome = outcome,
                                            exposures = exposures[i])))
}

# 02 run models in parralell
# start cluster
#sfInit(cpus = 25, parallel = T)
#sfLibrary(logistf)
cl <- makeCluster(2, type = "FORK")
clusterEvalQ(cl, library(logistf))

# big different in processing with 1 risk factor and processing with all 5, so best to process 
# in 5 chunks with equal length
form_length <- rbind(c(1, 1),
                     cbind(unname(unlist(lapply(b_modform_options, length))) + 1,
                           2:(length(b_modform_options) + 1)))
varlength = length(b_modform_options) + 1

cand_mod_form_0R <- list()
cand_mod_form_1R <- list()
cand_mod_form_2R <- list()
cand_mod_form_3R <- list()
cand_mod_form_4R <- list()
cand_mod_form_5R <- list()

for(i in 1:(length(exposures) + 1)){
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
#system.time(HT.model.wrapper(cand_mod_form_5R[[1]], mdf2 = mdf))
# with 100k dataset
# no covs = 0.730, length = 8
# 1 cov = 1.212, length = 40
# 2 cov = 1.212, length = 80
# 3 cov = 4.714 , length = 80
# 4 cov = 8.101 , length = 40
# 5 cov = 17.822 , length = 8

# probably looking at 4 cores on this 8 core machine.
# for 100k rows with 4 cores takes:
# observed: 859.887s = 14 minutes
# for full dataset with 4 cores would take:
# 12.65 hours , observed time with 8 cores = 4 hours





# recompiling AICs and selecting the best model
master_list = list() # is n(outcomes) long
for(i in 1:(length(exposures) + 1)){
  col_opts_0R = length(modAICs_0R) / (length(exposures)+1)
  col_opts_1R = length(modAICs_1R) / (length(exposures)+1)
  col_opts_2R = length(modAICs_2R) / (length(exposures)+1)
  col_opts_3R = length(modAICs_3R) / (length(exposures)+1)
  col_opts_4R = length(modAICs_4R) / (length(exposures)+1)
  col_opts_5R = length(modAICs_5R) / (length(exposures)+1)
  
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
for(i in 1:(length(exposures) + 1)){
  # AICs of all models
  exp_eval[[i]] = master_list[[i]]
  # selected optimal formula
  selected_RFs[[i]] = cand_mod_form[(32 * (i - 1) + 1): (32 * i)][which.min(exp_eval[[i]])][[1]]
  # identified final risk factors
  iden_RFs[[i]] = c("", b_modform_options)[which.min(exp_eval[[i]])][[1]]
}

# Refit final models now we know which risk factors shoudl be included
final_mods = parLapply(cl, selected_RFs, HT.model.wrapper, mdf2 = mdf, modrtn = TRUE)
Sys.time()
stopCluster(cl)


# 05 rename and save
# save each as a list where first item is identified risk factors and secodn is the model object
mod_names <- c("basemodAdj",
               "zikmodAdj",
               "denmodAdj",
               "chimodAdj",
               "farmmodAdj",
               "watmodAdj",
               "zikdenmodAdj",
               "zikchimodAdj",
               "zikyfmodAdj")

for(i in 1:length(final_mods)){
  Adjmod <- list(iden_RFs[[i]],
                 final_mods[[i]])
  save(Adjmod, file = paste0("MODELS/", case_type, "/", mod_names[i], ".Rdata"))
}