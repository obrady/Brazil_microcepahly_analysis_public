# reloads hypothesis testing models and extracts relative risk predictions

rm(list = ls())

case_type = "confirmed_microcephaly" # "confirmed_microcephaly" "suspected_microcephaly"

setwd("/Users/eideobra/Desktop/Brazil_data/BEST")

# functions
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
        if(names(newdata)[j] %in% rownames(coefs)){
          val = val + as.numeric(coefs[match(names(newdata)[j], rownames(coefs)), ]) * newdata[i, j]
        }
      }else{
        if(names(newdata)[j] %in% rownames(coefs)){
          val = val + as.numeric(coefs[match(newdata[i, j], rownames(coefs)), ])
        }
      }
    }
    # exp to reverse logit transform
    pred_df[i, ] = exp(val)/(1+exp(val))
    #pred_df[i, ] = exp(val)
  }
  return(pred_df)
}

RR.pred2B <- function(model, exposure_name, zeroval, medval){
  # extract model coefficients
  coefs1 <- cbind(coefficients(model), model$ci.lower, model$ci.upper)
  # create new data frame with median other characteristics
  # (other characteristics don't affect results)
  newdat1 <- data.frame(expo = zeroval,
                        mother_age = 26,
                        baby_sex = "M",
                        week = 82,
                        Region = "SouthEast",
                        SDI = 0.280103)
  
  newdat2 <- newdat1
  newdat2$expo = medval
  
  colnames(newdat1)[1] <- colnames(newdat2)[1] <-  exposure_name
  
  return(predict.logisf(coefs1, newdat2) / predict.logisf(coefs1, newdat1))
}

RR.pred2B.interaction <- function(model, exposure_names, zeroval, medval, zeroval2, medval2){
  # extract model coefficients
  coefs1 <- cbind(coefficients(model), model$ci.lower, model$ci.upper)
  # create new data frame with median other characteristics
  # (other characteristics don't affect results)
  newdat1 <- data.frame(expo = zeroval * zeroval2,
                        mother_age = 26,
                        baby_sex = "M",
                        week = 82,
                        Region = "SouthEast",
                        SDI = 0.280103)
  
  newdat2 <- newdat1
  newdat2$expo = medval * medval2
  
  if(grepl(":", names(coefs1[2, 1]))){
    colnames(newdat1)[1] <- colnames(newdat2)[1] <-  paste(exposure_names, collapse = ":")
  }else{
    colnames(newdat1)[1] <- colnames(newdat2)[1] <-  paste(exposure_names, collapse = "x")
  }
  
  return(predict.logisf(coefs1, newdat2) / predict.logisf(coefs1, newdat1))
}

# loading datasets

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
if(case_type == "suspected_microcephaly"){
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_SumEXP.RData")
}else{
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_SumEXP.RData")
}




# ok start preparing table

load(paste("WORKSPACE/MODELS/", case_type, "/basemodAdj.RData", sep = ""))
basemodAIC <- extractAIC(Adjmod[[2]])[2]
rm(Adjmod)
# construct the master table
mastab <- data.frame(Pri_expo = c("Dengue_PREG",
                                  "Chikungunya_PREG",
                                  "Cattle",
                                  "Water",
                                  "Zika_PREG",
                                  "Zika_PREG",
                                  "Zika_PREG",
                                  "Zika_PREG"),
                     Cofactors = c("",
                                   "",
                                   "",
                                   "",
                                   "",
                                   "Dengue_PREG",
                                   "Chikungunya_PREG",
                                   "YellowFeverV"),
                     Confounders = rep("-", 8),
                     ORUnAdj = rep(NA, 8),
                     ORUnAdjLow = rep(NA, 8),
                     ORUnAdjHigh = rep(NA, 8),
                     ORAdj = rep(NA, 8),
                     ORAdjLow = rep(NA, 8),
                     ORAdjHigh = rep(NA, 8),
                     DeltaAIC = rep(NA, 8))
mastab$Confounders = as.character(mastab$Confounders)

## exposure minimum and median values
# Zika
zikzero <- min(model_df[, arbo_type$zika_PREG], na.rm = T)
zikmed <- median(model_df[, arbo_type$zika_PREG][model_df[, arbo_type$zika_PREG] > zikzero], na.rm = T)
# Dengue
denzero <- min(model_df[, arbo_type$dengue_PREG], na.rm = T)
denmed <- median(model_df[, arbo_type$dengue_PREG][model_df[, arbo_type$dengue_PREG] > denzero], na.rm = T)
# Chik
chikzero <- min(model_df[, arbo_type$chik_PREG], na.rm = T)
chikmed <- median(model_df[, arbo_type$chik_PREG][model_df[, arbo_type$chik_PREG] > chikzero], na.rm = T)
# Farms
farmzero <- min(model_df$mun_farm_den, na.rm = T)
farmmed <- median(model_df$mun_farm_den, na.rm = T)
# Water
watzero <- min(model_df$mun_water_vul, na.rm = T)
watmed <- median(model_df$mun_water_vul, na.rm = T)
# Yellow Fever
yfzero <- min(model_df$YF_Vcov, na.rm = T)
yfmed <- median(model_df$YF_Vcov, na.rm = T)

############################
# 01 Primary exposures
############################

# Dengue 
# crude model:
denmodCrude <- logistf(as.formula(paste(case_type,"~", arbo_type$dengue_PREG)), data = model_df)

# model with selected covariates:
load(paste("WORKSPACE/MODELS/", case_type, "/denmodAdj.RData", sep = ""))
denmodAdj <- Adjmod
mastab[1, 3] = paste(denmodAdj[[1]], collapse = ", ")
mastab[1, 4:6] = RR.pred2B(model = denmodCrude, 
                           exposure_name = arbo_type$dengue_PREG,
                           zeroval = denzero,
                           medval = denmed)
mastab[1, 7:9] = RR.pred2B(model = denmodAdj[[2]],
                           exposure_name = arbo_type$dengue_PREG,
                           zeroval = denzero,
                           medval = denmed)
mastab[1, 10] = basemodAIC - extractAIC(denmodAdj[[2]])[2]

rm(denmodAdj, denmodCrude, Adjmod)

# Chikungunya
chimodCrude <- logistf(as.formula(paste(case_type,"~", arbo_type$chik_PREG)), data = model_df)
# model with selected covariates:
load(paste("WORKSPACE/MODELS/", case_type, "/chimodAdj.RData", sep = ""))
chimodAdj <- Adjmod
mastab[2, 3] = paste(chimodAdj[[1]], collapse = ", ")
mastab[2, 4:6] = RR.pred2B(model = chimodCrude, 
                           exposure_name = arbo_type$chik_PREG,
                           zeroval = chikzero,
                           medval = chikmed)
mastab[2, 7:9] = RR.pred2B(model = chimodAdj[[2]],
                           exposure_name = arbo_type$chik_PREG,
                           zeroval = chikzero,
                           medval = chikmed)
mastab[2, 10] = basemodAIC - extractAIC(chimodAdj[[2]])[2]

rm(chimodAdj, chimodCrude, Adjmod)

# Farms
farmmodCrude <- logistf(as.formula(paste(case_type,"~ mun_farm_den")), data = model_df)
# model with selected covariates:
load(paste("WORKSPACE/MODELS/", case_type, "/farmmodAdj.RData", sep = ""))
farmmodAdj <- Adjmod
mastab[3, 3] = paste(farmmodAdj[[1]], collapse = ", ")
mastab[3, 4:6] = RR.pred2B(model = farmmodCrude, 
                           exposure_name = "mun_farm_den",
                           zeroval = farmzero,
                           medval = farmmed)
mastab[3, 7:9] = RR.pred2B(model = farmmodAdj[[2]], 
                           exposure_name = "mun_farm_den",
                           zeroval = farmzero,
                           medval = farmmed)
mastab[3, 10] = basemodAIC - extractAIC(farmmodAdj[[2]])[2]

rm(farmmodAdj, farmmodCrude, Adjmod)


# Water
watmodCrude <- logistf(as.formula(paste(case_type,"~ mun_water_vul")), data = model_df)
# model with selected covariates:
load(paste("WORKSPACE/MODELS/", case_type, "/watmodAdj.RData", sep = ""))
watmodAdj <- Adjmod
mastab[4, 3] = paste(watmodAdj[[1]], collapse = ", ")
mastab[4, 4:6] = RR.pred2B(model = watmodCrude, 
                           exposure_name = "mun_water_vul",
                           zeroval = watzero,
                           medval = watmed)
mastab[4, 7:9] = RR.pred2B(model = watmodAdj[[2]], 
                           exposure_name = "mun_water_vul",
                           zeroval = watzero,
                           medval = watmed)
mastab[4, 10] = basemodAIC - extractAIC(watmodAdj[[2]])[2]

rm(watmodAdj, watmodCrude, Adjmod)

# Zika
zikmodCrude <- logistf(as.formula(paste(case_type,"~", arbo_type$zika_PREG)), data = model_df)
# model with selected covariates:
load(paste("WORKSPACE/MODELS/", case_type, "/zikmodAdj.RData", sep = ""))
zikmodAdj <- Adjmod
mastab[5, 3] = paste(zikmodAdj[[1]], collapse = ", ")
mastab[5, 4:6] = RR.pred2B(model = zikmodCrude, 
                           exposure_name = arbo_type$zika_PREG,
                           zeroval = zikzero,
                           medval = zikmed)
mastab[5, 7:9] = RR.pred2B(model = zikmodAdj[[2]], 
                           exposure_name = arbo_type$zika_PREG,
                           zeroval = zikzero,
                           medval = zikmed)
mastab[5, 10] = basemodAIC - extractAIC(zikmodAdj[[2]])[2]

rm(zikmodAdj, zikmodCrude, Adjmod)









############################
# 02 Zika cofactors
############################

# Zika with dengue (coinfection)
zikdenmodCrude <- logistf(as.formula(paste(case_type,"~", arbo_type$zika_PREG, ":", arbo_type$dengue_PREG)), data = model_df)
# model with selected covariates:
load(paste("WORKSPACE/MODELS/", case_type, "/zikdenmodAdj.RData", sep = ""))
zikdenmodAdj <- Adjmod
mastab[6, 3] = paste(zikdenmodAdj[[1]], collapse = ", ")
mastab[6, 4:6] = RR.pred2B.interaction(model = zikdenmodCrude, 
                                       exposure_names = c(arbo_type$zika_PREG, 
                                                          arbo_type$dengue_PREG),
                                       zeroval = zikzero,
                                       medval = zikmed,
                                       zeroval2 = denzero,
                                       medval2 = denmed)
mastab[6, 7:9] = RR.pred2B.interaction(model = zikdenmodAdj[[2]], 
                                       exposure_names = c(arbo_type$zika_PREG, 
                                                          arbo_type$dengue_PREG),
                                       zeroval = zikzero,
                                       medval = zikmed,
                                       zeroval2 = denzero,
                                       medval2 = denmed)
mastab[6, 10] = basemodAIC - extractAIC(zikdenmodAdj[[2]])[2]

rm(zikdenmodAdj, zikdenmodCrude, Adjmod)


# Zika with chikungunya (coinfection)
zikchimodCrude <- logistf(as.formula(paste(case_type,"~", arbo_type$zika_PREG, ":", arbo_type$chik_PREG)), data = model_df)

# model with selected covariates:
load(paste("WORKSPACE/MODELS/", case_type, "/zikchimodAdj.RData", sep = ""))
zikchimodAdj <- Adjmod
mastab[7, 3] = paste(zikchimodAdj[[1]], collapse = ", ")
mastab[7, 4:6] = RR.pred2B.interaction(model = zikchimodCrude, 
                                       exposure_names = c(arbo_type$zika_PREG, 
                                                          arbo_type$chik_PREG),
                                       zeroval = zikzero,
                                       medval = zikmed,
                                       zeroval2 = chikzero,
                                       medval2 = chikmed)
mastab[7, 7:9] = RR.pred2B.interaction(model = zikchimodAdj[[2]], 
                                       exposure_names = c(arbo_type$zika_PREG, 
                                                          arbo_type$chik_PREG),
                                       zeroval = zikzero,
                                       medval = zikmed,
                                       zeroval2 = chikzero,
                                       medval2 = chikmed)
mastab[7, 10] = basemodAIC - extractAIC(zikchimodAdj[[2]])[2]

rm(zikchimodAdj, zikchimodCrude, Adjmod)


# Zika with yellow fever vaccination
zikyfmodCrude <- logistf(as.formula(paste(case_type,"~", arbo_type$zika_PREG, ":YF_Vcov")), data = model_df)
# model with selected covariates:
load(paste("WORKSPACE/MODELS/", case_type, "/zikyfmodAdj.RData", sep = ""))
zikyfmodAdj <- Adjmod
mastab[8, 3] = paste(zikyfmodAdj[[1]], collapse = ", ")
mastab[8, 4:6] = RR.pred2B.interaction(model = zikyfmodCrude, 
                                       exposure_names = c(arbo_type$zika_PREG, 
                                                          "YF_Vcov"),
                                       zeroval = zikzero,
                                       medval = zikmed,
                                       zeroval2 = yfzero,
                                       medval2 = yfmed)
mastab[8, 7:9] = RR.pred2B.interaction(model = zikyfmodAdj[[2]], 
                                       exposure_names = c(arbo_type$zika_PREG, 
                                                          "YF_Vcov"),
                                       zeroval = zikzero,
                                       medval = zikmed,
                                       zeroval2 = yfzero,
                                       medval2 = yfmed)
# !!!! special edit for YF vaccine as appears to be protective
# therefore RR is for not vaccinated for YF
#mastab[9, 4:6] = 1 / mastab[9, 4:6]
#mastab[9, 7:9] = 1 / mastab[9, 7:9]
mastab[8, 4:6] = mastab[8, 4:6]
mastab[8, 7:9] = mastab[8, 7:9]


mastab[8, 10] = basemodAIC - extractAIC(zikyfmodAdj[[2]])[2]

rm(zikyfmodAdj, zikyfmodCrude, Adjmod)

# save summary table
write.csv(mastab, file = paste("FIGURES/Primary_cofactors_table_", case_type, ".csv", sep = ""))






# relaod top models
load(paste("WORKSPACE/MODELS/", case_type, "/zikmodAdj.RData", sep = ""))
zikmodAdj <- Adjmod
load(paste("WORKSPACE/MODELS/", case_type, "/zikdenmodAdj.RData", sep = ""))
zikdenmodAdj <- Adjmod
load(paste("WORKSPACE/MODELS/", case_type, "/zikchimodAdj.RData", sep = ""))
zikchimodAdj <- Adjmod
load(paste("WORKSPACE/MODELS/", case_type, "/zikyfmodAdj.RData", sep = ""))
zikyfmodAdj <- Adjmod

# test top 2 models:
anova(zikmodAdj[[2]], zikdenmodAdj[[2]],test="Chisq", method = "PLR") # 0.947
anova(zikmodAdj[[2]], zikchimodAdj[[2]],test="Chisq", method = "PLR") # 0.203
anova(zikmodAdj[[2]], zikyfmodAdj[[2]],test="Chisq", method = "PLR") # 2.026214e-05 





















########################
# Mean Z score hypothesis
########################
load(paste("WORKSPACE/MODELS/", case_type, "/zikmodAdj.RData", sep = ""))
zikmodAdj <- Adjmod
munZscore = read.csv("REFERENCE/Mun_Zscores.csv")
model_df$munZscore = munZscore$Zscore_mean[match(model_df$case_mun, munZscore$IBGE_CODE)]

zikmodAdjZscore = logistf(confirmed_microcephaly ~ mun_con_zik_incid_PREG + baby_sex + Region + munZscore,
                          data = model_df)
AdjZscore = logistf(confirmed_microcephaly ~ baby_sex + Region + munZscore,
                    data = model_df)
basemodAIC = -1506.615
basemodAIC - extractAIC(zikmodAdj[[2]])[2]
basemodAIC - extractAIC(zikmodAdjZscore)[2]
basemodAIC - extractAIC(AdjZscore)[2]

# Z score adds marginally to the model, but does very poorly on its own
# may be a slight contributing factor
anova(zikmodAdj[[2]], zikmodAdjZscore,test="Chisq", method = "nested") # 0.947








########################
# Race hypothesis
########################
load(paste("WORKSPACE/MODELS/", case_type, "/zikmodAdj.RData", sep = ""))
zikmodAdj <- Adjmod
zikmodAdj[[1]]
load("INDIVIDUAL/Indiv_reg_df_JUL18_170718_processed_confirmed_RACE.RData")
model_df_race = model_df

# new data so need to refit original model
zikmodAdj = logistf(confirmed_microcephaly ~ mun_con_zik_incid_PREG + baby_sex + Region,
                    data = model_df_race)
zikmodAdjRace = logistf(confirmed_microcephaly ~ mun_con_zik_incid_PREG + baby_sex + Region + mother_race,
                        data = model_df_race)
zikmodAdjRaceInt = logistf(confirmed_microcephaly ~ mun_con_zik_incid_PREG:mother_race + baby_sex + Region,
                           data = model_df_race)
AdjRace = logistf(confirmed_microcephaly ~ baby_sex +  Region + mother_race,
                  data = model_df_race)

extractAIC(zikmodAdj)[2]
extractAIC(zikmodAdjRace)[2]
extractAIC(zikmodAdjRaceInt)[2]
extractAIC(AdjRace)[2]

# Z score adds marginally to the model, but does very poorly on its own
# may be a slight contributing factor
anova(zikmodAdj, zikmodAdjRace,test="Chisq", method = "nested") # 0.430









########################
# alternative poisson model struture
########################


mun_demo <- read.csv("REFERENCE/GAUL_IBGE_conversion_full.csv")
muntab <- t(table(model_df$confirmed_microcephaly, model_df$case_mun))
muntab = data.frame(mun = rownames(muntab),
                    Birth_Norm = muntab[, 1],
                    Birth_MC = muntab[, 2],
                    Birth_MCrate = muntab[, 2] / muntab[, 1],
                    Birth_tot = muntab[, 2] + muntab[, 1])
# load Zika timeseries
zik_con_mun <- read.csv("MUN_COVS/JUL18_timeseries/zik_con_mun.csv")
muntab$CUMZIK_I <- rowSums(zik_con_mun[match(muntab$mun, zik_con_mun$IBGE_Municipality), 5:160], na.rm = T)
muntab$pop = mun_demo$IBGE_MUN_POP[match(muntab$mun, mun_demo$IBGE_CODE)]
muntab$CUMZIK_I = log((muntab$CUMZIK_I + 1) / muntab$pop)

# dengue and chikungunya
den_con_mun <- read.csv("MUN_COVS/JUL18_timeseries/den_con_mun.csv")
muntab$CUMDEN_I <- rowSums(den_con_mun[match(muntab$mun, den_con_mun$IBGE_Municipality), 5:160], na.rm = T)
muntab$CUMDEN_I = log((muntab$CUMDEN_I + 1) / muntab$pop)

chi_con_mun <- read.csv("MUN_COVS/JUL18_timeseries/chik_con_mun.csv")
muntab$CUMCHI_I <- rowSums(chi_con_mun[match(muntab$mun, chi_con_mun$IBGE_Municipality), 5:160], na.rm = T)
muntab$CUMCHI_I = log((muntab$CUMCHI_I + 1) / muntab$pop)

# BVDV and water toxins
cattle <- read.csv("BVDV/Bovines_mun_2015_rawdata.csv")
water_vul <- read.csv("WATER/Water_supply_vulnerability_mun_2014_2016.csv")

cattle$farm_den <- (cattle$rural.property_number / as.numeric(as.character(cattle$municipality_area_Km2)))
muntab$farm_den <- cattle$farm_den[match(muntab$mun, cattle$municipality_code)]

muntab$water_vul <- water_vul$PV.AGUA_2016[match(muntab$mun, water_vul$COD_IBGE)]


# add region and average mother age and SDI
muntab$Region = mun_demo$Region[match(muntab$mun, mun_demo$IBGE_CODE)]

avg_age <- aggregate(model_df$mother_age, by = list(model_df$case_mun), FUN = mean)
muntab$mother_age = avg_age[match(muntab$mun, avg_age[, 1]), 2]

avg_SDI <- aggregate(model_df$SDI, by = list(model_df$case_mun), FUN = mean)
muntab$SDI = avg_SDI[match(muntab$mun, avg_age[, 1]), 2]

# add YF vaccine coverage
YF_Vcov <- read.csv("MUN_COVS/YF_vaccination_2016_V2.csv")
G_codevec <- mun_demo$GAUL_CODE[match(muntab$mun, mun_demo$IBGE_CODE)]
muntab$YF_Vcov <- 100 * YF_Vcov$coverage_untargeted_unbiased_2016[match(G_codevec, YF_Vcov$GAUL_CODE)]

# add state
muntab$sta = mun_demo$IBGE_STATE_CODE[match(muntab$mun, mun_demo$IBGE_CODE)]

# remove any NAs - just removes 1 municipality with missing farm info
muntab = muntab[!is.na(muntab$farm_den), ]


### ok model fitting and selection

step.custom.poisson <- function(response, covariates, selectable, glmdata){
  # top loop - covariate removal loop
  candidates = covariates[selectable]
  responsedat = glmdata[, (colnames(glmdata) %in% response)]
  bestmodelAIC <- rep(NA, length(covariates[selectable]))
  bestmodelcovs <- list()
  #sfInit(parallel=TRUE, cpus=4)
  #sfExport("responsedat", "glmdata", local = TRUE)
  # pre-reduction model AIC
  #covdat = glmdata[, colnames(glmdata) %in% c(response, covariates)]
  covdat = glmdata
  fform <- as.formula(paste(response, "~", paste(covariates, collapse = "+")))
  candmod <- glm(fform, data = covdat, family = "poisson")
  fullmodAIC = extractAIC(candmod)[2]
  rm(candmod)
  
  for(i in length(covariates[selectable]):2){
    print(paste("Model including", (i-1), "of", sum(selectable), "candidates"))
    covlist = list()
    for(k in 1:length(candidates)){covlist[[(k)]] = c(covariates[!selectable], candidates[((1:length(candidates)) != k)])}
    #scovlist <- split(covlist, ceiling((1:length(covlist)) / round(i / 6)))
    #AICvec <- sfLapply(covlist, glm.fit.para)
    #system.time({gh  = glm.fit.para(covlist[[1]])})
    #AICvec <- unlist(lapply(covlist, glm.fit.para, glmdata = glmdata, responsedat = responsedat))
    AICvec = vector()
    for(k in 1:length(covlist)){
      fform <- as.formula(paste(response, "~", paste(covlist[[k]], collapse = "+")))
      covdat = glmdata[, colSums(sapply(colnames(glmdata), grepl , c(covlist[[k]], response))) > 0]
      candmod <- glm(fform, data = covdat, family = "poisson")
      AICvec[k] = extractAIC(candmod)[2]
      rm(candmod)
    }
    
    #if(i == length(covariates[selectable])){
    #  curbest = list(min(AICvec), flist[[which.min(AICvec)]])
    #}
    #if(min(AICvec) <= curbest[[1]]){curbest = list(min(AICvec), flist[[which.min(AICvec)]])}
    bestmodelAIC[i] = min(AICvec)
    bestmodelcovs[[i]] = covlist[[which.min(AICvec)]]
    
    print(paste("Dropping", candidates[which.min(AICvec)]))
    candidates = candidates[((1:length(candidates)) != which.min(AICvec))]
  }
  # add original model
  bestmodelAIC = c(bestmodelAIC, fullmodAIC)
  bestmodelcovs[[(length(bestmodelcovs) + 1)]] = covariates
  # final model
  fform <- as.formula(paste(response, "~", paste(bestmodelcovs[[which.min(bestmodelAIC)]], collapse = "+")))
  finalmod = glm(fform, glmdata, family = "poisson")
  return(list(bestmodelcovs[[which.min(bestmodelAIC)]], finalmod, bestmodelAIC, bestmodelcovs))
}

#basevars <- c("sta", "mother_age", "SDI")
basevars <- c("Region", "mother_age", "SDI")

# baseline model
basemod <- step.custom.poisson(response <- "Birth_MC",
                               covariates <- c("Birth_tot", basevars),
                               selectable <- c(FALSE, TRUE, TRUE, TRUE),
                               glmdata = muntab)

# 01 Dengue
denmod_poi <- step.custom.poisson(response <- "Birth_MC",
                                  covariates <- c("Birth_tot", "CUMDEN_I", basevars),
                                  selectable <- c(FALSE, FALSE, TRUE, TRUE, TRUE),
                                  glmdata = muntab)

# 02 Chik
chimod_poi <- step.custom.poisson(response <- "Birth_MC",
                                  covariates <- c("Birth_tot", "CUMCHI_I", basevars),
                                  selectable <- c(FALSE, FALSE, TRUE, TRUE, TRUE),
                                  glmdata = muntab)

# 03 BVDV
bvdvmod_poi <- step.custom.poisson(response <- "Birth_MC",
                                   covariates <- c("Birth_tot", "farm_den", basevars),
                                   selectable <- c(FALSE, FALSE, TRUE, TRUE, TRUE),
                                   glmdata = muntab)

# 04 water vulnerability
water_vulmod_poi <- step.custom.poisson(response <- "Birth_MC",
                                        covariates <- c("Birth_tot", "water_vul", basevars),
                                        selectable <- c(FALSE, FALSE, TRUE, TRUE, TRUE),
                                        glmdata = muntab)

# 05 Zika
zikmod_poi <- step.custom.poisson(response <- "Birth_MC",
                                  covariates <- c("Birth_tot", "CUMZIK_I", basevars),
                                  selectable <- c(FALSE, FALSE, TRUE, TRUE, TRUE),
                                  glmdata = muntab)

# 06 Zika and Dengue
zik_denmod_poi <- step.custom.poisson(response <- "Birth_MC",
                                      covariates <- c("Birth_tot", "CUMZIK_I:CUMDEN_I", basevars),
                                      selectable <- c(FALSE, FALSE, TRUE, TRUE, TRUE),
                                      glmdata = muntab)

# 07 Zika and Chik
zik_chimod_poi <- step.custom.poisson(response <- "Birth_MC",
                                      covariates <- c("Birth_tot", "CUMZIK_I:CUMCHI_I", basevars),
                                      selectable <- c(FALSE, FALSE, TRUE, TRUE, TRUE),
                                      glmdata = muntab)

# 07 Zika and YF vaccine
zik_yfvmod_poi <- step.custom.poisson(response <- "Birth_MC",
                                      covariates <- c("Birth_tot", "CUMZIK_I:YF_Vcov", basevars),
                                      selectable <- c(FALSE, FALSE, TRUE, TRUE, TRUE),
                                      glmdata = muntab)

# basemod AIC
AIC(basemod[[2]])

# reductions in AIC in candidate exposures:
# dengue
AIC(basemod[[2]]) - AIC(denmod_poi[[2]])
# chik
AIC(basemod[[2]]) - AIC(chimod_poi[[2]])
# BVDV
AIC(basemod[[2]]) - AIC(bvdvmod_poi[[2]])
# water vulnerability
AIC(basemod[[2]]) - AIC(water_vulmod_poi[[2]])
# Zika
AIC(basemod[[2]]) - AIC(zikmod_poi[[2]])

# Zika and dengue
AIC(basemod[[2]]) - AIC(zik_denmod_poi[[2]])
# Zika and chik
AIC(basemod[[2]]) - AIC(zik_chimod_poi[[2]])
# Zika and YF vaccine
AIC(basemod[[2]]) - AIC(zik_yfvmod_poi[[2]])









