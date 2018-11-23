rm(list = ls())
gc()

require(parallel)
require(logistf)
require(lubridate)

setwd("/Users/eideobra/Desktop/Brazil_data/BEST")
setwd("/home/admin/OBRADY/Brazil_data")

case_type <- "confirmed_microcephaly" #"confirmed_microcephaly" # "suspected_microcephaly"

testrun = FALSE

if(case_type == "confirmed_microcephaly"){
  # Full exposure history
  load("INDIVIDUAL/NOV18_issue/Zika_FES_processed_confirmed_FullEXP.RData")
  # accompanying detail
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_confirmed_FullEXP.RData")
}else{
  # Full exposure history
  load("INDIVIDUAL/NOV18_issue/Zika_FES_processed_suspected_FullEXP.RData")
  # accompanying detail
  load("INDIVIDUAL/NOV18_issue/Indiv_reg_df_NOV18_131118_processed_suspected_FullEXP.RData")
}

# renaming
model_df = complete_dat
rm(complete_dat)

# automatically fills in the type for arbovirus darta
if(case_type == "suspected_microcephaly"){
  zik = "mun_sus_zik_incid_PREG"
  micro = "suspected_microcephaly"
}else{
  zik = "mun_con_zik_incid_PREG"
  micro = "confirmed_microcephaly"
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
  FES = FES[c(micro_sam, normal_sam), ]
}



#########################
# 01 initial data processing
#########################



### adding finer spatiotemporal variables
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
model_df$Region_plus[model_df$case_sta == "BA"] = "Northeast_BA"
model_df$Region_plus[model_df$case_sta == "PE"] = "Northeast_PE"
model_df$Region_plus[model_df$case_sta == "RR"] = "North_RR"
model_df$Region_plus[model_df$case_sta == "PI"] = "Northeast_PI"
model_df$Region_plus = as.factor(model_df$Region_plus)

# time period modifier for finer temporal resolution
# Key period is October 2-15 - July 2016 (when national incidence first exceeded and finally
# descended below 4 MWSD per 10k live births)
# use simple year effect for the rest of the points
require(lubridate)
md_month <- as.Date("2015-01-01") + model_df$week * 7
md_month2 <- year(md_month) * 100 + month(md_month)
model_df$TF = year(md_month)
model_df$TF[md_month2 %in% c(201510:201512, 201601:201607)] = md_month2[md_month2 %in% c(201510:201512, 201601:201607)]
model_df$TF = as.factor(model_df$TF)
rm(md_month, md_month2)


#########################
# 02 aggregating exposure into thre trimesters and preconception
#########################

# week sequence
wkseq = seq(1, 55, 5)

# sero inflation transform
ZIT <- function(x){
  log(x + 0.5 * min(x[x > 0], na.rm = T))
}
FESvec = as.vector(FES[, 2:56])
FESvec = ZIT(FESvec)
FES[, 2:56] = FESvec


# collect prediction values for later RR predictions
#minval = min(FESvec, na.rm = T)
if(case_type == "suspected_microcephaly"){
  minval = -18.82506
  medianval = -12.64266
}else{
  minval = -20.64291 # confirmed
  medianval = -14.20591 # confirmed
}

#medianval = median(FESvec[FESvec != min(FESvec, na.rm = T)], na.rm = T)
medianMA = median(model_df$mother_age, na.rm = T)
medianweek = median(model_df$week, na.rm = T)
rm(FESvec)

# cols 2-11 are pre birth
# cols 12-51 are pregnancy - 12-24 = T1, 25-37 = T2, 38-51 = T3
# cols 52-56 are late deliveries / post birth
# bit of preprocessing to work out "post-birth"
deldate <- cbind(model_df$weeks_gestation, FES[, 12:56])
PBE <- rep(0, nrow(deldate))
T3exp <- rep(0, nrow(deldate))
for(i in 1:nrow(deldate)){
  PBE[i] = sum(deldate[i, (deldate[i, 1] + 1):ncol(deldate)], na.rm = T) / length((deldate[i, 1] + 1):ncol(deldate))
  if(deldate[i, 1] > 27){
    T3exp[i] <- sum(deldate[i, 28:deldate[i, 1]], na.rm = T) / 13
  }
}


# filter to only contain relevant variables
dat = data.frame(MC = model_df[, case_type],
                 Precon = rowSums(FES[, 2:11], na.rm = T) / length(2:11),
                 #T1 = rowSums(FES[, 12:24], na.rm = T),
                 #T2 = rowSums(FES[, 25:37], na.rm = T),
                 T3 = T3exp,
                 Postbirth = PBE,
                 T1T2 = rowSums(FES[, 12:37], na.rm = T) / 26,
                 mother_age = model_df$mother_age,
                 week = model_df$week,
                 Region = model_df$Region,
                 baby_sex = model_df$baby_sex,
                 case_sta = model_df$case_sta,
                 SDI = model_df$SDI,
                 Region_plus = model_df$Region_plus,
                 TF = model_df$TF)
# clear up
rm(deldate, PBE, T3exp, FES, model_df)

# trim to just those with data
# exclusions
sum(dat$Precon == 0)
sum(dat$T1 == 0)
sum(dat$T2 == 0)
sum(dat$T3 == 0)
dat = dat[!apply(dat[, 2:5], 1, function(x) any(x == 0)), ]

save(dat, file = paste0("INDIVIDUAL/NOV18_issue/FES_compiled_", case_type, ".RData"))

# assess correlation
#covcombs <- combn(c("Precon", "T1", "T2", "T3", "Postbirth", "T1T2"), 2)
#covcombs <- combn(c("Precon", "T1", "T2", "T3", "T1T2"), 2)
#cormat <- rep(NA, ncol(covcombs))

#for(i in 1:ncol(covcombs)){
#  cormat[i] = cor(dat[, covcombs[1, i]], dat[, covcombs[2, i]])
#}
#sum(sqrt(cormat^2) > 0.7) # hopefully equal to 0, but not
#covcombs[, (sqrt(cormat^2) > 0.7)]
# T1 and T2 colinear for both suspected and confirmed microcephaly
# - run all combinations of models and assess differences
# but model selection on default formula first

base_risk <- c("mother_age", "baby_sex","TF", "Region_plus", "SDI")
exposures <- c("Precon", "T1T2", "T3")

# variable selection options (31 total)
b_modform_options = c(as.list(data.frame(combn(base_risk, 5), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 4), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 3), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 2), stringsAsFactors = F)),
                      as.list(data.frame(combn(base_risk, 1), stringsAsFactors = F)))



# wrapper funciton for model running and extracting AIC
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

# now formulae with different outcomes
cand_mod_form <- list()
no_covs = as.formula(paste("MC", "~", paste(exposures, collapse = "+")))
cand_mod_form = c(cand_mod_form, c(no_covs,
                                   lapply(b_modform_options, 
                                          function(x, outcome, exposures){as.formula(paste(outcome, "~", paste(c(exposures, x), collapse = "+")))},
                                          outcome = "MC",
                                          exposures = exposures)))
# 02 run models in parralell
# start cluster
cl <- makeCluster(8, type = "FORK")
clusterEvalQ(cl, library(logistf))
# evaluate model combinations
Sys.time()
modAICs = parLapply(cl, cand_mod_form, HT.model.wrapper, mdf2 = dat)
Sys.time()
stopCluster(cl)

# now fit final model
final_form = cand_mod_form[which.min(modAICs)][[1]]
finalmod = HT.model.wrapper(final_form, mdf2 = dat, modrtn = T)

# save finalmodel
save(finalmod, file = paste0("MODELS/", case_type, "/FESmod.RData"))

















# reload model with selected covariates
load(file = paste0("WORKSPACE/MODELS/", case_type, "/FESmod.RData"))
mod <- finalmod

predict.logisf.FES <- function(coefs, newdata){
  # factors recognised in this analysis
  recog_factors <- c("Region_plus", "TF", "baby_sex")
  
  # start val processing
  val = coefs[1]
  for(i in 1:nrow(newdata)){
    for(j in 1:ncol(newdata)){
      # check if focus variable is a factor 
      if(colnames(newdata)[j] %in% recog_factors){
        # check if it is the one factor that is the referene value
        fac_detailed_name = paste0(colnames(newdata)[j], newdata[i, j])
        if(fac_detailed_name %in% rownames(coefs)){
          val = val + coefs[rownames(coefs) == fac_detailed_name]
        }
      }else{
        val = val + newdata[i, j] * coefs[rownames(coefs) == colnames(newdata)[j]]
      }
    }
  }
  # retransform to response scale
  val = exp(val)/(1+exp(val))
  return(val)
}

coefs1 <- cbind(coefficients(mod), mod$ci.lower, mod$ci.upper)
zikminval <- minval
zikmedianval = medianval

# predict relative risks
#EXPopts <- c("Precon", "T1T2", "T3", "Postbirth")
EXPopts <- c("Precon", "T1T2", "T3")
RRcol = matrix(NA, nrow = length(EXPopts), ncol = 3)
for(i in 1:length(EXPopts)){
  newdat1 = data.frame(Precon = zikminval,
                       #T1 = zikminval,
                       #T2 = zikminval,
                       T1T2 = zikminval,
                       T3 = zikminval,
                       #Postbirth = zikminval,
                       #mother_age = medianMA,
                       baby_sex = "M",
                       TF = "2016",
                       Region_plus = "Northeast")
  newdat2 = newdat1
  newdat1[, EXPopts[i]] = zikmedianval
  
  p1 <- predict.logisf.FES(coefs1, newdat1)
  p2 <- predict.logisf.FES(coefs1, newdat2)
  RRcol[i, 1:3] = as.numeric(p1 / p2)
}



RRcol = as.data.frame(RRcol)
RRcol$time = c("Pre-conception",
               #"Trimester 1",
               #"Trimester 2",
               "Trimester1-2",
               "Trimester 3")

pdf(width = 10, height = 5, file = paste("FIGURES/ZIKV_gestational_risk", case_type, ".pdf", sep = ""))
if(case_type == "suspected_microcephaly"){
  plot((1:3) + 0.5, RRcol[, 1], ylim = c(0, 7), xlab = "", xlim = c(1, 4),
       ylab = "Relative risk of microcephaly", pch = 16, cex = 1.5,
       cex.lab = 1.25, cex.axis = 1.5, axes = F)
}else{
  plot((1:3) + 0.5, RRcol[, 1], ylim = c(0, 7), xlab = "", xlim = c(1, 4),
       ylab = "Relative risk of MWSD", pch = 16, cex = 1.5,
       cex.lab = 1.25, cex.axis = 1.5, axes = F)
}


axis(2, seq(0, 7, 1), cex.axis = 1.25)
axis(1, (1:3) + 0.5, RRcol$time, cex.axis = 1.25)

#grid(nx = 0, ny = 15)
for(i in seq(0, 7, 1)){lines(c(-15, 45), c(i, i), col = "light gray", lty = 2, lwd = 0.75)}

abline(h = 1, col = "red", lwd = 1.5)

for(i in 1:nrow(RRcol)){
  polygon(x = c(i +0.5 - 0.45,  i +0.5 + 0.45, i +0.5 + 0.45,  i +0.5 - 0.45),
          y = c(RRcol[i, 2], RRcol[i, 2], RRcol[i, 3], RRcol[i, 3]),
          col = rgb(49 / 255, 130 / 255, 189 / 255, 0.5))
}
dev.off()






# probing other features of the model for the paper
newdat1 = data.frame(Precon = zikminval,
                     T1T2 = zikminval,
                     T3 = zikminval,
                     baby_sex = "M",
                     TF = "2016",
                     Region_plus = "Northeast")
newdat2 = newdat1
newdat2$baby_sex = "F"

p1 <- predict.logisf.FES(coefs1, newdat1)
p2 <- predict.logisf.FES(coefs1, newdat2)
as.numeric(p2 / p1)

# plotting spatiotemporal residuals of the model
# temporal
resids = dat$MC - mod$predict
resids_week = aggregate(resids, by = list(dat$week), FUN = sum)
require(mgcv)
trend = gam(x ~ s(Group.1), data = resids_week)
pdf(file = "FIGURES/Residual_TAC.pdf", height = 5, width = 8)
plot(resids_week$Group.1, resids_week$x, xlab = "Week", ylab = "Observed - predicted",
     ylim = c(-10, 20), main = "Model residuals over time")
#trend_preds <- predict(trend, se.fit = T)
#lines(resids_week$Group.1, trend_preds$fit, col = "red")
#lines(resids_week$Group.1, trend_preds$fit - 1.96 * trend_preds$se.fit, col = "black", lty = 2)
#lines(resids_week$Group.1, trend_preds$fit + 1.96 * trend_preds$se.fit, col = "black", lty = 2)
#abline(h = 0)
dev.off()

# space
resids_sta = aggregate(resids, by = list(dat$case_sta), FUN = sum)
pdf(file = "FIGURES/Residual_SAC.pdf", height = 5, width = 14)
plot(resids_sta, xlab = "state", ylab = "Observed - predicted")
dev.off()
