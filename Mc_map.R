# Microcephaly map
rm(list = ls())

require(logistf)
require(mgcv)
setwd("/Users/eideobra/Desktop/Brazil_data/BEST")

load("INDIVIDUAL/Indiv_reg_df_NOV18_170718.RData")


# fix errors in mother age
model_df$mother_age[model_df$mother_age > 55] = NA
model_df$mother_age[model_df$mother_age < 10] = NA

case_type = "confirmed_microcephaly"#"confirmed_microcephaly" # "suspected_microcephaly"

baseMc <- read.csv(paste0("MUN_COVS/Baseline_Mc_", case_type, ".csv"))
Zik_occ_messina <- read.csv("/Users/eideobra/Desktop/Brazil_data/BEST/REFERENCE/mun_zik_preds_messina.csv")


# 01 add a Socio Development index
sociocomps <- model_df[, c("mun_Income_pc",
                      "mun_Mean_years_edu",
                      "mun_Fert_rate")]
# standardize NAs
sociocomps[is.na(rowSums(sociocomps)), ] = NA
SDInaind <- is.na(sociocomps[, 1])
#Scale
sociocomps[!SDInaind, ] = apply(sociocomps[!SDInaind, ], 2, function(x) (x - min(x)) / (max(x) - min(x)))
SDI = rowSums(sociocomps) / 3
model_df$SDI = SDI


# prevalence of microcephaly at municipality level
mun_Mc = t(table(model_df[, case_type], model_df$case_mun))

# add region
sta_rs <- read.csv("STA_COVS/STA_regions.csv")
mun_demo <- read.csv("REFERENCE/GAUL_IBGE_conversion_full.csv")

mun_rs <- data.frame(mun = unique(model_df$case_mun))
model_df$Region <- sta_rs[match(model_df$case_sta, sta_rs$State), "Region"]
mun_rs$State = model_df$case_sta[match(mun_rs$mun, model_df$case_mun)]
mun_rs$Region = model_df$Region[match(mun_rs$mun, model_df$case_mun)]

# average mother age
mun_age = aggregate(model_df$mother_age, by = list(model_df$case_mun), FUN = mean, na.rm = T)

# average SDI
mun_SDI = aggregate(model_df$SDI, by = list(model_df$case_mun), FUN = mean, na.rm = T)

# Zika average probability of occurrence
mun_ZIKocc = Zik_occ_messina$Zik_pred[match(model_df$case_gaul_mun, Zik_occ_messina$MUN_GAUL)]
mun_ZIKocc = aggregate(mun_ZIKocc, by = list(model_df$case_mun), FUN = mean, na.rm = T)


mun_Mc = data.frame(mun = rownames(mun_Mc),
                    Norm_births = mun_Mc[, 1],
                    Mc_births = mun_Mc[, 2],
                    Mc_rate = mun_Mc[, 2] / (mun_Mc[, 1] + mun_Mc[, 2]),
                    #week = mun_MCtime$x,
                    SDI = mun_SDI$x,
                    Zik_occ = mun_ZIKocc$x,
                    mother_age = mun_age$x)
#mun_Mc$pred_MC = log(baseMc$actual_mean)[match(mun_Mc$mun, baseMc$mun)]

mun_Mc$total_births <- mun_Mc$Norm_births + mun_Mc$Mc_births



# add lat longs
munLL <- read.csv("REFERENCE/MUN_lat_long.csv")

mun_Mc$x = munLL$X1[match(mun_Mc$mun, munLL$IBGE)]
mun_Mc$y = munLL$X2[match(mun_Mc$mun, munLL$IBGE)]
#munLL$pred_MC = log(baseMc$actual_mean[match(munLL$IBGE, baseMc$mun)])
#munLL$pred_MC[is.na(munLL$pred_MC)] = median(munLL$pred_MC, na.rm = T)


# fit model
#Mcmod <- gam(cbind(Mc_births, Norm_births) ~ s(x, y, k = 150), mun_Mc, family = "binomial")

# differet options
mun_Mc = mun_Mc[apply(mun_Mc[, c("x",
                                 "y",
                                 "SDI",
                                 "mother_age",
                                 "Zik_occ")], 
                      1,
                      function(x) {!any(is.na(x))}), ]
Mcmod1 <- gam(cbind(Mc_births, Norm_births) ~ s(x, y, k = 1000, bs = "tp") + SDI + mother_age + Zik_occ, 
              mun_Mc, 
              family = "binomial",
              select = TRUE)
#Mcmod2 <- gam(cbind(Mc_births, Norm_births) ~ s(x, y, k = 1000) + Region, mun_Mc, family = "binomial")
#Mcmod3 <- gam(cbind(Mc_births, Norm_births) ~ s(x, y, k = 1000) + mother_age, mun_Mc, family = "binomial")
#Mcmod4 <- gam(cbind(Mc_births, Norm_births) ~ s(x, y, k = 1000) + Region + mother_age, mun_Mc, family = "binomial")
#Mcmod5 <- gam(cbind(Mc_births, Norm_births) ~ s(x, y, k = 1000) + pred_MC, mun_Mc, family = "binomial")

deviance(Mcmod1)
#deviance(Mcmod2)
#deviance(Mcmod4)
#deviance(Mcmod3)
#deviance(Mcmod5)
Mcmod <- Mcmod1
#Mcmod <- gam(cbind(Mc_births, Norm_births) ~ s(x, y, k = 1000) + pred_MC, mun_Mc, family = "binomial")
#AIC(Mcmod)
#AIC(Mcmod_alt1)

# make predictions to all municipalities
# assemble prediciton dataset
predset <- data.frame(x = munLL$X1, 
                      y = munLL$X2,
                      SDI = mun_SDI$x[match(munLL$IBGE, mun_SDI$Group.1)],
                      mother_age = mun_age$x[match(munLL$IBGE, mun_age$Group.1)],
                      Zik_occ = mun_ZIKocc$x[match(munLL$IBGE, mun_ZIKocc$Group.1)])
                      #pred_MC = munLL$pred_MC)
# replace NAs with median value
predset = apply(predset, 2, function(x){
  if(any(is.na(x))){
    gh <- x
    gh[is.na(gh)] = median(gh, na.rm = T)
    return(gh)
  }else(return(x))
})
predset = data.frame(predset)

Mcpreds <- predict(Mcmod, newdata = predset, se.fit = T)
ILOGIT <- function(x){exp(x) / (1 + exp(x))}

Mcpreds = data.frame(mun = munLL$GAUL, 
                     pred = ILOGIT(Mcpreds$fit), 
                     pred_low = ILOGIT(Mcpreds$fit - 1.96 * Mcpreds$se.fit),
                     pred_high = ILOGIT(Mcpreds$fit + 1.96 * Mcpreds$se.fit),
                     mun_IBGE = munLL$IBGE,
                     MC_obs = mun_Mc$Mc_rate[match(munLL$IBGE,mun_Mc$mun)])
# add region to Mcpreds
Mcpreds$Region = mun_demo$Region[match(Mcpreds$mun_IBGE, mun_demo$IBGE_CODE)]


# check stats of predictions are comparable to the data
datcheck <- data.frame(prediction = Mcpreds$pred, observed = Mcpreds$MC_obs,
                       ALLbirths = mun_Mc$total_births[match(Mcpreds$mun_IBGE, mun_Mc$mun)])
summary(log(datcheck$prediction + 0.000005))
summary(log(datcheck$observed + 0.000005))
# predictions have lower mean but are more dispersed
# what about areas with Mc cases
summary(log(datcheck$prediction[datcheck$observed > 0]))
summary(log(datcheck$observed[datcheck$observed > 0]))
hist(log(datcheck$prediction[datcheck$observed > 0]))
hist(log(datcheck$observed[datcheck$observed > 0]))
# still all slightly lower for prediction

# what do they mean in terms of number of Mc cases predicted?
sum(datcheck$prediction * datcheck$ALLbirths, na.rm = T)
sum(datcheck$observed * datcheck$ALLbirths, na.rm = T)
# almost identical, and not overdispersed, same distribution in the tail

# straight plot
plot(log10(datcheck[, 1]), log10(datcheck[, 2]), xlab = "Predicted by GAM smoother (log 10 scale)", 
     ylab = "Observed rate (log 10 scale)",
     cex = 10*datcheck$ALLbirths / max(datcheck$ALLbirths, na.rm = T),
     xlim = c(-5, -2.5), ylim = c(-5, -2))
lines(log(seq(0, 1, length.out = 100000)),
      log(seq(0, 1, length.out = 100000)))
# all seems to suggest that predictions are a good estiamte of rate



# calculate rates and RRs

Mcpreds$MCRR <- 10000 * Mcpreds$pred / baseMc$baseline_mean[match(Mcpreds$mun_IBGE, baseMc$mun)]
Mcpreds$MCRR_low <- 10000 * Mcpreds$pred_low / baseMc$baseline_low[match(Mcpreds$mun_IBGE, baseMc$mun)]
Mcpreds$MCRR_high <- 10000 * Mcpreds$pred_high / baseMc$baseline_high[match(Mcpreds$mun_IBGE, baseMc$mun)]

Mcpreds$MCrate <- 10000 * Mcpreds$pred - baseMc$baseline_mean[match(Mcpreds$mun_IBGE, baseMc$mun)]
Mcpreds$MCrate_low <- 10000 * Mcpreds$pred_low - baseMc$baseline_low[match(Mcpreds$mun_IBGE, baseMc$mun)]
Mcpreds$MCrate_high <- 10000 * Mcpreds$pred_high - baseMc$baseline_high[match(Mcpreds$mun_IBGE, baseMc$mun)]

# constrain to be positive
Mcpreds$MCrate[Mcpreds$MCrate < 0] = 0
Mcpreds$MCrate_low[Mcpreds$MCrate_low < 0] = 0
Mcpreds$MCrate_high[Mcpreds$MCrate_high < 0] = 0

# max rate
Mcpreds[as.numeric(which.max(Mcpreds$MCrate)), ]
# region wide rate in NE
mean(Mcpreds[Mcpreds$Region == "Northeast", "MCRR"], na.rm = T)
mean(Mcpreds[Mcpreds$Region == "Northeast", "MCRR_low"], na.rm = T)
mean(Mcpreds[Mcpreds$Region == "Northeast", "MCRR_high"], na.rm = T)

# max RR
Mcpreds[as.numeric(which.max(Mcpreds$MCRR)), ]

# process Mcpreds to be compatible with ArcGIS
require(rgdal)
ad2 <- readOGR("SHAPE/Admin2_BRA.shp", "Admin2_BRA")
UNIGAULs <- unique(ad2$GAUL_CODE)
all(UNIGAULs %in% Mcpreds$mun) # good all there
# fix NAs
Mcpreds$MC_obs[is.na(Mcpreds$MC_obs)] = 0
Mcpreds$MCRR[is.na(Mcpreds$MCRR)] = 1
Mcpreds$MCRR_low[is.na(Mcpreds$MCRR_high)] = 1
Mcpreds$MCRR_high[is.na(Mcpreds$MCRR_high)] = 1

Mcpreds$MCrate[is.na(Mcpreds$MCrate)] = 0
Mcpreds$MCrate_low[is.na(Mcpreds$MCrate_low)] = 0
Mcpreds$MCrate_high[is.na(Mcpreds$MCrate_high)] = 0

# ad an FID column
Mcpreds = data.frame(FID = 1:nrow(Mcpreds),
                     Mcpreds)


# save Mcpreds
write.csv(Mcpreds, paste("/Users/eideobra/Desktop/Brazil_data/BEST/FIGURES/Microcephaly_map_",
                         case_type,
                         ".csv",
                         sep = ""))






# Salvador extraction for Table 1 numbers:
# salvador GAUL code = 6865
Mcpreds2 = Mcpreds[Mcpreds$mun == 6865, ]

# absolute risk- community-level
Mcpreds2$pred * 10000
Mcpreds2$pred_low * 10000
Mcpreds2$pred_high * 10000

# absolute risk- individual-level
(Mcpreds2$pred * 10000) / 0.633
(Mcpreds2$pred_low * 10000) / 0.668
(Mcpreds2$pred_high * 10000) / 0.594

# relative risk community level
Mcpreds2[, 9:11]
# relative risk at an individual level
Mcpreds2[, 9:11] / c(0.633, 0.594, 0.668)


# National absolute rate
model_df$case_gaul_mun <- mun_demo$GAUL_CODE[match(model_df$case_mun, mun_demo$IBGE_CODE)]
mun_births <- aggregate(model_df$confirmed_microcephaly >=0, by = list(model_df$case_gaul_mun),
                        FUN = sum)
Mcpreds$births <- mun_births$x[match(Mcpreds$mun, mun_births$Group.1)]

# total Mc cases divided my total births
totB <- sum(mun_births$x)
10000 * sum(Mcpreds$births * Mcpreds$pred, na.rm = T) / totB
10000 * sum(Mcpreds$births * Mcpreds$pred_low, na.rm = T) / totB
10000 * sum(Mcpreds$births * Mcpreds$pred_high, na.rm = T) / totB

# national RR
Mcpreds$births[is.na(Mcpreds$births)] = 344
weighted.mean(x = Mcpreds$MCRR, w = Mcpreds$births)
weighted.mean(x = Mcpreds$MCRR_low, w = Mcpreds$births)
weighted.mean(x = Mcpreds$MCRR_high, w = Mcpreds$births)


# 
# # plotting
# require(rgdal)
# require(sp)
# library(RColorBrewer)
# require(ggplot2)
# ad2 <- readOGR("SHAPE/Admin2_BRA.shp", "Admin2_BRA")
# 
# ad2$predMCRR <- Mcpreds$MCRR[match(ad2$GAUL_CODE, Mcpreds$mun)]
# ad2$predMCrate <- Mcpreds$MCrate[match(ad2$GAUL_CODE, Mcpreds$mun)]
# tempval = 10000 * Mcpreds$MC_obs[match(ad2$GAUL_CODE, Mcpreds$mun)] - baseMc$baseline_mean[match(Mcpreds$mun_IBGE, baseMc$mun)]
# tempval[tempval < 0] = 0
# ad2$obsMCrate <- tempval
# 
# 
# # # Zika RR and rate of Zika- associated Mc map
# # my.palette2 <- rev(brewer.pal(n = 7, name = "Spectral"))
# pdf("FIGURES/MCRR_riskmap_RR.pdf", width = 11, height = 10)
# spplot(ad2, "predMCRR",  col.regions = my.palette2, cuts = 6, col = "transparent",
#         colorkey = list(labels = list( labels = rev(c("74", 
#                                                       "61", 
#                                                       "50", 
#                                                       "37",
#                                                       "25",
#                                                       "12",
#                                                       "0")),
#                                        width = 2.5, cex = 2.5)))
#  dev.off()
 # 
 # pdf("FIGURES/MCRR_riskmap_rate.pdf", width = 11, height = 10)
 # spplot(ad2, "predMCrate",  col.regions = my.palette2, cuts = 6, col = "transparent",
 #        colorkey = list(labels = list( labels = rev(c("42", 
 #                                                      "35", 
 #                                                      "28", 
 #                                                      "21", 
 #                                                      "14",
 #                                                      "7",
 #                                                      "0")),
 #                                       width = 2.5, cex = 2.5)))
 # dev.off()
# # 
# # pdf("FIGURES/MCRR_riskmap_OBS_rate.pdf", width = 11, height = 10)
# # spplot(ad2, "obsMCrate",  col.regions = my.palette2, cuts = 6, col = "transparent",
# #        colorkey = list(labels = list( labels = rev(c("42", 
# #                                                      "35", 
# #                                                      "28", 
# #                                                      "21", 
# #                                                      "14",
# #                                                      "7",
# #                                                      "0")),
# #                                       width = 2.5, cex = 2.5)))
# # dev.off()
# 
# 
# 
# ad2$predMCRR[is.na(ad2$predMCRR)] = 0
# ad2$predMCrate[is.na(ad2$predMCrate)] = 0
# ad2$obsMCrate[is.na(ad2$obsMCrate)] = 0
# ad2$predMCRR = as.numeric(ad2$predMCRR)
# ad2$predMCrate = as.numeric(ad2$predMCrate)
# ad2$obsMCrate = as.numeric(ad2$obsMCrate)
# 
# # now making maps in ArcGIS so export
# writeOGR(ad2, dsn = "FIGURES/Mcmaps_shapefile.shp", layer = "Mcmaps_shapefile",
#          driver="ESRI Shapefile")
