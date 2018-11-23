# Links together data on births (SINASC), microcepahly (RESP) and Zika infection (GAL)
# AUG17 version with greater focus on infection in 1st trimester

# Oliver Brady 13/09/17
rm(list = ls())
require(plyr)
require(foreign)
require(data.table)
require(pryr)
require(geosphere)
setwd("/Users/Admin/Desktop/Brazil_data/BEST")
setwd("/Users/eideobra/Desktop/Brazil_data/BEST")

###################
# part 1 load in data
###################

# resp and sinasc
resp <-read.csv("RESP/PREDICTIONS/RESP_AUG17update_V2_280917.csv")
sinasc <- data.frame(fread("SINASC/Sinasc_2015_2017_AUG17.csv"))
# early trim of sinasc
varlist <- c("dtnasc", "sexo", "racacor", "codmunres", "idademae", "codanomal", "semagestac",
             "peso", "escmae")
sinasc = sinasc[, varlist]

# admin unit lookup tables
mun_demo <- read.csv("REFERENCE/GAUL_IBGE_conversion_full.csv")
sta_gaul = read.csv("REFERENCE/STATE_GAULS.csv")
# compile population for each state
sta_list <- unique(mun_demo$IBGE_STATE_CODE)
sta_demo = aggregate(mun_demo$IBGE_MUN_POP, by = list(mun_demo$IBGE_STATE_CODE), FUN = sum)
colnames(sta_demo) = c("state", "pop")

# arbovirus timeseries
zik_sus_mun <- read.csv("MUN_COVS/JUL18_timeseries/zik_sus_mun.csv")
zik_con_mun <- read.csv("MUN_COVS/JUL18_timeseries/zik_con_mun.csv")

chik_sus_mun <- read.csv("MUN_COVS/JUL18_timeseries/chik_sus_mun.csv")
chik_con_mun <- read.csv("MUN_COVS/JUL18_timeseries/chik_con_mun.csv")

den_sus_mun <- read.csv("MUN_COVS/JUL18_timeseries/den_sus_mun.csv")
den_con_mun <- read.csv("MUN_COVS/JUL18_timeseries/den_con_mun.csv")

#zik_true_mun <- read.csv("MUN_COVS/Zika_true_mun_ts_AUG17_V2.csv")
#zik_con_mun <- read.csv("MUN_COVS/Zika_con_mun_ts_AUG17_V2.csv")
#zik_true_sp2 <- read.csv("MUN_COVS/Zika_true_mun_ts_sp2_AUG17_V2.csv")
#den_not_sta <- read.csv("STA_COVS/Dengue_not_sta_ts_AUG17_V2.csv")
#den_not_mun <- read.csv("MUN_COVS/Dengue_not_mun_ts_AUG17_V2.csv")
#chi_not_sta <- read.csv("STA_COVS/Chikungunya_not_sta_ts_AUG17_V2.csv")
#chi_not_mun <- read.csv("MUN_COVS/Chikungunya_not_mun_ts_AUG17_V2.csv")

# YF vaccine coverage
#YF_Vcov <- read.csv("MUN_COVS/YF_vaccination_2015.csv")
YF_Vcov <- read.csv("MUN_COVS/YF_vaccination_2016_V2.csv")
YF_Vcov$IBGE_mun = mun_demo$IBGE_CODE[match(YF_Vcov$GAUL_CODE, mun_demo$GAUL_CODE)]
YF_Vcov$YF_VC_2016 = 100 * YF_Vcov$coverage_untargeted_unbiased_2016

# !! check no extra columns have been inserted
Acols <- c("FID", "State", "IBGE_Municipality","Week_first_intro", 
           paste("Wk_", 1:156, sep = ""))
zik_con_mun = zik_con_mun[, colnames(zik_con_mun) %in% Acols]
zik_sus_mun = zik_sus_mun[, colnames(zik_sus_mun) %in% Acols]
chik_con_mun = chik_con_mun[, colnames(chik_con_mun) %in% Acols]
chik_sus_mun = chik_sus_mun[, colnames(chik_sus_mun) %in% Acols]
den_con_mun = den_con_mun[, colnames(den_con_mun) %in% Acols]
den_sus_mun = den_sus_mun[, colnames(den_sus_mun) %in% Acols]

# load in sociodemographic variables
# may want to not use household sewerage as many (2000) NAs
sociodemo <- read.csv("SOCIODEMO/Municipalities_sociodemographics_chars_UNICODE.csv", stringsAsFactors = FALSE)
colnames(sociodemo) = gsub(patter = "\\.", replacement = "_", x = colnames(sociodemo))
# summarise mean years of mother education for each municiaplity
munedu <- aggregate(as.numeric(sinasc$escmae), list(sinasc$codmunres), FUN = mean, na.rm = T)
trim.trailing <- function (x) sub("\\s+$", "", x)
munedu$Group.1 = as.numeric(trim.trailing(munedu$Group.1))

sociodemo_useful <- data.frame(IBGE_mun = floor(sociodemo$Codigo / 10),
                               HH_piped_w = sociodemo$Domic_lios___com__gua_encanada___pessoas__2000_,
                               HH_garbage_col = sociodemo$Domic_lios___com_servi_o_de_coleta_de_lixo___pessoas__2000_,
                               HH_sewage = sociodemo$Domic_lios___com_instala__o_adequada_de_esgoto___pessoas__1991_,
                               Fert_rate = sociodemo$Taxa_de_fecundidade__2000_,
                               Income_pc = sociodemo$Renda_per_capita___2000,
                               Perc_rural = sociodemo$Popula__o_rural__2000_ / sociodemo$Popula__o__2000_,
                               Perc_employed = sociodemo$Populacao_Ocupada__2000_ / sociodemo$Popula__o__2000_,
                               Hospitalisations_pc = sociodemo$Produ__o_Ambulatorial_do_SUS___2016 / sociodemo$Popula__o__2000_,
                               Years_edu = munedu$x[match(floor(sociodemo$Codigo / 10),
                                                          munedu$Group.1)])

sociodemo_useful$State = mun_demo$IBGE_STATE_CODE[match(sociodemo_useful$IBGE_mun, mun_demo$IBGE_CODE)]
sociodemo_useful$pop = mun_demo$IBGE_MUN_POP[match(sociodemo_useful$IBGE_mun, mun_demo$IBGE_CODE)]
sociodemo_useful$pop[is.na(sociodemo_useful$pop)] = mean(sociodemo_useful$pop, na.rm = T)
sociodemo_useful_state = apply(sociodemo_useful[, 2:9], 2, function(y){
  ddply(data.frame(var = y, State = sociodemo_useful$State, pop = sociodemo_useful$pop), .(State),
        function(x) data.frame(wv=weighted.mean(x$var, 
                                                x$pop, na.rm = T)))[, 2]
})
sociodemo_useful_state = data.frame(State = sort(unique(sociodemo_useful$State)), sociodemo_useful_state[1:27, ])

# cattle
cattle <- read.csv("BVDV/Bovines_mun_2015_rawdata.csv")
cattle$municipality_area_Km2 = as.numeric(cattle$municipality_area_Km2)
cattle$pop = mun_demo$IBGE_MUN_POP[match(cattle$municipality_code, mun_demo$IBGE_CODE)]
cattle$pop[is.na(cattle$pop)] = mean(cattle$pop, na.rm = T)
cattle_state = data.frame(tcattle = aggregate(cattle$Total_animals, 
                                              by = list(cattle$federation_unity), 
                                              FUN = sum),
                          tfarms = aggregate(cattle$rural.property_number, 
                                             by = list(cattle$federation_unity), 
                                             FUN = sum)[, 2],
                          pop = aggregate(mun_demo$IBGE_MUN_POP, 
                                          by = list(mun_demo$IBGE_STATE_CODE), 
                                          FUN = sum)[, 2])
# maternal characteristics
matchar <- read.csv("SOCIODEMO/Maternal_characteristics.csv")
matchar$State = mun_demo$IBGE_STATE_CODE[match(matchar$IBGE_mun, mun_demo$IBGE_CODE)]
matchar$pop = mun_demo$IBGE_MUN_POP[match(matchar$IBGE_mun, mun_demo$IBGE_CODE)]
matchar$pop[is.na(matchar$pop)] = mean(matchar$pop, na.rm = T)
matchar_state = apply(matchar[, 3:8], 2, function(y){
  ddply(data.frame(var = y, State = matchar$State, pop = matchar$pop), .(State),
        function(x) data.frame(wv=weighted.mean(x$var, 
                                                x$pop, na.rm = T)))[, 2]
})
matchar_state = data.frame(State = sort(unique(matchar$State)), matchar_state[1:27, ])

# water supply
water_vul <- read.csv("WATER/Water_supply_vulnerability_mun_2014_2016.csv")
water_vul$pop = mun_demo$IBGE_MUN_POP[match(water_vul$COD_IBGE, mun_demo$IBGE_CODE)]
water_vul$pop[is.na(water_vul$pop)] = mean(water_vul$pop, na.rm = T)
water_vul_state = data.frame(WV14 = ddply(water_vul, .(GO),
                                          function(x) data.frame(wv=weighted.mean(x$PV.AGUA_2014, 
                                                                                  x$pop))),
                             WV15 = ddply(water_vul, .(GO),
                                          function(x) data.frame(wv=weighted.mean(x$PV.AGUA_2015, 
                                                                                  x$pop)))[, 2],
                             WV16 = ddply(water_vul, .(GO),
                                          function(x) data.frame(wv=weighted.mean(x$PV.AGUA_2016, 
                                                                                  x$pop)))[, 2])

###################
# resp and sinasc processing
###################

# RESP (microcephaly) data processing

# remove those that were due to confirmed toxoplasmosis or sifilis
# not any more - tryign to capture this background rate as well
#resp = resp[!xor((resp$ST_RESULTADO_SIFILIS == "REAGENTE/POSITIVO"), 
#                  resp$ST_RESULTADO_TOXOPLASMOSE == "REAGENTE/POSITIVO"), ]

# weeks of gestation from RESP matched to SINASC
resp$weeks_gestation = resp$DW_SEMANA_NASC
resp$weeks_gestation[resp$weeks_gestation > 200] = NA

resp$TP_SEXO = as.character(resp$TP_SEXO)
resp$TP_SEXO[resp$TP_SEXO == "FEMININO"] = "F"
resp$TP_SEXO[resp$TP_SEXO == "MASCULINO"] = "M"

# convert string to 2 letter abbreviation for resp and give it a population column
OTS <- data.frame(sort(unique(resp$SG_UF_RESIDENCIA)),
                  sort(unique(mun_demo$STATE_NAME))[1:27])
OTS2 <- OTS[match(resp$SG_UF_RESIDENCIA, OTS[, 1]), 2]
resp$SG_UF_RESIDENCIA <- mun_demo$IBGE_STATE_CODE[match(OTS2, mun_demo$STATE_NAME)]

mm_sta_r = match(resp$SG_UF_RESIDENCIA, sta_demo$state)
mm_r = match(resp$CO_MUNICIPIO_IBGE, mun_demo$IBGE_CODE)
resp$mun_pop = mun_demo$IBGE_MUN_POP[mm_r]
resp$sta_pop = sta_demo$pop[mm_sta_r]


# sinasc processing
sinasc$codanomal = trim.trailing(sinasc$codanomal)

# sinasc race coding
sinasc$racacor = as.character(sinasc$racacor)
sinasc$racacor[sinasc$racacor == "1"] = "BRANCA"
sinasc$racacor[sinasc$racacor == "2"] = "NEGRA"
sinasc$racacor[sinasc$racacor == "3"] = "AMARELA"
sinasc$racacor[sinasc$racacor == "4"] = "PARDA"
sinasc$racacor[sinasc$racacor == "5"] = "INDIGENA"
sinasc$racacor[sinasc$racacor == " "] = "IGNORADO"

# reassign non standard date formats in SINASC

date.reassign <- function(x){
  if(nchar(x) == 7){
    return(paste(paste0("0", substr(x, 1, 1)), substr(x, 2, 3), substr(x, 4, 7), sep = "/"))
  }else{
    return(paste(substr(x, 1, 2), substr(x, 3, 4), substr(x, 5, 8), sep = "/"))
  }
}

sinasc$dtnasc = sapply(sinasc$dtnasc, date.reassign)
sinasc$dtnasc = as.Date(sinasc$dtnasc, format = "%d/%m/%Y")

# trim SINASC to dates of RESP data (01 Jan 2015 to 24 May 2017)
sinasc = sinasc[sinasc$dtnasc > "2014-12-31", ]
sinasc = sinasc[sinasc$dtnasc < "2017-05-24", ]

resp = resp[as.Date(resp$DT_NOTIFICACAO, format = "%d/%m/%Y") > "2014-12-31", ]
resp = resp[as.Date(resp$DT_NOTIFICACAO, format = "%d/%m/%Y") < "2017-05-24", ]


# give a state column to SINASC
sinasc$codmunres = trim.trailing(sinasc$codmunres)
mm <- match(sinasc$codmunres, mun_demo$IBGE_CODE)
sinasc$IBGE_STATE_CODE = mun_demo$IBGE_STATE_CODE[mm]
# IBGE codes for missing states

# compile population for each state
sta_list <- unique(mun_demo$IBGE_STATE_CODE)
sta_demo = aggregate(mun_demo$IBGE_MUN_POP, by = list(mun_demo$IBGE_STATE_CODE), FUN = sum)
colnames(sta_demo) = c("state", "pop")

# assigning municipalities and population characteristics to resp and sinasc
sinasc$mun_pop = mun_demo$IBGE_MUN_POP[mm]
mm_sta = match(sinasc$IBGE_STATE_CODE, sta_demo$state)
sinasc$sta_pop = sta_demo$pop[mm_sta]

# trimmign of resp adn sinasc
useful1 <- c("TP_CLASSIFICACAO_FINAL", "DT_NASCIMENTO_RN", "CO_RACA_COR", "TP_SEXO", 
             "NU_IDADE_GESTANTE", "weeks_gestation", "CO_MUNICIPIO_IBGE", "SG_UF_RESIDENCIA", 
             "sta_pop", "mun_pop", "NU_PESO", "Pred_result", "NU_PERIMETRO_CEFALICO")
resp_useful <- resp[, useful1]
names(resp_useful)[5] = "IDADEMAE"
rm(resp)

colnames(sinasc) <- c("DTNASC", "SEXO", "RACACORMAE", "CODMUNRES", "IDADEMAE", 
                      "CODANOMAL","SEMAGESTAC", "PESO","Edu_years", "IBGE_STATE_CODE", "mun_pop", 
                      "sta_pop")
sinasc_useful <- sinasc
rm(sinasc)

# quick RESP analysis
resp_useful_inv <- resp_useful
resp_useful_inv$TP_CLASSIFICACAO_FINAL[resp_useful_inv$TP_CLASSIFICACAO_FINAL %in% c("PROVÁVEL", "SEM CLASSIFICAÇÃO")] = "EM INVESTIGAÇÃO"
resptab <- t(table(resp_useful_inv$TP_CLASSIFICACAO_FINAL, resp_useful_inv$SG_UF_RESIDENCIA))
resptab = data.frame(state = rownames(resptab),
                     conf_rate = as.numeric((resptab[, 1] + resptab[, 2]) / (resptab[, 1] + resptab[, 2] + resptab[, 3])))
resptab$Region = mun_demo$Region[match(resptab$state, mun_demo$IBGE_STATE_CODE)]  
resptab = resptab[order(resptab$Region), ]

# plotting
require(RColorBrewer)
colPAL <- brewer.pal(5, "Dark2")
pdf("/Users/eideobra/Dropbox/06_Brazil_MoH/Writeup/Comments_versions/NewFIGS/Microcephaly_tested.pdf", height = 5, width = 13)
plot(resptab[, 2], col = colPAL[as.numeric(resptab$Region)], pch = 19,
     ylim = c(0, 1), axes = F, xlab = "State", ylab = "Proportion tested", cex = 2)
axis(1, 1:27, resptab[, 1])
axis(2, seq(0, 1, 0.1), seq(0, 1, 0.1))
upper = mean(resptab[, 2], na.rm = T) + 1.96 * sd(resptab[, 2], na.rm = T)
lower = mean(resptab[, 2], na.rm = T) - 1.96 * sd(resptab[, 2], na.rm = T)
if(upper > 1){upper = 1}
if(lower < 0){lower = 0}
abline(h = mean(resptab[, 2], na.rm = T), col = "grey")
abline(h = upper, col = "grey", lty = 2)
abline(h = lower, col = "grey", lty = 2)
legend(4, 1.25, legend = as.character(unique(resptab[, 3])), fill = colPAL,
       horiz = T, cex = 1, xpd = T)
dev.off()

#######################
# what to dow ith NAs?
#######################


# leavign as is for now (i.e. omitted from analysis)
# sex - 1318 missing
sinasc_useful$SEXO[sinasc_useful$SEXO == "0"] = NA
sinasc_useful$SEXO[sinasc_useful$SEXO == "1"] = "M"
sinasc_useful$SEXO[sinasc_useful$SEXO == "2"] = "F"
#sinasc_useful$SEXO[sinasc_useful$SEXO == "I"] = sample(sinasc_useful$SEXO[sinasc_useful$SEXO != "I"],
#                                                       sum(sinasc_useful$SEXO == "I"))
resp_useful$TP_SEXO[resp_useful$TP_SEXO == "NO INFORMADO"] = NA
#resp_useful$TP_SEXO[resp_useful$TP_SEXO == "NO INFORMADO"] = sample(resp_useful$TP_SEXO[resp_useful$TP_SEXO != "NO INFORMADO"],
#                                                                    sum(resp_useful$TP_SEXO == "NO INFORMADO"))
# race - 291,983 + 4872 missing
sinasc_useful$RACACORMAE[sinasc_useful$RACACORMAE == "IGNORADO"] = NA
#sinasc_useful$RACACORMAE[sinasc_useful$RACACORMAE == "IGNORADO"] = sample(sinasc_useful$RACACORMAE[sinasc_useful$RACACORMAE != "IGNORADO"], 
#                                                                   sum(sinasc_useful$RACACORMAE == "IGNORADO"))
resp_useful$CO_RACA_COR[resp_useful$CO_RACA_COR == "SEM INFORMACAO"] = NA
#resp_useful$CO_RACA_COR[resp_useful$CO_RACA_COR == "SEM INFORMACAO"] = sample(resp_useful$CO_RACA_COR[resp_useful$CO_RACA_COR != "SEM INFORMACAO"], 
#                                                                              sum(resp_useful$CO_RACA_COR == "SEM INFORMACAO"))
# mother age - 67 missing
sinasc_useful$IDADEMAE = as.numeric(sinasc_useful$IDADEMAE)
#sinasc_useful$IDADEMAE[is.na(sinasc_useful$IDADEMAE)] = sample(sinasc_useful$IDADEMAE[!is.na(sinasc_useful$IDADEMAE)], 
#                                                               sum(is.na(sinasc_useful$IDADEMAE)))
#resp_useful$IDADEMAE[is.na(resp_useful$IDADEMAE)] = as.numeric(sample(resp_useful$IDADEMAE[!is.na(resp_useful$IDADEMAE)], 
#                                                                      sum(is.na(resp_useful$IDADEMAE))))
# gestaional age - 165,222 + 557 missing
sinasc_useful$SEMAGESTAC = as.numeric(sinasc_useful$SEMAGESTAC)
#sinasc_useful$SEMAGESTAC[is.na(sinasc_useful$SEMAGESTAC)] = sample(sinasc_useful$SEMAGESTAC[!is.na(sinasc_useful$SEMAGESTAC)], 
#                                                                   sum(is.na(sinasc_useful$SEMAGESTAC)))
#resp_useful$weeks_gestation[is.na(resp_useful$weeks_gestation)] = as.numeric(sample(resp_useful$weeks_gestation[!is.na(resp_useful$weeks_gestation)], 
#                                                                                    sum(is.na(resp_useful$weeks_gestation))))
# newborn weight - 11752 missing
sinasc_useful$PESO = as.numeric(as.character(sinasc_useful$PESO))
sinasc_useful$PESO[sinasc_useful$PESO < 500] = NA
#sinasc_useful$PESO[is.na(sinasc_useful$PESO)] = sample(sinasc_useful$PESO[!is.na(sinasc_useful$PESO)], 
#                                                       sum(is.na(sinasc_useful$PESO)))
resp_useful$NU_PESO = as.numeric(as.character(resp_useful$NU_PESO))
resp_useful$NU_PESO[resp_useful$NU_PESO < 500] = NA
#resp_useful$NU_PESO[is.na(resp_useful$NU_PESO)] = as.numeric(sample(resp_useful$NU_PESO[!is.na(resp_useful$NU_PESO)], 
#                                                                    sum(is.na(resp_useful$NU_PESO))))

# converting factor variables to numeric
sinasc_useful$IDADEMAE = as.numeric(as.character(sinasc_useful$IDADEMAE))
sinasc_useful$SEMAGESTAC = as.numeric(as.character(sinasc_useful$SEMAGESTAC))
sinasc_useful$PESO = as.numeric(as.character(sinasc_useful$PESO))
sinasc_useful$CODMUNRES = as.numeric(as.character(sinasc_useful$CODMUNRES))

## birth defect classifications
# load in a file with these
source("REFERENCE/BD_ICD10s.R")
# Neural Tube defects
#NTDs <- c("Q000", "Q001", "Q002", "Q050", "Q051", "Q052", "Q053", "Q054", "Q055", "Q056", "Q057", "Q058", "Q059",
#          "Q010", "Q011", "Q012", "Q013", "Q014", "Q015", "Q016", "Q017", "Q018", "Q019")
# Cariovascular defects
#CVDs <- paste("Q", 200:289, sep = "")
# abdominal wall defects
#AWDs <- c("Q792", "Q793", "Q794")
# Hypospadias and epispadias
#HAEs <- c("Q540", "Q541", "Q542", "Q543", "Q544", "Q545", "Q546", "Q547", "Q548", "Q549", "Q640")

# cleft lip and cleft palate
#CLCP <- c("Q035", "Q351", "Q353", "Q355", "Q357", "Q359", "Q036", "Q360", "Q361", "Q369", "Q037",
#          "Q370", "Q371", "Q372", "Q373", "Q374", "Q375", "Q378", "Q379")
# abnormal limbs
#Limbs <- c("Q749")


#### RESP-SINASC linkage ####
# i.e. removing records from sinasc that are updated by RESP
#source("REFERENCE/resp_sinasc_matcher_algo.R")
#sinasc_delete = resp.sinasc.matcher(resp_useful, sinasc_useful)
#save(sinasc_delete, file = "REFERENCE/RESP_SINASC_LINK_JUL18.R")
load("REFERENCE/RESP_SINASC_LINK_JUL18.R")
munLL <- read.csv("REFERENCE/MUN_lat_long.csv")
# matching summary diagnostics
conf_cases = (resp_useful$TP_CLASSIFICACAO_FINAL == "CONFIRMADO")[resp_useful$CO_MUNICIPIO_IBGE %in% munLL$IBGE]
# % exact match
sum((rowSums(sinasc_delete[, 2:3]) == 0)[conf_cases]) / nrow(sinasc_delete[conf_cases, ]) # 63% exact match
# % exact match within 30 days
sum(((sinasc_delete[, 3] <= 30) & sinasc_delete[, 2] <= 0)[conf_cases]) / nrow(sinasc_delete[conf_cases, ]) # 86% within 30 days
# % exact match within 30 days and close municipality
sum(((sinasc_delete[, 3] <= 30) & sinasc_delete[, 2] <= 50000)[conf_cases]) / nrow(sinasc_delete[conf_cases, ]) # 96% within 30 days and close municipality

# select within a week and close municipality as the cutoff option
sinasc_delete = sinasc_delete[((sinasc_delete[, 3] <= 30) & (sinasc_delete[, 2] <= 50000) & conf_cases), ]
# trimming duplicated records
sinasc_useful = sinasc_useful[!((1:nrow(sinasc_useful)) %in% sinasc_delete[, 1]), ]


#######################
# setting up new data frame
#######################


nbirths <- nrow(sinasc_useful) + nrow(resp_useful)
model_df <- data.frame(suspected_microcephaly = rep(NA, nbirths), # response
                       confirmed_microcephaly = rep(NA, nbirths), # response1B
                       notified_birth_defect = rep(NA, nbirths), # response 2
                       BD_BRAIN = rep(NA, nbirths),
                       BD_EYE = rep(NA, nbirths),
                       BD_MSK = rep(NA, nbirths),
                       BD_HAN = rep(NA, nbirths),
                       mun_sus_zik_incid_PREG = rep(NA, nbirths),
                       mun_con_zik_incid_PREG = rep(NA, nbirths),
                       mun_sus_den_incid_PREG = rep(NA, nbirths),
                       mun_con_den_incid_PREG = rep(NA, nbirths),
                       mun_sus_chik_incid_PREG = rep(NA, nbirths),
                       mun_con_chik_incid_PREG = rep(NA, nbirths),
                       mun_sus_den_incid_PRE12 = rep(NA, nbirths),
                       mun_con_den_incid_PRE12 = rep(NA, nbirths),
                       mun_sus_chik_incid_PRE12 = rep(NA, nbirths),
                       mun_con_chik_incid_PRE12 = rep(NA, nbirths),
                       YF_Vcov = rep(NA, nbirths),
                       week = rep(NA, nbirths),
                       mother_race = rep(NA, nbirths),
                       mother_age = rep(NA, nbirths),
                       baby_sex = rep(NA, nbirths),
                       case_mun = rep(NA, nbirths),
                       case_gaul_mun = rep(NA, nbirths),
                       case_sta = rep(NA, nbirths),
                       case_gaul_sta = rep(NA, nbirths),
                       mun_pop = rep(NA, nbirths),
                       sta_pop = rep(NA, nbirths),
                       weeks_gestation = rep(NA, nbirths),
                       weight = rep(NA, nbirths),
                       mun_HH_piped_w = rep(NA, nbirths),
                       mun_HH_garbage_col = rep(NA, nbirths),
                       mun_HH_sewage = rep(NA, nbirths),
                       mun_Fert_rate = rep(NA, nbirths),
                       mun_Income_pc = rep(NA, nbirths),
                       mun_Perc_rural = rep(NA, nbirths),
                       mun_Perc_employed = rep(NA, nbirths),
                       mun_Hospitalisations_pc = rep(NA, nbirths),
                       mun_cattle_den = rep(NA, nbirths),
                       mun_farm_den = rep(NA, nbirths),
                       mun_water_vul = rep(NA, nbirths),
                       mun_Perc_vaginal = rep(NA, nbirths),
                       mun_Pre_natal_cons = rep(NA, nbirths),
                       mun_Perc_single = rep(NA, nbirths),
                       mun_Mean_mother_age = rep(NA, nbirths),
                       mun_Mean_mother_edu = rep(NA, nbirths),
                       mun_Mean_baby_apgar5 = rep(NA, nbirths),
                       mun_Mean_years_edu = rep(NA, nbirths))

# Fulll exposure time series - onyl do 1 at a time
#Zika_sus_FES <- matrix(NA, nrow = nrow(model_df), ncol = 55)
Zika_con_FES <- matrix(NA, nrow = nrow(model_df), ncol = 55)
#Dengue_true_FES <- matrix(NA, nrow = nrow(model_df), ncol = 55)


#######################
# Adding values from SINASC
#######################
model_df$notified_birth_defect[1:nrow(sinasc_useful)] = as.numeric(sinasc_useful$CODANOMAL != "")
model_df$suspected_microcephaly[1:nrow(sinasc_useful)] = 0
model_df$confirmed_microcephaly[1:nrow(sinasc_useful)] = 0

multi.BD.scan <- function(BDs, codestring){
  logmat = matrix(NA, nrow = length(codestring), ncol = length(BDs))
  for(i in 1:length(BDs)){
    logmat[, i] = grepl(BDs[i], codestring)
    #print(i)
  }
  logvec = rowSums(logmat)
  return(logvec)
}

model_df$BD_BRAIN[1:nrow(sinasc_useful)] = multi.BD.scan(BRAIN, sinasc_useful$CODANOMAL)
model_df$BD_EYE[1:nrow(sinasc_useful)] = multi.BD.scan(EYE, sinasc_useful$CODANOMAL)
model_df$BD_MSK[1:nrow(sinasc_useful)] = multi.BD.scan(MSK, sinasc_useful$CODANOMAL)
model_df$BD_HAN[1:nrow(sinasc_useful)] = multi.BD.scan(HAN, sinasc_useful$CODANOMAL)
#model_df$BD_MSK[1:nrow(sinasc_useful)] = as.numeric(sinasc_useful$CODANOMAL %in% MSK)
#model_df$BD_CNS[1:nrow(sinasc_useful)] = as.numeric(sinasc_useful$CODANOMAL %in% CNS)
#model_df$BD_PNS[1:nrow(sinasc_useful)] = as.numeric(sinasc_useful$CODANOMAL %in% PNS)
#model_df$BD_EYE[1:nrow(sinasc_useful)] = as.numeric(sinasc_useful$CODANOMAL %in% EYE)
#model_df$BD_CLCP[1:nrow(sinasc_useful)] = as.numeric(sinasc_useful$CODANOMAL %in% CLCP)
#model_df$BD_CVDs[1:nrow(sinasc_useful)] = as.numeric(sinasc_useful$CODANOMAL %in% CVDs)
#model_df$BD_AWDs[1:nrow(sinasc_useful)] = as.numeric(sinasc_useful$CODANOMAL %in% AWDs)
#model_df$BD_HAK[1:nrow(sinasc_useful)] = as.numeric(sinasc_useful$CODANOMAL %in% HAK)
#model_df$BD_LIMB[1:nrow(sinasc_useful)] = as.numeric(sinasc_useful$CODANOMAL %in% LIMB)

model_df$mother_race[1:nrow(sinasc_useful)] = as.character(sinasc_useful$RACACORMAE)
model_df$mother_age[1:nrow(sinasc_useful)] = sinasc_useful$IDADEMAE
model_df$baby_sex[1:nrow(sinasc_useful)] = as.character(sinasc_useful$SEXO)
model_df$week[1:nrow(sinasc_useful)] = floor(as.numeric(as.Date(sinasc_useful$DTNASC, format = "%d%m%Y") - as.Date("2015-01-01")) / 7)
model_df$case_mun[1:nrow(sinasc_useful)] = sinasc_useful$CODMUNRES
model_df$case_gaul_mun[1:nrow(sinasc_useful)] = mun_demo$GAUL_CODE[match(sinasc_useful$CODMUNRES, mun_demo$IBGE_CODE)]
model_df$case_sta[1:nrow(sinasc_useful)] = as.character(sinasc_useful$IBGE_STATE_CODE)
model_df$case_gaul_sta[1:nrow(sinasc_useful)] = sta_gaul$GAUL_CODE[match(sinasc_useful$IBGE_STATE_CODE, sta_gaul$State_IBGE_ID)]
model_df$mun_pop[1:nrow(sinasc_useful)] = sinasc_useful$mun_pop
model_df$sta_pop[1:nrow(sinasc_useful)] = sinasc_useful$sta_pop
model_df$weeks_gestation[1:nrow(sinasc_useful)] = sinasc_useful$SEMAGESTAC
model_df$weight[1:nrow(sinasc_useful)] = sinasc_useful$PESO

# adding household characteristics
HH_match = match(model_df$case_mun[1:nrow(sinasc_useful)], sociodemo_useful$IBGE_mun)
model_df$mun_HH_piped_w[1:nrow(sinasc_useful)] = sociodemo_useful[HH_match, 2]
model_df$mun_HH_garbage_col[1:nrow(sinasc_useful)] = sociodemo_useful[HH_match, 3]
model_df$mun_HH_sewage[1:nrow(sinasc_useful)] = sociodemo_useful[HH_match, 4]
model_df$mun_Fert_rate[1:nrow(sinasc_useful)] = sociodemo_useful[HH_match, 5]
model_df$mun_Income_pc[1:nrow(sinasc_useful)] = sociodemo_useful[HH_match, 6]
model_df$mun_Perc_rural[1:nrow(sinasc_useful)] = sociodemo_useful[HH_match, 7]
model_df$mun_Perc_employed[1:nrow(sinasc_useful)] = sociodemo_useful[HH_match, 8]
model_df$mun_Hospitalisations_pc[1:nrow(sinasc_useful)] = sociodemo_useful[HH_match, 9]
model_df$mun_Mean_years_edu[1:nrow(sinasc_useful)] = sociodemo_useful[HH_match, 10]

# cattle
cat_match = match(model_df$case_mun[1:nrow(sinasc_useful)], cattle$municipality_code)
model_df$mun_cattle_den[1:nrow(sinasc_useful)] = cattle[cat_match, 6] / cattle[cat_match, 7]
model_df$mun_farm_den[1:nrow(sinasc_useful)] = cattle[cat_match, 5] / cattle[cat_match, 7]

# adding water vulnerability
wat_match = match(model_df$case_mun[1:nrow(sinasc_useful)], water_vul$COD_IBGE)
wat_year_match = 2015 + floor(model_df$week[1:nrow(sinasc_useful)] / 52)
wat_year_match[wat_year_match == 2017] = 2016
model_df$mun_water_vul[1:nrow(sinasc_useful)] = apply(cbind(wat_year_match, water_vul[wat_match, 5:7]), 1,
                                                      function(x) if(is.na(x[1])){NA}else{x[(x[1] - 2012)]})

maternal_match = match(model_df$case_mun[1:nrow(sinasc_useful)], matchar$IBGE_mun)
model_df$mun_Perc_vaginal[1:nrow(sinasc_useful)] = matchar[maternal_match, 3]
model_df$mun_Pre_natal_cons[1:nrow(sinasc_useful)] = matchar[maternal_match, 4]
model_df$mun_Perc_single[1:nrow(sinasc_useful)] = matchar[maternal_match, 5]
model_df$mun_Mean_mother_age[1:nrow(sinasc_useful)] = matchar[maternal_match, 6]
model_df$mun_Mean_mother_edu[1:nrow(sinasc_useful)] = matchar[maternal_match, 7]
model_df$mun_Mean_baby_apgar5[1:nrow(sinasc_useful)] = matchar[maternal_match, 8]


# Yellow fever vaccine coverage - assume if not matched VC = 0% (outside of YF risk area)
model_df$YF_Vcov[1:nrow(sinasc_useful)] = YF_Vcov[match(model_df$case_mun[1:nrow(sinasc_useful)], YF_Vcov$IBGE_mun), "YF_VC_2016"]

# matching matrices for zika and arboviruses to speed up calculation
match_mat <- data.frame(zik_sus_mun_match = match(model_df$case_mun[1:nrow(sinasc_useful)], zik_sus_mun$IBGE_Municipality),
                        zik_con_mun_match = match(model_df$case_mun[1:nrow(sinasc_useful)], zik_con_mun$IBGE_Municipality))

match_mat_d <- data.frame(den_sus_mun_match = match(model_df$case_mun[1:nrow(sinasc_useful)], den_sus_mun$IBGE_Municipality),
                          den_con_mun_match = match(model_df$case_mun[1:nrow(sinasc_useful)], den_con_mun$IBGE_Municipality))
                          
match_mat_c <- data.frame(chik_sus_mun_match = match(model_df$case_mun[1:nrow(sinasc_useful)], chik_sus_mun$IBGE_Municipality),
                          chik_con_mun_match = match(model_df$case_mun[1:nrow(sinasc_useful)], chik_con_mun$IBGE_Municipality))

# convert character state names to numeric so the whole objects can be easily coerced to matrices
#zik_true_sta$State = as.numeric(zik_true_sta$State)
#den_not_sta$State = as.numeric(den_not_sta$State)
#chi_not_sta$State = as.numeric(chi_not_sta$State)


# trimester extractor function

TRI1.exp.match <- function(x, ts_mat){
  # x includes match place of municipality or state in ts_mat,
  # time of notification in weeks since the begining of 01 Jan 2015
  # gestational age at notification
  b <- as.numeric(x)
  # only proceed if there is a municipality match
  if(!is.na(b[1])){
    # if gestational age is missing assume 40 weeks
    if(is.na(b[3])){gesage = 40}else{gesage = b[3]}
    # work out date of becoming pregnant (in weeks) and end of 1st trimester
    pregweek = b[2] - gesage
    TRI1end = pregweek + 14
    # fill in data on exposure over the course of the 1st trimester
    expo = NA
    if(TRI1end > 0){
      ts = as.numeric(ts_mat[b[1], 4:159])
      if(pregweek < 1){pregweek = 1}
      ts = ts[pregweek:TRI1end]
      if(all(is.na(ts))){expo = NA}else{expo = mean(ts, na.rm = T)}
    }
    return(expo)
  }else{return(NA)}
}

# all pregnancy exposure function
ALLTRI.exp.match <- function(x, ts_mat){
  # x includes match place of municipality or state in ts_mat,
  # time of notification in weeks since the begining of 01 Jan 2015
  # gestational age at notification
  b <- as.numeric(x)
  # only proceed if there is a municipality match
  if(!is.na(b[1])){
    # if gestational age is missing assume 40 weeks
    if(is.na(b[3])){gesage = 40}else{gesage = b[3]}
    # work out date of becoming pregnant (in weeks) 
    pregweek = b[2] - gesage
    birthweek = b[2]
    #TRI1end = pregweek + 14
    # fill in data on exposure over the course of the pregnancy
    expo = NA
    if(birthweek > 0){
      ts = as.numeric(ts_mat[b[1], 4:159])
      if(pregweek < 1){pregweek = 1}
      ts = ts[pregweek:birthweek]
      if(all(is.na(ts))){expo = NA}else{expo = mean(ts, na.rm = T)}
    }
    return(expo)
  }else{return(NA)}
}
# and pre 1 months preganncy extractor function
PRE1.exp.match <- function(x, ts_mat){
  # x includes match place of municipality or state in ts_mat,
  # time of notification in weeks since the begining of 01 Jan 2015
  # gestational age at notification
  b <- as.numeric(x)
  # only proceed if there is a municipality match
  if(!is.na(b[1])){
    # if gestational age is missing assume 40 weeks
    if(is.na(b[3])){gesage = 40}else{gesage = b[3]}
    # work out date of becoming pregnant (in weeks) and end of 1st trimester
    pregweek = b[2] - gesage
    PRE1start = pregweek - 4
    # fill in data on exposure over the course of the 6 months before becoming pregnant
    expo = NA
    if(pregweek > 0){
      ts = as.numeric(ts_mat[b[1], 4:159])
      if(PRE1start < 1){PRE1start = 1}
      ts = ts[PRE1start:pregweek]
      if(all(is.na(ts))){expo = NA}else{expo = mean(ts, na.rm = T)}
    }
    return(expo)
  }else{return(NA)}
}

# and pre 12 months preganncy extractor function
PRE12.exp.match <- function(x, ts_mat){
  # x includes match place of municipality or state in ts_mat,
  # time of notification in weeks since the begining of 01 Jan 2015
  # gestational age at notification
  b <- as.numeric(x)
  # only proceed if there is a municipality match
  if(!is.na(b[1])){
    # if gestational age is missing assume 40 weeks
    if(is.na(b[3])){gesage = 40}else{gesage = b[3]}
    # work out date of becoming pregnant (in weeks) and end of 1st trimester
    pregweek = b[2] - gesage
    PRE6start = pregweek - 52
    # fill in data on exposure over the course of the 12 months before becoming pregnant
    expo = NA
    if(pregweek > 0){
      ts = as.numeric(ts_mat[b[1], 4:159])
      if(PRE6start < 1){PRE6start = 1}
      ts = ts[PRE6start:pregweek]
      if(all(is.na(ts))){expo = NA}else{expo = mean(ts, na.rm = T)}
    }
    return(expo)
  }else{return(NA)}
}

FES.exp.match <- function(x, ts_mat){
  # x includes match place of municipality or state in ts_mat,
  # time of notification in weeks since the begining of 01 Jan 2015
  # gestational age at notification
  b <- as.numeric(x)
  # template for exposure
  template = rep(NA, 55)
  # only proceed if there is a municipality match
  
  if(!is.na(b[1])){
    # if gestational age is missing assume 40 weeks
    if(is.na(b[3])){gesage = 40}else{gesage = b[3]}
    if(gesage > 45){
      gesage = 40
      b[2] = b[2] - (b[3] - 40)
     }
    if(b[2] > 0){
      # work out date of becoming pregnant (in weeks) and end of 1st trimester
      pregweek = b[2] - gesage
      birthweek = b[2]
      # 10 weeks before pregnancy week
      pweekm10 = pregweek - 10
      # date 45 weeks after getting pregnant
      pweek45 = pregweek + 45
      
      if(pregweek < 1){pregweek2 = 1}else{pregweek2 = pregweek}
      if(pweekm10 < 1){pweekm102 = 1}else{pweekm102 = pweekm10}
      if(pweek45 > 160){pweek452 = 160}else{pweek452 = pweek45}
      # fill in data on exposure over the course of the pregnancy
      ts = as.numeric(ts_mat[b[1], 4:159])
      
      if(birthweek > 1){
        #expo = ts[pweekm102:(birthweek - 1)]
        #template[(gesage + 10 - length(expo) + 1):(gesage + 10)] = expo
        expo = ts[pweekm102:(pweek452 - 1)]
        template[(pweekm102 - pweekm10 + 1):length(template)] = expo
        
        # chaning code for post birth
        #pbw = (gesage + 11)
        #if(pbw < length(template)){template[(gesage + 11):length(template)] = -888}
        # now picked up in post-birth analysis
      }
      return(template)
    }else{return(template)}
  }else{return(template)}
}


### ZIKA ###
# Zika exposure during pregnancy - Suspected
placeholder = apply(data.frame(match_place = match_mat$zik_sus_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(zik_sus_mun))
model_df$mun_sus_zik_incid_PREG[1:nrow(sinasc_useful)] = placeholder

# Zika exposure during pregnancy - Confirmed
placeholder = apply(data.frame(match_place = match_mat$zik_con_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(zik_con_mun))
model_df$mun_con_zik_incid_PREG[1:nrow(sinasc_useful)] = placeholder


### DENGUE ###

# dengue exposure in pregnancy - suspected
placeholder = apply(data.frame(match_place = match_mat_d$den_sus_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(den_sus_mun))
model_df$mun_sus_den_incid_PREG[1:nrow(sinasc_useful)] = placeholder

# dengue exposure in pregnancy - confirmed
placeholder = apply(data.frame(match_place = match_mat_d$den_con_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(den_con_mun))
model_df$mun_con_den_incid_PREG[1:nrow(sinasc_useful)] = placeholder

# dengue exposure BEFORE pregnancy - suspected
placeholder = apply(data.frame(match_place = match_mat_d$den_sus_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, PRE12.exp.match, ts_mat = as.matrix(den_sus_mun))
model_df$mun_sus_den_incid_PRE12[1:nrow(sinasc_useful)] = placeholder

# dengue exposure BEFORE pregnancy - confirmed
placeholder = apply(data.frame(match_place = match_mat_d$den_con_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, PRE12.exp.match, ts_mat = as.matrix(den_con_mun))
model_df$mun_con_den_incid_PRE12[1:nrow(sinasc_useful)] = placeholder



### CHIKUNGUNYA ###

# chikungunya exposure in pregnancy - suspected
placeholder = apply(data.frame(match_place = match_mat_c$chik_sus_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(chik_sus_mun))
model_df$mun_sus_chik_incid_PREG[1:nrow(sinasc_useful)] = placeholder

# chikungunya exposure in pregnancy - confirmed
placeholder = apply(data.frame(match_place = match_mat_c$chik_con_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(chik_con_mun))
model_df$mun_con_chik_incid_PREG[1:nrow(sinasc_useful)] = placeholder

# chikungunya exposure BEFORE pregnancy - suspected
placeholder = apply(data.frame(match_place = match_mat_c$chik_sus_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, PRE12.exp.match, ts_mat = as.matrix(chik_sus_mun))
model_df$mun_sus_chik_incid_PRE12[1:nrow(sinasc_useful)] = placeholder

# chikungunya exposure BEFORE pregnancy - confirmed
placeholder = apply(data.frame(match_place = match_mat_c$chik_con_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, PRE12.exp.match, ts_mat = as.matrix(chik_con_mun))
model_df$mun_con_chik_incid_PRE12[1:nrow(sinasc_useful)] = placeholder





###### Indiviudal full exposure histories
###########
# zika suspected and confirmed
placeholder = apply(data.frame(match_place = match_mat$zik_sus_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, FES.exp.match, ts_mat = as.matrix(zik_sus_mun))
Zika_sus_FES[1:nrow(sinasc_useful), ] = t(placeholder)

placeholder = apply(data.frame(match_place = match_mat$zik_con_mun_match, 
                               time = model_df$week[1:nrow(sinasc_useful)],
                               ges_age = model_df$weeks_gestation[1:nrow(sinasc_useful)]), 
                    1, FES.exp.match, ts_mat = as.matrix(zik_con_mun))
Zika_con_FES[1:nrow(sinasc_useful), ] = t(placeholder)







#######################
# PART B- RESP
#######################

# assign vectorable values
rindex = (nrow(sinasc_useful) + 1):nrow(model_df)
model_df$suspected_microcephaly[rindex] = 1
model_df$notified_birth_defect[rindex] = 0
model_df$confirmed_microcephaly[rindex] = as.numeric(resp_useful$TP_CLASSIFICACAO_FINAL == "CONFIRMADO")

# NA or 0 - RESP doesn't mention any other birth defects, but doesn;t mean they are not there

model_df$BD_BRAIN[rindex] = NA
model_df$BD_EYE[rindex] = NA
model_df$BD_MSK[rindex] = NA
model_df$BD_HAN[rindex] = NA
#model_df$BD_MSK[rindex] = NA
#model_df$BD_CNS[rindex] = NA
#model_df$BD_PNS[rindex] = NA
#model_df$BD_EYE[rindex] = NA
#model_df$BD_CLCP[rindex] = NA
#model_df$BD_CVDs[rindex] = NA
#model_df$BD_AWDs[rindex] = NA
#model_df$BD_HAK[rindex] = NA
#model_df$BD_LIMB[rindex] = NA

# intergrowth MC reclassification
#IGtable <- read.csv("REFERENCE/INTERGROWTH_head.csv")
#IG.micro.classifier <- function(ages, sexes, headCs, IGtable){
#  
#  Zscores <- rep(NA, length(ages))
#  for(i in 1:length(ages)){
#    if(!any(c(is.na(ages[i]), is.na(headCs[i])))){
#      # restrict to age range in INTERGROWTH - 24 to 42 weeks
#      if(ages[i] < 24){ages[i] = 24}
#      if(ages[i] > 42){ages[i] = 42}
#      # indexing
#      if(!(headCs[i] %in% IGtable$HeadCircumference)){
#        close <- which.min(sqrt((IGtable$HeadCircumference - headCs[i])^2))
#        headCs[i] = headCs[i] + (IGtable$HeadCircumference[close] - headCs[i])
#      }
#      temp = IGtable[IGtable$HeadCircumference == headCs[i], ]
#      temp = temp[temp$Age == (ages[i]*7 + 3), ]
#      if(sexes[i] == "FEMININO"){Zscores[i] = as.numeric(as.character(temp[temp$Sex == "Female", 11]))}else{
#        Zscores[i] = as.numeric(as.character(temp[temp$Sex == "Male", 11]))
#      }
#    }
#  }
#  # returns vector of microcepahly < -2, severe microcephaly < -3
#  return(cbind(as.numeric(Zscores < -2), as.numeric(Zscores < -3)))
#}

#placeholder = IG.micro.classifier(ages = resp_useful$weeks_gestation,
#                                  sexes = as.character(resp_useful$TP_SEXO),
#                                  headCs = resp_useful$NU_PERIMETRO_CEFALICO,
#                                  IGtable = IGtable)
#model_df$IG_microcephaly[rindex] = placeholder[, 1]
#model_df$IG_severe_microcephaly[rindex] = placeholder[, 2]

model_df$mother_race[rindex] = as.character(resp_useful$CO_RACA_COR)
model_df$mother_age[rindex] = resp_useful$IDADEMAE
model_df$baby_sex[rindex] = as.character(resp_useful$TP_SEXO)
model_df$week[rindex] = floor(as.numeric(as.Date(resp_useful$DT_NASCIMENTO_RN) - as.Date("2015-01-01")) / 7)
model_df$case_mun[rindex] = resp_useful$CO_MUNICIPIO_IBGE
model_df$case_gaul_mun[rindex] = mun_demo$GAUL_CODE[match(resp_useful$CO_MUNICIPIO_IBGE, mun_demo$IBGE_CODE)]
model_df$case_sta[rindex] = as.character(resp_useful$SG_UF_RESIDENCIA)
model_df$case_gaul_sta[rindex] = sta_gaul$GAUL_CODE[match(resp_useful$SG_UF_RESIDENCIA, sta_gaul$State_IBGE_ID)]
model_df$mun_pop[rindex] = resp_useful$mun_pop
model_df$sta_pop[rindex] = resp_useful$sta_pop
model_df$weeks_gestation[rindex] = resp_useful$weeks_gestation
model_df$weight[rindex] = resp_useful$NU_PESO

# adding household characteristics -RESP
HH_match = match(model_df$case_mun[rindex], sociodemo_useful$IBGE_mun)
model_df$mun_HH_piped_w[rindex] = sociodemo_useful[HH_match, 2]
model_df$mun_HH_garbage_col[rindex] = sociodemo_useful[HH_match, 3]
model_df$mun_HH_sewage[rindex] = sociodemo_useful[HH_match, 4]
model_df$mun_Fert_rate[rindex] = sociodemo_useful[HH_match, 5]
model_df$mun_Income_pc[rindex] = sociodemo_useful[HH_match, 6]
model_df$mun_Perc_rural[rindex] = sociodemo_useful[HH_match, 7]
model_df$mun_Perc_employed[rindex] = sociodemo_useful[HH_match, 8]
model_df$mun_Hospitalisations_pc[rindex] = sociodemo_useful[HH_match, 9]
model_df$mun_Mean_years_edu[rindex] = sociodemo_useful[HH_match, 10]

# Cattle
cat_match = match(model_df$case_mun[rindex], cattle$municipality_code)
model_df$mun_cattle_den[rindex] = cattle[cat_match, 6] / cattle[cat_match, 7]
model_df$mun_farm_den[rindex] = cattle[cat_match, 5] / cattle[cat_match, 7]

# adding water vulnerability
wat_match = match(model_df$case_mun[rindex], water_vul$COD_IBGE)
wat_year_match = 2015 + floor(model_df$week[rindex] / 52)
wat_year_match[wat_year_match == 2017] = 2016
model_df$mun_water_vul[rindex] = apply(cbind(wat_year_match, water_vul[wat_match, 5:7]), 1,
                                       function(x) if(is.na(x[1])){NA}else{x[(x[1] - 2012)]})

# adding maternal characteristics
maternal_match = match(model_df$case_mun[rindex], matchar$IBGE_mun)
model_df$mun_Perc_vaginal[rindex] = matchar[maternal_match, 3]
model_df$mun_Pre_natal_cons[rindex] = matchar[maternal_match, 4]
model_df$mun_Perc_single[rindex] = matchar[maternal_match, 5]
model_df$mun_Mean_mother_age[rindex] = matchar[maternal_match, 6]
model_df$mun_Mean_mother_edu[rindex] = matchar[maternal_match, 7]
model_df$mun_Mean_baby_apgar5[rindex] = matchar[maternal_match, 8]

# YF vaccine coverage - assume NA = 0% as outside of risk area
model_df$YF_Vcov[rindex] = YF_Vcov[match(model_df$case_mun[rindex], YF_Vcov$IBGE_mun), "YF_VC_2016"]
model_df$YF_Vcov[is.na(model_df$YF_Vcov)] = 0

# matching matrix to link RESP with timeseries from GAL
match_mat_r <- data.frame(zik_sus_mun_match = match(model_df$case_mun[rindex], zik_sus_mun$IBGE_Municipality),
                          zik_con_mun_match = match(model_df$case_mun[rindex], zik_con_mun$IBGE_Municipality))

match_mat_r_d <- data.frame(den_sus_mun_match = match(model_df$case_mun[rindex], den_sus_mun$IBGE_Municipality),
                            den_con_mun_match = match(model_df$case_mun[rindex], den_con_mun$IBGE_Municipality))

match_mat_r_c <- data.frame(chik_sus_mun_match = match(model_df$case_mun[rindex], chik_sus_mun$IBGE_Municipality),
                            chik_con_mun_match = match(model_df$case_mun[rindex], chik_con_mun$IBGE_Municipality))

#######################
# Exposure variables for microcephaly cases
#######################





### ZIKA ###
# Zika exposure during pregnancy - Suspected
placeholder = apply(data.frame(match_place = match_mat_r$zik_sus_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(zik_sus_mun))
model_df$mun_sus_zik_incid_PREG[rindex] = placeholder

# Zika exposure during pregnancy - Confirmed
placeholder = apply(data.frame(match_place = match_mat_r$zik_con_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(zik_con_mun))
model_df$mun_con_zik_incid_PREG[rindex] = placeholder


### DENGUE ###

# dengue exposure in pregnancy - suspected
placeholder = apply(data.frame(match_place = match_mat_r_d$den_sus_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(den_sus_mun))
model_df$mun_sus_den_incid_PREG[rindex] = placeholder

# dengue exposure in pregnancy - confirmed
placeholder = apply(data.frame(match_place = match_mat_r_d$den_con_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(den_con_mun))
model_df$mun_con_den_incid_PREG[rindex] = placeholder

# dengue exposure BEFORE pregnancy - suspected
placeholder = apply(data.frame(match_place = match_mat_r_d$den_sus_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, PRE12.exp.match, ts_mat = as.matrix(den_sus_mun))
model_df$mun_sus_den_incid_PRE12[rindex] = placeholder

# dengue exposure BEFORE pregnancy - confirmed
placeholder = apply(data.frame(match_place = match_mat_r_d$den_con_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, PRE12.exp.match, ts_mat = as.matrix(den_con_mun))
model_df$mun_con_den_incid_PRE12[rindex] = placeholder



### CHIKUNGUNYA ###

# chikungunya exposure in pregnancy - suspected
placeholder = apply(data.frame(match_place = match_mat_r_c$chik_sus_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(chik_sus_mun))
model_df$mun_sus_chik_incid_PREG[rindex] = placeholder

# chikungunya exposure in pregnancy - confirmed
placeholder = apply(data.frame(match_place = match_mat_r_c$chik_con_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, ALLTRI.exp.match, ts_mat = as.matrix(chik_con_mun))
model_df$mun_con_chik_incid_PREG[rindex] = placeholder

# chikungunya exposure BEFORE pregnancy - suspected
placeholder = apply(data.frame(match_place = match_mat_r_c$chik_sus_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, PRE12.exp.match, ts_mat = as.matrix(chik_sus_mun))
model_df$mun_sus_chik_incid_PRE12[rindex] = placeholder

# chikungunya exposure BEFORE pregnancy - confirmed
placeholder = apply(data.frame(match_place = match_mat_r_c$chik_con_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, PRE12.exp.match, ts_mat = as.matrix(chik_con_mun))
model_df$mun_con_chik_incid_PRE12[rindex] = placeholder



###### Indiviudal full exposure histories - RESP
##########
# zika suspected and confirmed
placeholder = apply(data.frame(match_place = match_mat_r$zik_sus_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, FES.exp.match, ts_mat = as.matrix(zik_sus_mun))
Zika_sus_FES[rindex, ] = t(placeholder)

placeholder = apply(data.frame(match_place = match_mat_r$zik_con_mun_match, 
                               time = model_df$week[rindex],
                               ges_age = model_df$weeks_gestation[rindex]), 
                    1, FES.exp.match, ts_mat = as.matrix(zik_con_mun))
Zika_con_FES[rindex, ] = t(placeholder)




#################
# final raw processing
#################

#model_df = model_df[, !(colnames(model_df) %in% c("weight", "weeks_gestation"))]

# cases to incidence
model_df$mun_sus_zik_incid_PREG = model_df$mun_sus_zik_incid_PREG / model_df$mun_pop
model_df$mun_con_zik_incid_PREG = model_df$mun_con_zik_incid_PREG / model_df$mun_pop

model_df$mun_sus_den_incid_PREG = model_df$mun_sus_den_incid_PREG / model_df$mun_pop
model_df$mun_con_den_incid_PREG = model_df$mun_con_den_incid_PREG / model_df$mun_pop
model_df$mun_sus_den_incid_PRE12 = model_df$mun_sus_den_incid_PRE12 / model_df$mun_pop
model_df$mun_con_den_incid_PRE12 = model_df$mun_con_den_incid_PRE12 / model_df$mun_pop

model_df$mun_sus_chik_incid_PREG = model_df$mun_sus_chik_incid_PREG / model_df$mun_pop
model_df$mun_con_chik_incid_PREG = model_df$mun_con_chik_incid_PREG / model_df$mun_pop
model_df$mun_sus_chik_incid_PRE12 = model_df$mun_sus_chik_incid_PRE12 / model_df$mun_pop
model_df$mun_con_chik_incid_PRE12 = model_df$mun_con_chik_incid_PRE12 / model_df$mun_pop


Zika_sus_FES = Zika_sus_FES / model_df$mun_pop
#Zika_sus_FES[Zika_sus_FES < 0] = -888

Zika_con_FES = Zika_con_FES / model_df$mun_pop
#Zika_con_FES[Zika_con_FES < 0] = -888


#################
# save main file
#################

save(model_df, file = "INDIVIDUAL/Indiv_reg_df_NOV18_170718.RData")


# Add FID column
Zika_sus_FES = cbind(1:nrow(Zika_sus_FES), Zika_sus_FES)
save(Zika_sus_FES, file = "INDIVIDUAL/FULL_TS/Zika_sus_FES_NOV18.RData")

Zika_con_FES = cbind(1:nrow(Zika_con_FES), Zika_con_FES)
save(Zika_con_FES, file = "INDIVIDUAL/FULL_TS/Zika_con_FES_NOV18.RData")












########################
# birth defect detailed summaries
########################
# Individual list of birth defects
bd_list = strsplit(sinasc_useful$CODANOMAL, "Q")
bd_list = c(bd_list, as.list(rep("", nrow(resp_useful))))
all_bds = unique(unlist(bd_list))
# get rid of X and D codes
bd_list = lapply(bd_list, function(x) x[!grepl("X", x)])
bd_list = lapply(bd_list, function(x) x[!grepl("D", x)])
all_bds = all_bds[!grepl("D", all_bds)]
all_bds = all_bds[!grepl("X", all_bds)]

# birth defect frequency
bd_freq = table(unlist(bd_list))
# which ones are in our main categories of BDs?
source("REFERENCE/BD_ICD10s.R")
# macro categories are:
# 1 Brain
# 2 Eye
# 3 musculoskeletal incl. limb defects
# 4 head and neck
# 
#All_catBD <- c(BRAIN, EYE, MSK, HAN)
#All_catBD = gsub("Q", "", All_catBD)
#bd_freq = bd_freq[names(bd_freq) %in% All_catBD]

# prelim trim - needs to have a frequency of at least 2 per 10,000 births
2 * nrow(model_df) / 10000 # 1372
#All_catBD = All_catBD[All_catBD %in% names(bd_freq)[bd_freq >= 343]]

#all_bds = all_bds[all_bds != ""]
#all_bds = all_bds[all_bds %in% All_catBD]

# or don't restrict to within category
all_bds = all_bds[all_bds != ""]
spec_bds = all_bds[all_bds %in% names(bd_freq)[bd_freq >= 1372]]
# all_bds = all_bds[all_bds %in% gsub("Q", "", BRAIN)]






# set up a matrix
BD_mat <- matrix(0, nrow = nbirths, ncol = length(spec_bds))

for(i in 1:length(spec_bds)){
  # focus birth defect
  fBD <- spec_bds[i]
  # which mothers is it in?
  mother_fBD = unlist(lapply(bd_list, function(x) fBD %in% x))
  # assign to matrix
  BD_mat[mother_fBD, i] = 1
  print(i)
}

colnames(BD_mat) = spec_bds
save(BD_mat, file = "INDIVIDUAL/Detailed_birth_defects_NOV18.RData")



#######################
# investigation only

# delving deper into brain and eye defects
# brain first
load("INDIVIDUAL/Detailed_birth_defects.RData")
brain_bds <- BRAIN
BD_mat <- matrix(0, nrow = nbirths, ncol = length(brain_bds))

for(i in 1:length(brain_bds)){
  # focus birth defect
  fBD <- brain_bds[i]
  fBD <- gsub("Q", "", fBD)
  # which mothers is it in?
  mother_fBD = unlist(lapply(bd_list, function(x) fBD %in% x))
  # assign to matrix
  BD_mat[mother_fBD, i] = 1
  print(i)
}
colnames(BD_mat) = brain_bds

load("INDIVIDUAL/Indiv_reg_df_JUL18_170718.RData")

BD_mat = BD_mat[!is.na(model_df$mun_con_zik_incid_PREG), ]
model_df = model_df[!is.na(model_df$mun_con_zik_incid_PREG), ]

model_df = data.frame(Zika = model_df$mun_con_zik_incid_PREG,
                      BD_mat)
expo = model_df[model_df$Zika != min(model_df$Zika, na.rm = T), ]
noexpo = model_df[model_df$Zika == min(model_df$Zika, na.rm = T), ]

BDcomp = data.frame(BD = colnames(model_df)[2:71],
                    Zikrate = apply(expo[, 2:71], 2, sum, na.rm = T) / nrow(expo),
                    noZikrate = apply(noexpo[, 2:71], 2, sum, na.rm = T)  / nrow(noexpo))

BDcomp$Zikrate = BDcomp$Zikrate * 10000
BDcomp$noZikrate = BDcomp$noZikrate * 10000
BDcomp$diff = BDcomp$Zikrate / BDcomp$noZikrate
BDcomp[rev(order(BDcomp$diff)), ]

# now eye
eye_bds <- EYE
BD_mat <- matrix(0, nrow = nbirths, ncol = length(eye_bds))

for(i in 1:length(eye_bds)){
  # focus birth defect
  fBD <- eye_bds[i]
  fBD <- gsub("Q", "", fBD)
  # which mothers is it in?
  mother_fBD = unlist(lapply(bd_list, function(x) fBD %in% x))
  # assign to matrix
  BD_mat[mother_fBD, i] = 1
  print(i)
}
colnames(BD_mat) = eye_bds

load("INDIVIDUAL/Indiv_reg_df_JUL18_170718.RData")

BD_mat = BD_mat[!is.na(model_df$mun_con_zik_incid_PREG), ]
model_df = model_df[!is.na(model_df$mun_con_zik_incid_PREG), ]

model_df = data.frame(Zika = model_df$mun_con_zik_incid_PREG,
                      BD_mat)
expo = model_df[model_df$Zika != min(model_df$Zika, na.rm = T), ]
noexpo = model_df[model_df$Zika == min(model_df$Zika, na.rm = T), ]

BDcomp = data.frame(BD = eye_bds,
                    Zikrate = apply(expo[, 2:52], 2, sum, na.rm = T) / nrow(expo),
                    noZikrate = apply(noexpo[, 2:52], 2, sum, na.rm = T)  / nrow(noexpo))

BDcomp$Zikrate = BDcomp$Zikrate * 10000
BDcomp$noZikrate = BDcomp$noZikrate * 10000
BDcomp$diff = BDcomp$Zikrate / BDcomp$noZikrate
BDcomp[rev(order(BDcomp$diff)), ]

# Q143, Congenital malformation of the choroid at crude RR of 8.05 only notable one
