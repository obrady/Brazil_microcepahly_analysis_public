# creates time series matrices for dengue, chik and Zika


rm(list = ls())

#require(data.table)
setwd("/Users/eideobra/Desktop/Brazil_data/BEST")


# arbovirs testing data
arbo <- read.csv("GAL/Combined_arbo_JUL18.csv") # 1,542,863 records
# municipality data
all_mun <- read.csv("REFERENCE/GAUL_IBGE_conversion_full.csv")
# extra state lookup table
d_state_lu <- read.csv("REFERENCE/DENV_STATE_IDs.csv", stringsAsFactors = FALSE)


### Part A: pre-processing

# remove any records that don't have a recognised municipality code (1421 records)
arbo <- arbo[arbo$mun %in% all_mun$IBGE_CODE, ]
# date pre-processing
arbo$date = as.Date(arbo$date)
# remove dates prior to 2014 (461 records)
arbo = arbo[arbo$date > as.Date("2015-01-01"), ]

# add a state code
arbo$state = all_mun$IBGE_STATE_CODE[match(arbo$mun, all_mun$IBGE_CODE)]



# list of all municipalities and states
all_muns = unique(all_mun$IBGE_CODE)
all_stas = as.character(unique(all_mun$IBGE_STATE_CODE))

# lit of municipalities with reports of suspected arbovirus cases
not_muns <- unique(arbo$mun) # ~ 800 don't report any arboviruses
not_stas <- as.character(unique(arbo$state)) # all states report
not_stas = not_stas[!is.na(not_stas)]

# function for creating the time series matrices
#arbodat <- arbo[arbo$zik_conf, ]
TSmat.fill <- function(arbodat, disease = "zik", admin = "mun", conf_type = "conf"){
  # 1: initialise matrix
  if(admin == "mun"){
    tsmat <- matrix(0, nrow = length(all_muns), ncol = 158)
    tsmat[, 1] = all_muns
    tsmat[, 2] = NA
  }else{
    tsmat <- as.data.frame(matrix(0, nrow = length(all_stas), ncol = 158))
    tsmat[, 1] = all_stas
    tsmat[, 2] = NA
  }
  
  # 2: Add case data
  for(i in 1:nrow(tsmat)){
    # looking for states or municipalities?
    colind_loc <- which.max(colnames(arbodat) %in% admin)
    if(tsmat[i, 1] %in% arbodat[, as.numeric(colind_loc)]){
      # subset to area of interest
      selection <- arbodat[arbodat[, as.numeric(colind_loc)] %in% tsmat[i, 1], ]
      
      # subset to disease of interest and type of confirmation method
      colind <- which.max(colnames(selection) %in% paste(disease, conf_type, sep = "_"))
      selection = selection[selection[, colind], ]
      
      if(nrow(selection) > 0){
        # convert sample date to epi week
        sel_dates <- floor(as.numeric(difftime(as.Date(selection$date), "2015-01-01", units = "weeks")))
        tsmat[i, 2] = min(sel_dates) + 1
        # table up number of cases per week then add them to the master time series matrix
        sel_dates_table = data.frame(as.numeric(names(table(sel_dates))), as.numeric(table(sel_dates)))
        for(j in 1:nrow(sel_dates_table)){
          tsmat[i, (sel_dates_table[j, 1] + 3)] = sel_dates_table[j, 2]
        }
        # checking
        #if(nrow(selection) != sum(tsmat[i, 3:158])){break}
      }
    }
  }
  
  # 3: make NA before first case reported at state level
  if(admin == "mun"){
    # first calculate state first arrival date
    mun_arr <- data.frame(mun = tsmat[, 1],
                          first_arr = tsmat[, 2],
                          state = all_mun$IBGE_STATE_CODE[match(tsmat[, 1], all_mun$IBGE_CODE)])
    sta_arr = aggregate(mun_arr$first_arr, by = list(mun_arr$state), FUN = min, na.rm = T)
    mun_arr$new_first_arr = sta_arr[, 2][match(mun_arr$state, sta_arr[, 1])]
    # now run through
    for(i in 1:nrow(tsmat)){
      tsmat[i, 3:(3 + mun_arr$new_first_arr[i] - 1)] = NA
    }
  }else{
    for(i in 1:nrow(tsmat)){
      tsmat[i, 3:(3 + tsmat[i, 2] - 1)] = NA
    }
  }
  
  # 4: add relevant headers
  if(admin == "mun"){
    colnames(tsmat) = c("IBGE_Municipality", "Week_first_intro", paste("Wk_", 1:156, sep = ""))
  }else{
    colnames(tsmat) = c("State", "Week_first_intro", paste("Wk_", 1:156, sep = ""))
  }
  
  # 5: add an FID column
  tsmat <- data.frame(FID = 1:nrow(tsmat), tsmat)
  
  # 6: return matrix
  return(tsmat)
}


### create time series matrices

# Zika
zik_con_mun <- TSmat.fill(arbodat = arbo[arbo$zik_conf, ],
                          disease = "zik",
                          admin = "mun",
                          conf_type = "conf")
zik_con_sta <- TSmat.fill(arbodat = arbo[arbo$zik_conf, ],
                          disease = "zik",
                          admin = "state",
                          conf_type = "conf")


zik_sus_mun <- TSmat.fill(arbodat = arbo[arbo$zik_sus, ],
                          disease = "zik",
                          admin = "mun",
                          conf_type = "sus")
zik_sus_sta <- TSmat.fill(arbodat = arbo[arbo$zik_sus, ],
                          disease = "zik",
                          admin = "state",
                          conf_type = "sus")

zik_tested_mun <- TSmat.fill(arbodat = arbo[arbo$zik_tested, ],
                             disease = "zik",
                             admin = "mun",
                             conf_type = "tested")
zik_tested_sta <- TSmat.fill(arbodat = arbo[arbo$zik_tested, ],
                             disease = "zik",
                             admin = "state",
                             conf_type = "tested")

# Chikungunya
chik_con_mun <- TSmat.fill(arbodat = arbo[arbo$chik_conf, ],
                           disease = "chik",
                           admin = "mun",
                           conf_type = "conf")
chik_con_sta <- TSmat.fill(arbodat = arbo[arbo$chik_conf, ],
                           disease = "chik",
                           admin = "state",
                           conf_type = "conf")

chik_sus_mun <- TSmat.fill(arbodat = arbo[arbo$chik_sus, ],
                           disease = "chik",
                           admin = "mun",
                           conf_type = "sus")
chik_sus_sta <- TSmat.fill(arbodat = arbo[arbo$chik_sus, ],
                           disease = "chik",
                           admin = "state",
                           conf_type = "sus")

chik_tested_mun <- TSmat.fill(arbodat = arbo[arbo$chik_tested, ],
                              disease = "chik",
                              admin = "mun",
                              conf_type = "tested")
chik_tested_sta <- TSmat.fill(arbodat = arbo[arbo$chik_tested, ],
                              disease = "chik",
                              admin = "state",
                              conf_type = "tested")

# Dengue
den_con_mun <- TSmat.fill(arbodat = arbo[arbo$den_conf, ],
                          disease = "den",
                          admin = "mun",
                          conf_type = "conf")
den_con_sta <- TSmat.fill(arbodat = arbo[arbo$den_conf, ],
                          disease = "den",
                          admin = "state",
                          conf_type = "conf")

den_sus_mun <- TSmat.fill(arbodat = arbo[arbo$den_sus, ],
                          disease = "den",
                          admin = "mun",
                          conf_type = "sus")
den_sus_sta <- TSmat.fill(arbodat = arbo[arbo$den_sus, ],
                          disease = "den",
                          admin = "state",
                          conf_type = "sus")

den_tested_mun <- TSmat.fill(arbodat = arbo[arbo$den_tested, ],
                             disease = "den",
                             admin = "mun",
                             conf_type = "tested")
den_tested_sta <- TSmat.fill(arbodat = arbo[arbo$den_tested, ],
                             disease = "den",
                             admin = "state",
                             conf_type = "tested")

# quick totals checks
sum(zik_con_mun[, 4:159], na.rm = T)
sum(zik_con_sta[, 4:159], na.rm = T)
sum(zik_sus_mun[, 4:159], na.rm = T)
sum(zik_sus_sta[, 4:159], na.rm = T)
sum(zik_tested_mun[, 4:159], na.rm = T)
sum(zik_tested_sta[, 4:159], na.rm = T)
# 17k confirmed, 111k tested, 195k suspected

sum(chik_con_mun[, 4:159], na.rm = T)
sum(chik_con_sta[, 4:159], na.rm = T)
sum(chik_sus_mun[, 4:159], na.rm = T)
sum(chik_sus_sta[, 4:159], na.rm = T)
sum(chik_tested_mun[, 4:159], na.rm = T)
sum(chik_tested_sta[, 4:159], na.rm = T)
# 110k confirmed, 221k tested, 436k suspected

sum(den_con_mun[, 4:159], na.rm = T)
sum(den_con_sta[, 4:159], na.rm = T)
sum(den_sus_mun[, 4:159], na.rm = T)
sum(den_sus_sta[, 4:159], na.rm = T)
sum(den_tested_mun[, 4:159], na.rm = T)
sum(den_tested_sta[, 4:159], na.rm = T)
# 161k confirmed, 566k tested, 908k suspected

# quick compare to previous time series
# Zika
sum(read.csv("MUN_COVS/Zika_con_mun_ts_AUG17_V2.csv")[, 4:159], na.rm = T)
sum(read.csv("MUN_COVS/Zika_not_mun_ts_AUG17_V2.csv")[, 4:159], na.rm = T)
# 800k confirmed, 1.4m suspected

# chik
sum(read.csv("MUN_COVS/Chikungunya_not_mun_ts_AUG17_V2.csv")[, 4:159], na.rm = T)
# 810k suspected

# dengue
sum(read.csv("MUN_COVS/Dengue_not_mun_ts_AUG17_V2.csv")[, 4:159], na.rm = T)
# 806k suspected



# checking co-linearity of suspected and confirmed cases
# Zika
plot(log(rowSums(zik_con_mun[, 4:159], na.rm = T)), log(rowSums(zik_sus_mun[, 4:159], na.rm = T)))
codf <- data.frame(conf = log(rowSums(zik_con_mun[, 4:159], na.rm = T)),
                   sus = log(rowSums(zik_sus_mun[, 4:159], na.rm = T)))
codf = codf[is.finite(rowSums(codf)), ]
comod = lm(conf ~ sus, data = codf)
comod = summary(comod)
comod$adj.r.squared # 0.59

# Chik
plot(log(rowSums(chik_con_mun[, 4:159], na.rm = T)), log(rowSums(chik_sus_mun[, 4:159], na.rm = T)))
codf <- data.frame(conf = log(rowSums(chik_con_mun[, 4:159], na.rm = T)),
                   sus = log(rowSums(chik_sus_mun[, 4:159], na.rm = T)))
codf = codf[is.finite(rowSums(codf)), ]
comod = lm(conf ~ sus, data = codf)
comod = summary(comod)
comod$adj.r.squared # 0.72

# Dengue
plot(log(rowSums(den_con_mun[, 4:159], na.rm = T)), log(rowSums(den_sus_mun[, 4:159], na.rm = T)))
codf <- data.frame(conf = log(rowSums(den_con_mun[, 4:159], na.rm = T)),
                   sus = log(rowSums(den_sus_mun[, 4:159], na.rm = T)))
codf = codf[is.finite(rowSums(codf)), ]
comod = lm(conf ~ sus, data = codf)
comod = summary(comod)
comod$adj.r.squared # 0.68

# conclusion: suspected and confirmed cases could give different 
# results for Zika and possibly dengue




# save new files
write.csv(zik_con_mun, file = "MUN_COVS/JUL18_timeseries/zik_con_mun.csv")
write.csv(zik_sus_mun, file = "MUN_COVS/JUL18_timeseries/zik_sus_mun.csv")
write.csv(zik_con_sta, file = "MUN_COVS/JUL18_timeseries/zik_con_sta.csv")
write.csv(zik_sus_sta, file = "MUN_COVS/JUL18_timeseries/zik_sus_sta.csv")
write.csv(zik_tested_mun, file = "MUN_COVS/JUL18_timeseries/zik_tested_mun.csv")
write.csv(zik_tested_sta, file = "MUN_COVS/JUL18_timeseries/zik_tested_sta.csv")

write.csv(chik_con_mun, file = "MUN_COVS/JUL18_timeseries/chik_con_mun.csv")
write.csv(chik_sus_mun, file = "MUN_COVS/JUL18_timeseries/chik_sus_mun.csv")
write.csv(chik_con_sta, file = "MUN_COVS/JUL18_timeseries/chik_con_sta.csv")
write.csv(chik_sus_sta, file = "MUN_COVS/JUL18_timeseries/chik_sus_sta.csv")
write.csv(chik_tested_mun, file = "MUN_COVS/JUL18_timeseries/chik_tested_mun.csv")
write.csv(chik_tested_sta, file = "MUN_COVS/JUL18_timeseries/chik_tested_sta.csv")

write.csv(den_con_mun, file = "MUN_COVS/JUL18_timeseries/den_con_mun.csv")
write.csv(den_sus_mun, file = "MUN_COVS/JUL18_timeseries/den_sus_mun.csv")
write.csv(den_con_sta, file = "MUN_COVS/JUL18_timeseries/den_con_sta.csv")
write.csv(den_sus_sta, file = "MUN_COVS/JUL18_timeseries/den_sus_sta.csv")
write.csv(den_tested_mun, file = "MUN_COVS/JUL18_timeseries/den_tested_mun.csv")
write.csv(den_tested_sta, file = "MUN_COVS/JUL18_timeseries/den_tested_sta.csv")




# figures for SI:

##################################
# 01 plots of suspected and confirmed cases over time
##################################
# suspected
pdf(file = "/Users/eideobra/Dropbox/06_Brazil_MoH/Writeup/Comments_versions/NewFIGS/Suspected_cases_time.pdf",
    height = 10, width = 5)
par(mfrow = c(3, 1), mar = c(5, 4, 4, 5) + 0.1)

# Zika
Z_S_TS <- colSums(zik_sus_sta[, 4:159], na.rm = T)
Z_C_TS <- 0.5 * max(Z_S_TS, na.rm = T) * colSums(zik_con_sta[, 4:159], na.rm = T) / max(colSums(zik_con_sta[, 4:159], na.rm = T), na.rm = T)

plot(as.Date("2015-01-01") + 7 * (1:156),
     Z_S_TS, type = "l",
     xlab = "Date",
     ylab = "Total suspected cases")
lines(as.Date("2015-01-01") + 7 * (1:156),
      Z_C_TS,
      col = "grey")
axis(4, seq(0, max(Z_S_TS), length.out = 5) * 0.5, round(seq(0, max(colSums(zik_con_sta[, 4:159], na.rm = T)), length.out = 5), 0))
text(as.Date("2018-07-12"), 4000, srt = 90, label = "PCR confirmed cases", xpd = T)
text(as.Date("2015-04-12"), max(colSums(zik_sus_sta[, 4:159]), na.rm = T) * 0.9, "Zika", cex = 2)


# chik
C_S_TS <- colSums(chik_sus_sta[, 4:159], na.rm = T)
C_C_TS <- 0.5 * max(C_S_TS, na.rm = T) * colSums(chik_con_sta[, 4:159], na.rm = T) / max(colSums(chik_con_sta[, 4:159], na.rm = T), na.rm = T)

plot(as.Date("2015-01-01") + 7 * (1:156),
     C_S_TS, type = "l",
     xlab = "Date",
     ylab = "Total suspected cases",
     col = "red")
lines(as.Date("2015-01-01") + 7 * (1:156),
      C_C_TS,
      col = "grey")
axis(4, seq(0, max(C_S_TS), length.out = 5) * 0.5, round(seq(0, max(colSums(chik_con_sta[, 4:159], na.rm = T)), length.out = 5), 0))
text(as.Date("2018-07-12"), 5000, srt = 90, label = "PCR confirmed cases", xpd = T)
text(as.Date("2015-06-12"), max(colSums(chik_sus_sta[, 4:159]), na.rm = T) * 0.9, "Chikungunya", cex = 2)


# dengue
D_S_TS <- colSums(den_sus_sta[, 4:159], na.rm = T)
D_C_TS <- 0.5 * max(D_S_TS, na.rm = T) * colSums(den_con_sta[, 4:159], na.rm = T) / max(colSums(den_con_sta[, 4:159], na.rm = T), na.rm = T)

plot(as.Date("2015-01-01") + 7 * (1:156),
     D_S_TS, type = "l",
     xlab = "Date",
     ylab = "Total suspected cases",
     col = "blue")
lines(as.Date("2015-01-01") + 7 * (1:156),
      D_C_TS,
      col = "grey")
axis(4, seq(0, max(D_S_TS), length.out = 5) * 0.5, round(seq(0, max(colSums(den_con_sta[, 4:159], na.rm = T)), length.out = 5), 0))
text(as.Date("2018-07-12"), 20000, srt = 90, label = "PCR confirmed cases", xpd = T)
text(as.Date("2015-04-12"), max(colSums(den_sus_sta[, 4:159]), na.rm = T) * 0.9, "Dengue", cex = 2)
dev.off()




##################################
# 02 plots of total reported ans suspected by state
##################################
require(RColorBrewer)
colPAL <- brewer.pal(5, "Dark2")

pdf(file = "/Users/eideobra/Dropbox/06_Brazil_MoH/Writeup/Comments_versions/NewFIGS/Suspected_cases_geography.pdf",
    height = 10, width = 10)
par(mfrow = c(3, 1))
# Zika
newdf = data.frame(state = zik_sus_sta[, 2],
                   tcase = rowSums(zik_sus_sta[, 4:159], na.rm = T))
newdf$Region =all_mun$Region[match(newdf$state, all_mun$IBGE_STATE_CODE)]
newdf = newdf[order(newdf$Region), ]

plot(newdf[, 2], col = colPAL[as.numeric(newdf$Region)], pch = 19, axes = F,
     xlab = "State", ylab = "Suspected cases (thousands)", ylim = c(0, 50000))
axis(1, 1:27, newdf[, 1])
axis(2, c(0, 10000, 20000, 30000, 40000, 50000), c(0, 10, 20, 30, 40, 50))
text(2, 50000, "Zika", cex = 2, xpd = T)
legend("top", fill = colPAL[as.numeric(unique(newdf$Region))], legend = unique(newdf$Region),
       horiz = T, cex = 1.5)

# Chik
newdf = data.frame(state = chik_sus_sta[, 2],
                   tcase = rowSums(chik_sus_sta[, 4:159], na.rm = T))
newdf$Region =all_mun$Region[match(newdf$state, all_mun$IBGE_STATE_CODE)]
newdf = newdf[order(newdf$Region), ]

plot(newdf[, 2], col = colPAL[as.numeric(newdf$Region)], pch = 19, axes = F,
     xlab = "State", ylab = "Suspected cases (thousands)", ylim = c(0, 150000))
text(2, 150000, "Chikungunya", cex = 2, xpd = T)
axis(1, 1:27, newdf[, 1])
axis(2, c(0, 50000, 100000, 150000), c(0, 50, 100, 150))


# dengue
newdf = data.frame(state = den_sus_sta[, 2],
                   tcase = rowSums(den_sus_sta[, 4:159], na.rm = T))
newdf$Region =all_mun$Region[match(newdf$state, all_mun$IBGE_STATE_CODE)]
newdf = newdf[order(newdf$Region), ]
plot(newdf[, 2], col = colPAL[as.numeric(newdf$Region)], pch = 19, axes = F,
     xlab = "State", ylab = "Reported cases (thousands)", ylim = c(0, 150000))
text(2, 150000, "Dengue", cex = 2, xpd = T)
axis(1, 1:27, newdf[, 1])
axis(2, c(0, 50000, 100000, 150000), c(0, 50, 100, 150))
dev.off()

##################################
# 03 municipality level stats
##################################
# - histogram of number of reported cases
# histogram of number of weeks of reported cases
# analysis of reporting decline

# part1: histogram of number of reported cases
pdf(file = "/Users/eideobra/Dropbox/06_Brazil_MoH/Writeup/Comments_versions/NewFIGS/Municipality_total_cases.pdf",
    height = 15, width = 10)
par(mfrow = c(3, 2))

# Zika
quant <- rowSums(zik_sus_mun[, 4:159], na.rm = T)
quant <- log(quant[quant > 0], 10)
hist(quant, axes = F,
     ylab = "Count of municipality", xlab = "Total cases",
     main = "Suspected zika")
abline(v = median(quant), col = "red")
text(median(quant), 1100, round(10^median(quant), 2), xpd = T)
axis(1, 0:5, 10^(0:5))
axis(2, seq(0, 1000, 200), seq(0, 1000, 200))

quant <- rowSums(zik_con_mun[, 4:159], na.rm = T)
quant <- log(quant[quant > 0], 10)
hist(quant, axes = F,
     ylab = "Count of municipality", xlab = "Total cases",
     main = "Confirmed zika")
abline(v = median(quant), col = "red")
text(median(quant), 700, round(10^median(quant), 2), xpd = T)
axis(1, 0:4, 10^(0:4))
axis(2, seq(0, 500, 100), seq(0, 500, 100))

# chikungunya
quant <- rowSums(chik_sus_mun[, 4:159], na.rm = T)
quant <- log(quant[quant > 0], 10)
hist(quant, axes = F,
     ylab = "Count of municipality", xlab = "Total cases",
     main = "Suspected chikungunya")
abline(v = median(quant), col = "red")
text(median(quant), 1000, round(10^median(quant), 2), xpd = T)
axis(1, 0:5, 10^(0:5))
axis(2, seq(0, 1000, 200), seq(0, 1000, 200))

quant <- rowSums(chik_con_mun[, 4:159], na.rm = T)
quant <- log(quant[quant > 0], 10)
hist(quant, axes = F,
     ylab = "Count of municipality", xlab = "Total cases",
     main = "Confirmed chikungunya")
abline(v = median(quant), col = "red")
text(median(quant), 1000, round(10^median(quant), 2), xpd = T)
axis(1, 0:4, 10^(0:4))
axis(2, seq(0, 1000, 200), seq(0, 1000, 200))


# dengue
quant <- rowSums(den_sus_mun[, 4:159], na.rm = T)
quant <- log(quant[quant > 0], 10)
hist(quant, axes = F,
     ylab = "Count of municipality", xlab = "Total cases",
     main = "Suspected dengue")
abline(v = median(quant), col = "red")
text(median(quant), 1100, round(10^median(quant), 2), xpd = T)
axis(1, 0:5, 10^(0:5))
axis(2, seq(0, 1000, 200), seq(0, 1000, 200))


quant <- rowSums(den_con_mun[, 4:159], na.rm = T)
quant <- log(quant[quant > 0], 10)
hist(quant, axes = F,
     ylab = "Count of municipality", xlab = "Total cases",
     main = "Confirmed dengue")
abline(v = median(quant), col = "red")
text(median(quant), 600, round(10^median(quant), 2), xpd = T)
axis(1, 0:4, 10^(0:4))
axis(2, seq(0, 500, 100), seq(0, 500, 100))

dev.off()

# part 2: histogram of outbreak duration
pdf(file = "/Users/eideobra/Dropbox/06_Brazil_MoH/Writeup/Comments_versions/NewFIGS/Municipality_outbreak_duration.pdf",
    height = 15, width = 10)
par(mfrow = c(3, 2))
# zika
quant <- apply(zik_sus_mun[, 4:159], 1, function(x) sum(x > 0, na.rm = T))
quant <- quant[quant > 0]
hist(quant, n = 100, xlab = "Total weeks of cases", 
     ylab = "Count of municipality", main = "Suspected zika")
abline(v = median(quant), col = "red")
text(median(quant) + 10, 800, round(median(quant), 2))

quant <- apply(zik_con_mun[, 4:159], 1, function(x) sum(x > 0, na.rm = T))
quant <- quant[quant > 0]
hist(quant, n = 100, xlab = "Total weeks of cases", 
     ylab = "Count of municipality", main = "Confirmed zika")
abline(v = median(quant), col = "red")
text(median(quant) + 10, 500, round(median(quant), 2))

# chik
quant <- apply(chik_sus_mun[, 4:159], 1, function(x) sum(x > 0, na.rm = T))
quant <- quant[quant > 0]
hist(quant, n = 100, xlab = "Total weeks of cases", 
     ylab = "Count of municipality", main = "Suspected chikungunya")
abline(v = median(quant), col = "red")
text(median(quant) + 10, 800, round(median(quant), 2))

quant <- apply(chik_con_mun[, 4:159], 1, function(x) sum(x > 0, na.rm = T))
quant <- quant[quant > 0]
hist(quant, n = 100, xlab = "Total weeks of cases", 
     ylab = "Count of municipality", main = "Confirmed chikungunya")
abline(v = median(quant), col = "red")
text(median(quant) + 10, 800, round(median(quant), 2))

# dengue
quant <- apply(den_sus_mun[, 4:159], 1, function(x) sum(x > 0, na.rm = T))
quant <- quant[quant > 0]
hist(quant, n = 100, xlab = "Total weeks of cases", 
     ylab = "Count of municipality", main = "Suspected dengue")
abline(v = median(quant), col = "red")
text(median(quant) + 10, 800, round(median(quant), 2))

quant <- apply(den_con_mun[, 4:159], 1, function(x) sum(x > 0, na.rm = T))
quant <- quant[quant > 0]
hist(quant, n = 100, xlab = "Total weeks of cases", 
     ylab = "Count of municipality", main = "Confirmed dengue")
abline(v = median(quant), col = "red")
text(median(quant) + 10, 800, round(median(quant), 2))
dev.off()

# 04 tables of totals suspected, tested and confirmed cases
tabdf <- data.frame(zik_suspected = sum(arbo$zik_sus),
                    zik_tested = sum(arbo$zik_tested),
                    zik_confirmed = sum(arbo$zik_conf),
                    chik_suspected = sum(arbo$chik_sus),
                    chik_tested = sum(arbo$chik_tested),
                    chik_confirmed = sum(arbo$chik_conf),
                    dengue_suspected = sum(arbo$den_sus),
                    dengue_tested = sum(arbo$den_tested),
                    dengue_confirmed = sum(arbo$den_conf))
# now manually enter
tabdf

# 05 regional analysis of testing and confirmation rate
require(RColorBrewer)
lookup = read.csv("REFERENCE/GAUL_IBGE_conversion_full.csv")
plot.func <- function(fdat, fcol, lab, yaxislab, plot_ledge = FALSE){
  # add region
  fdat$Region = lookup$Region[match(fdat$state, lookup$IBGE_STATE_CODE)]
  fdat = fdat[order(fdat$Region), ]
  
  # custom colour pallete
  colPAL <- brewer.pal(5, "Dark2")
  
  plot(fdat[, fcol], ylim = c(0, 1), xlab = "State", ylab = yaxislab,
       col = colPAL[as.numeric(fdat$Region)],
       axes = F, pch = 19, cex = 3)
  axis(1, 1:27, fdat$state)
  axis(2, seq(0, 1, 0.1), seq(0, 1, 0.1))
  abline(h = mean(fdat[, fcol]), col = "grey")
  upper = mean(fdat[, fcol]) + 1.95 * sd(fdat[, fcol])
  if(upper > 1){upper = 1}
  lower = mean(fdat[, fcol]) - 1.95 * sd(fdat[, fcol])
  if(lower < 0){lower = 0}
  
  abline(h = upper, col = "grey", lty = 2)
  abline(h = lower, col = "grey", lty = 2)
  if(plot_ledge){
    legend("top", legend = as.character(unique(fdat$Region)), fill = colPAL,
           horiz = T, cex = 1.5)
  }
  text(2, 1, lab, cex =1.75)
}

pdf("/Users/eideobra/Dropbox/06_Brazil_MoH/Writeup/Comments_versions/NewFIGS/Percentage_PCRtested.pdf", height = 15, width = 10)
par(mfrow = c(3, 1))
denobj <- data.frame(state = den_tested_sta$State,
                     conf_rate = rowSums(den_tested_sta[, 4:159], na.rm = T) / rowSums(den_sus_sta[, 4:159], na.rm = T))
chikobj <- data.frame(state = chik_tested_sta$State,
                     conf_rate = rowSums(chik_tested_sta[, 4:159], na.rm = T) / rowSums(chik_sus_sta[, 4:159], na.rm = T))
zikobj <- data.frame(state = ik_tested_sta$State,
                     conf_rate = rowSums(zik_tested_sta[, 4:159], na.rm = T) / rowSums(zik_sus_sta[, 4:159], na.rm = T))

plot.func(fdat = zikobj, fcol = 2, lab = "Zika", yaxislab = "Proportion PCR tested", plot_ledge = T)
plot.func(fdat = chikobj, fcol = 2, lab = "Chikungunya", yaxislab = "Proportion PCR tested")
plot.func(fdat = denobj, fcol = 2, lab = "Dengue", yaxislab = "Proportion PCR tested")
dev.off()





# 06 dropoff in confirmation rate over time- testing fatigue
# specifically for Zika i believe
zik_conf_rate = zik_tested_mun[, 4:159] / zik_sus_mun[, 4:159] 
coef_store = rep(NA, nrow(zik_conf_rate))
pval_store = rep(NA, nrow(zik_conf_rate))

for(i in 1:nrow(zik_conf_rate)){
  if(sum(is.finite(as.numeric(zik_conf_rate[i, ]))) > 1){
    fit_df <- data.frame(week = 1:length(zik_conf_rate[i, ]),
                         conf_rate = as.numeric(zik_conf_rate[i, ]))
    fit_df = fit_df[is.finite(fit_df[, 2]), ]
    # trim to first date of detection
    fit_df = fit_df[which.max(fit_df$conf_rate > 0):nrow(fit_df), ]
    if(nrow(fit_df) > 1){
      # now fit lm
      lm_mod <- lm(conf_rate ~ week, data = fit_df)
      coef_store[i] = lm_mod$coefficients[2]
      pval_store[i] = summary(lm_mod)$coefficients[2, 4]
    }
  }
}

# how many have statistically significant declines in 
# proportion reported after first identification?
sig0 <- ((coef_store < 0) + (pval_store < 0.05)) == 2
sum(sig0, na.rm = T)
# 152 municiaplities
# how much of a reduction are we talking?
hist(coef_store[sig0])
median(coef_store[sig0], na.rm = T) # -0.003158172
# over the median duration of 3 weeks:
median(coef_store[sig0], na.rm = T) * 2 # = a less than 1% reduction in PCR positivity rate over the three weeks.
