resp.sinasc.matcher <- function(resp_useful, sinasc_useful){
  # 01 data pre processing"
  sinasc_match = sinasc_useful[, c("DTNASC", "SEXO", "RACACORMAE", "CODMUNRES", "IDADEMAE", "IBGE_STATE_CODE", "CODANOMAL")]
  resp_match = resp_useful[, c("DT_NASCIMENTO_RN", "TP_SEXO", "CO_RACA_COR", "CO_MUNICIPIO_IBGE", "IDADEMAE", "SG_UF_RESIDENCIA")]
  
  colnames(sinasc_match) = c("DOB", "sex", "race", "mun", "age", "sta", "bds")
  colnames(resp_match) = c("DOB", "sex", "race", "mun", "age", "sta")
  
  # character to factor
  sinasc_match$sex = as.factor(sinasc_match$sex)
  resp_match$sex[resp_match$sex == "I"] = NA
  resp_match$sex = as.factor(resp_match$sex)
  sinasc_match$race = as.factor(sinasc_match$race)
  
  # plausibility constraints
  sinasc_match$age[((sinasc_match$age < 10) | (sinasc_match$age > 55))] = NA
  resp_match$age[((resp_match$age < 10) | (resp_match$age > 55))] = NA
  
  # add bds to resp
  resp_match$bds = TRUE
  
  # date transformation - convert to day month year
  #sinasc_match$DOB_day = mday(sinasc_match$DOB)
  sinasc_match$DOB_month = month(sinasc_match$DOB)
  sinasc_match$DOB_year = year(sinasc_match$DOB)
  
  #resp_match$DOB_day = mday(resp_match$DOB)
  resp_match$DOB_month = month(resp_match$DOB)
  resp_match$DOB_year = year(resp_match$DOB)
  
  sinasc_match$DOB = as.numeric(sinasc_match$DOB - as.Date("2015-01-01"))
  resp_match$DOB = as.numeric(as.Date(resp_match$DOB) - as.Date("2015-01-01"))
  
  #sinasc_match = sinasc_match[, colnames(sinasc_match) != "DOB"]
  #resp_match = resp_match[, colnames(resp_match) != "DOB"]
  
  # split sinasc into those with and without MC identified at birth
  sinasc_match$bds = grepl("Q02", sinasc_match$bds)
  #sinasc_match_MC <- sinasc_match[grepl("Q02", sinasc_match$bds), ]
  #sinasc_match_noMC <- sinasc_match[!grepl("Q02", sinasc_match$bds), ]
  #sinasc_match_MC = sinasc_match_MC[, colnames(sinasc_match_MC) != "bds"]
  #sinasc_match_noMC = sinasc_match_noMC[, colnames(sinasc_match_noMC) != "bds"]
  
  # add an UID to sinasc_match
  sinasc_match$UID = 1:nrow(sinasc_match)
  
  # trimming nonsense municipality records from sinasc - will get cut later in the analysis anyway
  munLL <- read.csv("REFERENCE/MUN_lat_long.csv")
  sinasc_match = sinasc_match[(sinasc_match$mun %in% munLL$IBGE), ]
  resp_match = resp_match[(resp_match$mun %in% munLL$IBGE), ]
  
  
  
  # linkage
  
  # which records to delete from sinasc list
  sin_del <- data.frame(rowmatch = rep(NA, nrow(resp_match)), 
                        mundist = rep(NA, nrow(resp_match)),
                        timedist = rep(NA, nrow(resp_match)))
  
  for(i in 1:nrow(resp_match)){
    
    # blocking stages
    # state
    sinasc_match_focus = sinasc_match[((sinasc_match$sta == resp_match$sta[i]) | is.na(sinasc_match$sta)), ]
    # sex of baby
    if(!is.na(resp_match$sex[i])){
      sinasc_match_focus = sinasc_match_focus[((sinasc_match_focus$sex == resp_match$sex[i]) | is.na(sinasc_match_focus$sex)), ]
    }
    # age of mother (+/- 1 year)
    if(!is.na(resp_match$age[i])){
      sinasc_match_focus = sinasc_match_focus[((sinasc_match_focus$age %in% c(resp_match$age[i] - 1,
                                                                              resp_match$age[i],
                                                                              resp_match$age[i] + 1)) | is.na(sinasc_match_focus$age)), ]
    }
    
    # create a distance matrix between resp record and candidate sinasc records
    # get lat longs
    resp_LLs = munLL[match(resp_match[i, "mun"], munLL$IBGE), c("X1", "X2")]
    sinasc_LLs = munLL[match(unique(sinasc_match_focus$mun), munLL$IBGE), c("X1", "X2")]
    
    #pull up a list of distances from the resp record
    dist_from_resp <- distm(rbind(resp_LLs, sinasc_LLs))[, 1]
    # now re-assign distance back to sinasc_match
    asign_table <- data.frame(mun = unique(sinasc_match_focus$mun),
                              distance = dist_from_resp[-1])
    sinasc_match_focus$dist = asign_table$distance[match(sinasc_match_focus$mun, asign_table$mun)]
    
    # calculate space time difference between records and time record
    mundist <- sinasc_match_focus$dist
    timedist <- sqrt((sinasc_match_focus$DOB - resp_match$DOB[i])^2)
    
    # scale
    mundist_s = (mundist - min(mundist)) / max((mundist - min(mundist)))
    timedist_s = (timedist - min(timedist)) / max((timedist - min(timedist)))
    
    # final result
    sinasc_match_focus_final = sinasc_match_focus[which.min(mundist_s + timedist_s), ]
    
    
    # record result
    sin_del[i, 1] = sinasc_match_focus_final$UID
    
    # record how close the space-time metrics were
    sin_del$mundist[i] = mundist[which.min(mundist_s + timedist_s)]
    sin_del$timedist[i] = timedist[which.min(mundist_s + timedist_s)]
  }
  
  # return matches
  return(sin_del)
}