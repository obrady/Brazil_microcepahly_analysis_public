
rm(list = ls())

###########################
### 01 - cofactors plot ###
###########################


# load in dataframe
setwd("/Users/eideobra/Desktop/Brazil_data/BEST")

#tab = read.csv("FIGURES/Primary_cofactors_table_confirmed_microcephaly.csv")
#outfile = "FIGURES/Cofactors_RR_plot_con.pdf"
tab = read.csv("FIGURES/Primary_cofactors_table_suspected_microcephaly.csv")
outfile = "FIGURES/Cofactors_RR_plot_sus.pdf"


require(RColorBrewer)

crudecol = brewer.pal(11, "Spectral")[1]
adjcol = brewer.pal(11, "Spectral")[11]

pdf(outfile, width = 14, height = 7)
par(mar = c(5, 14, 8, 8) + 0.1)
plot(tab$ORUnAdj, seq(8.33, 1.33, -1), axes = FALSE, pch = 17, ylab = "", xlab = "", xlim = c(0, 3), 
     col = crudecol, ylim = c(1, 10))
axis(1, at = seq(0, 3, length.out = 4), seq(0, 3, length.out = 4), xpd = T, cex.axis = 1.25)
text(1.5, -2.5, "Relative risk", xpd = T, cex = 1.5)

# shading and text
text_str = c("Dengue", 
             "Chikungunya", 
             "Bovine Viral Diarrhea Virus", 
             "Water toxins",
             "Zika",
             "Zika (coinfection with dengue)",
             "Zika (coinfection with Chikungunya)",
             "Zika (enhanced by prior yellow fever vaccination)")
for(i in seq(1, 8, 2)){
  #polygon(c(0, 7, 7, 0),
  #        c((i + 1), (i + 1), i, i),
  #        col = " light grey", border = NA)
  polygon(c(-7, 14, 14, -7),
          c((i + 1), (i + 1), i, i),
          col = " light grey", border = NA, xpd = T)
}
#for(i in 1:9){text(-2, i + 0.5, text_str[i], xpd = T, cex = 0.5)}


lines(c(1,1), c(1, 10), lty = 2, lwd = 1)
for(i in nrow(tab):1){
  lines(c(rev(tab$ORUnAdjLow)[i], rev(tab$ORUnAdjHigh)[i]), c(i + 0.33, i + 0.33))
}
points(tab$ORUnAdj, 8.33:1.33, col = crudecol, pch = 17, cex = 1.5)


# now add adjusted
for(i in nrow(tab):1){
  lines(c(rev(tab$ORAdjLow)[i], rev(tab$ORAdjHigh)[i]), c(i+0.66, i+0.66))
}
points(tab$ORAdj, 8.66:1.66, xpd = T, pch = 19, cex = 2, col = adjcol)

# legend
legend(1.5, 10, c("Unadjusted", "Adjusted"),
       col = c(crudecol, adjcol), pch = c(17, 19), bty = "n", horiz = T, xpd = T, cex = 1.5)
#grid()
dev.off()














###########################
### 02 - other causes ###
###########################

case_type = "confirmed_microcephaly"#"confirmed_microcephaly" # "suspected_microcephaly"


tab = read.csv(paste("FIGURES/Diff_outcomes_tab_", case_type, ".csv", sep = ""))
tab2 = read.csv(paste("FIGURES/Diff_outcomes_tab_INDIV_", case_type, ".csv", sep = ""))
# remove Q690 and Q699 which have been combines
tab2 = tab2[!(tab2$ICD10 %in% c("Q690", "Q699")), ]

# sdjust p-values for bonferroni correction with 12 hypothesis tests - now done directly so not needed
#tab$p_val = tab$p_val * 12
#tab$p_val[tab$p_val > 1] = 1
#tab2$Pval = tab2$Pval * 12
#tab2$Pval[tab2$Pval > 1] = 1

tab =data.frame(name = c(as.character(tab$outcome), as.character(tab2$ICD10)),
                ORUnAdj = c(tab$RR_unadjusted, tab2$CrudeRR),
                ORUnAdjLow = c(tab$RR_unadjusted, tab2$CrudeRR),
                ORUnAdjHigh = c(tab$RR_unadjusted, tab2$CrudeRR),
                ORAdj = c(tab$pred, tab2$AdjRR_mid),
                ORAdjLow = c(tab$pred_lower, tab2$AdjRR_low),
                ORAdjHigh = c(tab$pred_upper, tab2$AdjRR_high))

require(RColorBrewer)

crudecol = brewer.pal(11, "Spectral")[1]
adjcol = brewer.pal(11, "Spectral")[11]

pdf(paste("FIGURES/Other_BDs_RR_plot_", case_type, ".pdf", sep = ""), width = 3, height = 7)
#plot(tab$ORUnAdj, seq(12.33, 1.33, -1), axes = FALSE, pch = 17, ylab = "", xlab = "", 
plot(NA, NA, axes = FALSE, pch = 17, ylab = "", xlab = "", 
     xlim = c(0.5, 2), 
     #col = crudecol, ylim = c(1, 12), xpd = T)
     col = "white", ylim = c(1, 12), xpd = T)
axis(1, at = c(0.5, 1, 1.5, 2), c(0.5, 1, 1.5, 2), xpd = T, cex.axis = 1.25)
#text(1, -0.65, "Relative risk", xpd = T, cex = 1.5)

# shading
for(i in seq(1, 12, 2)){
  polygon(c(0, 7, 7, 0),
          c((i + 1), (i + 1), i, i),
          col = " light grey", border = NA)
}

lines(c(1,1), c(1, 13), lty = 2, lwd = 1, xpd = T)
#for(i in nrow(tab):1){
#  lines(c(rev(tab$ORUnAdjLow)[i], rev(tab$ORUnAdjHigh)[i]), c(i + 0.33, i + 0.33))
#}
#points(tab$ORUnAdj, 12.33:1.33, col = crudecol, pch = 17, cex = 1.5, xpd = T)


# now add adjusted
points(tab$ORAdj, 12.5:1.5, xpd = T, pch = 19, cex = 2, col = adjcol)
for(i in nrow(tab):1){
  lines(c(rev(tab$ORAdjLow)[i], rev(tab$ORAdjHigh)[i]), c(i+0.5, i+0.5))
}


# legend
#legend(0, 14.5, c("Crude", "Adjusted"),
#       col = c(crudecol, adjcol), pch = c(17, 19), bty = "n", horiz = T, xpd = T, cex = 1.75)
#grid()
dev.off()

