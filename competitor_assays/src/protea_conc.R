### COMPETITOR ASSAYS ###
# description:  determine optimal proteasome concentration for competitor assays
# input:        fluorescense measurements (2 hrs, vevery 30s)
#               LLE, VGR, LLVY + protea
# output:       time course plots
# author:       HR

library(dplyr)
library(readxl)

### INPUT ###
data = readxl::read_excel("data/20220314.xlsx", sheet = 1, range = c("B35:CU276"))

### MAIN PART ###
# ----- annotation -----
# found in pipetting_scheme.xlsx
# 20220314_1

annotation = data.frame(cell = paste0(rep(toupper(letters[1:8]), each = 12), rep(seq(1,12), 8)),
                        protea_conc = NA)

protea_conc = c(rep(10, 12), rep(20, 12), rep(50, 12), rep(0, 12), rep(50, 2), rep(0,2))
substrate = c(rep(rep(rep(c("LLE","VGR","LLVY"), each = 2),2),4), rep("none",4))
substrate_conc = c(rep(c(rep(50, 6), rep(100,6)), 4), rep(0,4))

annotation = annotation[c(1:52),]
annotation$protea_conc = protea_conc
annotation$substrate = substrate
annotation$substrate_conc = substrate_conc

# ----- subtract 1st measurement -----
# to account for pipetting delay
FLUORO = as.matrix(data[, c(3:c(2+nrow(annotation)))])
colnames(FLUORO) = annotation$cell

min1 = FLUORO[1,]
FLUORO = FLUORO - matrix(rep(min1, each=nrow(FLUORO)), nrow = nrow(FLUORO))


# ----- plot -----
data$Time[1:5]
time = seq(6, 6+3600*2, by=30) / 60

# two different substrate concentrations
# 3 different substrates

# su = "LLVY"
# suc = 50

x11()
par(mfrow = c(3,2))

for (su in c("LLE", "VGR", "LLVY")) {
  
  for (suc in c(50, 100)) {
    
    cnt = annotation[annotation$substrate == su & annotation$substrate_conc == suc, ]
    plot(x=time,
         y=rowMeans(FLUORO[,colnames(FLUORO) == cnt$cell[cnt$protea_conc == 50]]),
         type = "l", col = "palevioletred4",
         xlab = "time [min]", ylab = "RFU",
         main = paste0(su, ", concentration: ", suc, " uM"))
    lines(x=time,
          y=rowMeans(FLUORO[,colnames(FLUORO) == cnt$cell[cnt$protea_conc == 20]]),
          col = "orchid4")
    lines(x=time,
          y=rowMeans(FLUORO[,colnames(FLUORO) == cnt$cell[cnt$protea_conc == 10]]),
          col = "plum2")
    legend("topleft",
           legend = c("50 nM protea", "20 nM", "10 nM"),
           lty = rep("solid",3), col = c("palevioletred4","orchid4","plum2"),
           cex=.5)
    
  }
  
}

dev.copy2eps(file = "results/protea_conc.ps")



