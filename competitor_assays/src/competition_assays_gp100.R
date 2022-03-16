### COMPETITOR ASSAYS ###
# description:  plot inhibition of fluorescence by 20S
# input:        fluorescense measurements (2 hrs, vevery 30s)
#               K654, K669 K647 + VGR + 20S protea
# output:       time course plots
# author:       HR

library(dplyr)
library(readxl)
library(ggplot2)
library(viridis)

### INPUT ###
data = readxl::read_excel("data/20220315.xlsx", sheet = 1, range = c("B35:CU276"))


### MAIN PART ###
# ----- annotation -----
# found in pipetting_scheme.xlsx
# 2022031g

annotation = data.frame(cell = paste0(rep(toupper(letters[1:8]), each = 12), rep(seq(1,12), 8)),
                        protea_conc = NA)

protea_conc = c(rep(20,38), rep(0,2))
substrate = c(rep("VGR", 40))
substrate_conc = c(rep(50,38), rep(0,2))
inhibitor = c(rep("K654", 12), rep("K669", 12), rep("K647", 12), rep("none",4))
inhibitor_seq = c(rep("RTKAWNRQL", 12), rep("RTKAWNRQ", 12), rep("AWNRQLYPEW", 12), rep("none",4))
inhibitor_conc = c(rep(rep(c(300, 300/2, 300/2^2, 300/2^3, 300/2^4, 300/2^5),each=2), 3), rep(0,4))

annotation = annotation[c(which(annotation$cell == "A1"):which(annotation$cell == "D4")),]
annotation$protea_conc = protea_conc
annotation = annotation %>%
  mutate(substrate = substrate,
         substrate_conc = substrate_conc,
         inhibitor = inhibitor,
         inhibitor_seq = inhibitor_seq,
         inhibitor_conc = inhibitor_conc)

data$Time[1:5]
time = seq(6, 6+3600*2, by=30) / 60

# ----- subtract 1st measurement -----
# to account for pipetting delay
FLUORO = as.matrix(data[, colnames(data) %in% annotation$cell])
colnames(FLUORO) = annotation$cell

min1 = FLUORO[1,]
FLUORO = FLUORO - matrix(rep(min1, each=nrow(FLUORO)), nrow = nrow(FLUORO))

# ----- plot inhibition curves -----
intens_noinh =  FLUORO[,colnames(FLUORO) %in% c("D1","D2")]
mean_noinh = rowMeans(intens_noinh)
sd_noinh = apply(intens_noinh,1,sd)

cols = magma(n = 7)
plot_bars = seq(1,214, by = 10)

png(file = "results/competitor_gp100Peps.png",
    height = 12, width = 6, units = "in", res = 300)
par(mfrow = c(3,1))
for (inh in c("K654", "K669", "K647")) {
  
  cnt = annotation[annotation$inhibitor == inh, ]
  inh_conc = cnt$inhibitor_conc %>% unique() %>% sort()
  
  plot(x=time,
       y=mean_noinh,
       type = "l", col = "black",
       xlab = "time [min]", ylab = "RFU",
       main = paste0(inh, ", ",annotation$inhibitor_seq[annotation$inhibitor == inh][1]),
       ylim = c(0,35e03))
  arrows(x0 = time[plot_bars], y0 = mean_noinh[plot_bars]-sd_noinh[plot_bars],
         x1=time[plot_bars],y1 = mean_noinh[plot_bars]+sd_noinh[plot_bars],
         length=0.03, angle=90, code=3, lty = "solid", lwd = 1, col = "black")
  # polygon(c(time, rev(time)),
  #         c(mean_noinh+sd_noinh, rev(mean_noinh-sd_noinh)),
  #         col = "gray", border = NA, alpha=.3)
  
  for (i in 1:length(inh_conc)) {
    
    intens = FLUORO[,colnames(FLUORO) %in% cnt$cell[cnt$inhibitor_conc == inh_conc[i]]]
    mean_intens = rowMeans(intens)
    sd_intens = apply(intens,1,sd)
    
    lines(x=time, y=mean_intens, col = cols[i+1])
    arrows(x0 = time[plot_bars], y0 = mean_intens[plot_bars]-sd_intens[plot_bars],
           x1=time[plot_bars],y1 = mean_intens[plot_bars]+sd_intens[plot_bars],
           length=0.03, angle=90, code=3, lty = "solid", lwd = 1, col = cols[i+1])
    
  }
  
  legend("topleft",
         legend = c("no inhibitor", paste0(inh_conc, " µM")),
         lty = rep("solid", length(inh_conc)+1),
         col = c("black", cols[2:c(length(inh_conc)+1)]))
  
}

dev.off()


### OUTPUT ###


