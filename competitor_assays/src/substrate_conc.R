### COMPETITOR ASSAYS ###
# description:  determine Km values for the flurogenic substrates
# input:        fluorescense measurements (2 hrs, vevery 30s)
#               LLE, VGR, LLVY @ 6 different concentrations
# output:       Michaelis-Menten fits, Km values
# author:       HR

library(dplyr)
library(readxl)
library(latex2exp)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(drc)
library(matlib)
library(viridis)

### INPUT ###
data = readxl::read_excel("data/20220314_2.xlsx", sheet = 1, range = c("B35:CU276"))

### MAIN PART ###
# ----- annotation -----
# found in pipetting_scheme.xlsx
# 20220314_2
annotation = data.frame(cell = paste0(rep(toupper(letters[1:8]), each = 12), rep(seq(1,12), 8)),
                        protea_conc = NA)

protea_conc = c(rep(0, 6), rep(20, 3*12))
substrate = c(rep(c("LLE","VGR","LLVY"), each = 2), rep("LLE",12), rep("VGR",12), rep("LLVY",12))
substrate_conc = c(rep(400,6),
                   rep(rep(c(400, 400/3, 400/3^2, 400/3^3, 400/3^4, 400/3^5),each=2),3))


annotation = annotation[c(which(annotation$cell == "E7"):which(annotation$cell == "H12")),]
annotation$protea_conc = protea_conc
annotation$substrate = substrate
annotation$substrate_conc = substrate_conc

data$Time[1:5]
time = seq(6, 6+3600*2, by=30) / 60
plot_bars = seq(1,214, by = 10)

# ----- subtract 1st measurement -----
# to account for pipetting delay
FLUORO = as.matrix(data[, colnames(data) %in% annotation$cell])
colnames(FLUORO) = annotation$cell

min1 = FLUORO[1,]
FLUORO = FLUORO - matrix(rep(min1, each=nrow(FLUORO)), nrow = nrow(FLUORO))


# ----- plot all kinetics -----
png(file = "results/substrate_conc.png",
    height = 12, width = 6, units = "in", res = 300)
par(mfrow = c(3,1))
for (su in c("LLE", "VGR", "LLVY")) {
  
  cnt = annotation[annotation$substrate == su & annotation$protea_conc == 20, ]
  suc = cnt$substrate_conc %>% unique %>% sort(decreasing = T)
  cols = viridis(length(suc))
  
  intens = FLUORO[,colnames(FLUORO) == cnt$cell[cnt$substrate_conc == suc[1]]]
  mean_intens = rowMeans(intens)
  sd_intens = apply(intens,1,sd)
  
  plot(x=time,
       y=mean_intens,
       type = "l", col = cols[1],
       xlab = "time [min]", ylab = "RFU",
       main = su)
  arrows(x0 = time[plot_bars], y0 = mean_intens[plot_bars]-sd_intens[plot_bars],
         x1=time[plot_bars],y1 = mean_intens[plot_bars]+sd_intens[plot_bars],
         length=0.03, angle=90, code=3, lty = "solid", lwd = 1, col = cols[1])
  
  for (i in 2:length(suc)) {
    intens = FLUORO[,colnames(FLUORO) == cnt$cell[cnt$substrate_conc == suc[i]]]
    mean_intens = rowMeans(intens)
    sd_intens = apply(intens,1,sd)
    
    lines(x=time, y=mean_intens, col = cols[i])
    arrows(x0 = time[plot_bars], y0 = mean_intens[plot_bars]-sd_intens[plot_bars],
           x1=time[plot_bars],y1 = mean_intens[plot_bars]+sd_intens[plot_bars],
           length=0.03, angle=90, code=3, lty = "solid", lwd = 1, col = cols[i])
  }
  
  legend("topleft",
         legend = paste0(round(suc,4), " µM"),
         lty = rep("solid", length(suc)),
         col = cols)
  
  
}
dev.off()


# ----- get initial velocities -----
# take first x min to get initial velocities
# longer interval bc no MM behaviour
tlim = 20
V0 = data.frame(cell = colnames(FLUORO),
                v0 = (FLUORO[tlim*2+1, ] - FLUORO[2,])/(tlim))

kinetics = left_join(annotation, V0)

# remove F1 since there was not enough proteasome
kinetics = kinetics[-which(kinetics$cell == "F1"),]
# remove blanks
kinetics = kinetics[kinetics$protea_conc > 0, ]

# ----- plot Michaelis-Menten and Lineweaver-Burk -----
theme_set(theme_classic())
# Michaelis-Menten
MM = ggplot(kinetics, aes(x = substrate_conc, y = v0, color = factor(substrate))) +
  geom_point() +
  geom_smooth(aes(group = substrate), se=T) +
  scale_color_manual("[S] in µM",
                     values = c("midnightblue","mediumvioletred","mediumturquoise")) +
  ggtitle("", subtitle = "Michaelis-Menten") +
  xlab("[S] in µM") + ylab("v0")
MM

# Lineweaver-Burk
kinetics$v_rec = 1/kinetics$v0
kinetics$S_rec = 1/kinetics$substrate_conc

LB = ggplot(kinetics, aes(x = S_rec, y = v_rec, color = factor(substrate))) +
  geom_point() +
  geom_smooth(method = "glm", aes(group = substrate), fullrange=T, se=T) +
  scale_color_manual("[S] in µM",
                     values = c("midnightblue","mediumvioletred","mediumturquoise")) +
  ggtitle("", subtitle = "Lineweaver-Burk") +
  xlab("1/[S]") + ylab("1/v0") +
  xlim(-.1, .7) +
  stat_regline_equation()
LB

kinetic_plots = grid.arrange(MM, LB, nrow=1)
ggsave(plot = kinetic_plots, filename = "results/substrate_conc_kinetics.png",
       height = 5, width = 12, dpi = "retina")

# ----- determine Km -----
# just use the first 4 concentrations since the rest is a bit off
# fit hyerbolic model
png(file = "results/substrate_conc_KMvals.png",
    height = 8, width = 4, units = "in", res = 300)
par(mfrow = c(3,1))
for (su in c("LLE", "VGR", "LLVY")) {
  # hyperbolic fit
  mm.fit = drm(v0~substrate_conc,
               data = kinetics[kinetics$substrate == su, ],
               fct = MM.2())
  print(mm.fit)
  
  # Lineweaver-Burk fit
  lb.fit = lm(v0~substrate_conc,
              data = kinetics[kinetics$substrate == su, ])
  print(lb.fit)
  KM = lb.fit$coefficients[2]/lb.fit$coefficients[1]
  
  # plot the fit
  S0_fit = seq(0,500,0.5)
  v0_fit = (mm.fit$coefficients[1]*S0_fit) / (mm.fit$coefficients[2]+S0_fit)
  
  plot(v0_fit~S0_fit, type="l",
       xlab=TeX("$\\[S\\]_0$ in µM"),
       ylab=TeX("$v_0$"),
       main=su,
       # sub = paste0("KM = ", round(mm.fit$coefficients[2],4)), " (LB fit)",
       sub = paste0("KM = ",round(mm.fit$coefficients[2],4), " (MM fit)"))
  points(kinetics$v0[kinetics$substrate == su] ~ kinetics$substrate_conc[kinetics$substrate == su],
         pch = 16)
  abline(h = 0.5*mm.fit$coefficients[1], lty="dashed")
  abline(v = mm.fit$coefficients[2], lty="dashed")
  
}
dev.off()
