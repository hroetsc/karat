
dt = read.table("qiSPI/INPUT/quantitation_results/table_pept_int.txt", header = T)
k = which(dt$Sequence == "AHSSSAFTITDQVPFSVSVSQLRALDGGNK" & dt$Intensity.Ref. < dt$Intensity.C2.)
length(k) / length(which(dt$Sequence == "AHSSSAFTITDQVPFSVSVSQLRALDGGNK"))

# no inhibitor
kk = which(dt$Intensity.Ref. > dt$Intensity.C2.)
length(kk) / nrow(dt)

kk = which(dt$Intensity.C1. > dt$Intensity.C3.)
length(kk) / nrow(dt)

