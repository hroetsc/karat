### INHIBITOR KINETICS ###
# description:  merge quantification results of Distiller
# input:        batches - table_matches and table_pept_int
# output:       table_matches and table_pept_int
# author:       AM, adapted by HR

library(XML)
library(gtools)

# skip this many lines in table_pept_int.txt files --> modify!!!
skip_param = 103

### INPUT ###
# table matches
matches = list.files("data/Distiller/0+2", pattern = "[0-9].html", full.names = T)
# table peptide ints
ints = list.files("data/Distiller/0+2", pattern = "[0-9].txt", full.names = T)


### MAIN PART ###
# loading HTML table matches
d = list()
for (i in 1:length(matches)){
  
  print(i)
  i1 = readHTMLTable(matches[i], header = TRUE, as.data.frame = TRUE, which = 2)
  d[[i]] = i1
  
}
dtbl = plyr::ldply(d)
write.table(dtbl, file = "qiSPI/INPUT/quantitation_results/TSN5_0+2_table_matches.txt",
            sep = "\t", row.names = FALSE, append = FALSE)



# loading txt table pept ints
l = list()
for (i in 1:length(ints)){
  
  print(i)
  i1 = read.table(ints[i], header = TRUE, sep = "\t", skip = skip_param)
  l[[i]] = i1
  
}
ltbl = plyr::ldply(l)

write.table(ltbl, file = "qiSPI/INPUT/quantitation_results/TSN5_0+2_table_pept_int.txt",
            sep = "\t", row.names = FALSE, append = FALSE)
