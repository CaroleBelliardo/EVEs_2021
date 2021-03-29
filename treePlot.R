r=read.table('Bureau/clement_Fev2021/Plot_len_reads_FG_FS.txt',header=TRUE, sep='\t')
install.packages('dplyr')
library(ggplot2)
library(dplyr)

install.packages("ape")
install.packages("phangorn")
install.packages("phytools")
install.packages("geiger")
install.packages("devtools")
install_github("liamrevell/phytools")

library(devtools)
library(ape)

par(mfrow=c(1,3))
## read tree from string
text.string<-
  "((Limulus_polyphemus,(Centruroides_exilicauda,Parasteatoda_tepidariorum)),((Strigamia_maritima),(((Armadillidium_vulgare),((Locusta_migratoria),(((Oncopeltus_fasciatus), Acyrthosiphon_pisum),(((Bombus_terrestris), Apis_mellifera),((Nicrophorus_vespilloides, (Tribolium_castaneum, Diabrotica_virgifera)),
  ((Plutella_xylostella,(Heliconius_melpomene)),
  (Aedes_aegypti,(Musca_domestica,((Drosophila_virilis), Drosophila_melanogaster))))))))))));"

vert.tree<-read.tree(text=text.string)
plot(vert.tree,no.margin=TRUE)


## ---
labels=rev(vert.tree$tip.label)

#labels=c('Lpol','Cexi','Ptep', 'Smar', 'Avul','Lmig', 'Ofas', 'Apis','Bter', 'Nves', 'Tcas', 'Dvir','Pxyl', 'Hmel', 'Aaeg','Dviri', 'Dmel')



## plot nb eve

## plot type

