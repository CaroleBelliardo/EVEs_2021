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
library(RColorBrewer)
display.brewer.all()

## read tree from string
text.string<-
  "((Limulus_polyphemus,(Centruroides_exilicauda,Parasteatoda_tepidariorum)),((Strigamia_maritima),(((Armadillidium_vulgare),((Locusta_migratoria),(((Oncopeltus_fasciatus), Acyrthosiphon_pisum),(((Bombus_terrestris), Apis_mellifera),((Nicrophorus_vespilloides, (Tribolium_castaneum, Diabrotica_virgifera)),
  ((Plutella_xylostella,(Heliconius_melpomene)),
  (Aedes_aegypti,(Musca_domestica,((Drosophila_virilis), Drosophila_melanogaster))))))))))));"

vert.tree<-read.tree(text=text.string)

c=brewer.pal(n = 4, name = "PRGn")
## ---
vec=rev(vert.tree$tip.label)

d1=read.table('Bureau/clement_Avril2021/Len_eve_onlyG/sumEVE_exp-Noexp.txt2', sep='\t', stringsAsFactors=FALSE) #row.names = 1
names(d1)<-c( 'species', 'No expressed', 'Expressed' )
lab=c( 'No expressed', 'Expressed' )

data_new1 <- d1[match(vec, d1$species), ]

barplot(rbind(data_new1$`No expressed`, data_new1$Expressed), col=c( "grey", c[4]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend =lab, beside=FALSE)



## plot nb eve

## plot type
d2=read.table('Bureau/clement_Avril2021/Len_eve_onlyG/SynthesePlot1.txt', sep='\t',
              header = TRUE, stringsAsFactors=FALSE) #row.names = 1

data_new2 <- d2[match(vec, d2$species), ]

barplot(rbind(data_new2$pisRNA, data_new2$siRNA), col=c(c[1],c[2]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend = colnames(data_new2[,4:5]), beside=FALSE)


##--- DRAWING
par(mfrow=c(1,3))
par(mar = c(3, 1, 1, 1)) 

plot(vert.tree,no.margin=FALSE,cex=1.4)
barplot(rbind(data_new1$`No expressed`, data_new1$Expressed), col=c( "grey", c[4]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend =lab, beside=FALSE,  space=rep(0.8,19))


barplot(rbind(data_new2$pisRNA, data_new2$siRNA), col=c(c[1],c[2]),horiz=TRUE, xpd=TRUE,las=2,cex.names=0.5,
        legend = colnames(data_new2[,4:5]),, beside=FALSE, space=rep(0.8,19))

