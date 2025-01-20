## code to prepare `DATASET` dataset goes here
#usethis::use_data(DATASET, overwrite = TRUE)

library(readxl)
library(tidyverse)

# Yellow-bellied toad example from Prohl et al 2021
# data application 
toad = read_xlsx("data-raw/ProhlEtAl2021_Bombina variegata Microsat data.xlsx")


Toad_microsat_freq = vector(mode="list",length = 10)
#toad.diffl = matrix(NA,10,4)
#toad.diff.pairl = c()
for(i in 1:10){ # for each locus
  tmp = toad[,(2+i*2):(3+i*2)]
  tmp[tmp==-9] = NA
  alleles = sort(na.exclude(unique(unlist(tmp))))
  tmp[[1]] = factor(tmp[[1]],levels=alleles)
  tmp[[2]] = factor(tmp[[2]],levels=alleles)
  colnames(tmp)=c("A1","A2")
  tmp2 = bind_cols(Population=toad[[2]],tmp) %>% pivot_longer(A1:A2)
  Freq.Mat = table(tmp2$Population,tmp2$value)
  Freq.Mat = sweep(Freq.Mat,1,rowSums(Freq.Mat),"/")
  Toad_microsat_freq[[i]] = Freq.Mat
    #toad.diffl[i,]=Diff.M(Freq.Mat)
  
  #for(j in 1:46){
  #  for(k in 2:47){
  #    toad.diff.pairl = rbind(toad.diff.pairl ,c(i,j,k,Diff.M(Freq.Mat[c(j,k),])))
  #  }
  #}
}
names(Toad_microsat_freq) = paste0("locus_",colnames(toad)[seq(4,23,2)])

use_data(Toad_microsat_freq)

K=47

##ggplot version
library(ggpointdensity)

toad.tib = tibble(M=toad.diffl[,1],FST=toad.diffl[,2],GpST=toad.diffl[,3],D=toad.diffl[,4])


ggbounds = function(data=toad.tib,K=47){
  nudge = (mean(data$M)<0.5)*0.16-(mean(data$M)>=0.5)*0.25
  MFtmp = tibble(M= seq(0,1,0.00001), FST= Fup(K,seq(0,1,0.00001)) )
  resF = ggplot(data,aes(x=M,y=FST)) + geom_pointdensity() + 
    geom_point(data=tibble(M=mean(toad.diffl[,1]),FST=mean(toad.diffl[,2])),col="red",
               pch=16,size=3,stroke=2) + 
    geom_segment(data=tibble(M=mean(data$M),FST=mean(data$FST)), 
                 aes(x=M,xend=M,y=0,yend=Fup(K,M)), col="red",size=1 ) + 
    geom_label(data=tibble(M=mean(data$M),FST=mean(data$FST)),
               aes(x=M,y=FST,
                   label= paste0("mean FST=", format(FST,digits=2)," (",
                                 format(FST/Fup(K,M)*100,digits=2),"% of range)" )), 
               nudge_x = nudge,col="red") +
    geom_line(data=MFtmp,aes(x=M,y=FST)) + xlab(expression(italic(M))) +
    ylab(expression(italic(F[ST]))) + 
    coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) + scale_color_viridis_b() + theme_bw()
  
  resG=ggplot(data,aes(x=M,y=GpST)) + geom_pointdensity() + 
    geom_point(data=tibble(M=mean(toad.diffl[,1]),FST=mean(toad.diffl[,2])),col="red",
               pch=16,size=3,stroke=2) + 
    geom_segment(data=tibble(M=mean(data$M),GpST=mean(data$GpST)), 
                 aes(x=M,xend=M,y=0,yend=Gpup(K,M)), col="red",size=1 ) + 
    geom_label(data=tibble(M=mean(data$M),GpST=mean(data$GpST)),
               aes(x=M,y=FST,
                   label= paste0("mean FST=", format(FST,digits=2)," (",
                                 format(FST/Fup(K,M)*100,digits=2),"% of range)" )), 
               nudge_x = nudge,col="red") +
    geom_line(data=MFtmp,aes(x=M,y=FST)) + xlab(expression(italic(M))) +
    ylab(expression(italic(F[ST]))) + 
    coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) + scale_color_viridis_b() + theme_bw()
  
}

ggbounds(toad.tib)
