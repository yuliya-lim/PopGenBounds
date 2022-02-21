#functions
## compute statistics
Diff.M = function(freqs){
  K=nrow(freqs)
  M = max(colMeans(freqs))
  HS = mean(1-rowSums(freqs**2))
  HT = 1-sum(colMeans(freqs)**2)
  FST = 1-HS/HT
  GpST= FST*(K-1+HS)/(K-1)/(1-HS)
  D = K/(K-1)*(HT-HS)/(1-HS)
  return(c(M=M,FST = FST,GpST,D))
}

## bounds
Fup = function(K,M){
  s1 = K*M
  J = ceiling(1/s1)
  res = (K-1)*(1-s1*(J-1)*(2-J*s1))/(K-(1-s1*(J-1)*(2-J*s1)))*(s1<=1) +
    (K*(K-1)-s1**2+floor(s1)-2*(K-1)*(s1%%1)+(2*K-1)*(s1%%1)**2)/(K*(K-1) -s1**2-floor(s1)+2*s1-(s1%%1)**2)*(s1>1)
  return(res)
}

Dup = function(K,M){
  s1 = K*M
  res = 1*(s1<=1) + (1-(s1**2-(s1%%1)**2-floor(s1))/(K-1)/(K-1+(s1%%1)**2+(1-s1%%1)**2) )*(s1>1)
}

Gpup = function(K,M){
  s1 = K*M
  J = ceiling(1/s1)
  res = 1*(s1<=1) + #(K-1)*(1-s1*(J-1)*(2-J*s1))/(K-(1-s1*(J-1)*(2-J*s1)))*(s1<=1) +
    (K*(K-1)-s1**2+floor(s1)-2*(K-1)*(s1%%1)+(2*K-1)*(s1%%1)**2)/(K*(K-1) -s1**2-floor(s1)+2*s1-(s1%%1)**2)*(s1>1)/(K-1)/(K-1+(s1%%1)^2+(1-s1%%1)^2 )*(K^2 - K+1-(s1%%1)^2-(1-s1%%1)^2 )
  return(res)
}

# plotting functions
ggbounds = function(data=toad.tib,K=47){
  nudge = (mean(data$M)<0.5)*0.16-(mean(data$M)>=0.5)*0.25
  MFtmp = tibble(M= seq(0,1,0.00001), FST= Fup(K,seq(0,1,0.00001)) )
  ggplot(data,aes(x=M,y=FST)) + geom_pointdensity() +
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
}
