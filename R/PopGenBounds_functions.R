#' Computes differentiation statistics
#'
#' @param freqs An allele frequency matrix, where rows are populations and columns alleles
#'
#' @return A numerical vector
#' @export
#'
#' @examples
#' freqs <- matrix(c(1,0.5,0,0.5),nrow=2)
#' Diff(freqs)
Diff = function(freqs){
  K=nrow(freqs)
  M = max(colMeans(freqs))
  HS = mean(1-rowSums(freqs**2))
  HT = 1-sum(colMeans(freqs)**2)
  FST = 1-HS/HT
  GpST= FST*(K-1+HS)/(K-1)/(1-HS)
  D = K/(K-1)*(HT-HS)/(1-HS)
  return(c(M=M,FST = FST,Gp=GpST,D=D))
}


#' Computes the maximal FST in terms of the most frequent allele
#'
#' @param K Number of subpopulations
#' @param M Frequency of the most frequent allele
#'
#' @return A numerical value
#' @export
#'
#' @examples
#' #' Fup(2,0.5)
Fup = function(K,M){
  s1 = K*M
  J = ceiling(1/s1)
  res = (K-1)*(1-s1*(J-1)*(2-J*s1))/(K-(1-s1*(J-1)*(2-J*s1)))*(s1<=1) +
    (K*(K-1)-s1**2+floor(s1)-2*(K-1)*(s1%%1)+(2*K-1)*(s1%%1)**2)/(K*(K-1) -s1**2-floor(s1)+2*s1-(s1%%1)**2)*(s1>1)
  return(res)
}


#' Computes the maximal D in terms of the most frequent allele
#'
#' @param K Number of subpopulations
#' @param M Frequency of the most frequent allele
#'
#' @return A numerical value
#' @export
#'
#' @examples
#' Dup(2,0.5)
Dup = function(K,M){
  s1 = K*M
  res = 1*(s1<=1) + (1-(s1**2-(s1%%1)**2-floor(s1))/(K-1)/(K-1+(s1%%1)**2+(1-s1%%1)**2) )*(s1>1)
}

#' Computes the maximal G'ST in terms of the most frequent allele
#'
#' @param K Number of subpopulations
#' @param M Frequency of the most frequent allele
#'
#' @return A numerical value
#' @export
#'
#' @examples
#' Gpup(2,0.5)
Gpup = function(K,M){
  s1 = K*M
  J = ceiling(1/s1)
  res = 1*(s1<=1) +
    (K*(K-1)-s1**2+floor(s1)-2*(K-1)*(s1%%1)+(2*K-1)*(s1%%1)**2)/(K*(K-1) -s1**2-floor(s1)+2*s1-(s1%%1)**2)*(s1>1)/(K-1)/(K-1+(s1%%1)^2+(1-s1%%1)^2 )*(K^2 - K+1-(s1%%1)^2-(1-s1%%1)^2 )
  return(res)
}


#' Plots differentiation statistics along with the maximal and minimal values of the statistic
#'
#' @param freqs An allele frequency matrix, where rows are populations and columns alleles
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' freqs <- matrix(c(1,0.5,0,0.5),nrow=2)
#' ggbounds(freqs)
ggbounds = function(freqs){
  K=nrow(freqs)
  nudge = (mean(freqs$M)<0.5)*0.16-(mean(freqs$M)>=0.5)*0.25
  MFtmp = tibble(M= seq(0,1,0.00001), FST= Fup(K,seq(0,1,0.00001)) )
  ggplot(freqs,aes(x=M,y=FST)) + geom_pointdensity() +
    geom_point(freqs=tibble(M=mean(toad.diffl[,1]),FST=mean(toad.diffl[,2])),col="red",
               pch=16,size=3,stroke=2) +
    geom_segment(freqs=tibble(M=mean(freqs$M),FST=mean(freqs$FST)),
                 aes(x=M,xend=M,y=0,yend=Fup(K,M)), col="red",size=1 ) +
    geom_label(freqs=tibble(M=mean(freqs$M),FST=mean(freqs$FST)),
               aes(x=M,y=FST,
                   label= paste0("mean FST=", format(FST,digits=2)," (",
                                 format(FST/Fup(K,M)*100,digits=2),"% of range)" )),
               nudge_x = nudge,col="red") +
    geom_line(data=MFtmp,aes(x=M,y=FST)) + xlab(expression(italic(M))) +
    ylab(expression(italic(F[ST]))) +
    coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) + scale_color_viridis_b() + theme_bw()
}
