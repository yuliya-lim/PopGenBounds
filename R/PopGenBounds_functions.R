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
  if( !all(rowSums(freqs)==1) ){
    stop("Rows of freqs matrix do not sum to 1. Do rows correspond to populations and columns to alleles?")
  }
  
  K=nrow(freqs)
  M = max(colMeans(freqs))
  HS = mean(1-rowSums(freqs**2))
  HT = 1-sum(colMeans(freqs)**2)
  FST = 1-HS/HT
  GpST= FST*(K-1+HS)/(K-1)/(1-HS)
  D = K/(K-1)*(HT-HS)/(1-HS)
  
  stats = tibble(statistic=c("FST","G'ST","D"), value = c(FST,GpST,D) )
  
  # compute bounds
  bounds = tibble(statistic=c("FST","G'ST","D") , bound=c(Fup(K=K,M=M), Gpup(K=K,M=M), Dup(K=K,M=M) ) )
  # compute closeness to upper bound
  closeness = tibble( statistic=c("FST","G'ST","D"), bound_closeness = stats$value/bounds$bound)
  
  return(bind_rows(tibble(statistic="M",value=M),left_join(left_join(stats,bounds),closeness)))
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
  return(res)
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
#' @param FST Numerical vector of FST values for each locus
#' @param M Numerical vector of frequency of the most frequent allele for each locus
#' @param K Number of subpopulations
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' freqs_locus1 <- matrix(c(1,0.5,0,0.5),nrow=2)
#' freqs_locus2 <- matrix(c(1,0.8,0,0.2),nrow=2)
#' freqs_locus3 <- matrix(c(1,0.2,0,0.8),nrow=2)
#' data  <- rbind(Diff(freqs_locus1),Diff(freqs_locus2),Diff(freqs_locus2))
#' ggbounds(data[,1],data[,2])
ggbounds = function(M,FST,GpST=NULL,D=NULL,K=2){
  nudge = (mean(M)<0.5)*0.16-(mean(M)>=0.5)*0.25
  MFtmp = dplyr::tibble(M= seq(0.00001,1-0.00001,0.00001),
                        FST= Fup(K,seq(0.00001,1-0.00001,0.00001)) )
  MGptmp = dplyr::tibble(M= seq(0.00001,1-0.00001,0.00001),
                        GpST= Gpup(K,seq(0.00001,1-0.00001,0.00001)) )
  MDtmp = dplyr::tibble(M= seq(0.00001,1-0.00001,0.00001),
                        D= Dup(K,seq(0.00001,1-0.00001,0.00001)) )
  plotFST <- ggplot2::ggplot(dplyr::tibble(M=M,FST=FST),ggplot2::aes(x=M,y=FST)) + ggpointdensity::geom_pointdensity() +
    ggplot2::geom_line(data=MFtmp,ggplot2::aes(x=M,y=FST)) + ggplot2::xlab(expression(italic(M))) +
    ggplot2::geom_point(data=dplyr::tibble(M=mean(M),FST=mean(FST)),col="red",
               pch=16,size=3,stroke=2) +
    ggplot2::geom_segment(data=dplyr::tibble(M=mean(M),FST=mean(FST)),
                          ggplot2::aes(x=M,xend=M,y=0,yend=Fup(K,M)), col="red",size=1 ) +
    ggplot2::geom_label(data=dplyr::tibble(M=mean(M),FST=mean(FST)),
                        ggplot2::aes(x=M,y=FST,
                   label= paste0("mean FST=", format(FST,digits=2)," (",
                                 format(mean(FST)/Fup(K,mean(M))*100,digits=2),"% of range)" )),
               nudge_x = nudge,col="red") +
    ggplot2::ylab(expression(italic(F[ST]))) +
    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) +
    ggplot2::theme_bw()
  
  if(!is.null(GpST)){
    plotGpST <- ggplot2::ggplot(dplyr::tibble(M=M,GpST=GpST),ggplot2::aes(x=M,y=GpST)) + ggpointdensity::geom_pointdensity() +
      ggplot2::geom_line(data=MGptmp,ggplot2::aes(x=M,y=GpST)) + ggplot2::xlab(expression(italic(M))) +
    ggplot2::geom_point(data=dplyr::tibble(M=mean(M),GpST=mean(GpST)),col="red",
                        pch=16,size=3,stroke=2) +
    ggplot2::geom_segment(data=dplyr::tibble(M=mean(M),GpST=mean(GpST)),
                          ggplot2::aes(x=M,xend=M,y=0,yend=Gpup(K,M)), col="red",size=1 ) +
    ggplot2::geom_label(data=dplyr::tibble(M=mean(M),GpST=mean(GpST)),
                        ggplot2::aes(x=M,y=GpST,
                                     label= paste0("mean GpST=", format(GpST,digits=2)," (",
                                                   format(mean(GpST)/Gpup(K,mean(M))*100,digits=2),"% of range)" )),
                        nudge_x = nudge,col="red") +
    ggplot2::ylab(expression(italic(G[ST]))) +
    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) +
    ggplot2::theme_bw()
  }else{plotGpST = NULL}
  
  if(!is.null(D)){
  plotD <- ggplot2::ggplot(dplyr::tibble(M=M,D=D),ggplot2::aes(x=M,y=D)) + ggpointdensity::geom_pointdensity() +
    ggplot2::geom_line(data=MDtmp,ggplot2::aes(x=M,y=D)) + ggplot2::xlab(expression(italic(M))) +
    ggplot2::geom_point(data=dplyr::tibble(M=mean(M),D=mean(D)),col="red",
                        pch=16,size=3,stroke=2) +
    ggplot2::geom_segment(data=dplyr::tibble(M=mean(M),D=mean(D)),
                          ggplot2::aes(x=M,xend=M,y=0,yend=Dup(K,M)), col="red",size=1 ) +
    ggplot2::geom_label(data=dplyr::tibble(M=mean(M),D=mean(D)),
                        ggplot2::aes(x=M,y=D,
                                     label= paste0("mean D=", format(D,digits=2)," (",
                                                   format(mean(D)/Fup(K,mean(M))*100,digits=2),"% of range)" )),
                        nudge_x = nudge,col="red") +
    ggplot2::ylab(expression(italic(D))) +
    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) +
    ggplot2::theme_bw()
  }else{plotD=NULL}
  
  if(length(M)>2) plot <- plot + ggplot2::scale_color_viridis_b()
  return(list(plotFST,plotGpST,plotD) )
}

#' Plots differentiation statistics along with the maximal and minimal values of the statistic
#'
#' @param freqs a matrix of allele frequency
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' popa = matrix(c(0.3,0.3,0.3,0.1,rep(0,4*5),
#'   rep(0,4),0.3,0.3,0.3,0.1,rep(0,4*4),
#'   rep(0,4*2),0.3,0.3,0.3,0.1,rep(0,4*3),
#'   rep(0,4*3),0.3,0.3,0.3,0.1,rep(0,4*2),
#'   rep(0,4*4),0.3,0.3,0.3,0.1,rep(0,4),
#'   rep(0,4*5),0.3,0.3,0.3,0.1),ncol=24,byrow = T)
#' ggfreqtable(popa)
ggfreqtable = function(freqs){
  freqs.tmp = freqs
  colnames(freqs.tmp) = paste0(1:ncol(freqs.tmp))
  freqs.long = bind_cols(Subpopulation=paste0(1:6),as_tibble(freqs.tmp)) %>% pivot_longer(-Subpopulation,names_to = "Allele",values_to = "Frequency")

  freqs.long$Allele = factor(freqs.long$Allele,levels=unique(freqs.long$Allele))
  freqs.long$Subpopulation = factor(freqs.long$Subpopulation,levels=rev(unique(freqs.long$Subpopulation)) )

ggres = ggplot(freqs.long,aes(x=Allele,y=Subpopulation,fill=Frequency)) + geom_tile(color="black",size=0.25) + 
  geom_text(aes(label=Frequency))+
  theme_bw() + 
  scale_fill_gradientn(limits=c(0,1),colors=c("white","darkorange","darkmagenta"))+
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1)
  )

 return(ggres)
}