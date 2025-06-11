#' Plots differentiation statistics along with their ranges
#'
#' @param M Numerical vector of frequency of the most frequent allele for each locus
#' @param FST Numerical vector of FST values for each locus
#' @param GpST Numerical vector of G prime ST values for each locus
#' @param D Numerical vector of D values for each locus
#' @param K Number of subpopulations
#'
#' @returns A list of ggplot objects
#' @import ggplot2
#' @export
#'
#' @examples
#' library(tidyverse)
#' freqs_locus1 <- matrix(c(1,0.5,0,0.5),nrow=2)
#' freqs_locus2 <- matrix(c(1,0.8,0,0.2),nrow=2)
#' freqs_locus3 <- matrix(c(1,0.2,0,0.8),nrow=2)
#' data  <- rbind(Diff(freqs_locus1),Diff(freqs_locus2),Diff(freqs_locus2))
#' ggbounds_raw(M=data %>% filter(statistic=="M") %>% pull(value),
#'              FST=data %>% filter(statistic=="FST") %>% pull(value))
ggbounds_raw = function(M,FST,GpST=NULL,D=NULL,K=2){
  nudge = (mean(M,na.rm=T)<0.5)*0.16-(mean(M,na.rm=T)>=0.5)*0.25
  MFtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        FST= Fup(K,seq(0.001,1-0.001,0.001)) )
  MGptmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                         GpST= Gpup(K,seq(0.001,1-0.001,0.001)) )
  MDtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        D= Dup(K,seq(0.001,1-0.001,0.001)) )

  mean_M <- mean(M, na.rm = TRUE)
  mean_FST <- mean(FST, na.rm = TRUE)

  plotFST <-
    ggplot(dplyr::tibble(M=M,FST=FST),aes(x = M, y = FST)) +
    ggpointdensity::geom_pointdensity() +
    geom_line(data = MFtmp, aes(x = M, y = FST)) +
    geom_segment(data = dplyr::tibble(M = mean_M, FST = mean_FST),
                           aes(x = M, xend = M, y = 0, yend = 1),
                          col = "red", size = 0.8, linetype = "dashed") +
     geom_segment(data = dplyr::tibble(M = mean_M, FST_mean = mean_FST),
                           aes(x = 0, xend = 1, y = FST_mean, yend = FST_mean),
                          col = "red", size = 0.8, linetype = "dashed") +
     coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
     xlab(expression(italic(M))) +
     ylab(expression(italic(F[ST]))) +
     theme_bw()


  if(!is.null(GpST)){
    mean_GpST <- mean(GpST, na.rm = TRUE)

    plotGpST <-
       ggplot(dplyr::tibble(M=M,GpST=GpST), aes(x=M,y=GpST)) +
      ggpointdensity::geom_pointdensity() +
       geom_line(data=MGptmp, aes(x=M,y=GpST)) +
       geom_segment(data=dplyr::tibble(M=mean_M,GpST=mean_GpST),
                             aes(x=M,xend=M,y=0,yend=1),
                            col="red", size=0.8,linetype = "dashed") +
       geom_segment(data = dplyr::tibble(M = mean_M, GpST_mean = mean_GpST),
                             aes(x = 0, xend = 1, y = GpST_mean, yend = GpST_mean),
                            col = "red", size = 0.8, linetype = "dashed") +
       coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) +
       xlab(expression(italic(M))) +
       ylab(expression(italic(G[ST]))) +
       theme_bw()
  }
  else{plotGpST = NULL}

  if(!is.null(D)){
    mean_D <- mean(D, na.rm = TRUE)

    plotD <-
       ggplot(dplyr::tibble(M=M,D=D), aes(x=M,y=D)) +
      ggpointdensity::geom_pointdensity() +
       geom_line(data=MDtmp, aes(x=M,y=D)) +
       geom_segment(data=dplyr::tibble(M=mean_M,D=mean_D),
                             aes(x=M,xend=M,y=0,yend=1),
                            col="red", size=0.8,linetype = "dashed") +
       geom_segment(data = dplyr::tibble(M = mean_M, D_mean = mean_D),
                             aes(x = 0, xend = 1, y = D_mean, yend = D_mean),
                            col = "red", size = 0.8, linetype = "dashed") +
       coord_cartesian(xlim=c(0,1),ylim=c(0,1),expand = F) +
       xlab(expression(italic(M))) +
       ylab(expression(italic(D))) +
       theme_bw()
  }
  else{plotD=NULL}

  if(length(M)>2) plot <- plot +  scale_color_viridis_b()
  return(list(plotFST,plotGpST,plotD) )
}

#' Plots differentiation statistics normalized by their maximal values
#'
#' @inheritParams ggbounds_raw
#' @returns A list of ggplot objects
#' @export
#'
#' @examples
#' library(tidyverse)
#' freqs_locus1 <- matrix(c(1,0.5,0,0.5),nrow=2)
#' freqs_locus2 <- matrix(c(1,0.8,0,0.2),nrow=2)
#' freqs_locus3 <- matrix(c(1,0.2,0,0.8),nrow=2)
#' data  <- rbind(Diff(freqs_locus1),Diff(freqs_locus2),Diff(freqs_locus2))
#' ggbounds_norm(M=data %>% filter(statistic=="M") %>% pull(value),
#'               FST=data %>% filter(statistic=="FST") %>% pull(value))
ggbounds_norm = function(M,FST,GpST=NULL,D=NULL,K=2){
  nudge = (mean(M,na.rm=T)<0.75)*0.16-(mean(M,na.rm=T)>=0.75)*0.16
  MFtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        FST= Fup(K,seq(0.001,1-0.001,0.001)) )
  MGptmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                         GpST= Gpup(K,seq(0.001,1-0.001,0.001)) )
  MDtmp = dplyr::tibble(M= seq(0.001,1-0.001,0.001),
                        D= Dup(K,seq(0.001,1-0.001,0.001)) )
  # compute normalised FST
  FST_norm <- FST / sapply(M, function(m) Fup(K, m))
  MF_ST_tib <- dplyr::tibble(M=M,FST_n=FST_norm)

  mean_M <- mean(M, na.rm = TRUE)
  mean_FST_n <- mean(FST_norm, na.rm=T)

  plotFST_norm <-
     ggplot(MF_ST_tib,  aes(x = M, y = FST_n)) +
    ggpointdensity::geom_pointdensity() +
     geom_segment(data = dplyr::tibble(M = mean_M, FST = mean_FST_n),
                           aes(x = M, xend = M, y = 0, yend = 1),
                          col = "red", size = 0.8, linetype = "dashed") +
     geom_point(data=dplyr::tibble(M=mean_M,FST_n=mean_FST_n),
                        col="red", pch=16,size=3,stroke=2) +
     geom_label(data=dplyr::tibble(M=mean_M,FST_n=mean_FST_n),
                         aes(x=M,y=FST_n,
                                     label = paste0("mean FST=",
                                                    format(FST_n,digits=2))),
                        nudge_x = nudge,nudge_y=0,
                        col="red", size=3) +
     coord_cartesian(xlim = c(0.5, 1), ylim = c(0, 1), expand = F) +
     xlab(expression(italic(M))) +
     ylab(expression(italic(F[ST]))) +
     theme_bw()

  if(!is.null(GpST)){
    GST_norm <- GpST / sapply(M, function(m) Gpup(K, m))
    MG_ST_tib <- dplyr::tibble(M=M,GST_n=GST_norm)
    mean_GST_n <- mean(GST_norm, na.rm=T)

    plotGpST_norm <-
       ggplot(MG_ST_tib,  aes(x = M, y = GST_n)) +
      ggpointdensity::geom_pointdensity() +
       geom_segment(data = dplyr::tibble(M = mean_M, GST_n = mean_GST_n),
                             aes(x = M, xend = M, y = 0, yend = 1),
                            col = "red", size = 0.8, linetype = "dashed") +
       geom_point(data=dplyr::tibble(M=mean_M,GST_n=mean_GST_n),
                          col="red", pch=16,size=3,stroke=2) +
       geom_label(data=dplyr::tibble(M=mean_M,GST_n=mean_GST_n),
                           aes(x=M,y=GST_n,
                                       label = paste0("mean GpST=",
                                                      format(GST_n,digits=2))),
                          nudge_x = nudge,nudge_y=0,
                          col="red", size=3) +
       coord_cartesian(xlim = c(0.5, 1), ylim = c(0, 1), expand = F) +
       xlab(expression(italic(M))) +
       ylab(expression(italic(G[ST]))) +
       theme_bw()
  }
  else{
    plotGpST_norm = NULL
  }

  if(!is.null(D)){
    D_norm <- D / sapply(M, function(m) Dup(K, m))
    MD_tib <- dplyr::tibble(M=M,D_n=D_norm)
    mean_D_n <- mean(D_norm, na.rm=T)

    plotD_norm <-
       ggplot(MD_tib,  aes(x = M, y = D_n)) +
      ggpointdensity::geom_pointdensity() +
       geom_segment(data = dplyr::tibble(M = mean_M, D_n = mean_D_n),
                             aes(x = M, xend = M, y = 0, yend = 1),
                            col = "red", size = 0.8, linetype = "dashed") +
       geom_point(data=dplyr::tibble(M=mean_M,D_n=mean_D_n),
                          col="red", pch=16,size=3,stroke=2) +
       geom_label(data=dplyr::tibble(M=mean_M,D_n=mean_D_n),
                           aes(x=M,y=D_n,
                                       label = paste0("mean D=",
                                                      format(D_n,digits=2))),
                          nudge_x = nudge,nudge_y=0,
                          col="red", size=3) +
       coord_cartesian(xlim = c(0.5, 1), ylim = c(0, 1), expand = F) +
       xlab(expression(italic(M))) +
       ylab(expression(italic(D))) +
       theme_bw()
  }
  else{
    plotD_norm=NULL
  }

  if(length(M)>2) plot <- plot +  scale_color_viridis_b()
  return(list(plotFST_norm,plotGpST_norm,plotD_norm) )
}


#' Plot differentiation statistics values and their bounds (first column)
#' and the normalized values of statistics (second column)
#'
#'@inheritParams ggbounds_raw
#'
#' @returns A list of ggplot objects
#' @export
#'
#' @examples
#' library(tidyverse)
#' freqs_locus1 <- matrix(c(1,0.5,0,0.5),nrow=2)
#' freqs_locus2 <- matrix(c(1,0.8,0,0.2),nrow=2)
#' freqs_locus3 <- matrix(c(1,0.2,0,0.8),nrow=2)
#' data  <- rbind(Diff(freqs_locus1),Diff(freqs_locus2),Diff(freqs_locus2))
#' ggbounds_new(M=data %>% filter(statistic=="M") %>% pull(value),
#'              FST=data %>% filter(statistic=="FST") %>% pull(value))
ggbounds_new <- function(M,FST,GpST=NULL,D=NULL,K=2){
  gg1 <- ggbounds_raw(M, FST, GpST, D, K)
  gg2 <- ggbounds_norm(M, FST, GpST, D, K)
  return(list(gg1, gg2))
}
