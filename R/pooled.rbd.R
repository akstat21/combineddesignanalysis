
#' Pooled RCBD
#' @description Provides an ANOVA table for combined analysis of series of experiments conducted in RCBD over the locations of time
#' @param loc A vector of location (factor)
#' @param treat A vector of treatments in each experiment and over the locations (factor)
#' @param block A vector of blocks in each experiment and over the locations (factor)
#' @param Y A vector of output variable (numeric)
#'
#' @return An ANOVA table of combined analysis
#' @import agricolae
#' @export
#' @references
#' @examples
#' pooled.rcbd(loc=)

pooled.rcbd=function(loc,treat,block,Y){

  loc <- as.factor(loc)
  treat <- as.factor(treat)
  block <- as.factor(block)

  data=data.frame(loc,treat,block,Y)

  #Individual Analysis
  out=by(data,data[,"loc"], function(x) aov(Y~factor(block)+factor(treat),data=x))
  out=as.list(out)
  sapply(out,summary)
  lapply(out,TukeyHSD,"factor(treat)")
  MSE=sapply(1:6,function(i,out) sum(out[[i]]$residual^2)/out[[i]]$df.residual,out)

  residfs=sapply(1:6,function(i,out) out[[i]]$df.residual,out)

  #Check homogeneity of error variances across locations
  homogeneity.test=bartlett.test(out)

  homogeneity.test$p.value

  if(homogeneity.test$p.value<=0.05){
    tY=unlist(tapply(data$Y,data$loc,function(Y,MSE) Y/sqrt(MSE),MSE))
    lm2=aov(tY~loc+treat+loc/block-block+loc:treat,data=data)
    anova(lm2)
    A=anova(lm2)
    if (A$Pr[4]<=0.05) {

      A$"F value"=c(A$"F value"[1],A$"Mean Sq"[2]/A$"Mean Sq"[4],A$"F value"[3],A$"F value"[4],A$"F value"[5])
      A$"Pr(>F)"=c(A$"Pr(>F)"[1],1-pf(A$"F value"[2],A$Df[2],A$Df[4]),A$"Pr(>F)"[3],A$"Pr(>F)"[4],A$"Pr(>F)"[5])
      A
    } else {
      lm3=aov(tY~loc+treat+loc/block-block,data=data)
      anova(lm3)  }
  } else {
    lm1=aov(Y~loc+treat+loc/block-block+loc:treat,data=data)
    anova(lm1)
    if (anova(lm1)$Pr[4]<=0.05) {
      B=anova(lm1)
      B$"F value"=c(B$"F value"[1],B$"Mean Sq"[2]/B$"Mean Sq"[4],B$"F value"[3],B$"F value"[4],B$"F value"[5])
      B$"Pr(>F)"=c(B$"Pr(>F)"[1],1-pf(B$"F value"[2],B$Df[2],B$Df[4]),B$"Pr(>F)"[3],B$"Pr(>F)"[4],B$"Pr(>F)"[5])
      B
    } else {
      lm4=aov(Y~loc+treat+loc/block-block,data=data)
      anova(lm4)
    }
  }

}
