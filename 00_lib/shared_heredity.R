st=ls())
setwd('~/Dropbox/Dropbox/MV_PAIN/results/20181010/data/')
#setwd('D:/Yasha2018')
load('20181116_four_pains.RData')
CorPhenTr <- as.matrix(phe)
A0 <- as.matrix(rgs_cov)
h2 <- diag(A0)
CorGenTr <- cov2cor(A0)
Ntr <- length(h2)
CorPhenTr; CorGenTr; h2 #plot(hclust(as.dist(CorGenTr)))
library(Rsolnp)  #install.packages('Rsolnp')

#####################################
###              U %^% k          ###
#####################################
"%^%" <- function(U, k) {  
  UUU <- eigen(U, symmetric = TRUE)  # UUU = Uvec %*% diag(Uval) %*% t(Uvec)
  Uval <- UUU$val; Uvec <- UUU$vec
  Uvec <- Uvec[, Uval > 1.e-12];	Uval <- Uval[Uval > 1.e-12]
  Uvec %*% (t(Uvec) * (Uval ^ k))   #	Uvec %*% (diag(Uval ^ k) %*% t(Uvec))
}

ConVert <- function(k,nnn,Ntr){
  d<-c()
  for(i in 1:Ntr){d <- c(k %% nnn,d);	k <- k %/% nnn}
  return(d+1)
}

###############################################################
###         testing whether there Are  shared SNPs?         ###                   
###############################################################
Test.SharedSNP <- function(eps=0.1){ # significant threshold for correlation coefficients
  if (min(abs(CorGenTr)) > eps) {
    Z <- CorGenTr/abs(CorGenTr)
    if (qr(Z)$rank==1) {
      print('there are shared SNPs!!! ')
    }else {
      print('there are no shared SNPs!!! ')
    }} else print('there are no shared SNPs!!! ')
}

###############################################################
###  estimation of initial values for parameters            ###                   
###############################################################	
### Begin: parameter regulating Residuals Matrix	
INITIAL <- function(PosRes = TRUE,LossFun = 'L2'){
  ZZZ <- eigen(CorGenTr/abs(CorGenTr),sym=TRUE)$vec[,1]
  Znak_Weights <- ZZZ/abs(ZZZ)
  
  d <- p <- Loss.1 <- Loss.2 <- c()
  WW<-rep(NaN,Ntr)
  param <- sqrt(seq(0.005,1,by=0.1)); print(param)
  nnn <- length(param)
  Points <- (nnn^Ntr)-1
  if (PosRes==T) Begin <- 0 else Begin <- (-1)
  for (iii in 0:Points){
    dx <-ConVert(iii,nnn,Ntr)
    w <- as.vector(sapply(dx,function(x) param[x]) ) * Znak_Weights
    VVV <- w %*% t(w) + diag(1.- w^2)
    p.min <- min(CorGenTr - VVV + diag(Ntr)) # without diagonal elements
    p.max <- max(CorGenTr - VVV + diag(Ntr))
    if((Begin <= p.min) & (p.max <= 1)) {
      B <- VVV %*% InvCorGenTr
      Loss.1 <- c(Loss.1,sum(diag(B)) - log(abs(det(B))) - Ntr)
      Loss.2 <- c(Loss.2,sum(diag((B - diag(Ntr)) %^% 2)))
      d <- c(d,paste(round(w,2), collapse = ' '))
      WW <- cbind(WW,w)
    }	
  }
  names(Loss.1) <- names(Loss.2) <- d
  if(LossFun == 'L1') {
    qwe <- which(min(Loss.1) == Loss.1)
  } else {qwe <- which(min(Loss.2) == Loss.2)}
  initial <- WW[,qwe] 
  return(initial)
  #mat <- as.matrix(cbind(Loss.1,Loss.2));plot(mat[,2]);plot(mat[,1])
}

###########################################################################
###            minimization of loss function  L2                        ###
###########################################################################

fun <- function(w,InvCorGenTr = InvCorGenTr) {
  VVV <- w %*% t(w) + diag(1.- w^2)
  L2 <- sum(diag((VVV %*% InvCorGenTr - diag(Ntr)) %^% 2)); 	return(L2)
}

ineqfun <- function(w,InvCorGenTr = InvCorGenTr){
  VVV <- w %*% t(w) + diag(1.- w^2)
  min(CorGenTr - VVV + diag(Ntr))
}

OPTIM <- function(PosRes=TRUE){
  eps<-1e-7
  gr1 <- initial/abs(initial)
  gr2 <- rep(0,Ntr)
  gr12 <- cbind(gr1,gr2)
  gr12 <- cbind(apply(gr12,1,min),apply(gr12,1,max))
  if (PosRes==TRUE){
    res0 <- solnp(pars = initial, fun = fun, ineqfun = ineqfun, ineqLB = 0, ineqUB = 1-eps, 
                  LB = gr12[,1], UB = gr12[,2],InvCorGenTr = InvCorGenTr)
  }else{
    res0 <- solnp(pars = initial, fun = fun,
                  LB = gr12[,1], UB = gr12[,2],InvCorGenTr = InvCorGenTr)
  }
  w <- res0$pars
  return(list(w = w, fun_w = fun(w,InvCorGenTr),VVV = (w %*% t(w)) + diag(1.- w^2)))
}

#####################################################################
###            Maximization of shared heritability                ###
#####################################################################
MAXIM <- function(w){
  A0 <- CorGenTr   * sqrt(h2 %*% t(h2))             # total genotype corr.matrix component 
  A2 <- (w%*%t(w)) * sqrt(h2 %*% t(h2))			  # shared genotype corr.matrix component
  
  VarianceComp <- function(alphas){
    Phen = t(alphas) %*% CorPhenTr %*% alphas	
    total  = (t(alphas) %*% A0 %*% alphas)/Phen
    Shared = (t(alphas) %*% A2 %*% alphas)/Phen
    Unique = (t(alphas) %*% (A0-A2) %*% alphas)/Phen
    return(c(total,Shared,Unique,alphas))
  }
  
  a.1=eigen(((CorPhenTr - A0) %^% -1) %*% A0)$vec[,1] 				#max Gen/Phen
  a.2=eigen(((CorPhenTr - A2) %^% -1) %*% A2)$vec[,1] 				#max Shared/Phen
  a.3=eigen(((CorGenTr - (w%*%t(w))) %^% -1) %*% (w%*%t(w)))$vec[,1]	#max Shared/Gen
  
  study.1 <- VarianceComp(a.1);study.2 <- VarianceComp(a.2);study.3 <- VarianceComp(a.3)
  
  res <- rbind(study.1,study.2,study.3 )
  colnames(res) <- c('Gen/Phen','Shared/Phen','Unique/Phen',paste(' Alpha',1:length(h2),sep = '')) 
  rownames(res) <- c('max Gen/Phen','max Shared/Phen','max Shared/Gen')
  return(res)
}

InvCorGenTr  <- CorGenTr %^%-1
initial <- INITIAL(PosRes=TRUE,LossFun='L1')
test1 <- OPTIM(PosRes=TRUE)
resu1 <- MAXIM(test1$w)


test1$w; resu1


