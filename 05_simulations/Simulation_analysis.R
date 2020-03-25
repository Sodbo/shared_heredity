############################## SIMULATIONS #######################################
rm(list=ls())

library(rsvd) ###
library(Rsolnp) ###
library(psych)
library(fBasics)###
library(ggplot2)
library(gridExtra)
library(Matrix)

###################################################################
###     generalization of the 'truncated' normal value          ###
###################################################################
rtnorm <- function(n=1, mean = 0, sd = 0.5, min = 0, max = 1) {
    bounds <- pnorm(c(min, max), mean, sd)
    u <- runif(n, bounds[1], bounds[2])
    qnorm(u, mean, sd) 
}

#####################################
###              U %^% k          ###
#####################################
"%^%" <- function(U, k) {  
	UUU <- eigen(U, symmetric = TRUE)  # UUU = Uvec %*% diag(Uval) %*% t(Uvec)
	Uval <- UUU$val; Uvec <- UUU$vec
	Uvec <- Uvec[, Uval > 1.e-7];	Uval <- Uval[Uval > 1.e-7]
	Uvec %*% (t(Uvec) * (Uval ^ k))   #	Uvec %*% (diag(Uval ^ k) %*% t(Uvec))
}

###############################################################
###         testing whether there Are  shared SNPs?         ###                   
###############################################################
Test.SharedSNP.Old <- function(AAA,threshold=threshold.cor){ # significant threshold for correlation coefficients
if (min(abs(AAA)) > threshold) {
	if (qr(sign(AAA))$rank == 1) {
		return(TRUE)   #		print('there are shared SNPs!!! ');
	} else {
		return(FALSE)  #		print('there are no shared SNPs!!! ');
	}
	} else {
		return(FALSE)  #		print('there is no significant correlation coefficient')
		}
}

###############################################################
###         testing whether there Are  shared SNPs? (New)   ###                   
###############################################################
Test.SharedSNP <- function(AAA,threshold = threshold.cor){ # significant threshold for correlation coefficients
AAA <- ifelse(abs(AAA)<threshold,0,AAA)
if (qr(sign(AAA))$rank == 1) return(TRUE) else return(FALSE)
}


#####################################################################################
### minimization of loss function  L2:max.Lh, (or L1:min.sum.squared.residuals)   ###
#####################################################################################
fun <- function(w) {
	D <- (w %*% t(w) + diag(1.- w^2)) %*% abs.InvCorGenTr
	detD <- det(D)
	if (detD < 0) {
		return(1e+16)
		} else {
		return( sum(diag(D)) - log(detD) - Ntr)
	}
}

OPTIM <- function(initial=initial,Upper.bound=0.9){
	eps <- 1e-5
	ctrl <- list(tol=1e-8, trace=0)
	if((eps <= min(initial)) & (max(initial)<=Upper.bound)) {pars <- initial} else { pars <- rep(Upper.bound*0.5,Ntr) }
	w <- solnp(pars = pars, fun = fun,LB = rep(eps,Ntr), UB = rep(Upper.bound,Ntr),control=ctrl)$pars
	return(list(w = w, fun_w = fun(w)))
}

############################################################
###            generation of random matrix Upriv         ###
############################################################
generation.random.matrix <- function(Ntr,Sd,prob.zero=prob.zero){
		#mysamp <- rtnorm(Ntr^2, mean = 0, sd = Sd, min = -1, max = 1)  ### normal distribution, hist(mysamp)
		mysamp   <- runif( n = Ntr^2, min = -1*Sd, max = Sd) 			### uniform distribution
		Num0 <- round( (Ntr^2) * prob.zero )                 			### number of 0s in Upriv
		x    <- c(rep(0,Num0),rep(1,Ntr^2 - Num0))           			### vector from 0s and 1s
        mysamp   <- mysamp * sample(x) 
		
		Upriv <- matrix( mysamp, Ntr, Ntr)
		Upriv[lower.tri(Upriv)] <- t(Upriv)[lower.tri(Upriv)]  # making the symmetry for matrix
		diag(Upriv) <- rep(1,Ntr)
		return(Upriv)
}


###################################################################
###   MAIN: test of accuracy of W-matrix based on simulations   ###
###################################################################

Ntr           <- 4       ### number of traits
NNN           <- 100   ### number of simulations for each scenario
threshold.cor <- 1.e-1   ### threshold for significant correlations
prob.zero     <- 0.8     ### fraction of 0s in Upriv
Upper.bound   <- 0.95     ### upper bound of W
Traits        <- paste0('Trait.',c(1:Ntr))			 
weigth_sim    <- seq(from = 0.1, to = 0.9, by = 0.1)              ### trait loadings
sd_cor        <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8)  ### SD for private corelations

#path <- "/home/lima/Gulya/MAXHERED/"
path <- 'F:/ICG2019/PAPER_SUMMARY_2018/Yasha2018/MAXHERED/'
path <- paste0(path,'Ntr_',Ntr,'_NNN_',NNN,'_threshold_',threshold.cor,'_prob.zero_',prob.zero,'_UB_',Upper.bound,'_W_equal/')
if (!dir.exists(path)) dir.create(path)


Res   <- c()
I     <- diag(rep(1,Ntr))       ### identity matrix
Ushar <- matrix(1,Ntr,Ntr)      ### shared correlation matrix
	
for (wei in weigth_sim) {
for (Sd in sd_cor) {
	simW <- c()
	vecW <- rep(wei,Ntr);  ##vecW[1] <- 1-wei ########################################
	vecB <- sqrt(1-vecW^2)
	Matrix1 <- t(Ushar * vecW)*vecW      # sqrt(sqW2) %*% Ushar %*% sqrt(sqW2)
	Num_sim <- NNN
	for (iii in (1:Num_sim)){
		Upriv <- generation.random.matrix(Ntr,Sd,prob.zero)

		if (Test.SharedSNP(Upriv,threshold = threshold.cor) == FALSE) {   # checking for the existence of shared SNPs
		    if (min(eigen(Upriv, sym = TRUE, only.values = TRUE)$val) >= 0) { #  checking for positive determinant of matrix

				Ugen <- Matrix1 + t(Upriv * vecB)*vecB  # B %*% Upriv %*% B
				
				if (Test.SharedSNP(Ugen,threshold=threshold.cor) == TRUE) {   # checking for the existence of shared SNPs
					if (min(eigen(Ugen,sym=TRUE,only.values=TRUE)$val) >= 0){ #  checking for positive determinant of matrix

						CorGenTr <- Ugen
						abs.InvCorGenTr  <- (abs(CorGenTr)) %^% -1
						
						### OPTIM decomposition
						PCA <- rpca(CorGenTr%^%.5, k=1, center = FALSE, scale = FALSE) #  only PC1
						weig = as.vector(PCA$x)
						test <- OPTIM(initial = abs(weig),Upper.bound=Upper.bound)
						W <- test$w * sign(weig) 
						names(W) <- Traits
						simW <- cbind(simW,as.vector(W))
						
					} ### end testing Ugen
				} ### end testing Ugen
			} ### end testing Upriv
		} ### end testing Upriv
	
	} ### end iii
	
	if (is.null(simW)==FALSE) {
		simW <- as.matrix(simW)
		if (dim(simW)[2]>2) {
			Res <- rbind(Res, c(wei,Sd,rowMeans(abs(simW))-vecW, rowSds(abs(simW))/sqrt(dim(simW)[2]),dim(simW)[2],Num_sim  ) )
		} 	
	}
}
}

colnames(Res) <- c('weigt_nominal', 'SD_noise', paste0('mean.',names(W)), paste0('se.',names(W)),'Nsim','Num_sim')

save(Res, file = paste0(path,"simdata_Ntr=",Ntr,".RData"))
write.table(Res,file=paste0('simdata.Ntr=',Ntr,'.out'),quo=F,row=F,app=T,col=F)

Res000 <- Res
Res <- as.data.frame(Res)
Res$weigt_nominal <- as.factor(Res$weigt_nominal)
Res$SD_noise <- as.factor(Res$SD_noise)


#################
### Figure  1 ###
################# 
pdf(paste0(path,'numbers_of_events_Ntr_',Ntr,'.pdf'),width = 12, height = 8)
	p <- ggplot(Res, aes(x= weigt_nominal, y=Nsim, fill=SD_noise))  
	p <- p + geom_bar(stat="identity", color="black",position=position_dodge()) 
	p1 <- p + theme_bw() + ggtitle("Nsim") 
	print(p1)
dev.off()


#################
### Figure  2 ###
################# 
pdf(paste0(path,'W_Weigths_noise_Ntr_',Ntr,'.pdf'),width = 24, height = 16)
par(mfrow=c(2,2))
for (kk in 1:Ntr){
	ME = paste0('mean.',names(W)[kk])
	SE = paste0('se.',names(W)[kk])

	min0=min(Res[,ME]) - 0.1
	max0=max(Res[,ME]) + 0.1
	
	col=rainbow(length(sd_cor))
	xxx=barplot(Res[,ME],main=ME,col=col,xlab='nominal_weights',ylab = 'est_Weights',ylim = c(min0,max0))
	segments(xxx, Res[,ME] - Res[,SE], xxx, Res[,ME] + Res[,SE], lwd=1, cex = 1.5 )
	legend("topright", legend=paste0('Sd=',sd_cor),col=col, lwd=5, cex=0.9)
	axis(side=1,at=(length(sd_cor)+2)*c(0:(length(weigth_sim)-1)),labels=paste0('w=',weigth_sim))
}
dev.off()


#################
### Figure  3 ###
################# 
Bar_Plots <- function(Res,event_prob=0.5,Num_sim){

	Res <- Res[Res[,"Nsim"] > event_prob*Num_sim,]
	pdf(paste0(path,'W_Weigths_noise_Ntr_',Ntr,'_ggplot.pdf'),width = 12, height = 8)
	#par(mfrow=c(2,2))
	Res$weigt_nominal=as.factor(Res$weigt_nominal)
	Res$SD_noise=as.factor(Res$SD_noise)
	xxx <- ifelse(round(Res$Nsim/Res$Num_sim,2)==1,'',round(Res$Nsim/Res$Num_sim,2))
	### Trait 1
		p <- ggplot(Res, aes(x= weigt_nominal, y=mean.Trait.1, fill=SD_noise))  
		p <- p + geom_bar(stat="identity", color="black",position=position_dodge()) 
		p <- p + geom_errorbar(aes(ymin=mean.Trait.1-se.Trait.1, ymax=mean.Trait.1+se.Trait.1), col='grey70',width=.5,position=position_dodge(.9)) 
		p <- p + theme_bw() + ggtitle("mean.Trait.1") 
		p1 <- p + geom_text(aes(label=xxx),position=position_dodge(width=0.9), vjust=-0.5,size=2.5, col='grey10',fontface = "bold")
	### Trait 2
		p <- ggplot(Res, aes(x= weigt_nominal, y=mean.Trait.2, fill=SD_noise))  
		p <- p + geom_bar(stat="identity", color="black",position=position_dodge()) 
		p <- p + geom_errorbar(aes(ymin=mean.Trait.2-se.Trait.2, ymax=mean.Trait.2+se.Trait.2), col='grey70',width=.5,position=position_dodge(.9)) 
		p <- p + theme_bw() + ggtitle("mean.Trait.2") 
		p2 <- p + geom_text(aes(label=xxx),position=position_dodge(width=0.9), vjust=-0.5,size=2.5, col='grey10',fontface = "bold")

	
grid.arrange(p1,p2,nrow=2)
dev.off()
}

Bar_Plots(Res,0,Num_sim)

