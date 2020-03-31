
###########################################################################
### Figure  1: Number of adequate simulations under different scenarios ###
###########################################################################

 Nsim.Plots <- function(path,Res){
	pdf(paste0(path,'numbers_of_adequate_simulations.pdf'),width = 12, height = 8)
	p <- ggplot(Res, aes(x = weigt_nominal, y = Nsim, fill = SD_noise))  
	p <- p + geom_bar(stat = "identity", color = "black",position = position_dodge()) 
	p <- p + theme_bw() + ggtitle("Nsim")
	print(p)
	dev.off()
}

########################################################################
### Figure  2: barplots for mean W-weights under different scenarios ###
######################################################################## 

Bar.Plots <- function(path,Res,Ntr,sd_cor,weigth_sim){
Traits <- paste0('Trait.',c(1:Ntr))			 
pdf(paste0(path,'W_Weigths_noise.pdf'),width = 24, height = 16)
par(mfrow=c(2,2))
for (kk in 1:Ntr){
	ME = paste0('mean.',Traits[kk])
	SE = paste0('se.',Traits[kk])

	min0=min(Res[,ME]) - 0.1
	max0=max(Res[,ME]) + 0.1
	
	col=rainbow(length(sd_cor))
	xxx=barplot(Res[,ME],main=ME,col=col,xlab='nominal_weights',ylab = 'est_Weights',ylim = c(min0,max0))
	segments(xxx, Res[,ME] - Res[,SE], xxx, Res[,ME] + Res[,SE], lwd=1, cex = 1.5 )
	legend("topright", legend=paste0('Sd=',sd_cor),col=col, lwd=5, cex=0.9)
	axis(side=1,at=(length(sd_cor)+2)*c(0:(length(weigth_sim)-1)),labels=paste0('w=',weigth_sim))
}
dev.off()
}

###########################################################################
### Figure  3: bar.ggplots for mean W-weights under different scenarios ###
###########################################################################  

Bar.ggPlots <- function(path=path,Res=Res,event_prob=event_prob,Num_sim=Num_sim){

	Res <- Res[Res[,"Nsim"] > event_prob*Num_sim,]
	pdf(paste0(path,'W_Weigths_noise_ggplot.pdf'),width = 12, height = 8)
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

Test.SharedSNP <- function(AAA,threshold=threshold.cor){ # significant threshold for correlation coefficients
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

Test.SharedSNP.New <- function(AAA,threshold = threshold.cor){ # significant threshold for correlation coefficients
	AAA <- ifelse(abs(AAA)<threshold,0,AAA)
	if (qr(sign(AAA))$rank == 1) return(TRUE) else return(FALSE)
}


#####################################################################################
### minimization of loss function  L2:max.Lh, (or L1:min.sum.squared.residuals)   ###
#####################################################################################


OPTIM <- function(initial=initial,Upper.bound=0.9,Ntr=Ntr,abs.InvCorGenTr=abs.InvCorGenTr){
	fun <- function(w) {
		D <- (w %*% t(w) + diag(1.- w^2)) %*% abs.InvCorGenTr
		detD <- det(D)
		if (detD < 0) { return(1e+16) } else { return( sum(diag(D)) - log(detD) - Ntr)}
	}	
	
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

Main.Simulation <- function(path=path0,
							Ntr=Ntr,                     ### number of traits
							NNN=NNN,                     ### number of simulations for each scenario
							threshold.cor=threshold.cor, ### threshold for significant correlations
							prob.zero=prob.zero,         ### fraction of 0s in Upriv
							Upper.bound=Upper.bound,     ### upper bound of W
							all.equal.W=all.equal.W){    ### all W-weights are equal

	library(rsvd) 
	library(Rsolnp) 
	library(psych)
	library(fBasics)
	library(ggplot2)
	library(gridExtra)
	library(Matrix)

   
Traits        <- paste0('Trait.',c(1:Ntr))			 
weigth_sim    <- seq(from = 0.1, to = 0.9, by = 0.1)              ### trait loadings
sd_cor        <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8)  ### SD for private corelations


path <- paste0(path,'Ntr_',Ntr,'_NNN_',NNN,'_threshold_',threshold.cor,
	    '_prob.zero_',prob.zero,'_UB_',Upper.bound,'_W_equal_',all.equal.W,'/')
if (!dir.exists(path)) dir.create(path)

Res   <- c()
I     <- diag(rep(1,Ntr))       ### identity matrix
Ushar <- matrix(1,Ntr,Ntr)      ### shared correlation matrix
	
for (wei in weigth_sim) {
for (Sd in sd_cor) {
	simW <- c()
	vecW <- rep(wei,Ntr)  
	if(all.equal.W == FALSE) vecW[1] <- 1-wei ########################################
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
						test <- OPTIM(initial = abs(weig),Upper.bound=Upper.bound,Ntr=Ntr,abs.InvCorGenTr=abs.InvCorGenTr)
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
print(paste0("Results are written in directory: ", path))
print(Res)

save(Res, file = paste0(path,'results.RData'))
write.table(Res,file = paste0(path,'results.out'),quo=F,row=F,app=T,col=F)

Res <- as.data.frame(Res)
Res$weigt_nominal <- as.factor(Res$weigt_nominal)
Res$SD_noise <- as.factor(Res$SD_noise)

Nsim.Plots(path,Res)
Bar.Plots(path,Res,Ntr,sd_cor,weigth_sim)
Bar.ggPlots(path,Res,0.05,Num_sim)
}

