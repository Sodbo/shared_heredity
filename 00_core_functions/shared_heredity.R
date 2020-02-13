shared_heredity <- function(CovGenTr = NULL, CorPhenTr = NULL, CorGenTr=NULL, h2=NULL){
	if(is.null(CorPhenTr)){
		stop('Error: The phenotype correlation matrix is not loaded.')
	}else{
		if(!is.null(CovGenTr) && is.null(CorGenTr) && is.null(h2)){
			h2 <- diag(CovGenTr)
			CorGenTr <- cov2cor(CovGenTr)
			print('The input gene covariation matrix is converted into the correlation matrix and the heritability vector')
		} else if (is.null(CovGenTr) && !is.null(CorGenTr) && !is.null(h2)){
			print('The gene correlation matrix and the heritability vector is loaded')
		} else {
			stop('Error: Incorrect input data. Please, input gene covariation matrix or both gene correlation matrix and heritability vector')
		}
	}

	if(any(colnames(CorGenTr)!=colnames(CorPhenTr))){
		write('Error: The names of traits are not identical for phenotype correlation and genotype covariance matrices: ',stderr())
		write.table(cbind(c('pheno_corr_matrix:','gene_cov_matrix:'),rbind(colnames(CorPhenTr),colnames(CorGenTr))), stderr(), col.names=F, row.names=F, quote=F)
		stop()
	}
	'input data:'
	CorPhenTr; CorGenTr; h2 #plot(hclust(as.dist(CorGenTr)))

	############################## DATA SET 2 #######################################
	library(rsvd)
	library(Rsolnp) 
	library(psych)

	#####################################
	###              U %^% k          ###
	#####################################
	"%^%" <- function(U, k) {  
		UUU <- eigen(U, symmetric = TRUE)  # UUU = Uvec %*% diag(Uval) %*% t(Uvec)
		Uval <- UUU$val; Uvec <- UUU$vec
		Uvec <- Uvec[, Uval > 1.e-12];	Uval <- Uval[Uval > 1.e-12]
		Uvec %*% (t(Uvec) * (Uval ^ k))   #	Uvec %*% (diag(Uval ^ k) %*% t(Uvec))
	}

	###############################################################
	###         testing whether there Are  shared SNPs?         ###                   
	###############################################################
	Test.SharedSNP <- function(CorGenTr,eps=0.1){ #eps: significant threshold for correlation coefficients
	if (min(abs(CorGenTr)) > eps) {
		Z <- sign(CorGenTr)
		if (qr(Z)$rank==1) {
			print('there are shared SNPs!!! ');return(TRUE)
		}else {
			print('there are no shared SNPs!!! ');return(FALSE)
		}} else {
			print('there is no significant correlation coefficient')
			return(FALSE)
			}
	}

	################################################################################
	### minimization of loss function  L1:max.Lh, L2:min.sum.squared.residuals   ###
	################################################################################
	fun <- function(w) {
		D <- (w %*% t(w) + diag(1.- w^2)) %*% abs.InvCorGenTr
		detD <- det(D)
		if (detD < 0) {
			return(1.e+16)
			}else{
			return(sum(diag(D)) - log(detD) - Ntr)
		}
	}

	OPTIM <- function(initial=initial){
		eps <- 1.e-7
		w <- solnp(pars = initial, fun = fun,LB = rep(0+eps,Ntr), UB = rep(1-eps,Ntr))$pars
		return(list(w = w, fun_w = fun(w)))
	}

	#####################################################################
	###            Maximization of shared heritability                ###
	#####################################################################

	MAXIM <- function(w,method){
		A0 <- CorGenTr   * sqrt(h2 %*% t(h2))             # total genotype corr.matrix component 
		A2 <- (w %*% t(w)) * sqrt(h2 %*% t(h2))			  # shared genotype corr.matrix component

		VarianceComp <- function(alphas){
			Phen = t(alphas) %*% CorPhenTr %*% alphas	
			total  = (t(alphas) %*% A0 %*% alphas)/Phen
			Shared = (t(alphas) %*% A2 %*% alphas)/Phen
			Unique = (t(alphas) %*% (A0-A2) %*% alphas)/Phen
			return(c(total,Shared,Unique,alphas))
		}

		S1 = (CorPhenTr - A0) %^% -0.5
		S2 = (CorPhenTr - A2) %^% -0.5
		
		a.1=eigen((S1 %*% A0) %*% S1,symmetric=TRUE)$vec[,1] 				#max Gen/Phen
		a.2=eigen((S2 %*% A2) %*% S2,symmetric=TRUE)$vec[,1] 				#max Shared/Phen

		study.1 <- VarianceComp(a.1)
		study.2 <- VarianceComp(a.2)

		#Checking for positive Alpha1 values and alpha vectors sign inversion if necessary 
    if(a.1[1]<0){
      a.1<-a.1*(-1)
    }
    if(a.2[1]<0){
      a.2<-a.2*(-1)
    }
		#alpha normalization
		phen=function(a,phem){
		  h=sum(phem*(a%o%a))
		  return(h)
		}    
		
		var_sh1=phen(a=a.1,phem=CorPhenTr)
		a.1=a.1/sqrt(var_sh1)
		var_sh2=phen(a=a.2,phem=CorPhenTr)
		a.2=a.2/sqrt(var_sh2)
		
		alphas <-rbind(a.1,a.2)
		res <- rbind(study.1,study.2)
		colnames(res) <- c('h2(Gen/Phen)','h2(Shared/Phen)','h2(Unique/Phen)',
						 paste(' Alpha.',rownames(A0),sep = ''))
		colnames(alphas) <- c(paste0(' Alpha.',rownames(A0)))
		rownames(res) <- paste0(method,c('.max(Gen/Phen)','.max(Shared/Phen)'))
		rownames(alphas) <- paste0(method,c('.max(Gen/Phen)','.max(Shared/Phen)'))
		return(list(res, alphas))
	}

	####################################################################
	###                   BEGIN                                      ###
	####################################################################



	### testing: sort by trait id
	#CorPhenTr <- CorPhenTr[,order(colnames(CorPhenTr))]
	#CorPhenTr <- CorPhenTr[order(rownames(CorPhenTr)),]
	#A0 <- A0[,order(colnames(A0))]
	#A0 <- A0[order(rownames(A0)),]

	### testing:  change of signs in cor.matrices
	#znak <- c(1,4)
	#CorPhenTr[,znak] <- CorPhenTr[,znak]*-1
	#CorPhenTr[znak,] <- CorPhenTr[znak,]*-1
	#A0[,znak] <- A0[,znak]*-1
	#A0[znak,] <- A0[znak,]*-1

	#h2 <- diag(A0)           # vector of heritabilities
	#CorGenTr <- cov2cor(A0)  # convert Cov to Cor
	Ntr <- length(h2)
	#CorPhenTr; CorGenTr; h2 
	Test.SharedSNP(CorGenTr,0.1)

	abs.InvCorGenTr  <- (abs(CorGenTr)) %^% -1

	### OPTIM decomposition
	PCA <- rpca(CorGenTr%^%.5, k=1, center = FALSE, scale = FALSE) #  only PC1
	weig = as.vector(PCA$x)
	test <- OPTIM(initial = abs(weig))
	W <- test$w * sign(weig); names(W) <- names(h2)
	res <- MAXIM(W,method='OPTIM')

	#############################################################################
	### Measures of quality of  reconstructed matrix                         ###
	### nrmse    : The normalized root mean squared error                     ###
	### Loss_fun : normal loss function value                                 ###
	#############################################################################
		CorGenTr.re <- W %*% t(W) + diag(1.- W^2)
		nrmse <- sqrt(sum((CorGenTr - CorGenTr.re) ^ 2 ) / sum(CorGenTr ^ 2)) 
		Mean.Weigs.2 <- mean(W^2)            # contribution of Shared SNPs to the genetic component
		est <- list(weights=W,nrmse=nrmse,Loss_fun=test$fun_w ,Mean.Weigs.2=Mean.Weigs.2,res=res[[1]],alphas=res[[2]])
		return(est)
}
