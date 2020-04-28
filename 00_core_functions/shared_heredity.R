shared_heredity <- function(CovGenTr = NULL, CorPhenTr = NULL, CorGenTr=NULL, h2=NULL, UpperW=0.995){
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
	if (!require('rsvd')) install.packages('rsvd'); library('rsvd')
	if (!require('Rsolnp')) install.packages('Rsolnp'); library('Rsolnp')
	if (!require('psych')) install.packages('psych'); library('psych')

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

	OPTIM <- function(initial=initial,UpperW=UpperW){
		eps <- 1.e-5
		ctrl <- list(tol=1e-8, trace=0)
		LB <- rep(eps,Ntr)
		UB <- rep(UpperW,Ntr)
		if ((eps <= min(initial)) & (max(initial) <= UpperW)) pars <- initial else pars <- UB*0.5 
		w <- solnp(pars = pars, fun = fun,LB = LB, UB = UB, ,control=ctrl)$pars
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
	
  #UpperW             ### upper bound of W

	
	#CorPhenTr; CorGenTr; h2 
	Test.SharedSNP(CorGenTr,0.1)

	abs.InvCorGenTr  <- (abs(CorGenTr)) %^% -1

	### OPTIM decomposition
	PCA <- rpca(CorGenTr%^%.5, k=1, center = FALSE, scale = FALSE) #  only PC1
	weig = as.vector(PCA$x)
	test <- OPTIM(initial = abs(weig),UpperW = UpperW)
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
		
		Cov_to_GIP=CorGenTr*sqrt(h2 %*% t(h2))
		gips=add_gips(gcovm=as.matrix(Cov_to_GIP),phem=as.matrix(CorPhenTr),l=0)
	
		est <- list(weights=W,nrmse=nrmse,Loss_fun=test$fun_w ,Mean.Weigs.2=Mean.Weigs.2,res=res[[1]],alphas=res[[2]],GIPs=gips)
		return(est)
}


add_gips=function(gcovm,phem,ind=1,l=0){
  
	#loading functions
	
	H2=function(a,covm=covm,phem=phem){
	  h=sum(covm*(a%o%a))/sum(phem*(a%o%a))
	  return(h)
	}

	gen=function(a,covm=covm){
	  h=sum(covm*(a%o%a))
	  return(h)
	}

	phen=function(a,phem=phem){
	  h=sum(phem*(a%o%a))
	  return(h)
	}



	# phenotipic correlation of linear combination with a coeeficient with i trait
	cor_yi_alfa=function(a,i,phem=phem){
	  cov_gi_sum_giai=sum(phem[i,]*a)
	  var_gi=phem[i,i]
	  var_sum_giai=sum(phem*(a%o%a)) #gen(a)
	  
	  cor_g_a=cov_gi_sum_giai/sqrt(var_gi*var_sum_giai)
	  return(cor_g_a)
	}
	# phenotipic covariance of linear combination with a coeeficient with i trait
	cov_yi_alfa=function(a,i,phem=phem){
	  cov_gi_sum_giai=sum(phem[i,]*a)
	  var_gi=phem[i,i]
	  var_sum_giai=phen(a,phem=phem) #gen(a)
	  
	  cor_g_a=cov_gi_sum_giai
	  return(cor_g_a)
	}

	# phenotipic covariance a1 a2
	cov_yi_a1_a2=function(a1,a2,phem=phem){
	  cov_gi_sum_g_a1_a2=sum(phem*(a1%o%a2))
	  cor_g_a=cov_gi_sum_g_a1_a2
	  return(cor_g_a)
	}



	# genetic correlation of linear combination with a coeeficient a1 and a coeeficient a2
	cor_gi_a1_a2=function(a1,a2,covm=covm){
	  cov_gi_sum_g_a1_a2=sum(covm*(a1%o%a2))
	  var_g_a1=sum(covm*(a1%o%a1))
	  var_g_a2=sum(covm*(a2%o%a2))
	  
	  cor_g_a=cov_gi_sum_g_a1_a2/sqrt(var_g_a1*var_g_a2)
	  return(cor_g_a)
	}
	# genetic correlation of linear combination with a coeeficient with i trait
	cor_gi_alfa=function(a,i,covm=covm){
	  cov_gi_sum_giai=sum(covm[i,]*a)
	  var_gi=covm[i,i]
	  var_sum_giai=sum(covm*(a%o%a)) #gen(a)
	  
	  cor_g_a=cov_gi_sum_giai/sqrt(var_gi*var_sum_giai)
	  return(cor_g_a)
	}
	# genetic covaiance of linear combination with a coeeficient a1 and a coeeficient a2
	cov_gi_a1_a2=function(a1,a2,covm=covm){
	  cov_gi_sum_g_a1_a2=sum(covm*(a1%o%a2))
	  var_g_a1=sum(covm*(a1%o%a1))
	  var_g_a2=sum(covm*(a2%o%a2))
	  
	  cor_g_a=cov_gi_sum_g_a1_a2
	  return(cor_g_a)
	}

	# genetic covariacne of linear combination with a coeeficient with i trait
	cov_gi_alfa=function(a,i,covm=covm){
	  cov_gi_sum_giai=sum(covm[i,]*a)
	  var_gi=covm[i,i]
	  var_sum_giai=sum(covm*(a%o%a)) #gen(a)
	  
	  cor_g_a=cov_gi_sum_giai
	  return(cor_g_a)
	}

	H2_cov=function(a,covm=covm){
	  h=sum(covm*(a%o%a))
	  return(h)
	}

  
  
  gcovm_orig=gcovm
  
  #lambda=diag(nrow(gcovm))*(1-l)
  #gcovm=gcovm%*%lambda
  
  
  lambda=rep((1-l),nrow(gcovm))
  lambda=sqrt(lambda)%o%sqrt(lambda)
  diag(lambda)=1
  gcovm=gcovm*lambda
  
  
  
  eigens=eigen(gcovm)$vectors
  eigens=apply(eigens,MAR=2,FUN=function(x){if (x[ind]<0) {x=-x}; x})
  
  gcovm=gcovm_orig
  
  vars=apply(eigens,MAR=2,phen,phem=phem)
  sds=sqrt(vars)
  i=1
  for (i in 1:ncol(eigens)){
    eigens[,i]=eigens[,i]/sds[i]
  }
  
  #gen covs
  h2=diag(gcovm)
  nl=nrow(gcovm)
  
  out=array(NA,c(nl*2,nl*2))
  colnames(out)=rep("tr",nl*2)
  colnames(out)[1:nl]=colnames(gcovm)
  colnames(out)[(1+nl):(nl*2)]=paste0("GIP",1:nl)
  rownames(out)=colnames(out)
  out[1:nl,1:nl]=gcovm
  
  h2_cor=NULL
  i=1
  for (i in 1:nl){
    av=rep(0,n=nl)
    av=eigens[,i]
    out[nl+i,nl+i]=H2_cov(av,covm=gcovm)
    h2_cor=c(h2_cor,H2(av,covm=gcovm,phem=phem))
    j=1
    for (j in 1:nl){
      out[nl+i,j]=out[j,nl+i]=cov_gi_alfa(av,j,covm=gcovm)
    }
    j=1
    for (j in 1:nl){
      out[nl+i,j+nl]=out[j+nl,nl+i]=cov_gi_a1_a2(a1=av,a2=eigens[,j],covm=gcovm)
    }
  }
  
  #h2_cov_tmp=diag(out)
  #out[(1+nl):(nl*2),(1+nl):(nl*2)]=0
  #diag(out)=h2_cov_tmp
  out_cg=out
  
  #phen covs
  out=array(NA,c(nl*2,nl*2))
  colnames(out)=rep("tr",nl*2)
  colnames(out)[1:nl]=colnames(phem)
  colnames(out)[(1+nl):(nl*2)]=paste0("GIP",1:nl)
  rownames(out)=colnames(out)
  out[1:nl,1:nl]=as.matrix(phem)
  
  i=1
  for (i in 1:nl){
    av=rep(0,n=nl)
    av=eigens[,i]
    out[nl+i,nl+i]=phen(av,phem=phem)
    j=1
    for (j in 1:nl){
      out[nl+i,j]=out[j,nl+i]=cov_yi_alfa(av,j,phem=phem)
    }
    j=1
    for (j in 1:nl){
      out[nl+i,j+nl]=out[j+nl,nl+i]=cov_yi_a1_a2(av,a2=eigens[,j],phem=phem)
    }
  }
  out_cy=out
  
  vare=diag(out_cy)
  cor_y=out_cy*(sqrt(1/vare)%o%sqrt(1/vare))
  
  vare=diag(out_cg)
  cor_g=out_cg*(sqrt(1/vare)%o%sqrt(1/vare))
  
  rownames(eigens)=colnames(gcovm)
  colnames(eigens)=paste0("GIP",1:ncol(eigens))
  
  big_out=list(cov_g=out_cg,cov_y=out_cy,cor_g=cor_g,cor_y=cor_y,H2=diag(out_cg)/diag(out_cy),GIP_coeff=eigens)
  return(big_out)
}


