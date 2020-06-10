
#rm(list=ls())
setwd("/mnt/polyomica/projects/shared_heredity/gulnara_src/shared_heredity/05_simulations/")
source("Simulation_analysis.R")

# path to output directory
path0 <- "/mnt/polyomica/projects/shared_heredity/gulnara_src/Results/"
#path0 <- "/home/lima/Gulya/MAXHERED/"   
cor = 0.1
n = 10000
pr.zero <- c(0.3,0.8)
Ntrait <- c(3,4,5) #,6,10,25)
for (kk in Ntrait) {
	for (jj in pr.zero) {
		Main.Simulation(path=path0,Ntr=kk,NNN=n,threshold.cor=cor,prob.zero=jj,Upper.bound=0.995,all.equal.W=TRUE)
		Main.Simulation(path=path0,Ntr=kk,NNN=n,threshold.cor=cor,prob.zero=jj,Upper.bound=0.995,all.equal.W=FALSE)
	}
}



