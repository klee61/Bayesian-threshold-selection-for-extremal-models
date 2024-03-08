rm(list = ls(all.names = TRUE))

library(SpatialExtremes)

#Posterior in log
lnPost <- function(y,u,xi,sigma){
	f = sum(dgpd(y[y>u],loc=u,scale=sigma,shape=xi,log=TRUE))-log(sigma)-log(1+xi)-0.5*log(1+2*xi); return(f)
}

#Metropolis-Hastings updates 
MHSim <- function(y,u,T,batch,lnPost){
	xi0=0.02; sig0 =8; # initial parameter values; xi0(shape) and sig0(scale)
	Wsig=2; Wxi=1; # initial proposal walk size for sig and xi
	lnP0=lnPost(y,u,xi0,sig0) # log-posterior for the initial parameter values
	accept=0; # number of accepted proposals 
	Para=matrix(,T,2); # MCMC outputs for paramer values [shape,scale]
	
	### Metropolis-Hastings update ####
	for (t in 1:T){  
	sig=exp(rnorm(1,log(sig0),Wsig)); xi=rnorm(1,xi0, Wxi); xi=xi*(xi>(-0.5))+xi0*(xi<(-0.5))
	lnP=lnPost(y,u,xi,sig)
	if(runif(1)<exp(lnP-lnP0)){sig0=sig; xi0=xi; lnP0 = lnP; accept=accept+1; }
	Para[t,]=c(sig0,xi0); 
	
	# readjust proposal walk size 
	if((t/batch)==floor(t/batch)){ 
		Wsig=max(0.1,2.38*sd(log(Para[(t-batch+1):t,1])))
		Wxi=max(0.1,2.38*sd(Para[(t-batch+1):t,2])) }		}
	
	out=list(Para=Para,accept); return(out)
}


#######  main script  ####### 

U = seq(0,40,by=5); # threshold cadidates
Ni0=30; # number replicates
N=5e3; # number of posterior samples used for the p-value estimate
MS=matrix(,length(U),Ni0) #posterior predictive p-value for Ni0 replicates

T=1e4; # Total number of itereation (length of Markov chain)
batch=1e3; # the proposal walk size is adjusted at every "batch"-iteration

for(i0 in 1:Ni0){

# simulating Ny-number of data points
Ny=300; rn=runif(Ny); 
y=(rn<0.3)*runif(Ny,0,20)+(rn>0.3)*rgpd(Ny,loc=20,scale=8,shape=0.2) #mixture of uniform and GPD (first example)


for(i in 1:length(U)){
	# Simulate T - posterior samples 
	out=MHSim(y,U[i],T,batch,lnPost)
	
	# Posterior predictive p-value calculation
	PosP=0; 
	for(j in 1:N){ ind=sample(c(2e3:T),1) 
		# simulated data using a posterior sample
	ysim=rgpd(sum(y>U[i]),loc=U[i],scale=out$Para[ind,1],shape=out$Para[ind,2])
	# Here the test statisitc is a likelihood
	Tsim=mean(dgpd(ysim,loc=U[i],scale=out$Para[ind,1],shape=out$Para[ind,2]))
	Tobs=mean(dgpd(y[y>U[i]],loc=U[i],scale=out$Para[ind,1],shape=out$Para[ind,2]))
	PosP= PosP+as.numeric(Tsim<Tobs)
}
MS[i,i0]=PosP/N #p-value 
}
}

boxplot(t(MS))
