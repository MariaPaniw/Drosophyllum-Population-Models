# Scrip for Paniw et al. XXXX  

# Code to construct an IPM based on 1000 parameter values obtained via MCMC sampling (result of BayModel.R)

# Note that we create IPMs for each posterior parameter only for three vital rates (i.e., goSB, staySB, outSB);
# The remaining parameter are kept at their mean values 
# The extention to construct IPMs for the entire set of posterior samples is trivial:
# (a) Replace the mean function mean() in the vital-rate functions below with p (inside the empty spaces that index the parameter data frames [])
# And (b) Remove the outer for loop (defining par.choice) in the last part of the code  

rm(list=ls(all=TRUE))


############################################################################

# set working directory
setwd("YOUR_PATH")

### FIRST step: 

# load data

dataCont= read.csv("dataDroso.csv")
dataSB= read.csv("dataDrosoSB.csv")

# Order the levels of TSF categories

dataCont$TSF=factor(dataCont$TSF,levels = c("one", "two", "three",">three"))
dataSB$block=factor(dataSB$block)


# Load the posterior paramter values and reorder them to facilitate running the model functions

mcmc.out=read.csv("mcmcOUT.csv")[,-1]


# replace some characters in the column names
colnames(mcmc.out) = gsub("four",">three",colnames(mcmc.out))

attach(mcmc.out)

######################################################################################################
# Create density or probability distributions of above-ground vital rates

# SURVIVAL:

S.fun <- function(z,cov,site) {
 
    mu.surv=exp(mean(a0.surv)+mean(get(paste("a1.surv.",cov,sep="")))+mean(get(paste("aS.surv.",site,sep="")))+mean(bc.surv[])*z)
  
    return(mu.surv/(1+mu.surv))
}

# GROWTH (for simplicity we assume a constant variance)

GR.fun <- function(z,zz,cov,site){

    growth.mu=(mean(a0.gr)+mean(get(paste("a1.gr.",cov,sep="")))+mean(get(paste("aS.gr.",site,sep="")))+(mean(bc.gr[])+mean(get(paste("bcTSF.gr.",cov,sep=""))))*z)
    
    var.res=0.3986543 # constant variance based on data
    # Density distribution function of the normal distribution
    gr1 = sqrt(2*pi*var.res)
    gr2 = ((zz-growth.mu)^2)/(2*var.res)
  
  return(exp(-gr2)/gr1)
  
}

## SEEDLING SIZES (same approach as in growth function)

SDS.fun <- function(z,zz,cov,site){
  
  sds.mu=(mean(a0.sds)+mean(get(paste("a1.sds.",cov,sep="")))+mean(get(paste("aS.sds.",site,sep=""))))
  

  var.res=0.3329295# constant variance based on data

  # Density distribution function of the normal distribution
  sds1 = sqrt(2*pi*var.res)
  sds2 = ((zz-sds.mu)^2)/(2*var.res)
  
  return(exp(-sds2)/sds1)
  
}

# PROBABILITY OF FLOWERING 

FL.fun <- function(z,cov,site) {
  
  mu.fl=exp(mean(a0.fl)+mean(get(paste("a1.fl.",cov,sep="")))+mean(get(paste("aS.fl.",site,sep="")))+(mean(bc.fl[])+mean(get(paste("bcTSF.fl.",cov,sep=""))))*z)
 
  return(mu.fl/(1+mu.fl))
}

# NUMBER OF FLOWERING STALKS 

FS.fun <- function(z,cov,site) {
  
  mu.fs=exp(mean(a0.fs)+mean(get(paste("a1.fs.",cov,sep="")))+mean(get(paste("aS.fs.",site,sep="")))+mean(bc.fs[])*z)
  
  return(mu.fs)
}

# NUMBER OF FLOWERS PER STALK

FPS.fun <- function(z,cov,site) {
  
  mu.fps=exp(mean(a0.fps)+mean(get(paste("a1.fps.",cov,sep="")))+mean(get(paste("aS.fps.",site,sep="")))+mean(bc.fps[])*z)
  
  return(mu.fps)
}

##########################################################################
##creating IPMs for different parameter values
# The argument 'sb.choice' let's use which of the three vital rates to sample 
 

# IMMEDIATE GERMINATION (goCont):

goCont.fun <- function(p,pfs,sb.choice) {
  
  if(sb.choice=="goCont"){
    
    mu.goCont <- exp(a0.goCont[p]+get(paste("a1.goCont.",pfs,sep=""))[p]+0)}else{
      
      mu.goCont = exp(mean(a0.goCont[])+mean(get(paste("a1.goCont.",pfs,sep="")))+0)
    }

  return(mu.goCont/(1+mu.goCont))
}



# PERCENTAGE STAYING IN THE SEED BANK (staySB):

staySB.fun <- function(p,pfs,sb.choice) {
  
  if(sb.choice=="staySB"){
    
    mu.staySB <- exp(a0.staySB[p]+get(paste("a1.staySB.",pfs,sep=""))[p]+0)}else{
      
      mu.staySB <- exp(mean(a0.staySB[])+mean(get(paste("a1.staySB.",pfs,sep="")))+0)
    }
  
  return(mu.staySB/(1+mu.staySB))
}

# PERCENTAGE GERMINATING OUT OF THE SEED BANK (outSB):

outSB.fun <- function(p,pfs,sb.choice) {
  
  if(sb.choice=="outSB"){
    
    mu.outSB <- exp(a0.outSB[p]+get(paste("a1.outSB.",pfs,sep=""))[p]+0)}else{
      
      mu.outSB <- exp(mean(a0.outSB[])+mean(get(paste("a1.outSB.",pfs,sep="")))+0)
    }

  return(mu.outSB/(1+mu.outSB))
}

###############################################################################

# Make and save IPMs for (ideally all) posterior parameters 

# data frame with the correction factor sigma.s, see main text

corr=data.frame(sigma.s=c(0.84,0.84,0.82,0.8,0.82),
                TSF=c("zero","one","two","three",">three"))

c=0.18 # correction for germination (see Appendix 2 and main text)

minsize=0.9*min(dataCont$size[!is.na(dataCont$size)]) # minimum size
maxsize=1.1*max(dataCont$size[!is.na(dataCont$size)]) # maximum size

IPMkernel<-function(n) { # n defines the size (rows and columns) of the kernel resulting from integrating the IPM 

  # Function to define the IPM kernel (the midpoints for the integration)
  b <- minsize+c(0:n)*(maxsize-minsize)/n # interval that each cell of the matrix covers 
  h <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  
  # Define the entries for each TSF
  
  # The zero (right after fire) matrix
  if(cov=="zero"){
    outSB <- 0.81*(sigma.s) #first value obtained from greenhouse germination trial, second ontained from in-situ experiment
    
    staySB <- 0.1 # constant assumed
    
    goCont=0.01 # so that we don´t devide by zero! (see below)
    
    goSB=0
    
    S <-matrix(0,n,n) # survival, growth (G), and fecundity (FecALL) are all zero
    
    G <- matrix(0,n,n)
    
    FecALL=matrix(0,n,n)
    
    R <- (t(outer(h,h,SDS.fun,cov="one",site)))# the relevant non-0 transition is the size distribution of seedlings
    
    # TSF 1
    }else if(cov=="one"){
      
    outSB <- outSB.fun(p,pfs="burned",sb.choice)*sigma.s
   
    staySB <-0.05
    
    goCont=0.01 # so that we don´t devide by zero! (see below)
    
    goSB=0
    
    S <- diag(S.fun(h,cov,site)) # Survival Matrix 
    
    G <- t(outer(h,h,GR.fun,cov,site)) # Growth Matrix
    #Recruits distribution
    
    R <- (t(outer(h,h,SDS.fun,cov2,site)))
    
    FecALL=matrix(0,n,n)# no seeds produced, therefore 0 fecundity
    
    # scale G and R below so columns sum to 1
    G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    
    
  }else{
    outSB=outSB.fun(p,pfs="unburned",sb.choice)*c
    
    staySB=staySB.fun(p,pfs,sb.choice)

    goCont=goCont.fun(p,pfs=pfs.goCont,sb.choice)*c
    
    goSB=1-(goCont/c)-0.03
    
    S <- diag(S.fun(h,cov,site)) # Survival Matrix 
    
    G <- t(outer(h,h,GR.fun,cov,site)) # Growth Matrix
    
    #Recruits distribution
    R <- (t(outer(h,h,SDS.fun,cov2,site)))
  
    #Probability of flowering
    Fec01 = (diag(FL.fun(h,cov,site)))
    
    #Number of Flowering Stalks 
    Fec02 = (diag(FS.fun(h,cov,site)))
    
    #Number of flowers per Stalk
    
    Fec03= (diag(FPS.fun(h,cov,site)))
    
    #Number of seeds per flower that survive to become offspring 
    Fec04 = (diag(rep(9.1,n)))
    
    
    FecALL= Fec01*Fec02*Fec03*Fec04*goCont*sigma.s # add goCont and survSeedl to fecundity
    
    # scale D and G so columns sum to 1
    G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        
  }
 
  R <- R/matrix(as.vector(apply(R,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  
  return(list(S=S,G=G,FecALL=FecALL,R=R,meshpts=h,goCont=goCont,sigma.s=sigma.s,
              goSB=goSB,outSB=outSB,staySB=staySB))
  
}

# Decide which seed bank parameter to loop through

par.sb=c("goCont","staySB","outSB")

# Take a sample of the parameters (here we use just 10 posterior samples)
set.seed(12)

par.subset=sample(1:1000,10,replace=F)


#Loop though the parameters, TSF categories, and sites, and save IPMs
TSF.name=c("zero","one","two","three",">three")

# The two follwoing vectors link the seed-bank related vital-rate values obtained experimentally
# to TSF categories required for stochastic simulations

pfs.name=c("NA","NA","burned","unburned","unburned")
pfs.name.goCont=c("NA","NA","unburned","unburned","unburned")

site.name=as.character(unique(dataCont$site))

#decide how many bins you want to use (here only 10 to speed things up)

bins=10

#and how many discrete stages you have

disc=1

IPMs.pu=array(0,c(bins+disc,bins+disc,length(site.name),length(TSF.name),
                   length(par.sb),length(par.subset))) # array to hold IPMs in


for(sb in 1:length(par.sb)){ # loop through seed-bank vital rates
  sb.choice = par.sb[sb] 
    
  for(i in 1:length(par.subset)){ # loop through parameters
    p=par.subset[i]
    
    temp =NULL
    
    count=0
    
    for(t in 1:5){ # loop through TSF
      cov=TSF.name[t]
      
      #Associate IPM in TSF[t] with seedling size measured at TSF[t+1]
      # except for TSF >3
      
      if(t==5){cov2<-TSF.name[t]}else{cov2<-TSF.name[t+1]}
      
      sigma.s=corr[corr$TSF==cov,"sigma.s"]
      
      pfs=pfs.name[t]
      
      pfs.goCont=pfs.name.goCont[t]
      
      for(s in 1:5){ # loop through sites 
        site=site.name[s]
        
        M=IPMkernel(bins)
        
        ### Create P and F matrices including discrete stages
        
        Pmat.cont <- M$G%*%M$S
        Pmat.discr = c(M$staySB,M$outSB*M$R[,2])
        Pmat = cbind(Pmat.discr,rbind(rep(0,length(M$meshpts)),Pmat.cont))
        
        Fmat.cont <- M$R%*%M$FecALL
        Fmat.discr =rep(0,length(M$meshpts)+1)
        Fmat=cbind(Fmat.discr,rbind(diag(M$FecALL)*M$goSB/(M$goCont),Fmat.cont))
        
        mat <-Pmat+Fmat
        
        IPMs.pu[,,s,t,sb,i]=mat
        
      }
    }
    
    
    
      # Save output (recommended if creating IPMs for all parameters)
#       write.csv(temp,paste("PATH_TO_FOLDER/IPM_Bay/IPMlist",par.sb[sb],i,".csv",sep=""))
 
  }
  
}
