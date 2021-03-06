# Scrip for Paniw et al. XXXX  

# Function to create IPMs based on mean parameter values including site variation

# This function is very similar to the code provided in makeIPM.R

############################################################################

makeIPM.mu=function(mcmc.out,bins=10,disc=1){
  
  disc=disc

  # replace some characters in the column names
  colnames(mcmc.out) = gsub("four",">three",colnames(mcmc.out))
  
  attach(mcmc.out)

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
  

  # goCont
  
  goCont.fun <- function(pfs) {
    
    
    mu.goCont = exp(mean(a0.goCont)+mean(get(paste("a1.goCont.",pfs,sep=""))))
    
    return(mu.goCont/(1+mu.goCont))
  }
  
  
  
  # PERCENTAGE STAYING IN THE SEED BANK (staySB):
  
  staySB.fun <- function(pfs) {
    
    
    mu.staySB <- exp(mean(a0.staySB)+mean(get(paste("a1.staySB.",pfs,sep=""))))
    
    return(mu.staySB/(1+mu.staySB))
  }
  
  # PERCENTAGE GERMINATING OUT OF THE SEED BANK (outSB):
  
  outSB.fun <- function(pfs) {
    
 
    mu.outSB <- exp(mean(a0.outSB)+mean(get(paste("a1.outSB.",pfs,sep=""))))

    
    return(mu.outSB/(1+mu.outSB))
  }
  
  ###############################################################################
  
  # Make and save IPMs
  
  # data frame with the correction factor sigma.s, see main text
  
  corr=data.frame(sigma.s=c(0.84,0.84,0.82,0.8,0.82),
                  TSF=c("zero","one","two","three",">three"))
  
  c=0.18 # correction factor for germination (see Appendix 2 and main text)
  
  minsize=0 # minimum size
  maxsize=9.6 # maximum size
  
  IPMkernel<-function(n) { # n defines the size (rows and columns) of the kernel resulting from integrating the IPM 
    
    # Function to define the IPM kernel (the midpoints for the integration)
    b <- minsize+c(0:n)*(maxsize-minsize)/n # interval that each cell of the matrix covers 
    h <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
    
    # Define the entries for each TSF
    
    # The zero (right after fire) matrix
    if(cov=="zero"){
      outSB <- 0.81*(sigma.s) #first value obtained from greenhouse germination trial, second ontained from in-situ experiment
      
      staySB <- 0.1 # constant assumed
      
      goCont=0.01 # so that we don�t devide by zero! (see below)
      
      goSB=0
      
      S <-matrix(0,n,n) # survival, growth (G), and fecundity (FecALL) are all zero
      
      G <- matrix(0,n,n)
      
      Fec01=matrix(0,n,n)
      Fec02=matrix(0,n,n)
      Fec03=matrix(0,n,n)
      Fec04=matrix(0,n,n)
      
      R <- (t(outer(h,h,SDS.fun,cov="one",site)))# the relevant non-0 transition is the size distribution of seedlings
      
      # TSF 1
    }else if(cov=="one"){
      
      outSB <- outSB.fun(pfs="burned")*sigma.s
      
      staySB <-0.05
      
      goCont=0.01 # so that we don�t devide by zero! (see below)
      
      goSB=0
      
      S <- diag(S.fun(h,cov,site)) # Survival Matrix 
      
      G <- t(outer(h,h,GR.fun,cov,site)) # Growth Matrix
      #Recruits distribution
      
      R <- (t(outer(h,h,SDS.fun,cov2,site)))
      
      Fec01=matrix(0,n,n)
      Fec02=matrix(0,n,n)
      Fec03=matrix(0,n,n)
      Fec04=matrix(0,n,n)# no seeds produced, therefore 0 fecundity
      
      # scale G and R below so columns sum to 1
      G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
      
      
    }else{
      outSB=outSB.fun(pfs="unburned")*c
      
      staySB=staySB.fun(pfs)
      
      goCont=goCont.fun(pfs=pfs.goCont)*c
      
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
      
      # scale D and G so columns sum to 1
      G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
      
    }
    
    R <- R/matrix(as.vector(apply(R,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    
    return(list(S=S,G=G,Fec01=Fec01,Fec02=Fec02,Fec03=Fec03,Fec04=Fec04,R=R,meshpts=h,goCont=goCont,sigma.s=sigma.s,
                goSB=goSB,outSB=outSB,staySB=staySB))
    
  }
  
  # Decide which seed bank parameter to loop through
  
  par.sb=c("goCont","staySB","outSB")

  #Loop though the parameters, TSF categories, and sites, and save IPMs
  TSF.name=c("zero","one","two","three",">three")
  
  # The two follwoing vectors link the seed-bank related vital-rate values obtained experimentally
  # to TSF categories required for stochastic simulations
  
  pfs.name=c("NA","NA","burned","unburned","unburned")
  pfs.name.goCont=c("NA","NA","unburned","unburned","unburned")
  
  site.name=c("siteA","siteB","siteC","siteD","siteE")
  
  
  IPMs.mu.site=array(0,c(bins+disc,bins+disc,length(site.name),length(TSF.name))) # array to hold IPMs in
  
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
      
      Fmat.cont <- M$R%*%M$Fec01*M$Fec02*M$Fec03*M$Fec04*M$goCont*M$sigma.s
      Fmat.discr =rep(0,length(M$meshpts)+1)
      Fmat=cbind(Fmat.discr,rbind(diag(M$Fec01*M$Fec02*M$Fec03*M$Fec04*M$goCont*M$sigma.s)*M$goSB/(M$goCont),Fmat.cont))
      
      mat <-Pmat+Fmat
      IPMs.mu.site[,,s,t] =mat
      
    }
  }
 
  return(IPMs.mu.site)
}
