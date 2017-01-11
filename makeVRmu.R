# Scrip for Paniw et al. XXXX  

# Function to create vital-rate kernels based on mean parameter values 

############################################################################

makeVR.mu=function(mcmc.out,bins=10,discr=1, fire.freq=c(1/10,1/30,1/50,1/80)){

  # replace some characters in the column names
  colnames(mcmc.out) = gsub("four",">three",colnames(mcmc.out))
  
  attach(mcmc.out)

  # SURVIVAL:
  
  S.fun <- function(z,cov) {
    
    mu.surv=exp(mean(a0.surv)+mean(get(paste("a1.surv.",cov,sep="")))+mean(bc.surv[])*z)
    
    return(mu.surv/(1+mu.surv))
  }
  
  # GROWTH (for simplicity we assume a constant variance)
  
  GR.fun <- function(z,zz,cov){
    
    growth.mu=(mean(a0.gr)+mean(get(paste("a1.gr.",cov,sep="")))+(mean(bc.gr[])+mean(get(paste("bcTSF.gr.",cov,sep=""))))*z)
    
    var.res=0.3986543 # constant variance based on data
    # Density distribution function of the normal distribution
    gr1 = sqrt(2*pi*var.res)
    gr2 = ((zz-growth.mu)^2)/(2*var.res)
    
    return(exp(-gr2)/gr1)
    
  }
  
  ## SEEDLING SIZES (same approach as in growth function)
  
  SDS.fun <- function(z,zz,cov){
    
    sds.mu=(mean(a0.sds)+mean(get(paste("a1.sds.",cov,sep=""))))
    
    
    var.res=0.3329295# constant variance based on data
    
    # Density distribution function of the normal distribution
    sds1 = sqrt(2*pi*var.res)
    sds2 = ((zz-sds.mu)^2)/(2*var.res)
    
    return(exp(-sds2)/sds1)
    
  }
  
  # PROBABILITY OF FLOWERING 
  
  FL.fun <- function(z,cov) {
    
    mu.fl=exp(mean(a0.fl)+mean(get(paste("a1.fl.",cov,sep="")))+(mean(bc.fl[])+mean(get(paste("bcTSF.fl.",cov,sep=""))))*z)
    
    return(mu.fl/(1+mu.fl))
  }
  
  # NUMBER OF FLOWERING STALKS 
  
  FS.fun <- function(z,cov) {
    
    mu.fs=exp(mean(a0.fs)+mean(get(paste("a1.fs.",cov,sep="")))+mean(bc.fs[])*z)
    
    return(mu.fs)
  }
  
  # NUMBER OF FLOWERS PER STALK
  
  FPS.fun <- function(z,cov) {
    
    mu.fps=exp(mean(a0.fps)+mean(get(paste("a1.fps.",cov,sep="")))+mean(bc.fps[])*z)
    
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
      
      goCont=0.01 # so that we don´t devide by zero! (see below)
      
      goSB=0
      
      S <-matrix(0,n,n) # survival, growth (G), and fecundity are all zero
      
      G <- matrix(0,n,n)
      
      Fec01=matrix(0,n,n)
      Fec02=matrix(0,n,n)
      Fec03=matrix(0,n,n)
      Fec04=matrix(0,n,n)
      
      R <- (t(outer(h,h,SDS.fun,cov="one")))# the relevant non-0 transition is the size distribution of seedlings
      
      # TSF 1
    }else if(cov=="one"){
      
      outSB <- outSB.fun(pfs="burned")*sigma.s
      
      staySB <-0.05
      
      goCont=0.01 # so that we don´t devide by zero! (see below)
      
      goSB=0
      
      S <- diag(S.fun(h,cov)) # Survival Matrix 
      
      G <- t(outer(h,h,GR.fun,cov)) # Growth Matrix
      #Recruits distribution
      
      R <- (t(outer(h,h,SDS.fun,cov2)))
      
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
      
      S <- diag(S.fun(h,cov)) # Survival Matrix 
      
      G <- t(outer(h,h,GR.fun,cov)) # Growth Matrix
      
      #Recruits distribution
      R <- (t(outer(h,h,SDS.fun,cov2)))
      
      #Probability of flowering
      Fec01 = (diag(FL.fun(h,cov)))
      
      #Number of Flowering Stalks 
      Fec02 = (diag(FS.fun(h,cov)))
      
      #Number of flowers per Stalk
      
      Fec03= (diag(FPS.fun(h,cov)))
      
      #Number of seeds per flower that survive to become offspring 
      Fec04 = (diag(rep(9.1,n)))
      
      
      FecALL= Fec01*Fec02*Fec03*Fec04*goCont*sigma.s # add goCont and survSeedl to fecundity
      
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
  

  surv.TSF=array(0,c(bins,bins,length(TSF.name))) # array to hold survival kernel in
  gr.TSF=array(0,c(bins,bins,length(TSF.name))) # array to hold growth kernel in
  sds.TSF=array(0,c(bins,1,length(TSF.name))) # array to hold seedling size kernel in
  fl.TSF=array(0,c(bins,bins,length(TSF.name))) # array to hold flowering kernel in
  fs.TSF=array(0,c(bins,bins,length(TSF.name))) # array to hold # stalks kernel in
  fps.TSF=array(0,c(bins,bins,length(TSF.name))) # array to hold # flowers/stalks kernel in
  goCont.TSF=array(0,c(1,1,length(TSF.name))) # array to hold # immediate germination kernel in
  goSB.TSF=array(0,c(1,1,length(TSF.name))) # array to hold # seed-bank ingression kernel in
  staySB.TSF=array(0,c(1,1,length(TSF.name))) # array to hold # seed-bank stasis kernel in
  outSB.TSF=array(0,c(1,1,length(TSF.name))) # array to hold # seed-bank egression kernel in
  
  
  for(t in 1:5){ # loop through TSF

    
    cov=TSF.name[t]
    
    #Associate IPM in TSF[t] with seedling size measured at TSF[t+1]
    # except for TSF >3
    
    if(t==5){cov2<-TSF.name[t]}else{cov2<-TSF.name[t+1]}
    
    sigma.s=corr[corr$TSF==cov,"sigma.s"]
    
    pfs=pfs.name[t]
    
    pfs.goCont=pfs.name.goCont[t]

      
      M=IPMkernel(bins) 
      
      ### Create P and F matrices including discrete stages
      surv.TSF[,,t]=M$S
      gr.TSF[,,t]=M$G
      sds.TSF[,,t]=M$R[,1]
      fl.TSF[,,t]=M$Fec01
      fs.TSF[,,t]=M$Fec02
      fps.TSF[,,t]=M$Fec03
      goCont.TSF[,,t]=M$goCont
      goSB.TSF[,,t]=M$goSB
      staySB.TSF[,,t]=M$staySB
      outSB.TSF[,,t]=M$outSB
      
    
  }
  
  ### For each fire frequency in simulations calculate the stable distibutions of environmental trnastions,
  # i.e., the right eigenvector of the Markov chains
  # and multiply each TSF specific vital rates by the proportion the environment will be in that state 
  
  surv.mu=array(0,c(bins,bins,length(fire.freq))) 
  gr.mu=array(0,c(bins,bins,length(fire.freq))) 
  sds.mu=array(0,c(bins,1,length(fire.freq))) 
  fl.mu=array(0,c(bins,bins,length(fire.freq))) 
  fs.mu=array(0,c(bins,bins,length(fire.freq))) 
  fps.mu=array(0,c(bins,bins,length(fire.freq))) 
  goCont.mu=array(0,c(1,1,length(fire.freq))) 
  goSB.mu=array(0,c(1,1,length(fire.freq))) 
  staySB.mu=array(0,c(1,1,length(fire.freq))) 
  outSB.mu=array(0,c(1,1,length(fire.freq))) 
  
  # FUNCTION TOSIMULATE SEQUENCE OF ENVIRONMENTS
  
  simula <- function(trans,N) {
    # This function tells you sample the rownames (i.e., environment at t+1), each having a probability = column name that you got at previous iteration(environment at t) 
    state.at.N <- function(char,trans) {
      sample(rownames(trans),1,prob=trans[,char])
    }
    
    sim <- character(N)
    sim[1] <- "1" # start with 0 matrix and populate 1st character with "1"
    
    # Populate characters 2-N with the row names samples obtained from (state.at.N)
    for (b in 2:N) {
      sim[b] <- state.at.N(sim[b-1],trans)
    }
    
    sim
  }
  
  for(f in 1: length(fire.freq)){
    
    
    fire=fire.freq[f]; 
    
    P=matrix(c(fire,1-fire,rep(0,3),
               fire,0,1-fire,rep(0,2),
               fire,0,0,1-fire,0,
               fire,0,0,0,1-fire,
               fire,0,0,0,1-fire),nrow=5,ncol=5,byrow=F)
    colnames(P) <- c("1","2","3", "4","5")
    row.names(P) <- c("1","2","3", "4","5")
    
    vec_c <- eigen(P) 
    cmax  <- as.numeric(max(abs(vec_c$values)))
    loc   <- which(vec_c$values==max(abs(vec_c$values)))
    
    cstar <- as.numeric(vec_c$vec[,loc])  
    cstar <- cstar/sum(cstar)
    
    ####
    
    #for each time-since -fire environment
    
    for (i in 1:length(cstar)) {
      
      
      surv.mu[,,f]=surv.mu[,,f]+cstar[i]*surv.TSF[,,i]
      gr.mu[,,f]=gr.mu[,,f]+cstar[i]*gr.TSF[,,i]
      sds.mu[,,f]=sds.mu[,,f]+cstar[i]*sds.TSF[,,i]
      fl.mu[,,f]=fl.mu[,,f]+cstar[i]*fl.TSF[,,i]
      fs.mu[,,f]=fs.mu[,,f]+cstar[i]*fs.TSF[,,i]
      fps.mu[,,f]=fps.mu[,,f]+cstar[i]*fps.TSF[,,i]
      goCont.mu[,,f]=goCont.mu[,,f]+cstar[i]*goCont.TSF[,,i]
      goSB.mu[,,f]=goSB.mu[,,f]+cstar[i]*goSB.TSF[,,i]
      staySB.mu[,,f]=staySB.mu[,,f]+cstar[i]*staySB.TSF[,,i]
      outSB.mu[,,f]=outSB.mu[,,f]+cstar[i]*outSB.TSF[,,i]
      
    }  
    
    
  }
  
  
 
  return(list(surv=surv.mu,gr=gr.mu,sds=sds.mu,fl=fl.mu,
              fs=fs.mu,fps=fps.mu,goCont=goCont.mu,goSB=goSB.mu,
              staySB=staySB.mu,outSB=outSB.mu))
}
