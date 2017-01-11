# Scrip for Paniw et al. XXXX 

# Code to simulations of the stochstic population growth rate, its elasticities,
# and the probability of quasi-extinction.

# This code runs simulations for a subset of the fire return intervals discussed in the main text
# These intervals are: 10, 30, 50, and 80 years 


# This script consistst of several parts:

# PART A - load necessary fucntions to simulate Markov chains and link IPMs to each environmental state 
# PART B - prepare the data and parameters necessary to run simulations
# PART C - do the actual simulations over 2000 years including parameter uncertainty
# PART D - run simulations under environmental stochasticity only and calculate elasticities to changes in mean and variance of vital rates
# PART E - plot results

# !!!!!!!!!!!!!  NOTE  !!!!!!!!!!!!!!!!!!!!!

# In PART D, only the elasticities for goSB, staySB, and outSB are calculated in order to keep the coding to a minimum
# To calculate the relative elasticties (comparing all vital rates), the reader can follow the example given for the three
# seed-bank vital rates and calculate the elasticities for the remaining vital rates

#####################################################################################
# PART A - FUNCTIONS 

library(plyr)
library(ggplot2)
library(fields)
library(reshape2)
library(lme4)


### Add simulations 

# Simulation function (just one of many possibilities)
# This function takes a environmental transition matrix (trans) and the number of simulation (N years)
# and returns a character vector of the sequence of environments in N years 

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


#############################################################################

### PART B - Prepare parameters


num.of.sim = 10 # number of simulations to run (here fewer to speed process up); in the ms, we used 100

# Define frequency of fire (= 1/fire return interval) 

freqarray=c(1/10,1/30,1/50,1/80) # in the ms, we used 10 different frequencies 

# Define simulation time
tr=500 #the discard time
ts=2000 

# run each simulation for ts+tr years (in main text, we use 4500 years, here only 2500 to save computation time)  
trun=ts+tr

nsite = 5 # five sites used 

#IPMs
par.sb=c("goCont","staySB","outSB")

# How many parameters were sampled?
par.subset=dim(IPMs.pu)[6] # if you are working with the array otherwise type manually

# The different fire environments 
env=c("zero","one","two","three",">three")
env.n=length(env)


lambda.s.pu=array(0,c(length(freqarray),num.of.sim,par.subset,length(par.sb))) # save final stochstic lambda
  
  ######################################################################################
  
  ##  PART C - RUN SIMULATIONS (GO OVER T=RUN YEARS FOR EACH PAR.SUBSET PICK) 
  
for(sb in 1:length(par.sb)){ # loop over vital rates
  
  for(par in 1:par.subset){ # loop over parameters
    
    
    for(si in 1:num.of.sim){ # loop over simulations
      
      for (f in 1:length(freqarray)) {
        #create environmental transition matrices
        
        fire=freqarray[f]; 
        
        P=matrix(c(fire,1-fire,rep(0,3),
                   fire,0,1-fire,rep(0,2),
                   fire,0,0,1-fire,0,
                   fire,0,0,0,1-fire,
                   fire,0,0,0,1-fire),nrow=5,ncol=5,byrow=F)
        colnames(P) <- c("1","2","3", "4","5")
        row.names(P) <- c("1","2","3", "4","5")
        
        
        ## Simulations of TSF states
        simulate=simula(P,trun+1)
        simul.n=sample(nsite,length(simulate),replace=T) # pick a site for each state
        # Calculate stochastic lambda
        ########################################################
        nstage <- dim(IPMs.pu)[1]
        states <- as.numeric(simulate)
        growth <- array(0,ts)   
        
        # Initialize population vectors (n0)
        vec1=c(1000,rep(1,nstage-1)) # start with 1000 seeds 
        vec1 <- vec1/sum(vec1)
        vec1 <- t(vec1) 
        
        # ITERATION TO CALCULATE LAMBDA FOR EACH TIME STEP
        for (i  in 1:trun){
          i2 <- states[i]
          i3 <- simul.n[i]
          mat1 <-  IPMs.pu[,,i3,i2,sb,par]
          vec1 <- mat1%*%as.numeric(vec1)
          growth1 <- sum(vec1) # population growth at one time step 
          vec1 <- vec1/growth1
          if( i > tr){ # after the burn-in, save the growth rate for each time step      
            i1 <- i - tr
            growth[i1] <- growth1
          }
          
        }
        
        a1=sum(log(growth[1:ts]))
        
        lambda.s.pu[f,si,par,sb]= a1  
        
      }
      
      
    }
    
  }
}

## Calculate Contribution of parameter uncertainty to total variation 
ls.var=adply(lambda.s.pu,c(1,2,3,4))
colnames(ls.var)=c("fire","si","param","sb.rate","lambda")
ls.var$fire=factor(ls.var$fire)
ls.var$sb.rate=factor(ls.var$sb.rate)
levels(ls.var$sb.rate)=c("goSB","staySB","outSB")
ls.var$lambda=ls.var$lambda/ts

par.var=data.frame(var=rep(NA,12),sb.rate=rep(c("goSB","staySB","outSB"),4),
                   fire=rep(levels(ls.var$fire),each=3))

count=0
for(m in 1:length(levels(ls.var$fire))){
  
  for(i in 1:length(levels(ls.var$sb.rate))){
    count=count+1
    sub=ls.var[ls.var$fire==ls.var$fire[m]&ls.var$sb.rate==levels(ls.var$sb.rate)[i],]
    
    #fit a GLMM
    mod=lmer(lambda~ 1 + (1|param),data=sub)
    
    # exract random variance components of the model
    r.var=as.data.frame(VarCorr(mod))
    
    # define proportion of variance explained by parameter uncertainty:
    par.var$var[count]= r.var$vcov[1]/sum(r.var$vcov)
  }
}

par.var$var=round(par.var$var,2)

############################
# PART D - SIMULATIONS ON MEAN PARAMETERS 

# For this part, we need to load and used the function created in makeIPMmu.R

# set working directory (where you have the function and the needed mcmc samples)
setwd("YOUR_PATH")

# Load the posterior paramter values needed as argument in the function below

mcmc.out=read.csv("mcmcOUT.csv")[,-1]

# This is to load the IPMs created with mean parameters 
source("makeIPMmu.R")


IPM.site=makeIPM.mu(mcmc.out,bins=10,disc=1)

# This is to load the permuted IPMs 
source("perturbVR.R")


IPM.pert=perturbVR(mcmc.out)

# Needed for elasticity analysis 
t2 <- ts + 2 * tr # this is the time vector for the backward iteration to calculate the left eigenvectors

n2=(dim(IPM.site)[1])^2
n=dim(IPM.site)[1]

lambda.s=array(0,c(length(freqarray),num.of.sim)) # save final stochastic lambda

# Define arrays to hold elastiticties in 
sensE=array(0,c(n,n,length(freqarray),num.of.sim,2,3))# 2 - mu and sd elasticities; 3 - number of vital rates considered (here only a subset) 

######################################################################################

##  PART C - RUN SIMULATIONS (GO OVER T=RUN YEARS FOR EACH PAR.SUBSET PICK) 


for(si in 1:num.of.sim){ # loop over simulations
  
  for (f in 1:length(freqarray)) {
    #create environmental transition matrices
    
    fire=freqarray[f]; 
    
    P=matrix(c(fire,1-fire,rep(0,3),
               fire,0,1-fire,rep(0,2),
               fire,0,0,1-fire,0,
               fire,0,0,0,1-fire,
               fire,0,0,0,1-fire),nrow=5,ncol=5,byrow=F)
    colnames(P) <- c("1","2","3", "4","5")
    row.names(P) <- c("1","2","3", "4","5")
    
    
    ## Simulations of TSF states
    simulate=simula(P,t2+1)
    # Pick any of 5 sites for each state at each time step of the simulations
    simul.n=sample(nsite,length(simulate),replace=T) 
    
    # Calculate stochastic lambda
    ########################################################
    nstage <- dim(IPM.site)[1] 
    states <- as.numeric(simulate)
    growth <- array(0,ts)  
    uvecs  <- array(0,c(n,ts)) # right eigenvectors
    vvecs  <- array(0,c(n,ts)) # left eigenvectors
    
    # Initialize population vectors (n0)
    vec1=c(1000,rep(1,nstage-1)) # start with 1000 seeds 
    vec1 <- vec1/sum(vec1)
    vec1 <- t(vec1) 
    vec2 <- vec1
    
    # # ITERATION TO CALCULATE LAMBDA AND V AND W FOR EACH TIME STEP
    for (i  in 1:trun){
      i2 <- states[i]
      i3 <- simul.n[i]
      mat1 <-  IPM.site[,,i3,i2]
      vec1 <- mat1%*%as.numeric(vec1)
      growth1 <- sum(vec1) # population growth at one time step 
      vec1 <- vec1/growth1
      
      # forward iteration
      if( i > tr){      
        i1 <- i - tr
        uvecs[,i1] <- vec1
        growth[i1] <- growth1 # stachastic growth rate
      }
      
      # backward iteration
      j <- (t2 - i+1)
      i2 <- states[j+1]
      i3 <- simul.n[j+1]
      mat1 <- IPM.site[,,i3,i2]
      vec2 <- t(as.numeric(vec2)%*%mat1)
      vec2 <- vec2/(sum(vec2))
      if (i > tr)  {     
        vvecs[, j-tr] <- vec2
      }
    }
    a1=sum(log(growth[1:ts]))
    
    lambda.s[f,si]= a1  
    
    ### CALCULATE ELASTICITIES TO ELEMENTS (sensE)
    
    for (i in (tr+1):(tr+ts-1)){ #start after a time lag (tr)
      
      i2 <- states[i+1]
      i3 <- simul.n[i+1]
      
      mat1 <- IPM.site[,,i3,i2]       
      itime <- i+1-tr
      i1 <- i-tr;      
      
      # goSB perturbed 
      mat2A.goSB <- IPM.pert$IPMs.goSB.mu[,,i3,i2,f]
      mat2V.goSB <- IPM.pert$IPMs.goSB.mu[,,i3,i2,f]
      
      # staySB perturbed
      mat2A.staySB <- IPM.pert$IPMs.staySB.mu[,,i3,i2,f]
      mat2V.staySB <- IPM.pert$IPMs.staySB.mu[,,i3,i2,f]
      
      # outSB perturbed
      mat2A.outSB <- IPM.pert$IPMs.outSB.mu[,,i3,i2,f]
      mat2V.outSB <- IPM.pert$IPMs.outSB.mu[,,i3,i2,f]
      
      # FORMULA FOR STOCHASTIC ELASTICITY (Translation Tuljapurkar et al. 2003)
      
      # Scalar eigenvector product (denominator) 
      scale1 <-(t(vvecs[,itime]))%*%(uvecs[,itime])
      a=as.numeric((1/(scale1*growth[itime])))
      
      # Elasticities
      mat2A.goSB <- a*(diag(vvecs[,itime])%*%mat2A.goSB%*%diag(uvecs[,i1]))
      mat2V.goSB <- a*(diag(vvecs[,itime])%*%mat2V.goSB%*%diag(uvecs[,i1])) 
      
      sensE[,,f,si,1,1] <- sensE[,,f,si,1,1] + mat2A.goSB
      sensE[,,f,si,2,1] <- sensE[,,f,si,2,1] + mat2V.goSB
      
      mat2A.staySB <- a*(diag(vvecs[,itime])%*%mat2A.staySB%*%diag(uvecs[,i1]))
      mat2V.staySB <- a*(diag(vvecs[,itime])%*%mat2V.staySB%*%diag(uvecs[,i1]))
      
      sensE[,,f,si,1,2] <- sensE[,,f,si,1,2] + mat2A.staySB
      sensE[,,f,si,2,2] <- sensE[,,f,si,2,2] + mat2V.staySB
      
      mat2A.outSB <- a*(diag(vvecs[,itime])%*%mat2A.outSB%*%diag(uvecs[,i1]))
      mat2V.outSB <- a*(diag(vvecs[,itime])%*%mat2V.outSB%*%diag(uvecs[,i1])) 
     
      sensE[,,f,si,1,3] <- sensE[,,f,si,1,3] + mat2A.outSB
      sensE[,,f,si,2,3] <- sensE[,,f,si,2,3] + mat2V.outSB
    }
    
    sensE[,,f,si,,] <- sensE[,,f,si,,]/(ts-1)
  }
  
  
}

###################
# PROBABILITY OF QUASI-EXTINCTION


# Calculate probability of QUASI-EXTINCTION at t=ext.t

ext.t=c(50,100)
q=0.01 # extinction threshold

# Normal integral
int.norm=function(x){(2*pi)^(-.5)*exp(-x^2/2)}


m.lambda.s.pu=apply(lambda.s.pu,c(1,3,4),sum) 
m.lambda.s.pu=m.lambda.s.pu/c(num.of.sim*ts)

#Calculate the varaibility of lambda.s for each parameter sample 
m.lambda.v.pu=array(0,c(length(freqarray),par.subset,length(par.sb))) 


for(ps in 1:length(par.sb)){
  
  for(f in 1:length(freqarray)){
    
    for(p in 1:par.subset){
      
    for(i in 1:num.of.sim){
            
      m.lambda.v.pu[f,p,ps] = m.lambda.v.pu[f,p,ps] + (lambda.s.pu[f,i,p,ps] - m.lambda.s.pu[f,p,ps]*ts)^2 
    }
    m.lambda.v.pu[f,p,ps]=m.lambda.v.pu[f,p,ps]/(num.of.sim*ts)
    }
    
    
    
  }
  
}


# Calculate probability of quasi-extinction at t=ext.t

Pqt.pu=array(0,c(length(freqarray),par.subset,length(par.sb),length(ext.t))) 

for(t in 1:length(ext.t)){
  for(p in 1:par.subset){
    for(ps in 1:length(par.sb)){
      for(f in 1:length(freqarray)){
        
        z=(log(q)-ext.t[t]*m.lambda.s.pu[f,p,ps])/sqrt(m.lambda.v.pu[f,p,ps]*ext.t[t])
        z2=(log(q)+ext.t[t]*m.lambda.s.pu[f,p,ps])/sqrt(m.lambda.v.pu[f,p,ps]*ext.t[t])
        
        Pqt.pu[f,p,ps,t]=integrate(int.norm,-Inf,z)$value+exp(2*m.lambda.s.pu[f,p,ps]*log(q)/m.lambda.v.pu[f,p,ps])*integrate(int.norm,-Inf,z2)$value
        
        
      }
      
      
    }
  }

}
# PART E - PLOTS 

# STOCHASTIC LAMBDA

# simulations based on mean parameters
ls=adply(lambda.s,c(1,2))
colnames(ls)=c("fire","si","lambda")
ls$fire=factor(ls$fire)
ls$lambda=ls$lambda/ts
ls$sbR="mean"

# simulations based on parameter samples (here only for outSB)
ls.pu=adply(lambda.s.pu[,,,3],c(1,2,3))
colnames(ls.pu)=c("fire","si","sample","lambda")
ls.pu$fire=factor(ls.pu$fire)
ls.pu$lambda=ls.pu$lambda/ts
ls.pu$sbR="pu"


data=rbind(ls[,c(1,3,4)],ls.pu[,c(1,4,5)])
data$sbR=factor(data$sbR,levels=c("pu","mean"))

pd <- position_dodge(width = 0.8)
ggplot(data, aes(x=fire, y=lambda,colour=sbR))  +
  ylim(-0.1,.2)+
  geom_boxplot(fill="white",position=pd)+
  geom_hline(aes(yintercept=0),color="black", linetype="dashed", size=1)+
  
  guides(fill=FALSE)+ 
  theme_bw()+
  scale_color_manual(values=c("black", "grey"))+
  theme(legend.position="none")+
  xlab("Fire return") +
  ylab(expression(paste("log ",lambda[s],sep="")))+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=32))+
  theme(axis.title = element_text(size=38))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9))

# PROBABILITY OF QUASI EXTINCTION (Pqt)


Pqt.U=apply(Pqt.pu,c(1,3,4),quantile,0.975,na.rm=T)
Pqt.L=apply(Pqt.pu,c(1,3,4),quantile,0.025,na.rm=T)
Pqt.m=apply(Pqt.pu,c(1,3,4),mean,na.rm=T)

Pqt=adply(Pqt.m,c(1,2,3))
colnames(Pqt)=c("fire","sbR","time","Pqt")
Pqt$U=adply(Pqt.U,c(1,2,3))[,4]
Pqt$L=adply(Pqt.L,c(1,2,3))[,4]

levels(Pqt$sbR)=c("goSB","staySB","outSB")
levels(Pqt$time)=c("50 years","100 years")

Pqt$fire=as.numeric(Pqt$fire)

pd <- position_dodge(width = 0.6)
ggplot(Pqt, aes(x=fire, y=Pqt))  +
  ylim(0,0.9)+
  facet_wrap(~time)+
  geom_line(color="black",size=1) +
  geom_point(aes(color=sbR),size=4,position=pd)+
  geom_errorbar(data=Pqt,aes(ymin=L, ymax=U,color=sbR),width=.6,size=1.1,position=pd)+
  guides(fill=FALSE)+ 
  theme_bw()+
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9"),name="",
                     labels=c("goSB","staySB","outSB"))+
  theme(legend.justification=c(0,0.9), legend.position=c(0,1),
        legend.title = element_text(size=16, face="bold"),
        legend.text = element_text(size=25),
        legend.key.size = unit(2, "lines"),
        legend.background  = element_rect(fill = "transparent", colour = "transparent"))+
  xlab("Fire return (years)") +
  ylab(expression(paste(P[q],"(t)",sep="")))+
  theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size=32))+
  theme(axis.title = element_text(size=38))+
  theme(axis.title.x = element_text(vjust=-0.4),
        axis.title.y = element_text(vjust=0.9))+
  scale_x_continuous(breaks=seq(1,4),labels=c("10","30","50","80"))
                     

# RELATIVE ELASTICITIES 

sensE.mu=apply(sensE,c(1,2,3,5,6),mean)

sensE.sum=apply(sensE.mu,c(3,5),sum)

sensE.sumALL=apply(sensE.mu,c(3),sum)

# Relative elasticities (NOTE: these are just comparing three vital rates)

sensE.rel=sensE.sum
for(i in 1:nrow(sensE.sum)){
  
  sensE.rel[i,]=sensE.rel[i,]/sensE.sumALL[i]
}

colnames(sensE.rel)=c("goSB","staySB","outSB")

sensE.rel
