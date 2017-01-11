# Scrip for Paniw et al. XXXX 

# Code to run a Bayesian models to estimate Drosophyllum vital rates

# All priors used here are uninformative 

# For details on parameterizations of the models see Appendix S3

# The code is devided into 3 parts: 

# Part A - write, initialize, and run global Bayesian model that includes all vital rates; 
# Part B - save results 
# Part C - run diagnostics 

# For the likelihood functions describing number of flowering stalks (fs) and flowers per stalk (fps),
# the negative binomial distribution is expressed as a Poisson-Gamma mixture distribution (this is commonly done in Bayesian models)

# To make the model identifiable, we used the sum-to-zero constraint 
# on all categorical variables TSF and site (see Appendix S3 for detail) 
##############################################
# PART A 

# set working directory
setwd("...")

# load data (assumed they are found in the working directory)

dataCont= read.csv("dataDroso.csv")
dataSB= read.csv("dataDrosoSB.csv")

# Order the levels of TSF categories

dataCont$TSF=factor(dataCont$TSF,levels = c("one", "two", "three",">three"))
dataSB$block=factor(dataSB$block)

# load BRUGS

library(BRugs)  
library(boot)

# BUGS model must be written enclosed by "" (in the case of BRugs; other programs require a slightly different input but the model core is always similar)

modelstring = "

model {

#SURVIVAL (surv) LIKELIHOOD 

for ( j in 1:NtotalSURV ) {
surv[j] ~ dbern( mu.surv[j] )
mu.surv[j] <- 1/(1+exp(-( a0.surv + a1.surv[TSF.surv[j]] + aS.surv[site.surv[j]] + bc.surv * size.surv[j] )))
}

##################################################################
# SURVIVAL PRIORS

a0.surv ~ dnorm( 0 , 1.0E-06 ) 
bc.surv ~ dnorm( 0 , 1.0E-06) 
#
#STZ constaints:
a1.surv[1] <- -sum(a1.surv[2:NTSF.Lvl.P])

aS.surv[1] <- -sum(aS.surv[2:Nsite.Lvl])

for ( j1 in 2:NTSF.Lvl.P ) { 
a1.surv[j1] ~ dnorm( 0.0 , 1.0E-06 ) 
}

for ( jS in 2:Nsite.Lvl ) { aS.surv[jS] ~ dnorm( 0.0 , aS.surv.tau ) }
aS.surv.tau <- 1 / pow( aS.surv.SDunabs , 2 )

aS.surv.SDunabs ~ dunif(0,20)


##################################################################
# MEAN OF GROWTH (gr) LIKELIHOOD

for ( i in 1:NtotalGR ) {

gr[i] ~ dnorm( mu.gr[i] , tau.gr )
mu.gr[i] <- a0.gr + a1.gr[TSF.gr[i]] + aS.gr[site.gr[i]] + ( bc.gr + bcTSF.gr[TSF.gr[i]] ) * size.gr[i]
# residual.gr[i] <- gr[i] - mu.gr[i]
}

##################################################################
# GROWTH PRIORS
#
tau.gr <- pow( sigma.gr , -2 )
sigma.gr ~ dunif(0,20) 

a0.gr ~ dnorm(0,1.0E-06) 
bc.gr ~ dnorm(0,1.0E-06) 

#STZ constaints:
a1.gr[1] <- -sum(a1.gr[2:NTSF.Lvl.P])
bcTSF.gr[1] <- -sum(bcTSF.gr[2:NTSF.Lvl.P])
aS.gr[1] <- -sum(aS.gr[2:Nsite.Lvl])


for ( i1 in 2:NTSF.Lvl.P ) { 
a1.gr[i1] ~ dnorm( 0.0 ,1.0E-06 )
bcTSF.gr[i1] ~ dnorm( 0.0 , 1.0E-06 )
}

for ( iS in 2:Nsite.Lvl ) { aS.gr[iS] ~ dnorm( 0.0 , aS.gr.tau ) }
aS.gr.tau <- 1 / pow( aS.gr.SDunabs , 2 )

aS.gr.SDunabs ~ dunif(0,100)

##################################################################

# MEAN OF SEEDLING SIZE (sds) LIKELIHOOD

for ( s in 1:NtotalSDS ) {
sds[s] ~ dnorm( mu.sds[s] , tau.sds )
mu.sds[s] <- a0.sds + a1.sds[TSF.sds[s]] + aS.sds[site.sds[s]]
# residual.sds[s] <- sds[s] - mu.sds[s]
}

##################################################################
# SEEDLING SIZE PRIORS
#
tau.sds <- pow( sigma.sds , -2 )
sigma.sds ~ dunif(0,20) 

a0.sds ~ dnorm(0,1.0E-06) 

#STZ constaints:
a1.sds[1] <- -sum(a1.sds[2:NTSF.Lvl.P])
aS.sds[1] <- -sum(aS.sds[2:Nsite.Lvl])


for ( s1 in 2:NTSF.Lvl.P ) { a1.sds[s1] ~ dnorm( 0.0 , 1.0E-06 ) }

for ( sS in 2:Nsite.Lvl ) { aS.sds[sS] ~ dnorm( 0.0 , aS.sds.tau ) }
aS.sds.tau <- 1 / pow( aS.sds.SDunabs , 2 )

aS.sds.SDunabs ~ dunif(0,100)

##################################################################
# PROBABILITY OF FLOWERING (fl) LIKEHOOD

for ( f in 1:NtotalFL ) {
fl[f] ~ dbern( mu.fl[f] )

mu.fl[f] <- 1/(1+exp(-lpsi.lim[f]))

lpsi.lim[f] <- min(999, max(-999,lpsi[f]))
lpsi[f] <-a0.fl + a1.fl[TSF.fl[f]] + aS.fl[site.fl[f]] + ( bc.fl + bcTSF.fl[TSF.fl[f]] ) * size.fl[f]
}

##################################################################
# PROBABILITY OF FLOWERING PRIORS

a0.fl ~ dnorm( 0 , 1.0E-06 ) 
bc.fl ~ dnorm( 0 , 1.0E-06) 
#
#STZ constaints:
a1.fl[1] <- -sum(a1.fl[2:NTSF.Lvl.F])
bcTSF.fl[1] <- -sum(bcTSF.fl[2:NTSF.Lvl.F])
aS.fl[1] <- -sum(aS.fl[2:Nsite.Lvl])

for ( f1 in 2:NTSF.Lvl.F ) { 
a1.fl[f1] ~ dnorm( 0.0 , 1.0E-06 )
bcTSF.fl[f1] ~ dnorm( 0.0 , 1.0E-06 )
}


for ( fS in 2:Nsite.Lvl ) { aS.fl[fS] ~ dnorm( 0.0 , aS.fl.tau ) }
aS.fl.tau <- 1 / pow( aS.fl.SDunabs , 2 )

aS.fl.SDunabs ~ dunif(0,20)

##################################################################
# NUMBER OF FLOWERING STALKS (fs) LIKELIHOOD

for ( t in 1:NtotalFS ) {
fs[t] ~ dpois( mustar.fs[t] )
mustar.fs[t] <- rho.fs[t]*mu.fs[t]
log(mu.fs[t]) <- a0.fs + bc.fs*size.fs[t] + a1.fs[TSF.fs[t]] + aS.fs[site.fs[t]] 
rho.fs[t] ~ dgamma( alpha.fs , alpha.fs )

}

##################################################################
# NUMBER OF FLOWERING STALKS PRIORS
a0.fs ~ dnorm( 0 , 1.0E-06 ) 
bc.fs ~ dnorm( 0 , 1.0E-06) 
alpha.fs <- exp( logalpha.fs )
logalpha.fs ~ dnorm( 0 , 0.00001 )
#
#STZ constaints:
a1.fs[1] <- -sum(a1.fs[2:NTSF.Lvl.F])
aS.fs[1] <- -sum(aS.fs[2:Nsite.Lvl])

for ( t1 in 2:NTSF.Lvl.F ) { 
a1.fs[t1] ~ dnorm( 0.0 , 1.0E-06 )
}


for ( tS in 2:Nsite.Lvl ) { aS.fs[tS] ~ dnorm( 0.0 , aS.fs.tau ) }
aS.fs.tau <- 1 / pow( aS.fs.SDunabs , 2 )

aS.fs.SDunabs ~ dunif( 0 , 100 )

##################################################################
# NUMBER OF FLOWERS PER STALK (fps) LIKELIHOOD

for ( x in 1:NtotalFPS ) {
fps[x] ~ dpois( mustar.fps[x] )
mustar.fps[x] <- rho.fps[x]*mu.fps[x]
log(mu.fps[x]) <- a0.fps + bc.fps*size.fps[x] + a1.fps[TSF.fps[x]] + aS.fps[site.fps[x]]
rho.fps[x] ~ dgamma( alpha.fps , alpha.fps )

}
##################################################################

# NUMBER OF FLOWERS PER STALK PRIORS

a0.fps ~ dnorm( 0 , 1.0E-06 ) 
bc.fps ~ dnorm( 0 , 1.0E-06) 
alpha.fps <- exp( logalpha.fps )
logalpha.fps ~ dnorm( 0 , 0.00001 )
#
#STZ constaints:
a1.fps[1] <- -sum(a1.fps[2:NTSF.Lvl.F])
aS.fps[1] <- -sum(aS.fps[2:Nsite.Lvl])

for ( x1 in 2:NTSF.Lvl.F ) { a1.fps[x1] ~ dnorm( 0.0 , 1.0E-06 ) }


for ( xS in 2:Nsite.Lvl ) { aS.fps[xS] ~ dnorm( 0.0 , aS.fps.tau ) }
aS.fps.tau <- 1 / pow( aS.fps.SDunabs , 2 )

aS.fps.SDunabs ~ dunif( 0 , 100 )
# 
##################################################################
# IMMEDIATE GERMINATION (goCont) LIKELIHOOD
for ( g in 1:NtotalGOCONT ) {
goCont[g] ~ dbern( mu.goCont[g] )
mu.goCont[g] <- 1/(1+exp(-( a0.goCont + a1.goCont[PFS.goCont[g]] + aS.goCont[block.goCont[g]])))
}

##################################################################
# IMMEDIATE GERMINATION PRIORS  

a0.goCont ~ dnorm( 0 , 1.0E-06 ) 

#
#STZ constaints:
a1.goCont[1] <- -sum(a1.goCont[2:NPFS.Lvl])
aS.goCont[1] <- -sum(aS.goCont[2:Nblock.Lvl])

for ( g1 in 2:NPFS.Lvl ) { a1.goCont[g1] ~ dnorm( 0.0 , 1.0E-06 ) }
#
for ( gS in 2:Nblock.Lvl ) { aS.goCont[gS] ~ dnorm( 0.0 , aS.goCont.tau ) }
aS.goCont.tau <- 1 / pow( aS.goCont.SDunabs , 2 )

aS.goCont.SDunabs ~ dunif(0,20)


##################################################################
# GERMINATION OUT OF THE SEED BANK (outSB) LIKELIHOOD

for ( e in 1:NtotalOUTSB ) {
outSB[e] ~ dbern( mu.outSB[e] )
mu.outSB[e] <- 1/(1+exp(-( a0.outSB + a1.outSB[PFS.outSB[e]] + aS.outSB[block.outSB[e]])))
}

##################################################################
# OUT OF THE SEED BANK PRIORS 

a0.outSB ~ dnorm( 0 , 1.0E-06 ) 
#
#STZ constaints:
a1.outSB[1] <- -sum(a1.outSB[2:NPFS.Lvl])
aS.outSB[1] <- -sum(aS.outSB[2:Nblock.Lvl])

for ( e1 in 2:NPFS.Lvl ) { a1.outSB[e1] ~ dnorm( 0.0 , 1.0E-06 ) }

for ( eS in 2:Nblock.Lvl ) { aS.outSB[eS] ~ dnorm( 0.0 , aS.outSB.tau ) }
aS.outSB.tau <- 1 / pow( aS.outSB.SDunabs , 2 )

aS.outSB.SDunabs ~ dunif(0,20)

##################################################################
# STASIS IN SEED BANK (staySB) LIKELIHOOD
for ( h in 1:NtotalSTAYSB ) {
staySB[h] ~ dbern( mu.staySB[h] )
mu.staySB[h] <- 1/(1+exp(-( a0.staySB + a1.staySB[PFS.staySB[h]] + aS.staySB[block.staySB[h]])))
}

##################################################################
# STASIS IN SEED BANK PRIORS 

a0.staySB ~ dnorm( 0 , 1.0E-06) 
#
a1.staySB[1] <- -sum(a1.staySB[2:NPFS.Lvl])
aS.staySB[1] <- -sum(aS.staySB[2:Nblock.Lvl])

for ( h1 in 2:NPFS.Lvl ) { a1.staySB[h1] ~ dnorm( 0.0 , 1.0E-06 ) }

for ( hS in 2:Nblock.Lvl ) { aS.staySB[hS] ~ dnorm( 0.0 , aS.staySB.tau ) }
aS.staySB.tau <- 1 / pow( aS.staySB.SDunabs , 2 )

aS.staySB.SDunabs ~ dunif(0,20)

##################################################################
}
# ... end BUGS model
" 
# Create a model text file and send to BUGS for syntax check:
writeLines(modelstring,con="C:/Users/.../Bay_model.txt")
modelCheck( "C:/Users/../Bay_model.txt" )


# TO RUN MODEL: 

# Specify the data in a form that is compatible with BUGS, as a list 
# Avoid NA entries here

datalist = list(
  # a.) Define number of total oservations
  NTSF.Lvl.P = length(unique(dataCont$TSF)) ,
  NTSF.Lvl.F = length(unique(dataCont$TSF))-1 , # Since plants don4t reproduce in TSF1
  NPFS.Lvl = length(unique(dataSB$FireStatus[!is.na(dataSB$FireStatus)])) ,
  Nblock.Lvl = length(unique(dataSB$block[!is.na(dataSB$block)])) ,
  Nsite.Lvl = length(unique(dataCont$site)) ,
  #SURV
  NtotalSURV =length(dataCont$surv[!is.na(dataCont$surv)]) ,
  #GR
  NtotalGR = length(dataCont$size[!is.na(dataCont$size)&!is.na(dataCont$sizeNext)]) ,
  #SDS
  NtotalSDS = length(dataCont$sds[!is.na(dataCont$sds)]) ,
  #FL
  NtotalFL = length(dataCont$fl[!is.na(dataCont$fl)]) ,
  #FS
  NtotalFS = length(dataCont$fs[!is.na(dataCont$fs)]) ,
  #FPS
  NtotalFPS = length(dataCont$fps[!is.na(dataCont$fps)]) ,
  #goCont
  NtotalGOCONT = length(dataSB$response[dataSB$vitalRate=="goCont"]) ,
  #outSB
  NtotalOUTSB = length(dataSB$response[dataSB$vitalRate=="outSB"]) ,
  #staySB
  NtotalSTAYSB = length(dataSB$response[dataSB$vitalRate=="staySB"]) ,
  # b.) Define response and independent variables 
  #SURV
  surv = dataCont$surv[!is.na(dataCont$surv)] ,
  size.surv = dataCont$size[!is.na(dataCont$surv)] ,
  TSF.surv =as.numeric(dataCont$TSF[!is.na(dataCont$surv)]) ,
  site.surv =as.numeric(dataCont$site[!is.na(dataCont$surv)]), 
  #GR
  gr = dataCont$sizeNext[!is.na(dataCont$size)&!is.na(dataCont$sizeNext)] ,
  size.gr =dataCont$size[!is.na(dataCont$size)&!is.na(dataCont$sizeNext)] ,
  TSF.gr =as.numeric(dataCont$TSF[!is.na(dataCont$size)&!is.na(dataCont$sizeNext)]) ,
  site.gr =as.numeric(dataCont$site[!is.na(dataCont$size)&!is.na(dataCont$sizeNext)]), 
  # SDS
  sds =dataCont$sds[!is.na(dataCont$sds)],
  TSF.sds =as.numeric(dataCont$TSF[!is.na(dataCont$sds)]) ,
  site.sds =as.numeric(dataCont$site[!is.na(dataCont$sds)]) ,
  #FL
  fl = as.numeric(dataCont$fl[!is.na(dataCont$fl)]) ,
  size.fl = dataCont$size[!is.na(dataCont$fl)],
  TSF.fl =as.numeric(droplevels(dataCont$TSF[!is.na(dataCont$fl)])) ,
  site.fl =as.numeric(dataCont$site) ,
  #FS
  fs = dataCont$fs[!is.na(dataCont$fs)] , 
  size.fs =dataCont$size[!is.na(dataCont$fs)] ,
  TSF.fs =as.numeric(droplevels(dataCont$TSF[!is.na(dataCont$fs)])) ,
  site.fs =as.numeric(dataCont$site[!is.na(dataCont$fs)]) , 
  #FPS
  fps = dataCont$fps[!is.na(dataCont$fps)] ,
  size.fps =dataCont$size[!is.na(dataCont$fps)] ,
  TSF.fps =as.numeric(droplevels(dataCont$TSF[!is.na(dataCont$fps)])) ,
  site.fps =as.numeric(dataCont$site[!is.na(dataCont$fps)]) ,
  #goCont
  goCont=dataSB$response[dataSB$vitalRate=="goCont"] ,
  PFS.goCont = as.numeric(dataSB$FireStatus[!is.na(dataSB$FireStatus)&dataSB$vitalRate=="goCont"]) ,
  block.goCont= as.numeric(dataSB$block[!is.na(dataSB$block)&dataSB$vitalRate=="goCont"]) ,
  #outSB
  outSB=dataSB$response[dataSB$vitalRate=="outSB"] ,
  PFS.outSB = as.numeric(dataSB$FireStatus[!is.na(dataSB$FireStatus)&dataSB$vitalRate=="outSB"]) ,
  block.outSB= as.numeric(dataSB$block[!is.na(dataSB$block)&dataSB$vitalRate=="outSB"]) ,
  #staySB
  staySB=dataSB$response[dataSB$vitalRate=="staySB"] ,
  PFS.staySB = as.numeric(dataSB$FireStatus[!is.na(dataSB$FireStatus)&dataSB$vitalRate=="staySB"]) ,
  block.staySB= as.numeric(dataSB$block[!is.na(dataSB$block)&dataSB$vitalRate=="staySB"])
  
)
# Get the data into BRugs:
modelData( bugsData( datalist ) )

# INTIALIZE THE CHAINS.

nchain = 4 # if convergence issues occur, you may choose more chains

modelCompile( numChains = nchain ) #compile the model


# automatic generation of priors (done via modelGenInits()) usually doesn't work with complex and diffuse priors (as is the case here)

# Hence, initialization based on data
# For all vital rate models, initial parameter values are either based on means of the observed data
# or on parameter values obtained by fitting mixed effect models to the data 
# The origin of the initial values is specified for the survival priors (but are similar for the rest of the priors)

# In order to check for convergence, chains should be intitated with different values, ideally outside of the range of the posteriors
# Therefore, initial parameter values of each chain were offset:

Inits=list(NULL)

# mean values obtained from data

mean.values <- 
  list(
    #SURV
    a0.surv = logit((mean(datalist$surv))) , # mean survival logit transformed!!!
    bc.surv = 7.08e-01, # size intercept (obtained from mixed effect model)
    a1.surv = c(NA,logit(aggregate( dataCont$surv , list( dataCont$TSF ) , mean, na.rm=T )[,2][-1]) - logit((mean(datalist$surv)))) , #changes from mean survival for each of the TSF levels
    aS.surv = c(NA, 6.57e-01, 1.72e-01, -6.12e-02, -1.3) , # obtained from model
    aS.surv.SDunabs = sd( c(NA, 6.57e-01, 1.72e-01, -6.12e-02, -1.3), na.rm=T) ,
    #GR
    a0.gr = mean(datalist$gr) ,
    bc.gr = 5.56e-01,
    a1.gr = c(NA,aggregate(datalist$gr , list(datalist$TSF.gr) , mean, na.rm=T )[,2][-1] - mean(datalist$gr)) ,
    bcTSF.gr = c(NA, 1.02e-01, 2.74e-01, 6.96e-02),
    aS.gr = c(NA,-7.92e-02,1.31e-01, -2.51e-01,4.93e-01) ,
    sigma.gr = sd(datalist$gr)/2 , # standard prior
    aS.gr.SDunabs = sd(c(NA,-7.92e-02,1.31e-01, -2.51e-01,4.93e-01),na.rm=T),
    #SDS
    a0.sds = mean(datalist$sds),
    a1.sds = c(NA,aggregate( datalist$sds , list( datalist$TSF.sds ) , mean )[,2][-1] - mean(datalist$sds)) ,
    aS.sds = c(NA,-2.29e-01,-3.37e-01,-3.43e-01, 8.81e-01),
    sigma.sds = sd(datalist$sds)/2 , # standard prior
    aS.sds.SDunabs = sd(c(NA,-2.29e-01,-3.37e-01,-3.43e-01, 8.81e-01),na.rm=T),
    #FL
    a0.fl = logit(mean(datalist$fl)) ,
    bc.fl = 2.04e+0,
    a1.fl = c(NA,logit(aggregate( datalist$fl , list( datalist$TSF.fl ) , mean)[,2][-1]) - logit(mean(datalist$fl))) ,
    bcTSF.fl=c(NA, -1.04,1.82e-01) ,
    aS.fl =  c(NA,-1.49,-6.7e-01,-8.52e-01,4.33e-01),
    aS.fl.SDunabs =sd(c(NA,-1.49,-6.7e-01,-8.52e-01,4.33e-01),na.rm=T),
    #FS
    logalpha.fs = rnorm(1) , #standard
    rho.fs = rep(0.001,datalist$NtotalFS) , #standard
    a0.fs = log(mean(datalist$fs)) ,
    bc.fs = -4.04e+0,
    a1.fs = c(NA,log(aggregate( datalist$fs , list( datalist$TSF.fs) , mean)[,2][-1]) - log(mean(datalist$fs))) ,
    aS.fs = c(NA,-2.5e-01,2.68e-01,-3.29e-02,-1.85e-01),
    aS.fs.SDunabs = sd(c(NA,-2.5e-01,2.68e-01,-3.29e-02,-1.85e-01), na.rm=T) ,
    #FPS
    logalpha.fps = rnorm(1) ,
    rho.fps = rep(0.001,datalist$NtotalFPS) , 
    a0.fps = log(mean(datalist$fps)) ,
    bc.fps = 2.53e-01,
    a1.fps = c(NA,log(aggregate( datalist$fps , list( datalist$TSF.fps ) , mean, na.rm=T)[,2][-1]) - log(mean(datalist$fps))) ,
    aS.fps = c(NA,-3.22e-02,-6.17e-02,-8.9e-03,1.76e-01),
    aS.fps.SDunabs = sd(c(NA,-3.22e-02,-6.17e-02,-8.9e-03,1.76e-01),na.rm=T) ,
    #goCont
    a0.goCont = -2.52699 ,
    a1.goCont = c(NA,-1.302 ),
    aS.goCont = c(NA,rnorm(6,0.001,0.001))  ,
    aS.goCont.SDunabs = sd(c(NA,rnorm(6,0.001,0.001)) ,na.rm=T) ,
    #staySB
    a0.staySB = 1.06225 ,
    a1.staySB = c(NA,0.6 ) ,
    aS.staySB = c(NA,-3.95e-01,-2.86e-02,-6.59e-02,-1.36e-01,4.96e-01,3.72e-01) ,
    aS.staySB.SDunabs = sd(c(NA,-3.95e-01,-2.86e-02,-6.59e-02,-1.36e-01,4.96e-01,3.72e-01),na.rm=T) ,
    #outSB
    a0.outSB = -3.108948 ,
    a1.outSB = c(NA,-0.5817575) ,
    aS.outSB = c(NA,1.1e-01,-3.44e-01,4.36e-01,4.99e-01,2.08e-02,-5.19e-01) ,
    aS.outSB.SDunabs = sd(c(NA,1.1e-01,-3.44e-01,4.36e-01,4.99e-01,2.08e-02,-5.19e-01),na.rm=T)
    
  )

# initals smaller than mean values above 
Inits[[1]]<- lapply(mean.values, "*", 0.6) 

# initals smaller than mean values above 
Inits[[2]] <- lapply(mean.values, "*", 0.8) 

# initals larger than mean values above 
Inits[[3]] <- lapply(mean.values, "*", 1.2) 

Inits[[4]] <- lapply(mean.values, "*", 1.5) 

# Specify initial values for each chain 
for ( i in 1 : nchain ) {
  modelInits( bugsInits( Inits[i] ) )
  
}

##############################################################################
# RUN THE CHAINS 
BurnIn = 100000 # you may need a larger burnin period if chains don4t converge
modelUpdate( BurnIn )

# which parameters to monitor
samplesSet( c( "a0.surv" , "bc.surv" ,  "a1.surv" , "aS.surv" ,
               "a0.gr" , "bc.gr" ,  "a1.gr" , "bcTSF.gr" , "aS.gr" ,
               "a0.sds" ,  "a1.sds" , "aS.sds" ,
               "a0.fl" , "bc.fl" ,  "a1.fl" , "bcTSF.fl" , "aS.fl" ,
               "a0.fs" , "bc.fs" , "a1.fs" , "aS.fs" ,
               "a0.fps" , "bc.fps" , "a1.fps" , "aS.fps" ,
               "a0.goCont" ,  "a1.goCont" , "aS.goCont" ,
               "a0.staySB" ,  "a1.staySB" , "aS.staySB" ,
               "a0.outSB" , "a1.outSB" ,"aS.outSB" ))

stepsPerChain = 250
thinStep = 400 # autocorrelation within chains may be high. Therefore, one needs to be careful not to sample posterior values to close to each other
modelUpdate( stepsPerChain , thin=thinStep ) # does a 250 * 400 MCMC run for each chain

######################################################################################33
# PART B 

# GET RESULTS - SURVIVAL 

# Extract a values:
a0.surv = samplesSample( "a0.surv" )
bc.surv=samplesSample( "bc.surv" )
chainLength = length(a0.surv)
a1.surv = array( 0 , dim=c( datalist$NTSF.Lvl.P , chainLength ) )
for ( i in 1:datalist$NTSF.Lvl.P ) {
  a1.surv[i,] = samplesSample( paste("a1.surv[",i,"]",sep="") )
}

aS.surv = array( 0 , dim=c( datalist$Nsite.Lvl , chainLength ) )
for ( i in 1:datalist$Nsite.Lvl ) {
  aS.surv[i,] = samplesSample( paste("aS.surv[",i,"]",sep="") )
}

##############################################################################################
#GROWTH
# Extract a values:
a0.gr = samplesSample( "a0.gr" )
bc.gr=samplesSample( "bc.gr" )
# residual.gr=samplesHistory( "residual.gr" ,plot=F)

a1.gr = array( 0 , dim=c( datalist$NTSF.Lvl.P , chainLength ) )
for ( i in 1:datalist$NTSF.Lvl.P ) {
  a1.gr[i,] = samplesSample( paste("a1.gr[",i,"]",sep="") )
}

bcTSF.gr = array( 0 , dim=c( datalist$NTSF.Lvl.P , chainLength ) )
for ( i in 1:datalist$NTSF.Lvl.P ) {
  bcTSF.gr[i,] = samplesSample( paste("bcTSF.gr[",i,"]",sep="") )
}

aS.gr = array( 0 , dim=c( datalist$Nsite.Lvl , chainLength ) )
for ( i in 1:datalist$Nsite.Lvl ) {
  aS.gr[i,] = samplesSample( paste("aS.gr[",i,"]",sep="") )
}
##############################################################################################

# SEEDLING SIZE DISTRIBUTION 


# Extract a values:
a0.sds = samplesSample( "a0.sds" )


a1.sds = array( 0 , dim=c( datalist$NTSF.Lvl.P , chainLength ) )
for ( i in 1:datalist$NTSF.Lvl.P ) {
  a1.sds[i,] = samplesSample( paste("a1.sds[",i,"]",sep="") )
}

aS.sds = array( 0 , dim=c( datalist$Nsite.Lvl , chainLength ) )
for ( i in 1:datalist$Nsite.Lvl ) {
  aS.sds[i,] = samplesSample( paste("aS.sds[",i,"]",sep="") )
}

###########################################################################################

# PROBABILITY OF FLOWERING  

# Extract a values:
a0.fl = samplesSample( "a0.fl" )
bc.fl=samplesSample( "bc.fl" )

a1.fl = array( 0 , dim=c( datalist$NTSF.Lvl.F , chainLength ) )
for ( i in 1:datalist$NTSF.Lvl.F ) {
  a1.fl[i,] = samplesSample( paste("a1.fl[",i,"]",sep="") )
}
bcTSF.fl = array( 0 , dim=c( datalist$NTSF.Lvl.F , chainLength ) )
for ( i in 1:datalist$NTSF.Lvl.F ) {
  bcTSF.fl[i,] = samplesSample( paste("bcTSF.fl[",i,"]",sep="") )
}

aS.fl = array( 0 , dim=c( datalist$Nsite.Lvl , chainLength ) )
for ( i in 1:datalist$Nsite.Lvl ) {
  aS.fl[i,] = samplesSample( paste("aS.fl[",i,"]",sep="") )
}

##############################################################################################
#NUMBER OF FLOWERING STALKS   

# Extract a values:
a0.fs = samplesSample( "a0.fs" )
bc.fs=samplesSample( "bc.fs" )

a1.fs = array( 0 , dim=c( datalist$NTSF.Lvl.F , chainLength ) )
for ( i in 1:datalist$NTSF.Lvl.F ) {
  a1.fs[i,] = samplesSample( paste("a1.fs[",i,"]",sep="") )
}

aS.fs = array( 0 , dim=c( datalist$Nsite.Lvl , chainLength ) )
for ( i in 1:datalist$Nsite.Lvl ) {
  aS.fs[i,] = samplesSample( paste("aS.fs[",i,"]",sep="") )
}

##############################################################################################

# SAVE RESULTS - NUMBER OF FLOWERS PER STALK   


# Extract a values:
a0.fps = samplesSample( "a0.fps" )
bc.fps=samplesSample( "bc.fps" )

a1.fps = array( 0 , dim=c( datalist$NTSF.Lvl.F , chainLength ) )
for ( i in 1:datalist$NTSF.Lvl.F ) {
  a1.fps[i,] = samplesSample( paste("a1.fps[",i,"]",sep="") )
}
aS.fps = array( 0 , dim=c( datalist$Nsite.Lvl , chainLength ) )
for ( i in 1:datalist$Nsite.Lvl ) {
  aS.fps[i,] = samplesSample( paste("aS.fps[",i,"]",sep="") )
}

##############################################################################################

# IMMEDIATE GERMINATION   

# Extract a values:
a0.goCont = samplesSample( "a0.goCont" )

a1.goCont = array( 0 , dim=c( datalist$NPFS.Lvl , chainLength ) )
for ( i in 1:datalist$NPFS.Lvl) {
  a1.goCont[i,] = samplesSample( paste("a1.goCont[",i,"]",sep="") )
}

aS.goCont = array( 0 , dim=c( datalist$Nblock.Lvl , chainLength ) )
for ( i in 1:datalist$Nblock.Lvl ) {
  aS.goCont[i,] = samplesSample( paste("aS.goCont[",i,"]",sep="") )
}
##############################################################################################

# GERMINATION FROM SEED BANK   

# Extract a values:
a0.outSB = samplesSample( "a0.outSB" )

a1.outSB = array( 0 , dim=c( datalist$NPFS.Lvl , chainLength ) )
for ( i in 1:datalist$NPFS.Lvl) {
  a1.outSB[i,] = samplesSample( paste("a1.outSB[",i,"]",sep="") )
}
aS.outSB = array( 0 , dim=c( datalist$Nblock.Lvl , chainLength ) )
for ( i in 1:datalist$Nblock.Lvl ) {
  aS.outSB[i,] = samplesSample( paste("aS.outSB[",i,"]",sep="") )
}
##############################################################################################

# SAVE RESULTS - STASIS IN THE SEED BANK  
# Extract a values:
a0.staySB = samplesSample( "a0.staySB" )

a1.staySB = array( 0 , dim=c( datalist$NPFS.Lvl , chainLength ) )
for ( i in 1:datalist$NPFS.Lvl ) {
  a1.staySB[i,] = samplesSample( paste("a1.staySB[",i,"]",sep="") )
}

aS.staySB = array( 0 , dim=c( datalist$Nblock.Lvl , chainLength ) )
for ( i in 1:datalist$Nblock.Lvl ) {
  aS.staySB[i,] = samplesSample( paste("aS.staySB[",i,"]",sep="") )
}

##############################################################################
#### SAVE RESULTS

param=cbind( a0.surv , bc.surv ,  t(a1.surv) , t(aS.surv) ,
             a0.gr , bc.gr ,  t(a1.gr) , t(bcTSF.gr ), t(aS.gr ),
             a0.sds ,  t(a1.sds ), t(aS.sds ),
             a0.fl , bc.fl ,  t(a1.fl ), t(bcTSF.fl ), t(aS.fl ),
             a0.fs , bc.fs , t(a1.fs ), t(aS.fs ),
             a0.fps , bc.fps , t(a1.fps ), t(aS.fps ),
             a0.goCont ,  t(a1.goCont ), t(aS.goCont ),
             a0.staySB ,  t(a1.staySB ), t(aS.staySB ),
             a0.outSB , t(a1.outSB), t(aS.outSB))

colnames(param)=c( "a0.surv" , "bc.surv" ,  "a1.surv.one" ,"a1.surv.two","a1.surv.three","a1.surv.four", "aS.surv.siteA","aS.surv.siteB", "aS.surv.siteC","aS.surv.siteD","aS.surv.siteE",
                   "a0.gr" , "bc.gr" ,  "a1.gr.one","a1.gr.two","a1.gr.three","a1.gr.four", "bcTSF.gr.one","bcTSF.gr.two","bcTSF.gr.three","bcTSF.gr.four", "aS.gr.siteA","aS.gr.siteB","aS.gr.siteC","aS.gr.siteD","aS.gr.siteE",
                   "a0.sds" ,  "a1.sds.one","a1.sds.two","a1.sds.three","a1.sds.four",  "aS.sds.siteA", "aS.sds.siteB", "aS.sds.siteC","aS.sds.siteD","aS.sds.siteE", 
                   "a0.fl" , "bc.fl" ,  "a1.fl.two", "a1.fl.three","a1.fl.four","bcTSF.fl.two","bcTSF.fl.three","bcTSF.fl.four", "aS.fl.siteA","aS.fl.siteB","aS.fl.siteC","aS.fl.siteD","aS.fl.siteE",
                   "a0.fs" , "bc.fs" , "a1.fs.two","a1.fs.three", "a1.fs.four", "aS.fs.siteA","aS.fs.siteB","aS.fs.siteC","aS.fs.siteD","aS.fs.siteE",
                   "a0.fps" , "bc.fps" , "a1.fps.two", "a1.fps.three", "a1.fps.four",  "aS.fps.siteA" ,"aS.fps.siteB","aS.fps.siteC","aS.fps.siteD","aS.fps.siteE",
                   "a0.goCont" ,  "a1.goCont.burned","a1.goCont.unburned", "aS.goCont.1","aS.goCont.2", "aS.goCont.3","aS.goCont.4","aS.goCont.5","aS.goCont.6","aS.goCont.7",
                   "a0.staySB" ,  "a1.staySB.burned","a1.staySB.unburned", "aS.staySB.1","aS.staySB.2","aS.staySB.3","aS.staySB.4","aS.staySB.5","aS.staySB.6","aS.staySB.7",
                   "a0.outSB" ,  "a1.outSB.burned", "a1.outSB.unburned", "aS.outSB.1","aS.outSB.2","aS.outSB.3","aS.outSB.4","aS.outSB.5","aS.outSB.6","aS.outSB.7" ) 

write.csv(as.data.frame(param),"mcmcOUT_U.csv")

### Plot some results 
library(coda)

# convert param from a dataframe to an MCMC object required by coda
param.coda=mcmc.list(list(mcmc(param[1:250,]),mcmc(param[251:500,]),mcmc(param[501:750,]),mcmc(param[751:1000,])))


par(mar=c(2,2,2,2))

#Trace plots (to check if chains are well mixed)
plot(param.coda,smooth=F) # The different colors indicate different chains

#Check for autocorrelation

autocorr.plot(param.coda)

# Calculate the GRB statistic

gelman.diag(param.coda,multivariate=F)# If the the difference in witin-chain and between-chain variance is small (as it should be), the medeian value will be around 1 
