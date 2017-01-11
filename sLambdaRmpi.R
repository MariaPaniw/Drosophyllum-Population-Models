# Scrip for Paniw et al. XXXX 

# Code to simulations of the stochstic population growth rate using parallel processing (Rmpi)

# This is particularly useful if you include parameter uncertainty by sampling several hundred parameters.


# Initialize MPI
options(echo=FALSE)

# load Rmpi
library(Rmpi)


###################################### Rmpi PART 

# Produce workers (or slaves in Rmpi language)
Nprocs <- mpi.universe.size() # gives the number of processors available (here 50) 

# One processor is reserved for the master, the remaining ones (here 49) are the workers (or slaves),
# which you create with the following code:
mpi.spawn.Rslaves(nslaves=Nprocs-1,needlog=FALSE) 

#Set up a random number generator for each slaves so that independece of the simulations is ensured 
mpi.setup.rngstream(iseed=NULL, comm = 1)

# Quit Rmpi if you do not have suffiencient processors assigned 
if (mpi.comm.size() != 50) {
  print("Please initialize an MPI cluster of at least 50 processors.")
  print("Then, try again")
  mpi.quit()
}

.Last <- function(){
  if (is.loaded("mpi_initialize")){
    if (mpi.comm.size(1) > 0){
      print("Please use mpi.close.Rslaves() to close slaves.")
      mpi.close.Rslaves()
    }
    print("Please use mpi.quit() to quit R")
    .Call("mpi_finalize")
  }
}


#####################################################################################

### wrap your simulations in a function!!!!


lambdaSimul <- function(par.subset){
  ### Add simulations 
  
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
  
  par.subset= parameters[sub==foldNumber]
  
  
  # The different fire environments 
  env=c("zero","one","two","three",">three")
  env.n=length(env)
  
  
  lambda.s.pu=array(0,c(length(freqarray),num.of.sim,length(par.subset),length(par.sb))) # save final stochstic lambda
  
  ######################################################################################
  
  ##  PART C - RUN SIMULATIONS (GO OVER T=RUN YEARS FOR EACH PAR.SUBSET PICK) 
  
  for(sb in 1:length(par.sb)){ # loop over vital rates
    
    for(par in 1:length(par.subset)){ # loop over parameters
      
      
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
            mat1 <-  IPMs.pu[,,i3,i2,sb,par.subset[par]]
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
  
  return(list(lambda.s.pu=lambda.s.pu))
}

####################################################################
####################################################################
# Back to Rmpi

# cut the posterior parameters into 39 parts 
sub=as.numeric(cut(1:1000,49))
parameters=1:1000

# Send both the cut and uncut parameters as input to the workers
mpi.bcast.Robj2slave(sub)
mpi.bcast.Robj2slave(parameters)
mpi.bcast.cmd(foldNumber <- mpi.comm.rank())

# Send the function that does the stochastic simulations to the workers
mpi.bcast.Robj2slave(lambdaSimul)

# Call the function in all the workers, and retrieve the results
lambdaS <- mpi.remote.exec(lambdaSimul())

# Save output
save(lambdaS,file="~/FOLDER/lambdaS")


# close slaves and exit
mpi.close.Rslaves()
mpi.quit(save = "no")


