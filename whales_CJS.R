# Workflow for CJS and jags using whale data as an example
# Can be applied to any species but may need some tuning
# in some spots such as length of burn-in for MCMC chain
#

# Sept 10

library(mcmcplots)
library(R2jags)


# read data file 
whales <- read.csv( "whales.csv")

# make encounter history matrix
EH <- as.matrix( whales[,2:11] )

nYears <- ncol( EH )
nAnimal <- nrow( EH )

# determine the occasion of the first capture of each individual
get.first <- function(x){
  # a function to identify the time period of the first capture from an encounter history matrix
  return( min(which( x != 0 ) ) )
}
fall <- apply( EH, 1, get.first )

# remove animals caught only in last year... no whales were caught in the last year so skip this step or it deletes the whole matrix
#EH <- EH[-which(fall==nYears),] # this just takes out animals that were only caught in period 7 (nYears is just a number set to sampling periods above)
f <- fall#[-which(fall==nYears)]

nAnimal <- nrow( EH )

# Define a list of data to be passed to JAGS
cjs.data <- list( y=EH, f=f, nind=nAnimal, nocc=nYears, female=whales$female )

# Define a function to generate the initial values and return them as a list
ch.init <- function( ch, f ){
  for( i in 1:dim(ch)[1] ){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

# build z matrix which is our initial starting value for when every individual is alive
# we initialize everyone to be alive from the time they're first caught to the end of the study
z <- ch.init( EH, f )    # z is partially observed, so has to agree with y
z <- ifelse( !is.na(z), 1, z )
cjs.inits <- function(){
  list(
    z = z,
    b0.phi = runif(1, -3, 3),
    b1.phi = runif(1, -3, 3),
    b0.p = runif(1, -3, 3)
  )
}

# set parameters to track in JAGS
cjs.parms <- c("b0.phi","b1.phi", "b0.p","mean.phi.male", "mean.phi.female","mean.p")

# set up for MCMC run
ni <- 10000 # iterations
nt <- 1    # thinning (use 1 if you want no thinning because it divides by this number)
nb <- 5000 # burn in... these iterations get thrown out
nc <- 3    # number of chains... you want to do a few, usually 3, so you can see if they converge on the same spot

# run the MCMC chain in JAGS
cjs.result <- jags( cjs.data, 
                    cjs.inits,
                    cjs.parms,
                    "cjs_phiFemale_pdot.txt",
                    n.chains=nc, 
                    n.iter=ni, 
                    n.burnin=nb,
                    n.thin=nt
)
cjs.result  # p is recapture probability, phi is survival, in this case, annual survival 
# because each capture period is 1 year apart 

mcmcplot(cjs.result) # note: when you look at this, b1.phi is the DIFFERENCE between male and female survival
# and this value is close to 0 indicating little difference 


# Note: cool thing Paul did in class Sept 10 lab, 13:26:00

# "One thing we can do with this bayesian method that we can't with likelihood 
# is look at the probability male and female survival is different
hist(x)
x <- cjs.result$BUGSoutput$sims.matrix[,"mean.phi.female"]
length(x) # 15,000
y <- x[x<0.5] # grab the values of x where x is less than 0.5 (because thats where female survival overlapped with male)
length(y) # so there was 1,092 out of 15,000 where female and male survival overlapped
length(y)/length(x) # so the probability that female survival is lower than males is 0.07 ie not great





