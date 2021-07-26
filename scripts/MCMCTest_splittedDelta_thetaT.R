rm(list=ls())
step.time.input <- function(current_time_input,current_delta,from, current_theta,current_time_std, f, q, HR_Multiplier, number_accepted) {
  #proposal_time_input <- q(current_time_input)
  proposal_time_std <- q(current_time_std)
  #HR_q= HR_Multiplier(proposal.time.input, current_time_input)
  HR_q= HR_Multiplier(proposal_time_std, current_time_std)
  #alpha <- min(1, f(proposal) * HR_q / f(current_value)  )
  #print(current_delta)
  #proposal_time_std<- current_time_input / current_theta
  proposal_time_input<- proposal_time_std * current_theta
  f_proposal = f(proposal_time_std, current_delta, 0 )
  f_current = f(current_time_std, current_delta, 0 )
  if (f_proposal >0 && f_current >0 && proposal_time_input >=from)
  {
    log_alpha <- min(0, log(f_proposal) + log(HR_q) -  log(f_current) )
    #log_alpha <- min(0, logf(proposal) + log(HR_q) -  logf(current_value))
    if (log(runif(1)) < log_alpha)
    {
      current_time_input <- proposal_time_input
      current_time_std = proposal_time_std
      number_accepted <- number_accepted +1
    }
  }
  return(list(current_time_input,current_time_std, number_accepted))
}
step <- function(current_value, f, q, HR_Multiplier, number_accepted) {
  proposal <- q(current_value)
  HR_q= HR_Multiplier(proposal, current_value)
  #alpha <- min(1, f(proposal) * HR_q / f(current_value)  )
  log_alpha <- min(0, log(f(proposal)) + log(HR_q) -  log(f(current_value)) )
  #log_alpha <- min(0, logf(proposal) + log(HR_q) -  logf(current_value))
  if (log(runif(1)) < log_alpha)
  {
    current_value <- proposal
    number_accepted <- number_accepted +1
  }
  return(list(current_value, number_accepted))
}
step.DeltaT.theta <- function(current_deltaT,current_theta, f_delta, f_theta, q_multiplier, HR_Multiplier, number_accepted) {
  proposal.delta <- q_multiplier(current_deltaT)
  proposal.theta <- q_multiplier(current_theta)
  HR_q_delta= HR_Multiplier(proposal.delta, current_deltaT)
  HR_q_theta= HR_Multiplier(proposal.theta, current_theta)
  #alpha <- min(1, f(proposal) * HR_q / f(current_value)  )
  log_alpha <- min(0, log(f_delta(proposal.delta)) +log(f_theta(proposal.theta))+log(HR_q_delta)+log(HR_q_theta) -  log(f_delta(current_deltaT)) -log(f_theta(current_theta)))
  #log_alpha <- min(0, logf(proposal) + log(HR_q) -  logf(current_value))
  if (log(runif(1)) < log_alpha)
  {
    current_deltaT <- proposal.delta
    current_theta <- proposal.theta
    number_accepted <- number_accepted +1
  }
  return(list(current_deltaT, current_theta, number_accepted))
}
step.DeltaT <- function(current_deltaT,current_theta, f_delta, f_theta, q_multiplier, HR_Multiplier, number_accepted) {
  proposal.delta <- q_multiplier(current_deltaT)
  HR_q_delta= HR_Multiplier(proposal.delta, current_deltaT)
  #alpha <- min(1, f(proposal) * HR_q / f(current_value)  )
  log_alpha <- min(0, log(f_delta(proposal.delta)) +log(HR_q_delta) -  log(f_delta(current_deltaT)) )
  #log_alpha <- min(0, logf(proposal) + log(HR_q) -  logf(current_value))
  if (log(runif(1)) < log_alpha)
  {
    current_deltaT <- proposal.delta
    number_accepted <- number_accepted +1
  }
  return(list(current_deltaT, number_accepted))
}
step.theta <- function(current_deltaT,current_theta, f_delta, f_theta, q_multiplier, HR_Multiplier, number_accepted) 
{
  proposal.theta <- q_multiplier(current_theta)
  HR_q_theta= HR_Multiplier(proposal.theta, current_theta)
  #alpha <- min(1, f(proposal) * HR_q / f(current_value)  )
  log_alpha <- min(0,  log(f_theta(proposal.theta))+log(HR_q_theta) -log(f_theta(current_theta)))
  #log_alpha <- min(0, logf(proposal) + log(HR_q) -  logf(current_value))
  if (log(runif(1)) < log_alpha)
  {
    current_theta <- proposal.theta
    number_accepted <- number_accepted +1
  }
  return(list( current_theta, number_accepted))
}
HR_multiplier_bounded=function(proposal, current_value) {
  #  result= proposal * (proposal -from)/ (current_value *(current_value -from))
  result= (to-proposal) * (proposal -from)/ ((to-current_value) *(current_value -from))
  return(result)
}

HR_multiplier_unbounded=function(proposal, current_value) {
  result= proposal / current_value
  return(result)
}

run <- function(time_input,delta,theta, from,   f, q, f_delta, f_theta, q_multiplier,  nsteps, HR_multiplier_bounded,HR_multiplier_unbounded, number_accepted_pair,number_accepted_time ) 
{
  res_time_input <- matrix(NA, nsteps, length(time_input))
  res_delta <- matrix(NA, nsteps, length(delta))
  res_theta <- matrix(NA, nsteps, length(theta))
  
  l=list()
  current_delta<-delta
  current_theta<-theta
  current_time_input<-time_input
  current_time_std<-current_time_input / current_theta
  res_time_std <- matrix(NA, nsteps, length(current_time_std))
  number_accepted_delta=0
  number_accepted_theta=0
  number_accepted_time=0
  for (i in seq_len(nsteps))
  {
    print(paste("iteration ",i, sep=""))
    #current_deltaT,current_theta, f_delta, q_delta,f_theta, q_theta, HR_Multiplier, number_accepted
    l <-step.DeltaT(current_delta,current_theta ,f_delta, f_theta, q_multiplier, HR_multiplier_unbounded, number_accepted_delta )
    current_delta<- l[[1]] 
    res_delta[i,] <- current_delta
    number_accepted_delta <- l[[2]]
    l <-step.theta(current_delta,current_theta ,f_delta, f_theta, q_multiplier, HR_multiplier_unbounded, number_accepted_theta )
    current_theta<- l[[1]]
    res_theta[i,] <- current_theta
    number_accepted_theta<- l[[2]]
    current_time_std<-current_time_input / current_theta
    #l <- step.time.input(current_time_input,current_delta,from, current_theta, current_time_std, f, q, HR_multiplier_bounded, number_accepted_time )
    l <- step.time.input(current_time_input,current_delta,from, current_theta, current_time_std, f, q, HR_multiplier_unbounded, number_accepted_time )
    #current_time_input,current_delta,from, current_theta,current_time_std, f, q, HR_Multiplier, number_accepted
    current_time_input <- l[[1]] 
    res_time_input[i,] <- current_time_input
    current_time_std<-l[[2]]
    res_time_std[i,] <- current_time_std
    number_accepted_time <- l[[3]] 
  }
  print(number_accepted_theta)
  print(number_accepted_delta)
  print(number_accepted_time)
  drop(res_delta)
  drop(res_theta)
  drop(res_time_input)
  drop(res_time_std)
  return(list(res_delta,res_theta,res_time_input,res_time_std,number_accepted_theta,number_accepted_delta,number_accepted_time ))
}
f <- function(x, delta, from)
{
  if (x >= from){
    term=delta*exp(-1.0*delta*(x))/(1-exp(-1.0*delta*(x)))
    result = delta*term * exp(-1*term) / (1-exp(-1.0*delta*(x)))
    return(result)
  }
  else {
    return(0)
  }
}
f_delta <- function(x)
{
  result= lambdaExponentialPriorDelta * exp(-1.0 * lambdaExponentialPriorDelta* x)
  return(result)
  
}
f_theta <- function(x)
{
  result= lambdaExponentialPriorTheta * exp(-1.0 * lambdaExponentialPriorTheta* x)
  return(result)
  
}
logf <- function(x, delta)
{
  term1 = exp(-1.0*delta*x);
  term2 = delta * term1;
  term3 = 1.0-term1;
  result = log(delta * term2 /(term3 * term3));
  result = result - term2/term3;
  return(result)
}

q_transformed_multiplier <- function(current.value)
{
  u=runif(1,0,1)
  m=exp(2*log(b)*(u-0.5))
  #currentTransformed = current.value /(current.value - from)
  currentTransformed = (to - current.value) /(current.value - from)
  newTransformed = m * currentTransformed
  
  #proposed.value= from * newTransformed /( newTransformed -1) 
  proposed.value= ( to+ (newTransformed * from))  /( 1+newTransformed) 
  proposed.value
}
q_multiplier <- function(current.value)
{
  u=runif(1,0,1)
  m=exp(2*log(b)*(u-0.5))
  proposed.value = m * current.value
  proposed.value
}

draw_random<-function(delta, fromTimeSTD)
{
  u= runif(1) 
  result=fromTimeSTD+ (1.0 /delta) *log(1-(delta)/log(u))
  return(result)
  
}
unconditional.distribution=function(x){
  
  result= lambdaExponentialPriorDelta / (lambdaExponentialPriorDelta +x)^2
  return(result)
  
}
rexpexp=function(n, lambda, from ){
  random.unif= runif(n, 0, 1)
  result=lapply(random.unif, FUN=function(x, from, lambda)
  {
    res=from
    return(res + (lambda*x)/(1-x))
    
  },
  from,
  lambda)
  return(unlist(result))
  
}
expected.time.origin=function(n, lambda, from ){
  random.exp= rexp(n, rate=1)
  result=lapply(random.exp, FUN=function(x, from, lambda)
  {
    res=from
    return(res + (1.0/lambda)*log(1+ (lambda/x)))
  },
  from,
  lambda)
  return(mean(unlist(result)))
  
}
get.thinning=function(chain_data,number.iterations, lags, threshold_corr){
  if (length(lags)>0 && threshold_corr>0 && threshold_corr <1)
  {
    mcmc_data<-mcmc(data= chain_data, start = 1, end = number.iterations, thin = 1)
    mcmc.autocorr=autocorr(mcmcTime, lags = lags, relative=TRUE)
    print(mcmc.autocorr)
    thinning=lags[which(mcmc.autocorr<threshold_corr)[1]]
    return(thinning)
  }
  
}
thin.chain.data=function(chain_data,number.iterations, percent_burn_in, thinning){
  if (length(chain_data)>0 && thinning >1 && percent_burn_in >=0 && percent_burn_in< 1 )
  {
    chain_data_after_burn_in=removeBurnIn(chain_data,number.iterations, percent_burn_in)
    indexes<-seq(0,length(chain_data_after_burn_in),thinning)
    thinned_data=chain_data_after_burn_in[indexes]
    return(thinned_data)
  }
  
}
removeBurnIn<-function(chain_data,number.iterations, percent_burn_in){
  mcmc_data<-mcmc(data= chain_data, start = 1, end = number.iterations, thin = 1)
  number_interations_burn_in= floor(percent_burn_in * length(chain_data))
  chain_data_after_burn_in=chain_data[(number_interations_burn_in+1):length(chain_data)]
  return(chain_data_after_burn_in)
}
runChain=function(from, lambdaExponentialPriorDelta,lambdaExponentialPriorTheta, number.iterations,  
                  f, q_transformed_multiplier, f_delta, f_theta, q_multiplier, 
                  HR_multiplier_bounded,HR_multiplier_unbounded){
  
  number_accepted_pair=0
  number_accepted_delta=0
  delta=rexp(1,rate=lambdaExponentialPriorDelta)
  theta =rexp(1,rate=lambdaExponentialPriorTheta)
  fromTimeSTD= from *delta
  time_input = draw_random(delta, fromTimeSTD)
  list_res <- run(time_input, delta,theta, from,  f, q_transformed_multiplier, f_delta, f_theta, q_multiplier, number.iterations, 
                  HR_multiplier_bounded,HR_multiplier_unbounded,  number_accepted_pair,number_accepted_time )
  return(list_res)
}
##############################

from= 0.000408897000000000025634 
to= 5
#from=0.035
lambdaExponentialPriorDelta= 0.01
lambdaExponentialPriorTheta=100

delta=rexp(1,rate=lambdaExponentialPriorDelta)
theta =rexp(1,rate=lambdaExponentialPriorTheta)
#delta=200
curve(f(x, 1.0/lambdaExponentialPriorDelta, from) , col="red", n=200, xlim=c(from,0.4), main ="Conditional density  Delta=100, lower bound 0.0004",xlab= "time of origin", ylab="density")

number.iterations=500000
list.values <-  vector("list", number.iterations)
list.values[[1]]<-time
length.interval.multiplier=0.008
length.interval.multiplier=0.7
b=  (length.interval.multiplier)/2.0 + sqrt(1 + (length.interval.multiplier * length.interval.multiplier)* 0.25);



list_res1 <- runChain( from, lambdaExponentialPriorDelta,lambdaExponentialPriorTheta, number.iterations,
                       f, q_transformed_multiplier, f_delta, f_theta, q_multiplier,  
                       HR_multiplier_bounded,HR_multiplier_unbounded)

list_res2 <- runChain( from, lambdaExponentialPriorDelta,lambdaExponentialPriorTheta, number.iterations,
                       f, q_transformed_multiplier, f_delta, f_theta, q_multiplier,  
                       HR_multiplier_bounded,HR_multiplier_unbounded )

list_res3 <- runChain( from, lambdaExponentialPriorDelta,lambdaExponentialPriorTheta, number.iterations,
                       f, q_transformed_multiplier, f_delta, f_theta, q_multiplier,  
                       HR_multiplier_bounded,HR_multiplier_unbounded )

list_res4 <- runChain( from, lambdaExponentialPriorDelta,lambdaExponentialPriorTheta, number.iterations,
                       f, q_transformed_multiplier, f_delta, f_theta, q_multiplier,  
                       HR_multiplier_bounded,HR_multiplier_unbounded )

res_delta_chain1 = list_res1[[1]]
res_theta_chain1 = list_res1[[2]]
res_time_input_chain1 = list_res1[[3]]
res_time_std_chain1 = list_res1[[4]]
number_accepted_theta1  =  list_res1[[5]]
number_accepted_delta1  =  list_res1[[6]]
number_accepted_time1  =  list_res1[[7]]

res_delta_chain2 = list_res2[[1]]
res_theta_chain2 = list_res2[[2]]
res_time_input_chain2 = list_res2[[3]]
res_time_std_chain2 = list_res2[[4]]
number_accepted_theta2  =  list_res2[[5]]
number_accepted_delta2  =  list_res2[[6]]
number_accepted_time2  =  list_res2[[7]]

res_delta_chain3 = list_res3[[1]]
res_theta_chain3 = list_res3[[2]]
res_time_input_chain3 = list_res3[[3]]
res_time_std_chain3 = list_res3[[4]]
number_accepted_theta3  =  list_res3[[5]]
number_accepted_delta3  =  list_res3[[6]]
number_accepted_time3  =  list_res3[[7]]

res_delta_chain4 = list_res4[[1]]
res_theta_chain4 = list_res4[[2]]
res_time_input_chain4 = list_res4[[3]]
res_time_std_chain4 = list_res4[[4]]
number_accepted_theta4  =  list_res4[[5]]
number_accepted_delta4  =  list_res4[[6]]
number_accepted_time4  =  list_res4[[7]]

layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
plot(res_delta_chain1, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)
usr <- par("usr")
xx <- seq(usr[3], usr[4], length=301)
plot(f_delta(xx), xx, type="l", yaxs="i", axes=FALSE, xlab="")

hist(res_delta_chain1, 50, freq=FALSE, main="", las=1,
     xlab="x", xlim=c(0, 1000), ylab="Probability density")
z <- integrate(f_delta, 0, Inf)$value
curve(f_delta(x) / z, add=TRUE, col="red", n=200)

hist(res_theta_chain1, 50, freq=FALSE, main="", las=1,
     xlab="x", xlim=c(0, 0.1), ylab="Probability density")
z <- integrate(f_theta, 0, Inf)$value
curve(f_theta(x) / z, add=TRUE, col="red", n=200)

hist(res_time_input_chain1, 50, freq=FALSE, main="", las=1,
     xlab="x", xlim=c(from, 0.1), ylab="Probability density")

hist(res_time_std_chain2, 50, freq=FALSE, main="", las=1,
     xlab="x", ylab="Probability density")

library(coda)
mcmcTime<-coda::mcmc(data= res_time_std_chain2, start = 1, end = number.iterations, thin = 1)
mcmc_list=as.mcmc.list(mcmcTime)
#gelman.diag(mcmc_list, confidence = 0.95,  autoburnin=TRUE,
#         multivariate=TRUE)
lags = c(100, 200, 500, 1000)
threshold_corr=0.1
thinning.delta_chain1<-get.thinning(res_delta_chain1,number.iterations, lags, threshold_corr)
thinning.time.std_chain1<-get.thinning(res_time_std_chain1,number.iterations, lags, threshold_corr)
thinning.time.std_chain2<-get.thinning(res_time_std_chain2,number.iterations, lags, threshold_corr)
thinning.time.std_chain3<-get.thinning(res_time_std_chain3,number.iterations, lags, threshold_corr)
thinning.time.std_chain4<-get.thinning(res_time_std_chain4,number.iterations, lags, threshold_corr)

thinning.time.input_chain1<-get.thinning(res_time_input_chain1,number.iterations, lags, threshold_corr)
thinning.time.input_chain2<-get.thinning(res_time_input_chain2,number.iterations, lags, threshold_corr)
thinning.time.input_chain3<-get.thinning(res_time_input_chain3,number.iterations, lags, threshold_corr)
thinning.time.input_chain4<-get.thinning(res_time_input_chain4,number.iterations, lags, threshold_corr)



percent_burn_in=0.1

remove_burnIn_res_delta_chain1=removeBurnIn(res_delta_chain1,number.iterations, percent_burn_in)
remove_burnIn_res_time_chain1=removeBurnIn(res_time_std_chain1,number.iterations, percent_burn_in)
remove_burnIn_res_time_chain2=removeBurnIn(res_time_std_chain2,number.iterations, percent_burn_in)
remove_burnIn_res_time_chain3=removeBurnIn(res_time_std_chain3,number.iterations, percent_burn_in)
remove_burnIn_res_time_chain4=removeBurnIn(res_time_std_chain4,number.iterations, percent_burn_in)

remove_burnIn_res_time_input_chain1=removeBurnIn(res_time_input_chain1,number.iterations, percent_burn_in, thinning.time.input_chain1)
remove_burnIn_res_time_input_chain2=removeBurnIn(res_time_input_chain2,number.iterations, percent_burn_in, thinning.time.input_chain2)
remove_burnIn_res_time_input_chain3=removeBurnIn(res_time_input_chain3,number.iterations, percent_burn_in, thinning.time.input_chain3)
remove_burnIn_res_time_input_chain4=removeBurnIn(res_time_input_chain4,number.iterations, percent_burn_in, thinning.time.input_chain4)

remove_burnIn_res_time_std=c(remove_burnIn_res_time_chain1, remove_burnIn_res_time_chain2, remove_burnIn_res_time_chain3, remove_burnIn_res_time_chain4)

remove_burnIn_res_time_input=c(remove_burnIn_res_time_input_chain1, remove_burnIn_res_time_input_chain2, remove_burnIn_res_time_input_chain3, remove_burnIn_res_time_input_chain4)



thinned_res_time_chain1=thin.chain.data(res_time_std_chain1,number.iterations, percent_burn_in, thinning.time.std_chain1)
thinned_res_time_chain2=thin.chain.data(res_time_std_chain2,number.iterations, percent_burn_in, thinning.time.std_chain2)
thinned_res_time_chain3=thin.chain.data(res_time_std_chain3,number.iterations, percent_burn_in, thinning.time.std_chain3)
thinned_res_time_chain4=thin.chain.data(res_time_std_chain4,number.iterations, percent_burn_in, thinning.time.std_chain4)

thinned_res_time_input_chain1=thin.chain.data(res_time_input_chain1,number.iterations, percent_burn_in, thinning.time.input_chain1)
thinned_res_time_input_chain2=thin.chain.data(res_time_input_chain2,number.iterations, percent_burn_in, thinning.time.input_chain2)
thinned_res_time_input_chain3=thin.chain.data(res_time_input_chain3,number.iterations, percent_burn_in, thinning.time.input_chain3)
thinned_res_time_input_chain4=thin.chain.data(res_time_input_chain4,number.iterations, percent_burn_in, thinning.time.input_chain4)

max_thinning_time_std = max(thinning.time.std_chain1, thinning.time.std_chain2, thinning.time.std_chain3, thinning.time.std_chain4)
max_thinning_time_input=max(thinning.time.input_chain1, thinning.time.input_chain2, thinning.time.input_chain3, thinning.time.input_chain4)

thinned_res_time_std=c(thinned_res_time_chain1, thinned_res_time_chain2, thinned_res_time_chain3, thinned_res_time_chain4)
thinned_res_time_input=c(thinned_res_time_input_chain1, thinned_res_time_input_chain2, thinned_res_time_input_chain3, thinned_res_time_input_chain4)


hist(thinned_res_time_std, 100000, freq=FALSE, main="", las=1,
     xlab="x", xlim=c(0, 100), ylab="Probability density")

print(summary(thinned_res_time_std))
print(summary(thinned_res_time_input))



library(ggplot2)
require(data.table)
require(dplyr)
require(ggplot2)
require(reshape2)
require(ggpubr)
require(grid)
require(coda)


######################################
source("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/ComparePriorPosterior.R")
number.iterations =500000
lambda1=100
path="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/NoDataFixedT=0.5/0AfixedTimeOriginExpPriors_2priors/500000_not_thinned/ThetaDeltaTTimeOrigin"
chain1path<-paste(path, "log00_new.log", sep="/")
chain2path<-paste(path, "log01_new.log", sep="/")
chai3path<-paste(path, "log02_new.log", sep="/")
chain4path<-paste(path, "log03_new.log", sep="/")
#list_paths=c(chain1path,chain2path, chain3path, chain2path )
list_paths=c(chain1path, chain2path, chai3path, chain4path )
compareExponentialPriorsPosterior(list_paths,number.iterations, 0.1, "theta", lambda1, 0, path)

lambda2=0.01
compareExponentialPriorsPosterior(list_paths,number.iterations, 0.1, "deltaT", lambda2, 0, path)

plotPosterior(list_paths,number.iterations, 0.1,remove_burnIn_res_time_std, max_thinning_time_std,  "time_origin_std", path)

plotPosterior(list_paths,number.iterations, 0.1,remove_burnIn_res_time_input,max_thinning_time_input,  "time_origin_input", path)
