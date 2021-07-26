rm(list=ls())
step.time.std <- function(current_time_input,current_delta,from, current_theta,current_time_std,
                          f, q, HR_Multiplier, number_accepted, current_coalescent_events) 
  
{
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
  log_f_proposal_coal_events= likelihood_coalescent_times_std(current_delta, proposal_time_std, unlist(current_coalescent_events))
  log_f_current_coal_events= likelihood_coalescent_times_std(current_delta, current_time_std, unlist(current_coalescent_events))
  if (f_proposal >0 && f_current >0 && proposal_time_input >=from)
  {
    log_alpha <- min(0, log(f_proposal) +log_f_proposal_coal_events+  log(HR_q) - log_f_current_coal_events- log(f_current) )
    
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
step.DeltaT <- function(current_deltaT,current_theta,current_time_std,current_coalescent_events, 
                        f, f_delta, f_theta, q_multiplier, HR_Multiplier, number_accepted)
{
  proposal.delta <- q_multiplier(current_deltaT)
  HR_q_delta= HR_Multiplier(proposal.delta, current_deltaT)
  #alpha <- min(1, f(proposal) * HR_q / f(current_value)  )
  f_proposal = f(current_time_std, proposal.delta, 0 )
  f_current = f(current_time_std, current_deltaT, 0 )
  log_f_proposal_coal_events= likelihood_coalescent_times_std(proposal.delta, current_time_std, unlist(current_coalescent_events))
  log_f_current_coal_events= likelihood_coalescent_times_std(current_deltaT, current_time_std, unlist(current_coalescent_events))
  f_delta_proposal_delta = f_delta(proposal.delta)
  f_delta_current_deltaT=f_delta(current_deltaT)

  if (f_delta_proposal_delta >0 && f_proposal >0 && HR_q_delta>0 && f_delta_current_deltaT>0 && f_current>0 )
    {
     log_alpha <- min(0, log(f_delta_proposal_delta) + log(f_proposal) + log_f_proposal_coal_events+log(HR_q_delta) 
                      -log(f_current)- log_f_current_coal_events- log(f_delta_current_deltaT) )
     }
  else{
     log_alpha = -Inf
     }
  #log_alpha <- min(0, logf(proposal) + log(HR_q) -  logf(current_value))
  
  if (log(runif(1)) < log_alpha)
  {
    current_deltaT <- proposal.delta
    number_accepted <- number_accepted +1
  }
  return(list(current_deltaT, number_accepted))
}
step.theta <- function(current_deltaT,current_theta,current_time_input, current_time_std,coalescent_events_time_input,
                       f, f_delta, f_theta, q_multiplier, HR_Multiplier, number_accepted, current_coalescent_events) 
  
{
  proposal.theta <- q_multiplier(current_theta)
  HR_q_theta= HR_Multiplier(proposal.theta, current_theta)
  proposal_time_std= current_time_input / proposal.theta
  f_proposal = f(proposal_time_std, current_deltaT, 0 )
  f_current = f(current_time_std, current_deltaT, 0 )
  
  proposal_coalescent_events = lapply(coalescent_events_time_input, FUN=function(x, proposal.theta) x/proposal.theta, proposal.theta )
  
  log_f_proposal_coal_events= likelihood_coalescent_times_std(current_deltaT, proposal_time_std, unlist(proposal_coalescent_events))
  log_f_current_coal_events= likelihood_coalescent_times_std(current_deltaT, current_time_std, unlist(current_coalescent_events))
  
  #alpha <- min(1, f(proposal) * HR_q / f(current_value)  )
  f_theta_proposal_theta=f_theta(proposal.theta)
  if (f_theta_proposal_theta>0 && f_proposal>0 )
     log_alpha <- min(0,  log(f_theta_proposal_theta)+log(f_proposal)+log_f_proposal_coal_events+log(HR_q_theta) -log(f_current)-log_f_current_coal_events -log(f_theta(current_theta)))
  else
    log_alpha = -Inf
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
likelihood_coalescent_times_std=function(current_delta, current_time_std, current_coalescent_events)
{
  sample.size=length(current_coalescent_events)+1
  num.alive.cells=sample.size
  result=0
  for(i in 1:length(current_coalescent_events))
  {
    timeCurrentEvent= current_coalescent_events[i]
    if (num.alive.cells >1){
      temp=log(num.alive.cells*(num.alive.cells-1)/2)
      result= result +temp
      temp=-1.0* logH(timeCurrentEvent,current_time_std, current_delta)
      result= result +temp
      temp= (num.alive.cells / 2.0)* (num.alive.cells - 1.0) * FmodelSTD(timeCurrentEvent, current_time_std, current_delta)
      result= result -temp
      num.alive.cells= num.alive.cells-1
    }
    
  }
  return(result)
  
  
}
FmodelSTD=function(t, current_time_std, current_delta){
  a = exp(current_delta * t) - 1.0;
  b = 1.0 - exp(-1.0 * current_delta * current_time_std);
  c = 1.0 - exp(-1.0 * current_delta * (current_time_std - t));
  ModelTimeF = a * b / (current_delta * c);
  return(ModelTimeF)
}
logH=function(t, current_time_std, current_delta){
  a = 1.0 - exp(-1.0 * current_delta * (current_time_std - t));
  firstTerm = 2.0 * log(a);
  secondTerm = -1.0 * current_delta * t;
  AboveTerm = firstTerm + secondTerm;
  b = 1.0 - exp(-1.0 * current_delta * current_time_std);
  BelowTerm = 2.0 * log(b);
  logH = AboveTerm - BelowTerm;
  return(logH)
}


run <- function(time_input,coalescent_events_time_input, delta,theta, from,   f, q, f_delta, f_theta, q_multiplier,  
                nsteps, HR_multiplier_bounded,HR_multiplier_unbounded, 
                number_accepted_pair,number_accepted_time ) 
{
  res_time_input <- matrix(NA, nsteps, 1)
  res_delta <- matrix(NA, nsteps, 1)
  res_theta <- matrix(NA, nsteps, 1)
  res_time_std <- matrix(NA, nsteps, 1)
  
  l=list()
  current_delta<-delta
  current_theta<-theta
  current_time_input<-time_input
  current_time_std<-current_time_input / current_theta
  current_coalescent_events = lapply(coalescent_events_time_input, FUN=function(x, current_theta) x/current_theta, current_theta )
  
  number_accepted_delta=0
  number_accepted_theta=0
  number_accepted_time=0
  for (i in seq_len(nsteps))
  {
    print(paste("iteration ",i, sep=""))
    #current_deltaT,current_theta, f_delta, q_delta,f_theta, q_theta, HR_Multiplier, number_accepted
    l <-step.DeltaT(current_delta,current_theta ,current_time_std,current_coalescent_events, 
                    f,f_delta, f_theta, q_multiplier, HR_multiplier_unbounded, number_accepted_delta )
    current_delta<- l[[1]] 
    res_delta[i,] <- current_delta
    number_accepted_delta <- l[[2]]
    l <-step.theta(current_delta,current_theta ,current_time_input, current_time_std,coalescent_events_time_input, 
                   f, f_delta, f_theta, q_multiplier, HR_multiplier_unbounded, number_accepted_theta, current_coalescent_events )
    current_theta<- l[[1]]
    res_theta[i,] <- current_theta
    number_accepted_theta<- l[[2]]
    current_time_std<-current_time_input / current_theta
    current_coalescent_events = lapply(coalescent_events_time_input, FUN=function(x, current_theta) x/current_theta, current_theta )
    #l <- step.time.input(current_time_input,current_delta,from, current_theta, current_time_std, f, q, HR_multiplier_bounded, number_accepted_time )
    l <- step.time.std(current_time_input,current_delta,from, current_theta, current_time_std, 
                       f, q_multiplier, HR_multiplier_unbounded, number_accepted_time, current_coalescent_events )
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
  require(coda)
  if (length(lags)>0 && threshold_corr>0 && threshold_corr <1)
  {
    mcmc_data<-coda::mcmc(data= chain_data, start = 1, end = number.iterations, thin = 1)
    mcmc.autocorr=coda::autocorr(mcmc_data, lags = lags, relative=TRUE)
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
  require(coda)
  mcmc_data<-coda::mcmc(data= chain_data, start = 1, end = number.iterations, thin = 1)
  number_interations_burn_in= floor(percent_burn_in * length(chain_data))
  chain_data_after_burn_in=chain_data[(number_interations_burn_in+1):length(chain_data)]
  return(chain_data_after_burn_in)
}
runChain=function(from,coalescent_events_time_input, lambdaExponentialPriorDelta,lambdaExponentialPriorTheta, number.iterations,  
                  f, q_transformed_multiplier, f_delta, f_theta, q_multiplier, 
                  HR_multiplier_bounded,HR_multiplier_unbounded){
  
  number_accepted_pair=0
  number_accepted_delta=0
  delta=rexp(1,rate=lambdaExponentialPriorDelta)
  theta =rexp(1,rate=lambdaExponentialPriorTheta)
  fromTimeSTD= from *theta
  time_input = draw_random(delta, fromTimeSTD)
  print(time_input)
  list_res <- run(time_input,coalescent_events_time_input, delta,theta, from,  
                  f, q_transformed_multiplier, 
                  f_delta, f_theta, 
                  q_multiplier, number.iterations, 
                  HR_multiplier_bounded,HR_multiplier_unbounded,  number_accepted_pair,number_accepted_time )
  
  
  return(list_res)
}
draw_coalescent_times= function(sample_size, delta, TOrigin)
{
  result_list=list()
  currentTime=0.0
  current_sample_size=sample_size
  i=1
  while(current_sample_size>=2){
      ThisRateCA = current_sample_size * (current_sample_size - 1) / 2.0;
      ThisTimeCA_W = rexp (1, ThisRateCA) ;
      ThisTimeCA_V1 = FmodelTstandard (currentTime, TOrigin, delta);
      ThisTimeCA_V1 = ThisTimeCA_V1 + ThisTimeCA_W;
      
      ThisTimeCA_V2 = GstandardTmodel(ThisTimeCA_V1, TOrigin, delta);
      currentTime = ThisTimeCA_V2
      result_list[i]=currentTime
      i=i+1
      current_sample_size = current_sample_size-1
    
  }
  return(unlist(result_list))
}
FmodelTstandard =function (t, TOrigin, delta)
{
  a = exp(delta * t) - 1.0;
  b = 1.0 - exp(-1.0 * delta * TOrigin);
  c = 1.0 - exp(-1.0 * delta * (TOrigin - t));
  ModelTimeF = a * b / (delta * c);
  return(ModelTimeF)
}
GstandardTmodel=function(V, TOrigin, delta)
{
  firstTerm = TOrigin;
  secondTerm = 1 / delta;
  a =  exp(-1.0 * delta * TOrigin);
  b = (1 - a) * (1 - a) * (1.0 / a);
  c = 1 - a;
  d = (1 - a) * (1.0 / a);
  e = V + d;
  thirdTerm = log(1 - b / e);
  thirdTerm = log(1 - ((1 - a) * (1 - a) * (1.0 / a)) / (V * delta + (1 - a) * (1.0 / a)));
  
  thirdTerm = log(1 + delta * V - a) - log(1 + (delta * V - 1) * a);
  
  StandardTimeG = secondTerm * thirdTerm;
  
  
}
##############################

from= 0.000408897000000000025634 
sample_size=2
#to= 5
#from=0.035
lambdaExponentialPriorDelta= 0.01
lambdaExponentialPriorTheta=100


initial_values_delta=list()
initial_values_theta=list()
initial_values_time_std=list()
initial_values_from=list()


delta=rexp(1,rate=lambdaExponentialPriorDelta)
initial_values_delta[1]=delta
theta =rexp(1,rate=lambdaExponentialPriorTheta)
initial_values_theta[1]=theta
time_std = draw_random(delta, 0)
initial_values_time_std[1]=time_std
time_input = time_std * theta
list_coal_times_std = draw_coalescent_times(sample_size, delta, time_std)
coalescent_events_time_input= list_coal_times_std *theta
from= coalescent_events_time_input[length(coalescent_events_time_input)]
initial_values_from[1]=from
#delta=200
#curve(f(x, 1.0/lambdaExponentialPriorDelta, from) , col="red", n=200, xlim=c(from,0.4), main ="Conditional density  Delta=100, lower bound 0.0004",xlab= "time of origin", ylab="density")

number.iterations=500000
list.values <-  vector("list", number.iterations)
list.values[[1]]<-time
length.interval.multiplier=0.008
length.interval.multiplier=1.5
b=  (length.interval.multiplier)/2.0 + sqrt(1 + (length.interval.multiplier * length.interval.multiplier)* 0.25);



list_res1 <- runChain( time_input,coalescent_events_time_input,   lambdaExponentialPriorDelta,lambdaExponentialPriorTheta, number.iterations,
                       f, q_transformed_multiplier, f_delta, f_theta, q_multiplier,  
                       HR_multiplier_bounded,HR_multiplier_unbounded)

delta=rexp(1,rate=lambdaExponentialPriorDelta)
initial_values_delta[2]=delta
theta =rexp(1,rate=lambdaExponentialPriorTheta)
initial_values_theta[2]=theta
time_std = draw_random(delta, 0)
initial_values_time_std[2]=time_std
time_input = time_std * theta
list_coal_times_std = draw_coalescent_times(sample_size, delta, time_std)
coalescent_events_time_input= list_coal_times_std *theta
from= coalescent_events_time_input[length(coalescent_events_time_input)]
initial_values_from[2]=from

list_res2 <- runChain( time_input,coalescent_events_time_input,  lambdaExponentialPriorDelta,lambdaExponentialPriorTheta, number.iterations,
                       f, q_transformed_multiplier, f_delta, f_theta, q_multiplier,  
                       HR_multiplier_bounded,HR_multiplier_unbounded )

delta=rexp(1,rate=lambdaExponentialPriorDelta)
initial_values_delta[3]=delta
theta =rexp(1,rate=lambdaExponentialPriorTheta)
initial_values_theta[3]=theta
time_std = draw_random(delta, 0)
initial_values_time_std[3]=time_std
time_input = time_std * theta
list_coal_times_std = draw_coalescent_times(sample_size, delta, time_std)
coalescent_events_time_input= list_coal_times_std *theta
from= coalescent_events_time_input[length(coalescent_events_time_input)]
initial_values_from[3]=from

list_res3 <- runChain( time_input,coalescent_events_time_input,  lambdaExponentialPriorDelta,lambdaExponentialPriorTheta, number.iterations,
                       f, q_transformed_multiplier, f_delta, f_theta, q_multiplier,  
                       HR_multiplier_bounded,HR_multiplier_unbounded )

delta=rexp(1,rate=lambdaExponentialPriorDelta)
initial_values_delta[4]=delta
theta =rexp(1,rate=lambdaExponentialPriorTheta)
initial_values_theta[4]=theta
time_std = draw_random(delta, 0)
initial_values_time_std[4]=time_std
time_input = time_std * theta
list_coal_times_std = draw_coalescent_times(sample_size, delta, time_std)
coalescent_events_time_input= list_coal_times_std *theta
from= coalescent_events_time_input[length(coalescent_events_time_input)]
initial_values_from[4]=from

list_res4 <- runChain( time_input, list(from ), lambdaExponentialPriorDelta,lambdaExponentialPriorTheta, number.iterations,
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
plot(res_delta_chain1, type="s", xpd=NA, ylab="Delta", xlab="Sample", las=1)
usr <- par("usr")
xx <- seq(usr[3], usr[4], length=301)


layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
plot(res_theta_chain1, type="s", xpd=NA, ylab="theta", xlab="Sample", las=1)
usr <- par("usr")
xx <- seq(usr[3], usr[4], length=301)


layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
plot(res_time_std_chain1, type="s", xpd=NA, ylab="Time std", xlab="Sample", las=1)
usr <- par("usr")
xx <- seq(usr[3], usr[4], length=301)

hist(res_delta_chain1, 50, freq=FALSE, main="", las=1,
     xlab="x", xlim=c(0, 10000), ylab="Probability density")


hist(res_theta_chain1, 50, freq=FALSE, main="", las=1,
     xlab="x", xlim=c(0, 0.01), ylab="Probability density")


hist(res_time_input_chain1, 50, freq=FALSE, main="", las=1,
     xlab="x", xlim=c(from, 0.1), ylab="Probability density")

hist(res_time_std_chain2, 50, freq=FALSE, main="", las=1,
     xlab="x", ylab="Probability density")

library(coda)
mcmcTime<-coda::mcmc(data= res_time_std_chain2, start = 1, end = number.iterations, thin = 1)
mcmc_list=as.mcmc.list(mcmcTime)
#gelman.diag(mcmc_list, confidence = 0.95,  autoburnin=TRUE,
#         multivariate=TRUE)
lags = c(100, 200, 500, 1000, 2000, 3000, 5000)
threshold_corr=0.1
thinning.delta_chain1<-get.thinning(res_delta_chain1,number.iterations, lags, threshold_corr)
thinning.delta_chain2<-get.thinning(res_delta_chain2,number.iterations, lags, threshold_corr)
thinning.delta_chain3<-get.thinning(res_delta_chain3,number.iterations, lags, threshold_corr)
thinning.delta_chain4<-get.thinning(res_delta_chain4,number.iterations, lags, threshold_corr)

thinning.theta_chain1<-get.thinning(res_theta_chain1,number.iterations, lags, threshold_corr)
thinning.theta_chain2<-get.thinning(res_theta_chain2,number.iterations, lags, threshold_corr)
thinning.theta_chain3<-get.thinning(res_theta_chain3,number.iterations, lags, threshold_corr)
thinning.theta_chain4<-get.thinning(res_theta_chain4,number.iterations, lags, threshold_corr)

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
remove_burnIn_res_delta_chain2=removeBurnIn(res_delta_chain2,number.iterations, percent_burn_in)
remove_burnIn_res_delta_chain3=removeBurnIn(res_delta_chain3,number.iterations, percent_burn_in)
remove_burnIn_res_delta_chain4=removeBurnIn(res_delta_chain4,number.iterations, percent_burn_in)

remove_burnIn_res_theta_chain1=removeBurnIn(res_theta_chain1,number.iterations, percent_burn_in)
remove_burnIn_res_theta_chain2=removeBurnIn(res_theta_chain2,number.iterations, percent_burn_in)
remove_burnIn_res_theta_chain3=removeBurnIn(res_theta_chain3,number.iterations, percent_burn_in)
remove_burnIn_res_theta_chain4=removeBurnIn(res_theta_chain4,number.iterations, percent_burn_in)

remove_burnIn_res_time_chain1=removeBurnIn(res_time_std_chain1,number.iterations, percent_burn_in)
remove_burnIn_res_time_chain2=removeBurnIn(res_time_std_chain2,number.iterations, percent_burn_in)
remove_burnIn_res_time_chain3=removeBurnIn(res_time_std_chain3,number.iterations, percent_burn_in)
remove_burnIn_res_time_chain4=removeBurnIn(res_time_std_chain4,number.iterations, percent_burn_in)

remove_burnIn_res_time_input_chain1=removeBurnIn(res_time_input_chain1,number.iterations, percent_burn_in)
remove_burnIn_res_time_input_chain2=removeBurnIn(res_time_input_chain2,number.iterations, percent_burn_in)
remove_burnIn_res_time_input_chain3=removeBurnIn(res_time_input_chain3,number.iterations, percent_burn_in)
remove_burnIn_res_time_input_chain4=removeBurnIn(res_time_input_chain4,number.iterations, percent_burn_in)

remove_burnIn_res_delta=c(remove_burnIn_res_delta_chain1, remove_burnIn_res_delta_chain2, remove_burnIn_res_delta_chain3, remove_burnIn_res_delta_chain4)
remove_burnIn_res_theta=c(remove_burnIn_res_theta_chain1, remove_burnIn_res_theta_chain2, remove_burnIn_res_theta_chain3, remove_burnIn_res_theta_chain4)

remove_burnIn_res_time_std=c(remove_burnIn_res_time_chain1, remove_burnIn_res_time_chain2, remove_burnIn_res_time_chain3, remove_burnIn_res_time_chain4)

remove_burnIn_res_time_input=c(remove_burnIn_res_time_input_chain1, remove_burnIn_res_time_input_chain2, remove_burnIn_res_time_input_chain3, remove_burnIn_res_time_input_chain4)

thinned_res_delta_chain1=thin.chain.data(res_delta_chain1,number.iterations, percent_burn_in, thinning.delta_chain1)
thinned_res_delta_chain2=thin.chain.data(res_delta_chain2,number.iterations, percent_burn_in, thinning.delta_chain2)
thinned_res_delta_chain3=thin.chain.data(res_delta_chain3,number.iterations, percent_burn_in, thinning.delta_chain3)
thinned_res_delta_chain4=thin.chain.data(res_delta_chain4,number.iterations, percent_burn_in, thinning.delta_chain4)

thinned_res_theta_chain1=thin.chain.data(res_theta_chain1,number.iterations, percent_burn_in, thinning.theta_chain1)
thinned_res_theta_chain2=thin.chain.data(res_theta_chain2,number.iterations, percent_burn_in, thinning.theta_chain2)
thinned_res_theta_chain3=thin.chain.data(res_theta_chain3,number.iterations, percent_burn_in, thinning.theta_chain3)
thinned_res_theta_chain4=thin.chain.data(res_theta_chain4,number.iterations, percent_burn_in, thinning.theta_chain4)

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

max_thinning_delta = max(thinning.delta_chain1, thinning.delta_chain2, thinning.delta_chain3, thinning.delta_chain4)
max_thinning_theta=max(thinning.theta_chain1, thinning.theta_chain2, thinning.theta_chain3, thinning.theta_chain4)


thinned_res_delta=c(thinned_res_delta_chain1, thinned_res_delta_chain2, thinned_res_delta_chain3, thinned_res_delta_chain4)
thinned_res_theta=c(thinned_res_theta_chain1, thinned_res_theta_chain2, thinned_res_theta_chain3, thinned_res_theta_chain4)

thinned_res_time_std=c(thinned_res_time_chain1, thinned_res_time_chain2, thinned_res_time_chain3, thinned_res_time_chain4)
thinned_res_time_input=c(thinned_res_time_input_chain1, thinned_res_time_input_chain2, thinned_res_time_input_chain3, thinned_res_time_input_chain4)


hist(thinned_res_time_std, 100000, freq=FALSE, main="", las=1,
     xlab="x", xlim=c(0, 100), ylab="Probability density")

print(summary(thinned_res_delta_chain1))
print(summary(thinned_res_delta_chain2))  
print(summary(thinned_res_delta_chain3))
print(summary(thinned_res_delta_chain4))
print(summary(thinned_res_theta))
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
