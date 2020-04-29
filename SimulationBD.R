#simulation BD comparison 
#following Prof. Carsten idea
# author Fausto Fabian Crespo Fernandez
compute.gamma=function(lambda, mu, rho)
{
  output= abs(lambda-mu)/max(rho*lambda, mu-(1-rho)*lambda)
  return(output)
}
log.conditional.probability.number.ancestor.population=function(k, m, mprime)
  {
  #mprime is the number of ancestors of the population the first time there will be  k-1 ancestors in the sample
  #m is the current number of ancestors of the population and k current number of ancestorc in the sample size
  #so mprime >= k-1 and also we should have that mprime <= m
  #output = log(k) + log(k-1) + logfact(m+k-1) + logfact(m-k) - logfact(mprime+k) -logfact(mprime-k+1)
  output = log(k) + log(k-1) + log.prod.between(mprime-k+1,m+k-1)  + log.prod.between(mprime+k,m-k) 
  #print(output)
  #output=output + logfact(mprime) + logfact(mprime-1) - logfact(m) - logfact(m-1)
  output=output + log.prod.between(m,mprime) + log.prod.between(m-1,mprime-1)
  #print(output)
  return(output)
}

log.conditional.probability.number.ancestor.population.present.time=function(n, mprime)
{
  output=log(n)+log(n-1)+logfact(mprime)+logfact(mprime-1)- logfact(mprime-n+1)-logfact(mprime+n)
    #sum(log((mprime):(mprime-n+1)))+sum(log((mprime-1):(mprime+n)))
    #logfact(mprime)+logfact(mprime-1)- logfact(mprime-n+1)-logfact(mprime+n)
  return(output)
}
fact=function(n){
  return(prod(1:n))
}
logfact=function(n){
  return(sum(log(1:n)))
}
log.prod.between=function(from,to){
  if (to >=from){
    return(sum(log(max(from,1):to)))
  }
  else
    return(-1*sum(log(max(to,1):from)))
}
ToStandardTime<-function(t, Time1, Delta )
{
  a=exp(-1*Delta*Time1)
  c=exp(Delta*t)
  result<-(c-1)*(1-a)/(1-a*c)
  print(result)
  return(result)
}
ToModelTime<-function(u, Time1, Delta )
{
  a=exp(-1*Delta*Time1)
  b=1+u-a
  c=1+(u-1)*a
  if ((b== 0) || (c==0))
  {
    return((1.0/Delta)*(log(1+u)-(u*u*a)/(1+u)))
  } 
  else{
    return((1.0/Delta)*log(b)-log(c))
  }
}
sample.conditional.coalescent.time.distribution = function(Delta, alpha, k, m){
  require(stats)
     U= rbeta(1, k, m-k+1)
     T=(-1.0/Delta)* log(U/(alpha+(1-alpha)*U))
     return(T)
}
sample.first.coalescent.time.distribution = function(Delta, Gamma, k){
  require(stats)
  U= rgamma(1, shape=k, rate=1)
  T=(-1.0 / Delta) * log(1- (Gamma/(U+Gamma)))
  return(T)
}
density.time.origin = function(delta, gamma, n, t){
  exp.delta.t = exp(-1*delta* t)
  output =n* delta* gamma* exp.delta.t*(1- exp.delta.t)**(n-1)/ (1+ (gamma-1)*(exp.delta.t)**(n+1))
  return(output)
}
conditional.density.first.time.k.ancestors = function(delta, gamma, n, time.Origin, times){
  exp.delta.t = exp(-1*delta* time.Origin)
  output = ((delta* gamma)**(n-1)) * (1+ (gamma-1)*exp.delta.t)**(n-1) / (1- exp.delta.t)**(n-1)
  for(i in  2:n)
    {
    if (i<n){
      
      output = i* output *exp(-1*delta* time.Origin) /  (1+ (gamma-1)*(exp.delta.t)**2)
    }
     else{
       
       output = output *exp(-1*delta* time.Origin) /  (1+ (gamma-1)*(exp.delta.t)**2)
     }
   }
  return(output)
}
density.first.time.k.ancestors = function(delta, gamma, n, time.Origin, times){
  exp.delta.t = exp(-1*delta* time.Origin)
  output = ((delta* gamma)**(n)) 
  for(i in  1:n)
  {
      exp.delta.t = exp(-1*delta* times[i])
      output = i* output *exp(-1*delta* time.Origin) /  (1+ (gamma-1)*((exp.delta.t)**2))
 
  }
  return(output)
}

get.most.probable.number.of.ancestors.population.when.k.minus.1.ancestors.sample=function(k,m)
  #k and m current sample  and population size respectively
  {
    probs=list()
    mprime=k
    i=0
    u=runif(1,0,1)
    #print("u")
    #print(u)
    mprime =(k-1)
    while(sum(unlist(probs)) < u && mprime< m )
      {
      currentProb=log.conditional.probability.number.ancestor.population(k, m,mprime )
      probs[mprime-(k-1)+1]=exp(currentProb)
      
      mprime=mprime+1
      #print(sum(unlist(probs)))
    }
    return(mprime)
}

get.most.probable.number.of.ancestors.population.present.time=function(sample.size)
{
  probs=list()
  #while(sum(unlist(probs)) <= 0.8)
  i=0
  probs=list()
  mprime=k
  i=0
  u=runif(1,0,1)
  #print("u")
  #print(u)
  mprime =sample.size
  while(sum(unlist(probs)) < u )
  {
    currentProb=log.conditional.probability.number.ancestor.population.present.time(sample.size, mprime)
    probs[mprime-sample.size+1]=exp(currentProb)
    
    mprime=mprime+1
    #print(sum(unlist(probs)))
  }
  return(mprime)
 
}


simulate.list.number.ancestors.population=function(sample.size){
  population.present.time= get.most.probable.number.of.ancestors.population.present.time(sample.size)
  current.population.size=population.present.time
  print(current.population.size)
  list.populations.sizes=list(current.population.size)
  for(i in  (sample.size):2)
  {
    current.population.size= get.most.probable.number.of.ancestors.population.when.k.minus.1.ancestors.sample(i,current.population.size)
    print(current.population.size)
    list.populations.sizes[[ length(list.populations.sizes)+1]]=current.population.size
  }
  
  return( unlist(list.populations.sizes))
}

simulate.coalescent.times.A1=function(lambda, mu, rho,sample.size, list.number.ancestors.population)
  {
      gamma=compute.gamma(lambda, mu, rho)
      delta=lambda - mu
      Gamma= sample.size * gamma
      Delta= sample.size * delta
      list.coal.times=list()
      u=unlist(sample.first.coalescent.time.distribution(Delta, Gamma, sample.size))
      u=u[[1]]
      list.coal.times[1]=u
      for(i in  (sample.size-1):2)
        {
         m= list.number.ancestors.population[sample.size-i]
         alpha=1 - exp(-1.0 * Delta * u)
         u=unlist(sample.conditional.coalescent.time.distribution(Delta, alpha, i, m))
         list.coal.times[length(list.coal.times)+1]=u[[1]]
         u=u[[1]]
         
        }
      return(unlist(list.coal.times))
}
#############################################################################################
#Scenario A with stochastic population size
sample.size=20
list.number.ancestors.population=simulate.list.number.ancestors.population(sample.size)
list.number.ancestors.population

lambda=1
mu=0.99
rho= 0.8
list.coal.times= simulate.coalescent.times.A1(lambda, mu, rho,sample.size, list.number.ancestors.population)
list.coal.times


#Scenario A1 with deterministic  population size(expected population size)