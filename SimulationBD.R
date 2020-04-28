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
  #m is the current number of ancestors of the population and k current sample size
  #so mprime >= k-1 and also we should expect that mprime <= m
  output=log(k)+log(k-1)+logfact(m-k)+logfact(m+k-1)-logfact(mprime-k+1)-logfact(mprime+k)
  #print(output)
  output=output+logfact(mprime)+logfact(mprime-1)- logfact(m)-logfact(m-1)
  #print(output)
  return(output)
}

log.conditional.probability.number.ancestor.population.present.time=function(n, mprime)
{
  output=log(n)+log(n-1)+logfact(mprime)+logfact(mprime-1)- logfact(mprime-n+1)-logfact(mprime+n)
  return(output)
}
fact=function(n){
  if(n<=1) return(1)
  output=1
  for(i in 1:n) {
    output = output * i
  }
  return(output)
}
logfact=function(n){
  if(n<=1) return(0)
  output=0
  for(i in 1:n) {
    output = output +log(i)
  }
  return(output)
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
simulate.distribution = function(Delta, alpha, k, m){
  require(stats)
     U= rbeta(1, k, m-k+1)
     T=(-1.0/Delta)* log(U/(alpha+(1-alpha)*U))
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
    #while(sum(unlist(probs)) <= 0.8)
    i=0
    for(mprime in k:m)
      {
      #print(mprime)
      currentProb=log.conditional.probability.number.ancestor.population(k, m,mprime )
      #print(currentProb)
      
      probs[mprime-k+1]=exp(currentProb)
      #mprime=mprime+1
       
    }
   # print(sum(unlist(probs)))
    norm.probabilities= unlist(probs)[1:length(unlist(probs))] / sum(unlist(probs))
    output=sample(k:m, 1, prob= norm.probabilities)
    return(output)
}

get.most.probable.number.of.ancestors.population.present.time=function(sample.size)
{
  probs=list()
  #while(sum(unlist(probs)) <= 0.8)
  i=0
  for(mprime in sample.size:(1000*sample.size))
  {
    #print(mprime)
    currentProb=log.conditional.probability.number.ancestor.population.present.time(sample.size, mprime)
    #print(currentProb)
    probs[mprime-sample.size+1]=exp(currentProb)
    #mprime=mprime+1
    
  }
 # print(sum(unlist(probs)))
  norm.probabilities= unlist(probs)[1:length(unlist(probs))] / sum(unlist(probs))
  output=sample(sample.size:(1000*sample.size), 1, prob= norm.probabilities )
  return(output)
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
#############################################################################################
sample.size=20
list.number.ancestors.population=simulate.list.number.ancestors.population(sample.size)
list.number.ancestors.population