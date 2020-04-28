#simulation BD comparison 
#following Prof. Carsten idea
# author Fausto Fabian Crespo Fernandez
compute.gamma=function(lambda, mu, rho)
{
  output= abs(lambda-mu)/max(rho*lambda, mu-(1-rho)*lambda)
  return(output)
}
log.conditional.probability.number.ancestor.population=function(num.ancestors.samplekminus1, m, mprime)
  {
  k= num.ancestors.samplekminus1 +1
  output=log(k)*log(k-1)+logfact(m-k)+logfact(m+k+1)-logfact(mprime-k+1)-logfact(mprime+k)
  print(output)
  output=output+logfact(mprime)+logfact(mprime-1)- logfact(m)-logfact(m-1)
  print(output)
  return(output)
}
conditional.probability.number.ancestor.population.present.time=function(n, m, mprime)
{

  output=n*(n-1)*(m-k)*fact(mprime)*fact(mprime-1)/ (fact(mprime-n+1)*fact(mprime+n))
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

get.number.of.ancestors.population.when.k.ancestors.sample=function(k,m)
  {
    probs=list()
    for(i in k:1000*k) 
      {
      currentProb=log.conditional.probability.number.ancestor.population(k, m,i )
      print(currentProb)
      if (exp(currentProb) < 0.00000001)
      {
        break;
        }
      
      probs[i]=currentProb
    }
    output=sample(k:(k+length(probs)), 1, prob= unlist(exp(probs)))
    return(output)
}


get.number.of.ancestors.population.when.k.ancestors.sample(3,100)

