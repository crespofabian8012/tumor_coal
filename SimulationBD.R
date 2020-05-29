#simulation BD comparison 
#following Prof. Carsten Wiuf  idea
# author Fausto Fabian Crespo Fernandez
#############################################
#functions for the cancer case
# Relative population size

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
relative <- function(Torigin,Delta,t)
  exp(-Delta*t)*(1-exp(-Delta*(Torigin-t)))^2/(1-exp(-Delta*Torigin))^2

# random origin
torigin <- function(Delta) log(1+Delta/rexp(1,rate=1))/Delta

# coalescent times with random origin
# returns origin plus n-1 times until i=1,..,n-1 ancestors
sampletimes <- function(n,Delta) {
  t <- torigin(Delta)
  w <- modelCoal(n,t,Delta)
  c(t,unlist( lapply(1:(n-1),function(i) sum(w[i:(n-1)])) ))
}

# Standard time to model time
standard2model <- function(Torigin,Delta,u) 
  ( log(1+Delta*u-exp(-Delta*Torigin)) - log(1+(Delta*u-1)*exp(-Delta*Torigin)) )/Delta

# Model time to standard time
model2standard <- function(Torigin,Delta,t)
  (exp(Delta*t)-1)*(1-exp(-Delta*Torigin))/(1-exp(-Delta*(Torigin-t)))/Delta

# Check conversion of time is correct (result should be x)
x <- 2
model2standard(10,1,standard2model(10,1,x))

# Coalescence times with constant population size
standardCoal <- function(n) unlist(lapply(2:n,function(i) rexp(1,rate=choose(i,2))))

# Coalescence times in model time
modelCoal <- function(n,Torigin,Delta) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  u <- standard2model( Torigin,Delta,unlist( lapply(1:(n-1),function(i) sum(w[i:(n-1)])) ) )
  u <- c(u,0)
  unlist(lapply(1:(n-1),function(i) u[i]-u[i+1]))
}

# Coalescence times with exponentially decreasing popualtion size
expCoal <- function(n,beta) {
  w <- standardCoal(n)
  # time until i=n-1,...,1 ancestors
  u <- unlist( lapply(1:(n-1),function(i) sum(w[i:(n-1)])) )
  u <- log(1+beta*u)/beta
  u <- c(u,0)
  unlist(lapply(1:(n-1),function(i) u[i]-u[i+1]))
}

#######################################################################
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
  output=log(n)+log(n-1)+log.prod.between(mprime-n+1,mprime)+log.prod.between(mprime+n,mprime-1)
  #  logfact(mprime)+logfact(mprime-1)- logfact(mprime-n+1)-logfact(mprime+n)
  return(output)
}
fact=function(n){
  return(prod(1:n))
}
logfact=function(n){
  return(sum(log(1:n)))
}
log.prod.between=function(from,to){
  if (min(from, to)>0)
  {
    if (to >=from ){
      return(sum(log(max(from,1):to)))
    }
    else 
      return(-1*sum(log(max(to,1):from)))
    
  }
 else 
   return(0)
}
ToStandardTime<-function(t, Time1, Delta )
{
  a=exp(-1*Delta*Time1)
  c=exp(Delta*t)
  result<-(c-1)*(1-a)/(1-a*c)
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
  U= stats::rgamma(1, shape=k, rate=1)
  T=(-1.0 / Delta) * log(1.0 - (Gamma/(U+Gamma)))
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
  if ( mprime< m){
    return(mprime)
  }
  else if( mprime== m )
  {
    
    return(m)
  }
  else{
    
    return(-1)
  }
}
get.most.probable.number.of.ancestors.population.when.k.minus.1.ancestors.sample2=function(k,m)
{
  accepted=FALSE
  iteration=1
  while(!accepted){
    print(paste("iteration",iteration, sep=" "))
    u=runif(1,0,1)
    x= 2*(k-1)/(1-u) - 2*(k-1)
    mprime= k-2+ floor(x)
    iteration= iteration+1
    if (mprime>0 &&  mprime >=(k-1) && mprime <m )
    {
      
      log.prob= log.prod.between(mprime-k+1,mprime-1)+ log.prod.between(mprime+k, mprime) 
      log.prob = log.prob + 2* log(x+2*(k-1))
      prob= exp(log.prob) 
      v=runif(1,0,1)
  
      print(paste0(" log.prob = ", log.prob," prob = ", prob, " x = ", x, " v = ", v, " mprime =", mprime, " m = ", m, " k = ", k, sep=" "))
      if (is.nan(prob) || prob > 3 || iteration > 500 )
        {
        print("something went wrong!")
      }
      if (v<prob){
        accepted=TRUE
      }
    }
    else{
      accepted=FALSE
    }
  }
  return(mprime)
}
solve.cumulative.density.number.ancestors.population.present.time=function(u, sample.size)
{
  require(rootSolve)
  fun <- function (x)  exp(log.prod.between(x-sample.size+2,x)+log.prod.between(x+sample.size,x+2)) - u
  #x <- rootSolve::uniroot.all(fun, c(sample.size-2, 1000*sample.size))
  x = sample.size 
  while(fun(x)<0)
     x=x+1
  
  return(x)

}
solve.cumulative.density.number.ancestors.population=function(u, k,m)
{
   require(rootSolve)
   fun <- function (x)  exp(log.prod.between(x+1-k,x+1)+log.prod.between(x+k,x)+ log.prod.between(m-1,m+k-1) + log.prod.between(m,m-k)) - u
    #x <- rootSolve::uniroot.all(fun, c(sample.size-2, 1000*sample.size))
   x = sample.size -1
   while(fun(x)<0)
     x=x+1
   
   return(x)
}
get.most.probable.number.of.ancestors.population.present.time2=function(sample.size)
{

    u=runif(1,0,1)
    print(u)
    mprime=solve.cumulative.density.number.ancestors.population.present.time(u, sample.size)
    return(mprime)
  
}
get.most.probable.number.of.ancestors.population.when.k.minus.1.ancestors.sample3=function(k,m)
{
  accepted=FALSE
  while(!accepted)
  {
    u=runif(1,0,1)
    mprime=solve.cumulative.density.number.ancestors.population(u, k, m)
    
      if(mprime <m)
      {
        accepted=TRUE
      }
 
 
  }
  return(mprime)
}
get.most.probable.number.of.ancestors.population.present.time=function(sample.size)
{
  probs=list()
  #while(sum(unlist(probs)) <= 0.8)
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
  population.present.time= get.most.probable.number.of.ancestors.population.present.time2(sample.size)
  current.population.size=population.present.time
  list.populations.sizes=list(current.population.size)
  list.populations.sizes = c(list(current.population.size), lapply((sample.size):2, FUN =
    get.most.probable.number.of.ancestors.population.when.k.minus.1.ancestors.sample3,m=current.population.size))
  # for(i in  (sample.size):2)
  # {
  #   current.population.size= get.most.probable.number.of.ancestors.population.when.k.minus.1.ancestors.sample(i,current.population.size)
  #   list.populations.sizes[[ length(list.populations.sizes)+1]]=current.population.size
  # }
  
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
  #for(i in  (sample.size-1):1)
  for(i in  (sample.size-1):1)
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
sim=10
sample.size=20
lambda=1
mu=0.99
rho= 0.8
list.number.ancestors.population.sim <- array(0,dim=c(sim,sample.size))
coal.events.times.sim <- array(0,dim=c(sim,sample.size))
library(future.apply)
#Scenario A with stochastic population size
invisible(lapply(1:sim, FUN=function(i, sample.size, lambda, mu, rho, coal.events.times.sim, list.number.ancestors.population.sim)
  {
 
   list.number.ancestors.population= simulate.list.number.ancestors.population(sample.size)
   print(unlist(list.number.ancestors.population))
  list.number.ancestors.population.sim[i,] <<-list.number.ancestors.population
  list.coal.times= simulate.coalescent.times.A1(lambda, mu, rho,sample.size, list.number.ancestors.population)
  coal.events.times =cumsum(list.coal.times)
  print(unlist(coal.events.times))
  #coal.events.times=c(0,coal.events.times)
  coal.events.times.sim[i,] <<- coal.events.times
  print(paste0("finished sim",i, sep=" "))
      
    },
  sample.size, lambda, mu,rho, coal.events.times.sim=coal.events.times.sim,
  list.number.ancestors.population.sim=list.number.ancestors.population.sim))
  


mean.events.times <- colSums(coal.events.times.sim)/sim
mean.number.population.ancestors <- rev(colSums(list.number.ancestors.population.sim)/sim)
saveRDS(mean.events.times, "~/project3/test/mean.events.times.rds")
saveRDS(mean.number.population.ancestors, "~/project3/test/mean.events.times.rds")
population=colSums(list.number.ancestors.population.sim)/sim

padded.coal.events.times.sim=cbind(rep(0,sim) , coal.events.times.sim)
padded.list.number.ancestors.population.sim=cbind( list.number.ancestors.population.sim, rep(0,sim))
rev.padded.coal.events.times.sim<-t(apply(padded.coal.events.times.sim,1,function(x) rev(x)))
all.times.population=as.vector(t(padded.coal.events.times.sim))
sample.sizes=rep(0:sample.size,sim+1)
rev.padded.list.number.ancestors.population.sim<-t(apply(padded.list.number.ancestors.population.sim,1,function(x) rev(x)))
population.sizes=as.vector(t(rev.padded.list.number.ancestors.population.sim))
num.sim=rep(1:(sim), rep(sample.size+1, sim))
population.dataframe=data.frame(sample.sizes, event.time=c(all.times.population,0,mean.events.times) , 
                                population.size=c(population.sizes, 0,mean.number.population.ancestors), num.sim=c(num.sim, rep(sim+1,sample.size+1)))
population.dataframe$simulation= paste0("simulation", population.dataframe$num.sim)


pdf("~/project3/test/plotScenarioAStochastic5.pdf")
plot(mean.events.times,1:sample.size,ylab="Origin (i=1), Number of ancestors (i>1)",xlab="Time until i-1 ancestors",pch=19,col="orange")

dev.off()


pdf("~/project3/test/plotScenarioAStochastic2.pdf")
plot(0:(sample.size-1),mean.events.times,xlab="Origin (i=1), Number of ancestors (i>1)",ylab="Time until i-1 ancestors",pch=19,col="orange")

dev.off()

pdf("~/project3/test/plotScenarioAStochastic3.pdf")
plot(population,mean.events.times,xlab="Origin (i=1), Number of ancestors population(i>1)",ylab="Time until i-1 ancestors",pch=19,col="orange")

dev.off()
#Scenario A with deterministic  population size(expected population size)
# for (k in 1:length(Torigin)) {
#   for (i in 1:sim) times[i,] <- modelCoal(n,Torigin[k],Delta)
#   xx <- colSums(times)/sim
# 
# }

#Scenario B with deterministic  population size(expected population size)

library(ggplot2)

# This example uses the ChickWeight dataset, which comes with ggplot2
# First plot
dev.off()

plot.new()
p1 <- ggplot(population.dataframe, aes(x=event.time, y=population.size))+
  geom_line(aes(colour=simulation))+
  ylim(19,2000)+
  xlim(0,310)+
  xlab("Time(forward)") +
  ylab("Number of population ancestors") +
  theme(legend.position = "none") +
  ggtitle("Stochastic population growth")
  #geom_point(aes(x=c(0,mean.events.times),y=c(0,mean.number.population.ancestors)), 
  #     pch = 16,
  #      colour = "yellow")+
p1 + theme(legend.position = "none") 

# Second plot
p2 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet)) +
  geom_point(alpha=.3) +
  geom_smooth(alpha=.2, size=1) +
  ggtitle("Fitted growth curve per diet")

# Third plot
p3 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, colour=Diet)) +
  geom_density() +
  ggtitle("Final weight, by diet")

# Fourth plot
p4 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, fill=Diet)) +
  geom_histogram(colour="black", binwidth=50) +
  facet_grid(Diet ~ .) +
  ggtitle("Final weight, by diet") +
  theme(legend.position="none")        # No legend (redundant in this graph)    

multiplot(p1, p2, p3, p4, cols=2)