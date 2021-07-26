
estimateAutoCorrTime=function( values)
{
  sum = 0.0;
  autocorrelationFunction = estimateAutocorrelationFunction(values)
  for (i in  1:length(autocorrelationFunction))
    sum =sum + 2 * autocorrelationFunction[[i]]
  return(1.0 + sum)
}

estimateAutocorrelationFunction=function( series)
{
  sampleVar = var(series)
  data = auto =list(rep(0,length(sampleVar)))
  for (i in  1:length(series))
    data[[i]] = series[[i]]
  auto = autocovariance(data)
  result = list()
   n = length(series)
  factor = 1.0 / n / sampleVar
  result[[1]]=1.0
  for (i in  2:n)
  {
    if (i > 1 && auto[[i]]+ auto[[i-1]] < 0.0) 
      break
    result[[length(result)+1]]= factor * auto[[i]]
  }
  return(result)
}

autocovarianceShift=function( x, maxShift) 
{
  total = sum(x)
  mean = mean(x)
  stop = min(length(x), maxShift);
  auto = list(rep(0, stop))
  for (i in  1:stop)
    {
    for (j in  1:(length(x) - i)) 
      {
      auto[[i]] = auto[[i]] + (x[[j]]-mean) * (x[[j + i]]-mean);
      }
    }
  return(auto)
}
autocovariance=function( x) 
  {
  return(autocovarianceShift(x, length(x)))
}
step <- function(current_value, f, q, HastingsRatioFun, number_accepted) {
  proposal <- q(current_value)
  HR_q=HastingsRatioFun(proposal, current_value)
  alpha <- min(1, f(proposal) * HR_q / f(current_value)  )
  if (runif(1) < alpha)
    {
    current_value <- proposal
    number_accepted <- number_accepted +1
  }
  return(list(current_value, number_accepted))
}
HR_Multiplier=function(proposal, current_value) {proposal / current_value}

run <- function(x, f, q, nsteps, HastingsRatioFun, number_accepted) {
  res <- matrix(NA, nsteps, length(x))
  l=list()
  for (i in seq_len(nsteps))
  {
    l <- step(x, f, q, HastingsRatioFun, number_accepted )
    x <- l[[1]] 
    res[i,] <- l[[1]]  
    number_accepted <- l[[2]] 
  }
  print(number_accepted)
  drop(res)
}
f <- function(x)
  {dexp(x, lambda, log=FALSE )}

q <- function(current.value)
{
  u=runif(1,0,1)
  m=exp(2*log(b)*(u-0.5))
  proposed.value= m * current.value
  proposed.value
  }


number.iterations=500000
lambda=100
initial.value=rexp(1,lambda)
list.values <-  vector("list", number.iterations)
list.values[[1]]<-initial.value
length.interval.multiplier=0.2
b=  (length.interval.multiplier)/2.0 + sqrt(1 + (length.interval.multiplier * length.interval.multiplier)* 0.25);
number_accepted=0

res <- run(initial.value, f, q, number.iterations, HR_Multiplier, number_accepted)

layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
plot(res, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)
usr <- par("usr")
xx <- seq(usr[3], usr[4], length=301)
plot(f(xx), xx, type="l", yaxs="i", axes=FALSE, xlab="")



hist(res, 50, freq=FALSE, main="", ylim=c(0, 100), las=1,
     xlab="x", ylab="Probability density")
z <- integrate(f, 0, Inf)$value
curve(f(x) / z, add=TRUE, col="red", n=200)

###########################################################################
  
#lapply(2:number.iterations, FUN=function(i,  b,lambda, list.values)
for (i in 2:number.iterations)
{
    #print(paste("iteration ",i , sep =""))
    current.value = list.values[[i-1]]
    u=runif(1,0,1)
    m=exp(2*log(b)*(u-0.5))
    proposed.value= m * current.value
    logMH= log(dexp(proposed.value, lambda, log=FALSE ))-log(dexp(current.value, lambda, log=FALSE ))+log(m)
    if (logMH >= 0)
      {logMH=0}
    r=runif(1,0,1)
    if (log(r)<logMH)
    {
      list.values[[i]]<-proposed.value
    }
    else
    {
      list.values[[i]]<-current.value
    }
    
  }
#, b, lambda, list.values)

mean(unlist(list.values))
sd(unlist(list.values))
#autocovariance(unlist(list.values))
hist(unlist(list.values), 50, freq=FALSE, main="", ylim=c(0, 110), las=1,
     xlab="x", ylab="Probability density")
z <- integrate(f, 0, Inf)$value
curve(f(x) / z, add=TRUE, col="red", n=200)

prior.values <-  rexp(length(list.values), rate = lambda)
print(ks.test(prior.values, unlist(list.values),
              alternative = c("two.sided")))
