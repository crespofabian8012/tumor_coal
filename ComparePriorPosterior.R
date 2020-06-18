library(ff)
library(data.table)

compareLogUniformPriorsPosterior=function(chainPath_List, thinning, percent_burn_in, column_name, a, b, do_remove_last, path_to_save)
{
  require(data.table)
  require(dplyr)
  list_chains_data<- lapply(chainPath_List, FUN=function(x) data.table::fread(x, header=T,  sep = "\t"))
  filtered_list_chains_data= list_chains_data
  
  if(do_remove_last==1)
  {
    filtered_list_chains_data<- lapply(list_chains_data, FUN= function(x)
    {
      return(x[1:(nrow(x)-1),])
    })
  }
  filtered_list_chains_data2= filtered_list_chains_data
  if(percent_burn_in>0 && percent_burn_in <1)
  {
    filtered_list_chains_data2<- lapply(filtered_list_chains_data, FUN= function(x)
    {
      number_interations_burn_in= floor(percent_burn_in * nrow(x))
      return(x[x$state >number_interations_burn_in,])
    })
  }
  
  list_thinned_chains_data<-lapply(filtered_list_chains_data2, FUN=function(x, thinning)
  {
    indexes<-seq(0,nrow(x),thinning)
    return(x[indexes, ])
  },
  thinning)
  
  allchainsdat<- do.call("rbind", list_thinned_chains_data)
  names(allchainsdat)<-unlist(lapply(names(allchainsdat), function(x) {
    x1= gsub("\\(", "_", x)
    gsub("\\)", "", x1)
  }
  ))
  
  number.breaks= max(floor(nrow(allchainsdat) / 100), 50)
  x <- exp(runif(nrow(allchainsdat), a, b))
  c1 <- rgb(0,0,255,max = 255, alpha = 80, names = "blue")
  c2 <- rgb(255,255,0, max = 255, alpha = 80, names = "yellow")
  posterior_values= allchainsdat[[ grep(paste("^", column_name, sep=""), names(allchainsdat), value = TRUE)[1] ]] 
  histPrior<-hist(x,   breaks=number.breaks , plot = FALSE)
  histPosterior <- hist(posterior_values, breaks=number.breaks, plot = FALSE) 
  i <- which.max(histPosterior$density)

  pdf(paste(path_to_save, "/", column_name, ".pdf", sep=""))
  dev.cur<-dev.cur()
  png(paste(path_to_save, "/",column_name, ".png", sep=""))
  dev.control("enable")
  plot(histPrior, col = c2, main=paste("Prior vs posterior thinned",column_name, sep = " " ),  xlab=column_name,xlim=c(0, exp(b))) # Plot 1st histogram using a transparent color
  plot(histPosterior, col = c1, add = TRUE ) 
  legend("topright", legend=c(paste("Posterior thinned", thinning, sep=" "), "Prior"),
         col=c("blue", "yellow"), lty=1:2, cex=0.8)
  dev.copy(which=dev.cur)
  dev.off()
  dev.off()
  
}

compareExponentialPriorsPosterior=function(chainPath_List, thinning,percent_burn_in, column_name, lambda, do_remove_last, path_to_save)
{
  require(data.table)
  require(dplyr)
  list_chains_data<- lapply(chainPath_List, FUN=function(x) data.table::fread(x, header=T,  sep = "\t"))
  filtered_list_chains_data= list_chains_data
  if(do_remove_last==1)
  {
    filtered_list_chains_data<- lapply(list_chains_data, FUN= function(x)
    {
      return(x[1:(nrow(x)-1),])
    })
  }
  filtered_list_chains_data2= filtered_list_chains_data
  if(percent_burn_in>0 && percent_burn_in <1)
  {
    filtered_list_chains_data2<- lapply(filtered_list_chains_data, FUN= function(x)
    {
      number_interations_burn_in= floor(percent_burn_in * nrow(x))
      return(x[x$state >number_interations_burn_in,])
    })
  }
  
  list_thinned_chains_data<-lapply(filtered_list_chains_data2, FUN=function(x, thinning)
  {
    indexes<-seq(0,nrow(x),thinning)
    return(x[indexes, ])
  },
  thinning)
  
  allchainsdat<- do.call("rbind", list_thinned_chains_data)
  names(allchainsdat)<-unlist(lapply(names(allchainsdat), function(x) {
    x1= gsub("\\(", "_", x)
    gsub("\\)", "", x1)
  }
  ))
  
  number.breaks= max(floor(nrow(allchainsdat) / 100), 50)
  print(nrow(allchainsdat))
  print(number.breaks)
  x <-  rexp(nrow(allchainsdat), rate = lambda)
  c1 <- rgb(0,0,255,max = 255, alpha = 80, names = "blue")
  c2 <- rgb(255,255,0, max = 255, alpha = 80, names = "yellow")
  posterior_values= allchainsdat[[ grep(paste("^", column_name, sep=""), names(allchainsdat), value = TRUE)[1] ]] 
  histPrior<-hist(x, breaks=number.breaks , plot = FALSE)
  histPosterior <- hist(posterior_values, breaks=number.breaks, plot = FALSE) 
  i <- which.max(histPosterior$density)
  
  pdf(paste(path_to_save, "/", column_name, ".pdf", sep=""))
  dev.cur<-dev.cur()
  png(paste(path_to_save, "/",column_name, ".png", sep=""))
  dev.control("enable")
  plot(histPrior, col = c2, main=paste("Prior vs posterior thinned",column_name, sep = " " ),  xlab=column_name,xlim=c(0, max(max(x), max(posterior_values)))) # Plot 1st histogram using a transparent color
  plot(histPosterior, col = c1, add = TRUE ) 
  legend("topright", legend=c(paste("Posterior thinned", thinning, sep=" "), "Prior"),
         col=c("blue", "yellow"), lty=1:2, cex=0.8)
  dev.copy(which=dev.cur)
  dev.off()
  dev.off()
  
}



thinning=100
#path = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/NoDataFixedT=0.5/LogUniformPriorMultiplierKernel"
path = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/NoDataFixedT=0.5/LogUniformPriorMultiplierKernel"

chain1path<-paste(path, "log00.log", sep="/")
chain2path<-paste(path, "log01.log", sep="/")
############################################################
a=-20
b=-10
compareLogUniformPriorsPosterior(c(chain1path,chain2path ), thinning,0.10, "mut_rate", a, b, 1, path)
# hist(log(x), breaks=100, xlim=c(a,b))
###########################################################
a=7
b=11
compareLogUniformPriorsPosterior(c(chain1path,chain2path ), thinning,0.10, "effect_pop_size", a, b, 1, path)
#hist(log(x), breaks=100, xlim=c(a, b))
##################################################
a=-20
b=log(log(2))
compareLogUniformPriorsPosterior(c(chain1path,chain2path ), thinning, 0.10,"growth_rate", a, b, 1, path)
#hist(log(x), breaks=100, xlim=c(a, b))
#####################################################
#log uniform --2 priors
a=-20
b=-3
thinning=2000
path="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/NoDataFixedT=0.5/0AfixedTimeOriginLogUnifPrior_2priors"
chain1path<-paste(path, "log0_new.log", sep="/")
chain2path<-paste(path, "log1_new.log", sep="/")
compareLogUniformPriorsPosterior(c(chain1path,chain2path ), thinning, 0.10, "theta", a, b, 0, path)

a=-20
b=1
compareLogUniformPriorsPosterior(c(chain1path,chain2path ), thinning,  0.10, "r", a, b, 0, path)
#######################################################################################
#exponential  --2 priors
thinning=2500
lambda1=1000
path="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/NoDataFixedT=0.5/0AfixedTimeOriginExpPriors_2priors/500000_not_thinned"
chain1path<-paste(path, "log0.log", sep="/")
chain2path<-paste(path, "log1.log", sep="/")
compareExponentialPriorsPosterior(c(chain1path,chain2path ), thinning, 0.10, "theta", lambda1, 0, path)

lambda2=10
compareExponentialPriorsPosterior(c(chain1path,chain2path ), thinning, 0.10, "r", lambda2, 0, path)
