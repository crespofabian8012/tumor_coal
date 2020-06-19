library(ff)
library(data.table)

compareLogUniformPriorsPosterior=function(chainPath_List, thinning, percent_burn_in, column_name, a, b, do_remove_last, path_to_save)
{
  require(data.table)
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
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
  all.dat<-data.frame(prior=x, posterior=posterior_values)
  data<- melt(all.dat)
  
  
  g1<-ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)+ggtitle(paste("Prior vs thinned posterior of ",column_name, " by ", thinning,  sep = " " ))+ labs(y="density", x = column_name)
  g2<-ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
  g3<-ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot()
  ggpubr::ggarrange(g1, g2, g3, 
            labels = c("A", "B", "C"),
            ncol = 2, nrow = 2)
  ggsave(paste(path_to_save, "/", column_name, ".pdf", sep=""))
  ggsave(paste(path_to_save, "/",column_name, ".png", sep=""))
  #ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
  #ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot()
  
}

compareExponentialPriorsPosterior=function(chainPath_List, thinning,percent_burn_in, column_name, lambda, do_remove_last, path_to_save)
{
  require(data.table)
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
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
  
  print(nrow(allchainsdat))
  x <-  rexp(nrow(allchainsdat), rate = lambda)
  
  mean.prior <- mean(x)
  SE.prior   <- sd(x) / sqrt(length(x))
  
  c1 <- rgb(0,0,255,max = 255, alpha = 80, names = "blue1")
  c2 <- rgb(255,255,0, max = 255, alpha = 80, names = "yellow1")
  posterior_values= allchainsdat[[ grep(paste("^", column_name, sep=""), names(allchainsdat), value = TRUE)[1] ]] 

  mean.posterior <- mean(posterior_values)
  
  SE.posterior   <- sd(posterior_values) / sqrt(length(posterior_values))
  print(mean.posterior)
  print(SE.posterior)
  SE      <- sqrt( SE.prior^2 + mean.posterior^2) 
  z       <- (mean.prior - mean.posterior) / SE
  P.z     <- pnorm(z, lower.tail = TRUE)
  print(P.z)
  
  all.dat<-data.frame(prior=log(x), posterior=log(posterior_values))
  data<- melt(all.dat)
  

  g1<-ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)+
          ggtitle(paste("Prior vs posterior of log ", "by ", thinning,  sep = " " ))+ 
          labs(y="density", x =paste("log ", column_name, sep=""))
  g2<-ggplot(data,aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
  g3<-ggplot(data,aes(x=variable, y=value, fill=variable)) + geom_boxplot()
  ggpubr::ggarrange(g1, g2, g3, 
            labels = c("A", "B", "C"),
            ncol = 2, nrow = 2)
  ggsave(paste(path_to_save, "/", column_name, ".pdf", sep=""))
  ggsave(paste(path_to_save, "/",column_name, ".png", sep=""))

}



thinning=100
#path = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/NoDataFixedT=0.5/LogUniformPriorMultiplierKernel"
path = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/NoDataFixedT=0.5/LogUniformPriorMultiplierKernel"

chain1path<-paste(path, "log00.log", sep="/")
chain2path<-paste(path, "log01.log", sep="/")
############################################################
a=-20
b=-10
#compareLogUniformPriorsPosterior(c(chain1path,chain2path ), thinning,0.10, "mut_rate", a, b, 1, path)
# hist(log(x), breaks=100, xlim=c(a,b))
###########################################################
a=7
b=11
#compareLogUniformPriorsPosterior(c(chain1path,chain2path ), thinning,0.10, "effect_pop_size", a, b, 1, path)
#hist(log(x), breaks=100, xlim=c(a, b))
##################################################
a=-20
b=log(log(2))
#compareLogUniformPriorsPosterior(c(chain1path,chain2path ), thinning, 0.10,"growth_rate", a, b, 1, path)
#hist(log(x), breaks=100, xlim=c(a, b))
#####################################################
#log uniform --2 priors
a=-20
b=-3
thinning=2000
path="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/NoDataFixedT=0.5/0AfixedTimeOriginLogUnifPrior_2priors"
chain1path<-paste(path, "log0_new.log", sep="/")
chain2path<-paste(path, "log1_new.log", sep="/")
#compareLogUniformPriorsPosterior(c(chain1path,chain2path ), thinning, 0.10, "theta", a, b, 0, path)

a=-20
b=1
#compareLogUniformPriorsPosterior(c(chain1path,chain2path ), thinning,  0.10, "r", a, b, 0, path)
#######################################################################################
#exponential  --2 priors
thinning=1700
lambda1=100
path="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/NoDataFixedT=0.5/0AfixedTimeOriginExpPriors_2priors/500000_not_thinned"
chain1path<-paste(path, "log00_new.log", sep="/")
chain2path<-paste(path, "log01_new.log", sep="/")
chai3path<-paste(path, "log02_new.log", sep="/")
chain4path<-paste(path, "log03_new.log", sep="/")
compareExponentialPriorsPosterior(c(chain1path,chain2path, chai3path,  chain4path), thinning, 0.10, "theta", lambda1, 0, path)

lambda2=0.01
thinning=1600
compareExponentialPriorsPosterior(c(chain1path,chain2path ), thinning, 0.10, "deltaT", lambda2, 0, path)
