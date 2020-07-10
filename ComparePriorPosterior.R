library(ff)
library(data.table)

compareLogUniformPriorsPosterior=function(chainPath_List, thinning, percent_burn_in, column_name, a, b, do_remove_last, path_to_save)
{
  require(data.table)
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
  library(ggplot2)
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

compareExponentialPriorsPosterior=function(chainPath_List,number.iterations,percent_burn_in, column_name, lambda, do_remove_last, path_to_save)
{
  require(data.table)
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
  require(grid)
  require(coda)
  list_chains_data<- lapply(chainPath_List, FUN=function(x)
    as.data.frame(data.table::fread(x, header=T,  sep = "\t", colClasses=list(numeric=1:12))))
  list_chains_data<- lapply(list_chains_data, FUN = function(x)
           x[, colSums(is.na(x)) != nrow(x)])
  
  filtered_list_chains_data= list_chains_data
  
  mcmc_list=build_mcmc_list(list_chains_data,number.iterations)
  
  #gelman.diag.output<-gelman.diag(mcmc_list, confidence = 0.95,multivariate=FALSE)
  #print(gelman.diag.output)
  
  lags = c(100, 200, 500,1000)
  list.autocorr= lapply(mcmc_list, FUN=function(x, lags, column_name) 
  {
    column.index=which(names(as.data.frame(x)) == grep(paste("^", column_name, sep=""), names(as.data.frame(x)), value = TRUE)[1] )
    chain_values= x[,column.index]
    mcmc.autocorr=autocorr(chain_values, lags = lags, relative=TRUE)
    
    result=lags[which(mcmc.autocorr<0.1)[1]]
    return(result)
  },
  lags,column_name)
  
  thinning= max(unlist(list.autocorr))
  print(thinning)
  

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
  prior.values <-  rexp(nrow(allchainsdat), rate = lambda)
  
  mean.prior <- mean(prior.values)
  SE.prior   <- sd(prior.values) / sqrt(length(prior.values))
  
  c1 <- rgb(0,0,255,max = 255, alpha = 80, names = "blue1")
  c2 <- rgb(255,255,0, max = 255, alpha = 80, names = "yellow1")
  posterior_values= allchainsdat[[ grep(paste("^", column_name, sep=""), names(allchainsdat), value = TRUE)[1] ]] 

  mean.posterior <- mean(posterior_values)
  SE.posterior   <- sd(posterior_values) / sqrt(length(posterior_values))
 
  SE      <- sqrt( SE.prior^2 + mean.posterior^2) 
  z       <- (mean.prior - mean.posterior) / SE
  P.z     <- pnorm(z, lower.tail = TRUE)
  print(P.z)
  KStest=ks.test(prior.values, posterior_values,
          alternative = c("two.sided"))
  print(KStest)
  
  all.dat<-data.frame(prior=log(prior.values), posterior=log(posterior_values))
  data.melted<- melt(all.dat)
  log.density <- function(lambda,x) lambda*exp(x)*exp(-lambda*exp(x))
  grob1 <- grobTree(textGrob(paste("p-value Z test",signif(P.z, digits = 3),sep=" "), x=0.0,  y=0.3, hjust=0,
                            gp=gpar(col="black", fontsize=10, fontface="bold")))
  
  g1<-ggplot(data.melted,aes(x=value, fill=variable)) + geom_density(alpha=0.3)+
       scale_fill_manual( values=c('green','blue' ), labels=c(expression('Prior'),expression('Posterior')))+
        stat_function(fun=function(x) { log.density(lambda, x) } , geom='line',size = 1, aes(colour = "Theor. Prior") )+
          ggtitle(paste("Prior vs posterior(thin.:",  thinning, ")", sep = " " ))+ 
          labs(y="Density", x =paste("log ", column_name, sep=""))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=15),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=15),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.title = element_blank(),
          legend.position = c(.30,.70),
          legend.text.align = 0)+
           annotation_custom(grob1)
        #+
         # scale_colour_manual(values = c("red"), labels = c("Theor. Prior"))
 
  grob2 <- grobTree(textGrob(paste("p-value KS test", signif(KStest$p.value, digits = 3),sep=" "), x=0.0,  y=0.5, hjust=0,
                            gp=gpar(col="black", fontsize=10, fontface="bold")))
    
  g2<-ggplot(data.melted,aes(x=value,colour=variable)) + geom_step(aes(y=..y..),stat="ecdf")+     
     ggtitle(paste("Prior vs posterior(thin.:",  thinning, ")", sep = " " ))+ 
    labs(y="Cum Density", x =paste("log ", column_name, sep=""))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=15),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=15),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.title = element_blank(),
          legend.position = c(.30,.90),
          legend.text.align = 0)+
    scale_colour_manual( values=c('green','blue' ), labels=c(expression('Prior'),expression('Posterior')))+
    annotation_custom(grob2)
    
  g3<-ggplot(data.melted,aes(x=variable, y=value, fill=variable)) + geom_boxplot()+
    ggtitle(paste("Prior vs posterior(thin.: ",  thinning, ")", sep = " " ))+
    labs(y="", x =paste("log ", column_name, sep=""))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+
          scale_fill_manual( values=c('green','blue'), labels=c(expression('Prior'),expression('Posterior')))
  
  ggpubr::ggarrange(g1, g2, g3, 
            labels = c("A", "B", "C"),
            ncol = 2, nrow = 2)
  ggsave(paste(path_to_save, "/", column_name, ".pdf", sep=""))
  ggsave(paste(path_to_save, "/",column_name, ".png", sep=""))

}
plotPosterior=function(chainPath_List,number.iterations,percent_burn_in,values_to_compare,  column_name, path_to_save)
{
  require(data.table)
  require(dplyr)
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
  require(grid)
  require(coda)
  list_chains_data<- lapply(chainPath_List, FUN=function(x)
    as.data.frame(data.table::fread(x, header=T,  sep = "\t", colClasses=list(numeric=1:12))))
  list_chains_data<- lapply(list_chains_data, FUN = function(x)
    x[, colSums(is.na(x)) != nrow(x)])
  
  filtered_list_chains_data= list_chains_data
  
  mcmc_list=build_mcmc_list(list_chains_data,number.iterations)
  
  #gelman.diag.output<-gelman.diag(mcmc_list, confidence = 0.95,multivariate=FALSE)
  #print(gelman.diag.output)
  
  lags = c(100, 200, 500,1000)
  list.autocorr= lapply(mcmc_list, FUN=function(x, lags, column_name) 
  {
    column.index=which(names(as.data.frame(x)) == grep(paste("^", column_name, sep=""), names(as.data.frame(x)), value = TRUE)[1] )
    chain_values= x[,column.index]
    mcmc.autocorr=autocorr(chain_values, lags = lags, relative=TRUE)
    
    result=lags[which(mcmc.autocorr<0.1)[1]]
    return(result)
  },
  lags,column_name)
  
  thinning= max(unlist(list.autocorr))
  print(thinning)
  
  
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

  
  c1 <- rgb(0,0,255,max = 255, alpha = 80, names = "blue1")
  c2 <- rgb(255,255,0, max = 255, alpha = 80, names = "yellow1")
  posterior_values= allchainsdat[[ grep(paste("^", column_name, sep=""), names(allchainsdat), value = TRUE)[1] ]] 
  
 
  
  if (!is.null(values_to_compare) && length(values_to_compare)>0)
    {
    all.dat<-data.frame( Posterior=log(posterior_values), PosteriorR=log(values_to_compare))
     color.values=c('blue', 'green')
     distrib.labels=c(expression('Posterior'),expression('PosteriorR'))
     
     mean.prior <- mean(values_to_compare)
     SE.prior   <- sd(values_to_compare) / sqrt(length(values_to_compare))
     mean.posterior <- mean(posterior_values)
     SE.posterior   <- sd(posterior_values) / sqrt(length(posterior_values))
     
     SE      <- sqrt( SE.prior^2 + mean.posterior^2) 
     z       <- (mean.prior - mean.posterior) / SE
     P.z     <- pnorm(z, lower.tail = TRUE)
     print(P.z)
     KStest=ks.test(values_to_compare, posterior_values,
                    alternative = c("two.sided"))
     print(KStest)
     grob1 <- grobTree(textGrob(paste("p-value Z test",signif(P.z, digits = 3),sep=" "), x=0.0,  y=0.3, hjust=0,
                                gp=gpar(col="black", fontsize=10, fontface="bold")))
     
     grob2 <- grobTree(textGrob(paste("p-value KS test", signif(KStest$p.value, digits = 3),sep=" "), x=0.0,  y=0.5, hjust=0,
                                gp=gpar(col="black", fontsize=10, fontface="bold")))
  }
  else{
    
    all.dat<-data.frame( posterior=log(posterior_values), log(values_to_compare))
    color.values=c('blue')
    distrib.labels=c(expression('Posterior'))
    grob1 <- grobTree(textGrob(""))
    
    grob2 <- grobTree(textGrob(""))
  }
  data.melted<- melt(all.dat)

  
  g1<-ggplot(data.melted,aes(x=value, fill=variable)) + geom_density(alpha=0.3)+
    scale_fill_manual( values=color.values, labels=distrib.labels)+
    ggtitle(paste("Posterior(thin.:",  thinning, ")", sep = " " ))+ 
    labs(y="Density", x =paste("log ", column_name, sep=""))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=15),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=15),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.title = element_blank(),
          legend.position = "bottom",
          legend.text.align = 0)+
          annotation_custom(grob1)
 
  
  
  g2<-ggplot(data.melted,aes(x=value,colour=variable)) + geom_step(aes(y=..y..),stat="ecdf")+     
    ggtitle(paste("Posterior(thin.:",  thinning, ")", sep = " " ))+ 
    labs(y="Cum Density", x =paste("log ", column_name, sep=""))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=15),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=15),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.title = element_blank(),
          legend.position ="bottom",
          legend.text.align = 0)+
         annotation_custom(grob2)+
    scale_colour_manual( values=color.values, labels=distrib.labels)
  
  g3<-ggplot(data.melted,aes(x=variable, y=value, fill=variable)) + geom_boxplot()+
    ggtitle(paste("Posterior(thin.: ",  thinning, ")", sep = " " ))+
    labs(y="", x =paste("log ", column_name, sep=""))+
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+
    scale_fill_manual( values=color.values, labels=distrib.labels)
  
  ggpubr::ggarrange(g1, g2, g3, 
                    labels = c("A", "B", "C"),
                    ncol = 2, nrow = 2)
  ggsave(paste(path_to_save, "/", column_name, ".pdf", sep=""))
  ggsave(paste(path_to_save, "/",column_name, ".png", sep=""))
  
}
build_mcmc_list<-function(filtered_list_chains_data, number.iterations)
  {
  require(coda)
  mcmc_list<-lapply(filtered_list_chains_data, FUN=function(x)
    {
    res<-mcmc(data= x, start = 1, end = number.iterations, thin = 1)
    return(res)
   }
   )
  mcmc_list=as.mcmc.list(mcmc_list)
  return(mcmc_list)
}
get.thinning=function(chain_data,number.iterations, lags, threshold_corr){
  if (length(lags)>0 && threshold_corr>0 && threshold_corr <1)
  {
    mcmc_data<-mcmc(data= chain_data, start = 1, end = number.iterations, thin = 1)
    mcmc.autocorr=autocorr(mcmcTime, lags = lags, relative=TRUE)
    
    thinning=lags[which(mcmc.autocorr<threshold_corr)[1]]
    return(thinning)
  }
  
}
thin.chain.data=function(chain_data,number.iterations, percent_burn_in, thinning){
  if (length(chain_data)>0 && thinning >1 && percent_burn_in >=0 && percent_burn_in< 1 )
  {
    mcmc_data<-mcmc(data= chain_data, start = 1, end = number.iterations, thin = 1)
    number_interations_burn_in= floor(percent_burn_in * length(chain_data))
    chain_data_after_burn_in=chain_data[(number_interations_burn_in+1):length(chain_data)]
    indexes<-seq(0,length(chain_data_after_burn_in),thinning)
    thinned_data=chain_data_after_burn_in[indexes]
    return(thinned_data)
  }
  
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
