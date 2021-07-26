rm(list=ls())


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
analyzeChains<-function(path.Results, path_to_save,  number.iterations, column_name,
                        column_name_true_values,
                        percent_burn_in,
                        prior_mean,
                        dorethinning,
                        includeHiperPriorMean){
  
  all.files.paths <- list.files(path = path.Results, pattern = "\\.log$",recursive = TRUE, full.names=TRUE)
  
  list_chains_data<-  lapply(all.files.paths, FUN=function(x)
    as.data.frame(data.table::fread(x,  sep = "\t")))
  
  list_chains_data<- lapply(list_chains_data, FUN = function(x)
    x[, colSums(is.na(x)) != nrow(x)])
  
  list_chains_data= list_chains_data[sapply(list_chains_data, function(x) dim(x)[1]) > 0]
  
  
  first_4_lines<-lapply(all.files.paths, FUN=function(x)
    scan(file= x, what = "",  blank.lines.skip
         =TRUE, nlines = 4)[-1]
  )
  
  true_values_string =lapply(first_4_lines, FUN=function(x) {
    values= grep(paste("^", column_name_true_values,"[\\(0-9a-z\\)]*=", sep=""), x, value = TRUE)[1]  
    return(values[1])
  }
  )  
  
  true_values_first_line =lapply(true_values_string, FUN=function(x) {
    as.numeric(gsub("[^[:digit:].]", "",  x))
  }
  )  
  
 
  
  filtered_list_chains_data= list_chains_data
  
  #mcmc_list=build_mcmc_list2(list_chains_data,number.iterations)
  
  lags = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  list.autocorr= lapply(list_chains_data, FUN=function(x, lags, column_name, dorethinning) 
  {
    require(coda)
    res<-mcmc(data= x, start = x$state[1], end = x$state[nrow(x)], thin =  x$state[2]- x$state[1])
    column.index=which(names(as.data.frame(res)) == grep(paste("^", column_name, sep=""), names(as.data.frame(res)), value = TRUE)[1] )
    chain_values= res[,column.index]
    mcmc.autocorr=autocorr(chain_values, lags = lags, relative=TRUE)
    print(mcmc.autocorr) 
    if  (dorethinning){
      if (length(which(mcmc.autocorr<0.1))>1)
        {
        result=lags[which(mcmc.autocorr<0.1)[1]]* (x$state[2]- x$state[1])
        }  
      else
        {
        
          result=lags[length(lags)]*(x$state[2]- x$state[1])
        }
        
    }
    else{
      result=x$state[2]- x$state[1]
    }
    if (is.na(result))
      {
      print("oops")
    }
    print(result)
    return(result)
  },
  lags,column_name, dorethinning )
  
  print(unlist(list.autocorr))
   
  thinning= max(unlist(list.autocorr))
     
   
  thinning= max(unlist(list.autocorr))
  print(thinning)
  
  
  filtered_list_chains_data2= filtered_list_chains_data
  # if(percent_burn_in>0 && percent_burn_in <1)
  # {
  #   filtered_list_chains_data2<- lapply(filtered_list_chains_data, FUN= function(x)
  #   {
  #     number_interations_burn_in= floor(percent_burn_in * nrow(x))
  #     return(x[x$state >number_interations_burn_in,])
  #   })
  # }
  
  list_thinned_chains_data=filtered_list_chains_data2
  
  if (dorethinning){
    print("rethinnig")
    list_thinned_chains_data<-lapply(filtered_list_chains_data2, FUN=function(x, thinning)
    {
      rethinning.length<-thinning /( x$state[2]- x$state[1])
      indexes<-seq(0,nrow(x),rethinning.length)
      print(thinning /( x$state[2]- x$state[1]))
      if (  nrow(as.data.frame(x)) / rethinning.length >=100)
       {
        return(x[indexes, ])
        
      }
      else{
        
        return(x)
      }
      
    },
    thinning)
    
  }

  
  allchainsdat<- do.call("rbind", list_thinned_chains_data)
  names(allchainsdat)<-unlist(lapply(names(allchainsdat), function(x) {
    x1= gsub("\\(", "_", x)
    gsub("\\)", "", x1)
  }
  ))
  
  print(nrow(allchainsdat))
  posterior_values= allchainsdat[[ grep(paste("^", column_name, sep=""), names(allchainsdat), value = TRUE)[1] ]] 
  
  
  print (quantile(posterior_values, probs = c(0.05, 0.5, 0.95)))
  overall_mean =mean(posterior_values)
  overall_LowerCI= quantile(posterior_values, probs = c(0.05))
  overall_UpperCI= quantile(posterior_values, probs = c(0.95))
  
  list_posterior_values=lapply(list_thinned_chains_data, FUN=function(x) 
    x[[ grep(paste("^", column_name, sep=""), names(x), value = TRUE)[1] ]] 
  )
  
  vector_means  = lapply(list_posterior_values, FUN=function(x) mean(x))
  vector_LowerCI  = lapply(list_posterior_values, FUN=function(x) quantile(x, probs = c(0.05)))
  vector_UpperCI  = lapply(list_posterior_values, FUN=function(x) quantile(x, probs = c(0.95)))
  
  library(ggplot2)
  library(magrittr)
  library(tidyr)
  library(dplyr)
  library(ggfan)
  

  df <- data.frame(x =1:length(vector_means),
                   F =unlist(vector_means),
                   true=unlist(true_values_first_line),
                   L =unlist(vector_LowerCI) ,
                   U =unlist(vector_UpperCI) )
  

  print("mean normalized MAP-true value")
  df$stat<- abs(df$F-df$true)/(df$U-df$L)
  print(mean(df$stat))
  print("count hits")
  library(dplyr)
  print(nrow(filter(df, (true >=L & U>=true))))
  print(nrow(filter(df, (true >=L & U>=true)))/nrow(df))


  
  true_values <- data.frame(x =1:length(vector_means),
                            true_value=unlist(true_values_first_line))
  library(grid)
  text_prior_mean <- textGrob("Prior\nmean", gp=gpar(fontsize=13, fontface="bold"))
  text_overall_posterior_mean <- textGrob("Combined posterior\nmean", gp=gpar(fontsize=13, fontface="bold"))
  text_overall_posterior_lowerCI <- textGrob("Combined posterior\nlower CI", gp=gpar(fontsize=13, fontface="bold"))
  text_overall_posterior_upperCI <- textGrob("Combined posterior\nupper CI", gp=gpar(fontsize=13, fontface="bold"))
  
  Prior_mean <- data.frame( x = c(-Inf, Inf), y = prior_mean, Prior_mean = factor(prior_mean) )
  #Combined_CI <- data.frame(yintercept = c(prior_mean,overall_LowerCI,overall_mean, overall_UpperCI), Lines = c("Prior mean","Combined posterior lowerCI", "Combined posterior mean", "Combined posterior upperCI"))
  Combined_CI <- data.frame(yintercept = c(prior_mean), Lines = c("Prior mean"))
  
  require(ggplot2)
  p<-ggplot(df, aes(x = x, y = F)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymax = U, ymin = L))+
    geom_point(data=true_values, mapping=aes(x = x, y=true_value ),size = 1, col ="red") +
  #  geom_hline(aes(yintercept = yintercept, linetype = Lines), Combined_CI)+
    ggtitle(paste("Credible intervals for ", column_name, "(thin.: ",  thinning, ")", sep = " " ))+
    labs(y="", x ="simulation")+
    theme(plot.margin = unit(c(1,1,2,1), "lines")) +
    theme(text = element_text(size=15),axis.text.x=element_text(angle=0,hjust=1,vjust=.9,size=10),
          axis.text.y=element_text(angle=0,hjust=1,vjust=.8,size=10),
          plot.title = element_text(size=15, face="bold", vjust=2),
          axis.title.x=element_text(size=10,face="bold",vjust=-0.5,hjust=0.5),
          axis.title.y=element_text(size=10,face="bold",vjust=1.5,hjust=0.5),
          legend.text.align = 0)+ theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  if (includeHiperPriorMean){
    p <- p +  geom_hline(aes(yintercept = yintercept, linetype = Lines), Combined_CI)
    
  }
    
  ggsave(paste(path_to_save, "/", column_name, ".pdf", sep=""))
  ggsave(paste(path_to_save, "/",column_name, ".png", sep=""))
}
#########################################################################################
path.Results = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/xcode/inference/Debug/Results"
path.Results = "/Volumes/LACIE_SHARE/ResultadosMeanPriorDelta=100/InitialValuesFirstLine/ResultsAll2"
path.Results = "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/UsingData/DeltaT=100T=0.5/Results"
path_to_save="/Volumes/LACIE_SHARE/ResultadosMeanPriorDelta=100/InitialValuesFirstLine/"
path.Results ="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/UsingData/DeltaT=100T=0.5/Results"
path.Results ="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/UsingData/DeltaT=1n=20/Results3"
source("/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/ComparePriorPosterior.R")

path_to_save= path.Results
number.iterations =500000
column_name_true_values="deltaT"
column_name="deltaT"
percent_burn_in=0.0
prior_mean=1
dorethinning=FALSE


analyzeChains(path.Results, path_to_save,  number.iterations, column_name,
                        column_name_true_values,
                        percent_burn_in,
                         prior_mean, TRUE, TRUE )

column_name_true_values="theta"
column_name="theta"
prior_mean=0.01
analyzeChains(path.Results, path_to_save,  number.iterations, column_name,
              column_name_true_values,
              percent_burn_in,
              prior_mean, dorethinning=TRUE, TRUE )


column_name_true_values="TOriginSTD"
column_name="time_origin_std"

analyzeChains(path.Results, path_to_save,  number.iterations, column_name,
              column_name_true_values,
              percent_burn_in,
              prior_mean, dorethinning=TRUE, FALSE )
##########################################################
path.Results ="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking/UsingData/DeltaT=100n=20/Results"
path_to_save= path.Results
number.iterations =500000
percent_burn_in=0.0
dorethinning=FALSE

column_name_true_values="deltaT"
column_name="deltaT"
prior_mean=100
analyzeChains(path.Results, path_to_save,  number.iterations, column_name,
              column_name_true_values,
              percent_burn_in,
              prior_mean, TRUE, TRUE )

column_name_true_values="theta"
column_name="theta"
prior_mean=0.01
analyzeChains(path.Results, path_to_save,  number.iterations, column_name,
              column_name_true_values,
              percent_burn_in,
              prior_mean, dorethinning=TRUE, TRUE )

column_name_true_values="TOriginSTD"
column_name="time_origin_std"
column_name_true_values="time_origin_std"
analyzeChains(path.Results, path_to_save,  number.iterations, column_name,
              column_name_true_values,
              percent_burn_in,
              prior_mean, dorethinning=TRUE, FALSE )
# 
# N_time <- 50
# N_sims <- 1000 
# time <- 1:N_time
# mu <- time**2 * 0.03 + time * 0.3
# sds <- exp(time**2 * -0.001 + time * 0.1)
# 
# # simulate 1000 samples from each time point
# fake_data <- sapply(time, function(i) rnorm(N_sims, mu[i], sds[i]))
# 
# # gather into a long-form, tidy dataset
# fake_df <- data.frame(x=time, t(fake_data)) %>% gather(key=Sim, value=y, -x)
# head(fake_df)
# p <- ggplot(fake_df, aes(x=x,y=y)) + geom_interval()
# print(p)