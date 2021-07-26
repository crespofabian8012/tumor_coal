library(ape)
library(phytools)
library(tidytree)
#library(dplyr)
library(ggtree)
library(adephylo)
path= "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/benchmarking"
setwd(path)
path0.1= "ResultsDelta=0.1T=5/trees_dir/"
path1= "ResultsDelta=1T=5/trees_dir/"
path2= "ResultsDelta=2T=5/trees_dir/"
path3= "ResultsDelta=3T=5/trees_dir/"
path4= "ResultsDelta=4T=5/trees_dir/"
path5= "ResultsDelta=5T=5/trees_dir/"
path10= "ResultsDelta=10T=5/trees_dir/"
path11= "ResultsDelta=0.1T=0.5/trees_dir/"
path12= "ResultsDeltaT=100T=0.5/"

pathslist=list(  path12)
listFullPaths=unlist(lapply(pathslist, FUN=function(x, path) paste(path,"/" ,x, sep=""), path))
#lapply(listFullPaths, FUN=function(x) system(paste("cd ", x,"; cat trees_0* >> allTrees.tree;")))#for the trees in model time
#lapply(listFullPaths, FUN=function(x) system(paste("cd ", x," ; cat ", "trees_2", "*","  >> allTreesPhysicalTime.tre ",sep='', intern = T)))

#lapply(listFullPaths, FUN=function(x) paste("cd ", x," ; cat ", "trees_2", "*"," >> allTreesPhysicalTime.tre ",sep=''))


phylo.treeClosest.to = NULL
current.tree = NULL
mutation.rate=10^-7


saveMostRepresentativeTree<-function(path, AllTreesFileName, Delta,mutation.rate )
{
  sumMRCAPhysicalTime=0.0;
  bestTree= NULL
  list.Physical.timeMRCA = list()
  MinDistance = 1000
  allTrees <- ape::read.tree(file=paste(path, "/",AllTreesFileName, sep=''))
  lapply(1:length(allTrees), function(i, allTrees, sumMRCAPhysicalTime, list.Physical.timeMRCA)
  {
    current.tree = allTrees[[i]]
    current.tree<-ape::as.phylo(current.tree)
    x <- as_tibble(current.tree)
    tips.labels<-current.tree$tip.label
    tumor.tip.labels= tips.labels[-which(tips.labels=="healthycell")]
    mrca= getMRCA(current.tree, tumor.tip.labels)
    scaledtimeRoot= distRoot(current.tree, tips = tumor.tip.labels, method = "patristic")
    scaledtimeMRCA =scaledtimeRoot- x[mrca,]$branch.length
    physicalTimeMRCA = scaledtimeMRCA / mutation.rate
    
    stopifnot(sum(physicalTimeMRCA - rep(mean(physicalTimeMRCA), length(physicalTimeMRCA)))==0)
    
    avgPhysicalTimeMRCA=mean(physicalTimeMRCA)
    list.Physical.timeMRCA[i] <<- avgPhysicalTimeMRCA 
    sumMRCAPhysicalTime <<- sumMRCAPhysicalTime + avgPhysicalTimeMRCA
    
  },allTrees=allTrees, sumMRCAPhysicalTime= sumMRCAPhysicalTime, list.Physical.timeMRCA= list.Physical.timeMRCA
  )
  meanMRCAPhysicalTime= sum(unlist(list.Physical.timeMRCA)) / length(list.Physical.timeMRCA)
  distances= abs(unlist(list.Physical.timeMRCA) - meanMRCAPhysicalTime)
  minDistance= min(distances)
  print(minDistance)
  indexMinDistance= which(distances == minDistance)
  
  bestTree=allTrees[[indexMinDistance[1]]]
  
  return(bestTree)
}

library(ape)
 allTrees <- "allTreesPhysicalTime.tre"
Delta.list=c(100)
lapply( 1:length(listFullPaths), FUN= function(i,listFullPaths, Delta.list, allTrees, mutation.rate ){
  path=listFullPaths[i]
  print(path)
  delta=Delta.list[i]
  bestTree= saveMostRepresentativeTree(path, allTrees, delta,mutation.rate )
  write.tree(bestTree,file=paste(path, "/bestTreePhysicalTime.tre", sep=""))
},
listFullPaths=listFullPaths, Delta.list=Delta.list, allTrees=allTrees, mutation.rate=mutation.rate
)
