library(ape)
library(phytools)
library(tidytree)
#library(dplyr)
library(ggtree)
library(adephylo)
meanModelTimeMRCA=c(1.303748, 0.989492, 0.286708,  0.051329, 0.007415  )
path= "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/benchmarking"
setwd(path)
path0.1= "ResultsDelta=0.1T=5/trees_dir"
path1= "ResultsDelta=1T=5/trees_dir"
path2= "ResultsDelta=2T=5/trees_dir"
path3= "ResultsDelta=3T=5/trees_dir"
path4= "ResultsDelta=4T=5/trees_dir"
path5= "ResultsDelta=5T=5/trees_dir"
path10= "ResultsDelta=10T=5/trees_dir"

pathslist=list(path0.1,path1, path2, path3, path4, path5, path10)
listFullPaths=lapply(pathslist, FUN=function(x, path) paste(path,"/" ,x, sep=""), path)
lapply(listFullPaths, FUN=function(x) system(paste("cd ", x,"; cat trees_0* >> allTrees.tree;")))



phylo.treeClosest.to = NULL
current.tree = NULL
mutation.rate=10^-7
meanMRCAmodelTime= meanModelTimeMRCA[1]


saveMostRepresentativeTree<-function(path, AllTreesFileName, Delta,mutation.rate )
{
  require(ape)
  
  sumMRCAModeTime=0.0;
  bestTree= NULL
  list.timeMRCA = list()
  MinDistance = 1000
  allTrees <- ape::read.tree(file=paste(path, "/",AllTreesFileName, sep=''))
  for (i in 1:length(allTrees))
  {
    current.tree = allTrees[[i]]
    current.tree<-as.phylo(current.tree)
    x <- as_tibble(current.tree)
    tips.labels<-current.tree$tip.label
    tumor.tip.labels= tips.labels[-which(tips.labels=="healthycell")]
    mrca= getMRCA(current.tree, tumor.tip.labels)
    scaledtimeRoot= distRoot(current.tree, tips = tumor.tip.labels, method = "patristic")
    scaledtimeMRCA =scaledtimeRoot- x[mrca,]$branch.length
    timeMRCA = scaledtimeMRCA / mutation.rate
    
    avgtimeMRCA=mean(timeMRCA)
    list.timeMRCA[i]=avgtimeMRCA 
    sumMRCAModeTime = sumMRCAModeTime + avgtimeMRCA
    
  }
  meanMRCAmodelTime= sumMRCAModeTime / length(allTrees)
  distances= abs(unlist(list.timeMRCA) - meanMRCAmodelTime)
  minDistance= min(distances)
  print(minDistance)
  indexMinDistance= which(distances == minDistance)
  
  bestTree=allTrees[[indexMinDistance[1]]]
  
  return(bestTree)
}


allTrees <- "allTrees.tree"
Delta.list=c(0.1,1,2,3,4,5,10)
for (i in 1:length(listFullPaths)){
  path=listFullPaths[[i]]
  delta=Delta.list[i]
  bestTree= saveMostRepresentativeTree(path, allTrees, delta,mutation.rate )
  write.tree(bestTree,file=paste(path, "/bestTree2.tre", sep=''))
}

