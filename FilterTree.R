library(ape)
library(phytools)
library(tidytree)
library(dplyr)
library(ggtree)
library(adephylo)
ModelTimeMRCA=c(1.303748, 0.989492, 0.286708,  0.051329, 0.007415  )
path= "/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/xcode/simul/Debug/ResultsDelta=0.1T=5/trees_dir/"
allTrees <- ape::read.tree(file="/Users/faustofabiancrespofernandez/Downloads/tumor_coal_last_last/tumor_coal/xcode/simul/Debug/ResultsDelta=0.1T=5/trees_dir/alltreesModelTime.txt")

phylo.treeClosest.to = NULL
current.tree = NULL
mutation.rate=10^-7
meanMRCAmodelTime= ModelTimeMRCA[1]
MinDistance = 1000
bestTree= NULL
for (i in 1:length(allTrees))
{
  str(i)
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
  distance= abs(avgtimeMRCA - meanMRCAmodelTime)
  if (distance <MinDistance)
    {
     
    MinDistance = distance 
    bestTree=current.tree
   }
}

 str(bestTree) 
 x <- as_tibble(bestTree)
 tips.labels<-bestTree$tip.label
 tumor.tip.labels= tips.labels[-which(tips.labels=="healthycell")]
 mrca= getMRCA(bestTree, tumor.tip.labels)
 scaledtimeRoot= distRoot(bestTree, tips = tumor.tip.labels, method = "patristic")
 scaledtimeMRCA =scaledtimeRoot- x[mrca,]$branch.length
 timeMRCA = scaledtimeMRCA / mutation.rate
 avgtimeMRCA=mean(timeMRCA)
 distance= abs(avgtimeMRCA - meanMRCAmodelTime)
 
 print(distance)
 print(MinDistance)
 stopifnot(bestTree==MinDistance)
 write.tree(bestTree,file=paste(path, "bestTree.tre"))
 
