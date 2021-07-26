#@dposada 27/11/18
# I use this to calculate the RF distance of the true (generating) tree and the ML tree 
# for benchmarking CellCoal
# trees are in Newick format; alignments are in Phylip format

library(ape)
library(phangorn)
library(phytools)

#prepare shortcurts
true_tree_file = "/Users/dposada/Coding/C/*Simulation\ tools/CellCoal\ project/results/trees_dir/trees.0001"
#here 
ML_tree_file = "/Users/dposada/Coding/C/*Simulation\ tools/CellCoal\ project/results/snv_haplotypes_dir/snv_hap.0001.treefile"
snv_hap_file_comm = "/Users/dposada/Coding/C/*Simulation\\ tools/CellCoal\\ project/results/snv_haplotypes_dir/snv_hap.0001"
iqtree_bin_comm = "/Users/dposada/Coding/C/*Simulation\\ tools/CellCoal\\ project/stuff/iqtree/iqtree"

#get rid of SNV position line in the CellCoal snv haplotype file
system (paste("sed -i '' 2d", snv_hap_file_comm))

#run iqtree to generate an ML tree from the snv ML haplotypes
system (paste(iqtree_bin_comm, "-s", snv_hap_file_comm , "-mredo -redo -nt AUTO -m TESTONLY -mrate E"))

#read the trees
true_tree <- read.tree(file=true_tree_file)
ML_tree <- read.tree(file=ML_tree_file)

#reroot at healthy cell
ML_treeR <-root(ML_tree, "healthycell")

#get rid of the healthy tip just in case
true_treeT<-drop.tip(true_tree,"healthycell")
ML_treeT<-drop.tip(ML_treeR,"healthycell")

#calculate distances among true tree and ML tree
treedist(true_treeT, ML_treeT, check.labels = TRUE)

#mirror true and ML trees
co<-cophylo(true_treeT, ML_treeT)
plot(co,fsize=0.8)
