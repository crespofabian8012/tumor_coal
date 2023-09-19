


rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

list_packages <-
  c(
    "phangorn",
    #  "rlang",
    "parallel",
    "ape",
    "tibble",
    "foreach",
    "doParallel",
    "stats",
    "tcltk",
    "TreeTools",
    "dplyr",
    "stringr",
    "tidyr",
    "tidytree",
    "castor",
    "phylotools",
    "NLP",
    "doFuture",
    "foreach",
    "phytools",
    "dcGOR"
    
  )

if (!require("Rcpp", quietly = TRUE))
  install.packages("Rcpp", repos="https://rcppcore.github.io/drat",INSTALL_opts = '--no-lock', type="source")
#if (!require("pillar", quietly = TRUE))

if (!require("devtools", quietly = TRUE))
  devtools::install_github("hadley/devtools")
 
  
library(devtools)
if (!require("pillar", quietly = TRUE)){
  
  devtools::install_version("pillar", version = "1.8.1", upgrade= c("always"),repos = "http://cran.us.r-project.org",dependencies=TRUE)
}

if (!require("ape", quietly = TRUE))
 devtools::install_version("ape", version = "5.4-1", upgrade= c("always"),repos = "http://cran.us.r-project.org",dependencies=TRUE)



if (!require("TreeTools", quietly = TRUE))
  devtools::install_version("TreeTools", version = "1.4.2", upgrade= c("always"),repos = "http://cran.us.r-project.org",dependencies=TRUE)

if (!require("phytools", quietly = TRUE))
  devtools::install_version("phytools", version = "0.7-70", upgrade= c("always"),repos = "http://cran.us.r-project.org",dependencies=TRUE)


library("pillar", lib.loc="~/dev/pillar/v1.8.1")
#library("pillar")

if (!require("ellipsis", quietly = TRUE))
  install.packages("ellipsis", dependencies=TRUE, repos = "http://cran.us.r-project.org")


if (!require("vctrs", quietly = TRUE))
  devtools::install_version("vctrs", version = "0.5.2", upgrade= c("always"),repos = "http://cran.us.r-project.org",dependencies=TRUE)
  #install.packages("vctrs", dependencies=TRUE,INSTALL_opts = '--no-lock', repos = "http://cran.us.r-project.org")



install_required_packages = function(list.of.packages) {
  new.packages <-
    list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages) > 0)
    install.packages(new.packages, dependencies = TRUE, INSTALL_opts = '--no-lock', repos = "http://cran.us.r-project.org")
  
}


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", dependencies = TRUE,INSTALL_opts = '--no-lock', repos = "http://cran.us.r-project.org")

#if (!require("ggtree", quietly = TRUE))
#  BiocManager::install("ggtree")

if (!require("supraHex", quietly = TRUE))
  BiocManager::install("supraHex", ask=FALSE)

if (!require("Rgraphviz", quietly = TRUE))
  BiocManager::install("Rgraphviz",  ask=FALSE)



#install_required_packages(list_packages)

#if (!require("rstan", quietly = TRUE))
#  install.packages("rstan", dependencies=TRUE, INSTALL_opts = '--no-lock', repos = "http://cran.us.r-project.org")

#if (!require("ellipsis", quietly = TRUE))
#  install.packages("ellipsis", dependencies=TRUE, repos = "http://cran.us.r-project.org")


#if (!require("vctrs", quietly = TRUE))
#  install.packages("vctrs", dependencies=TRUE,INSTALL_opts = '--no-lock', repos = "http://cran.us.r-project.org")



if (!require("posterior", quietly = TRUE))
  install.packages("posterior", INSTALL_opts = '--no-lock', dependencies = TRUE)

if (!require("cmdstanr", quietly = TRUE))
   install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")),  dependencies = TRUE)

lapply(list_packages,
       library,
       character.only = TRUE)

#library(rstan)
library(posterior)


get_current_file_location <-  function()
{
  require(tibble)
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(
      col = value,
      into = c("key", "value"),
      sep = "=",
      fill = 'right'
    ) %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file) == 0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}
path = get_current_file_location()

IUPAC_to_numeric<-function(IUPAC_code){
  upper_case_code<- toupper(IUPAC_code)
  res = switch(  
    upper_case_code,  
    "A"= 1,  
    "C"= 2,  
    "G"= 3,  
    "T"= 4,
    "M"= 5,
    "R"= 6,
    "W"= 7,
    "S"= 8,
    "Y"= 9,
    "K"= 10,
    "M"= 11,
    "R"= 12,
    "W"= 13,
    "S"= 14,
    "Y"= 15,
    "K"= 16,
    "N"= 0,
    "?"= 0
  ); 
  
  return(res)
}
from_list_IUPAC_to_numeric_codes<-function(IUPACsequence){
  list<-strsplit(IUPACsequence, "")[[1]]
  result<-lapply(1:length(list), FUN= function(i){
    res = IUPAC_to_numeric(list[i])
    return(res);
  })
  return(unlist(result));
}
list_IUPAC_codes<-function(IUPACsequence){
  list<-strsplit(IUPACsequence, "")[[1]]
  result<-lapply(1:length(list), FUN= function(i){
     
    return(list[i]);
  })
  return(unlist(result));
}



#########################################################################
#########################################################################
print(args)
if (length(args)==0) {
  stop("At least supply the start index of the files to run", call.=FALSE)
} 

start_idx = as.integer(args[1])

print("start_idx=")
print(start_idx)



#########################################################################
#########################################################################
K = 0.8
path_to_save = get_current_file_location()


healthycell = "healthycell"

tip_pattern<-"*_C1_*"

path_all_data_folders="/data"


path_all_trees<-paste0(path_to_save,path_all_data_folders,"/trees_dir")
path_all_genotypes<-paste0(path_to_save,path_all_data_folders,"/snv_haplotypes_dir")
path_all_IUPACgenotypes<-paste0(path_to_save,path_all_data_folders,"/true_haplotypes_dir")#IUPAC codes
path_all_IUPACgenotypes<-path_all_genotypes
path_all_genotypes<-paste0(path_to_save,path_all_data_folders,"/full_genotypes_dir")

print(paste("path_all_trees=", path_all_trees))
print(paste("path_all_geenotypes=", path_all_genotypes))

path_all_genotypes<-paste0(path_to_save,path_all_data_folders,"/snv_genotypes_dir")




list_genotypes_filenames<-list.files(path_all_genotypes, all.files=FALSE,
                                     pattern = "*\\.[phy txt]*",
                                     full.names=TRUE)

list_IUPACgenotypes_filenames<-list.files(path_all_IUPACgenotypes, all.files=FALSE,
                                          pattern = "*\\.[phy txt]",
                                          full.names=TRUE)

list_IUPACgenotypes_only_filenames<-list.files(path_all_IUPACgenotypes, all.files=FALSE,
                                          pattern = "*\\.[phy txt]",
                                          full.names=FALSE)

list_genotypes_data<-lapply(list_genotypes_filenames, FUN=function(filename){
  phylotools::read.phylip(filename)
})

if (length(list_IUPACgenotypes_filenames)==0){
  
  stop("Error: The number of IUPAC files is 0!")
}

list_IUPACgenotypes_data<-lapply(list_IUPACgenotypes_filenames, FUN=function(filename){
  phylotools::read.phylip(filename)
})


list_IUPAC_numeric_genotypes_data<- lapply(list_IUPACgenotypes_data, FUN=function(list_data){
  list_data<-list_data[list_data$seq.name !=healthycell, ]
  ras<-list_data$seq.text
  list_numeric_seq<-lapply(1:(length(ras)-1), FUN=function(j){
    seq<-NLP::as.String(ras[j])
    from_list_IUPAC_to_numeric_codes(seq)
  })
  matrix(unlist(list_numeric_seq), nrow=(length(ras)-1),ncol=length(list_numeric_seq[[1]]),  byrow = T)
})



list_phyDat<- lapply(list_IUPACgenotypes_data, FUN=function(list_data,  healthycell, tip_pattern){
  ras<-list_data$seq.text
  names<-list_data$seq.name
  idxs=grepl(tip_pattern, names)
  thealthycell = gsub('-', '_', healthycell)
  
  idxs[names==thealthycell]=FALSE
  filtered<-ras[idxs]
  filtered_names=names[idxs]
  min_number_loci = min(unlist(lapply(filtered, FUN= function(p) nchar(p))))
  
  list_numeric_seq<-lapply(1:(length(filtered)), FUN=function(j){
    seq<-NLP::as.String(filtered[j])
    list_IUPAC_codes(seq)
  })
 
  res<-phangorn::phyDat(matrix(unlist(list_numeric_seq), nrow=(length(filtered)),ncol=length(list_numeric_seq[[1]]),  byrow = TRUE, 
                                 dimnames= 
                                   list(filtered_names,  1:length(list_numeric_seq[[1]]))), 
                              levels=c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "N"), compress=TRUE) 
  res
},  healthycell, tip_pattern)

#TODO :remove duplicate sequences
#mapping <- map_duplicates(dat)

list_IUPAC_compressed_phyDat<-lapply(1:length(list_phyDat), FUN=function(i){
  phy_dat<- list_phyDat[[i]]
  compress_dat<-subset(phy_dat, 1:length(phy_dat), compress=TRUE, site.pattern = TRUE)
  weight= attr(phy_dat, "weight")
  nr <- as.integer(attr(phy_dat, "nr"))
  nc <-as.integer(attr(phy_dat, "nc"))
  levels <- attr(phy_dat, "levels")
  index_all= attr(phy_dat, "index")
  index_unique_pos = which(!duplicated(index_all))
  phyDat_unique_patterns<- subset(phy_dat, select=index_unique_pos, site.pattern = FALSE) 
  
  char_matrix<-as.character(subset(phyDat_unique_patterns, select=1:nr, site.pattern = FALSE))
  rownames(char_matrix)=c()
  numeric_IUPAC_matrix = apply(char_matrix, c(1,2), IUPAC_to_numeric)
  
  numeric_IUPAC_matrix
})

list_site_weights_compressed_phyDat<-lapply(1:length(list_phyDat), FUN=function(i){
  phy_dat<- list_phyDat[[i]]
  compress_dat<-subset(phy_dat, 1:length(phy_dat), compress=TRUE, site.pattern = TRUE)
  weight= attr(phy_dat, "weight")
  weight
})


 
order_of_seq<- lapply(list_IUPACgenotypes_data, FUN=function(list_data, healthycell, tip_pattern){
  
  names<-list_data$seq.name
  filtered<-names[grepl(tip_pattern, names)]
  thealthycell = gsub('-', '_', healthycell)
  filtered = filtered[filtered != thealthycell]
  df=data.frame(label =filtered, order_id= 1:length(filtered))
  df
}, healthycell, tip_pattern)




true_deltas<-unlist(lapply(1:length(list_IUPACgenotypes_filenames), FUN=function(i, list_IUPACgenotypes_filenames){
  res <- str_match(list_IUPACgenotypes_filenames[[i]], "[Gamma|Delta|G]=\\s*(.*?)\\s*_")
  as.numeric(res[,2])
}, list_IUPACgenotypes_filenames))

true_torigins<-unlist(lapply(1:length(list_IUPACgenotypes_filenames), FUN=function(i, list_IUPACgenotypes_filenames){
  res <- str_match(list_IUPACgenotypes_filenames[[i]], "T=\\s*(.*?)\\s*_")
  as.numeric(res[,2])
}, list_IUPACgenotypes_filenames))

true_sample_sizes<-unlist(lapply(1:length(list_IUPACgenotypes_filenames), FUN=function(i, list_IUPACgenotypes_filenames){
  res <- str_match(list_IUPACgenotypes_filenames[[i]], "_n=\\s*(.*?)\\s*_")
  as.numeric(res[,2])
}, list_IUPACgenotypes_filenames))

list_number_var_sites<-unlist(lapply(1:length(list_IUPACgenotypes_filenames), FUN=function(i, list_IUPACgenotypes_filenames){
  res <- str_match(list_IUPACgenotypes_filenames[[i]], "_nS=\\s*(.*?)\\s*_")
  as.numeric(res[,2])
}, list_IUPACgenotypes_filenames))

true_thetas<-unlist(lapply(1:length(list_IUPACgenotypes_filenames), FUN=function(i, list_IUPACgenotypes_filenames){
  res <- str_match(list_IUPACgenotypes_filenames[[i]], "theta=\\s*(.*?)\\s*_")
  as.numeric(res[,2])
}, list_IUPACgenotypes_filenames))


true_nMUs<-unlist(lapply(1:length(list_IUPACgenotypes_filenames), FUN=function(i, list_IUPACgenotypes_filenames){
  res <- str_match(list_IUPACgenotypes_filenames[[i]], "nMU=\\s*(.*?)\\s*_")
  as.numeric(res[,2])
}, list_IUPACgenotypes_filenames))
######################################################################
#ML tree topologies


ML_tree_filenames <- list.files(paste0(path_to_save,path_all_data_folders, "/Cellphy_ML_trees/"), pattern="*raxml.bestTree", full.names=TRUE, all.files=FALSE)

ML_tree_only_filenames <- list.files(paste0(path_to_save,path_all_data_folders, "/Cellphy_ML_trees/"), pattern="*raxml.bestTree", full.names=FALSE, all.files=FALSE)

if (length(ML_tree_filenames)==0){

  stop("Error: The number of Cellphy trees is 0!")
}

list_ML_trees<-lapply(ML_tree_filenames, FUN=function(x) {
  castor::read_tree(file=x)
})


list_ML_trees<-lapply(list_ML_trees, FUN=function(x) {
  root(x, outgroup=healthycell)
})





tip_pattern<-"*_C1_*"
list_ML_all_heights_pop1<- lapply(list_ML_trees, FUN=function(tree, tip_pattern){
  tree_without_healthy<- drop.tip(tree, healthycell, rooted=TRUE)
  
  list_tips<-tree_without_healthy$tip.label
  list_tips_pop<-list_tips[!grepl(tip_pattern, list_tips)]
  idx_mrca<-ape::getMRCA(tree_without_healthy, list_tips_pop)
  sub_tree<-drop.tip(tree_without_healthy, list_tips_pop)
  number.tips<-length(sub_tree$tip.label)
  list_depths<-node.depth.edgelength(sub_tree)
  max_depth= max(list_depths)
  list_heights_internal_nodes<-rep(max_depth,
                                   length(list_depths)-number.tips)-list_depths[(number.tips+1):length(list_depths)]
  sort(list_heights_internal_nodes)
}, tip_pattern)
print(paste("true sample sizes="))
print(true_sample_sizes)
num_coal_events = unlist(lapply(list_ML_all_heights_pop1, FUN=function(x) length(x)))
print(paste("num coal events="))
print(num_coal_events)

#stopifnot((true_sample_sizes-1)==num_coal_events)

tips_ML_indices_pop1<- lapply(list_ML_trees, FUN=function(tree, tip_pattern){
  tree_without_healthy<- drop.tip(tree, healthycell, rooted=TRUE)
  list_tips<-tree_without_healthy$tip.label
  result<-which(list_tips %in% list_tips[grepl(tip_pattern, list_tips)] )
}, tip_pattern)
library(dplyr, warn.conflicts = FALSE)

print(length(list_ML_trees))
ML_association<-lapply(list_ML_trees, FUN = function(ML_tree){
  tree_without_healthy<- drop.tip(ML_tree, healthycell, rooted=TRUE)
  tree_tibble<-tidytree::as_tibble(tree_without_healthy)
  
  
  tree_tibble<-tree_tibble[!is.na(tree_tibble$label),]
  table<-as.data.frame(tree_tibble)
  
  table$tip_id <- with(table, as.numeric(str_match(table$label, "tip_i00\\s*(.*?)\\s*_C")[,2])+1)
  
  table<-table[, c(1,2,5)]
  #table<-table[order(table$coal_time),]
  table
})


# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

list_ML_tables <-lapply(list_ML_trees, FUN = function(tree){
  tree_without_healthy<- drop.tip(tree, healthycell, rooted=TRUE)
  tree_tibble<-tidytree::as_tibble(tree_without_healthy)
  tree_tibble<-tree_tibble[!is.na(tree_tibble$branch.length),]
  table<-tree_tibble %>% group_by(parent) %>% summarise(childs = paste0(node, collapse = ","))
  number.tips<-length(tree_without_healthy$tip.label)
  list_depths<-node.depth.edgelength(tree_without_healthy)
  max_depth= max(list_depths)
  list_heights_internal_nodes<-rep(max_depth,
                                   length(list_depths)-number.tips)-list_depths[(number.tips+1):length(list_depths)]
  table<-as.data.frame(table)
  table<-table %>% separate(childs, into=c("left", "right"), ",")
  table$coal_time <-  list_heights_internal_nodes
  #table<-table[order(table$coal_time),]
  table
})

list_ML_ordered_tables <-lapply(list_ML_trees, FUN = function(tree){
  tree_without_healthy<- drop.tip(tree, healthycell, rooted=TRUE)
  tree_tibble<-tidytree::as_tibble(tree_without_healthy)
  tree_tibble<-tree_tibble[!is.na(tree_tibble$branch.length),]
  table<-tree_tibble %>% group_by(parent) %>% summarise(childs = paste0(node, collapse = ","))
  number.tips<-length(tree_without_healthy$tip.label)
  list_depths<-node.depth.edgelength(tree_without_healthy)
  max_tip_depth= max(list_depths[1:number.tips])
  max_depth= max(list_depths)
  list_heights_internal_nodes<-rep(max_depth,
                                   length(list_depths)-number.tips)-list_depths[(number.tips+1):length(list_depths)]
  table<-as.data.frame(table)
  table<-table %>% separate(childs, into=c("left", "right"), ",")
  table$coal_time <-  list_heights_internal_nodes
  table$depth <-  list_depths[(number.tips+1):length(list_depths)]
  table<-table[order(-1*table$depth),]
  table
})


list_ML_postorder_tables <-lapply(list_ML_trees, FUN = function(tree){
  tree_without_healthy<- drop.tip(tree, healthycell, rooted=TRUE)
  #tree_without_healthy<-reorder(tree_without_healthy, order = "postorder")
  ntips<- length(tree_without_healthy$tip.label)
  mm <- matrix(nrow = ntips-1, ncol = 3)
  pos=1
  postorder_list<-postorder(tree_without_healthy)
  for (i in (1:length(postorder_list))){
    
    if (i %%2 ==1)
    {
      father <- tree_without_healthy$edge[postorder_list[i],1]
      left_child=tree_without_healthy$edge[postorder_list[i],2]
      right_child=tree_without_healthy$edge[postorder_list[i+1],2]
      mm[pos, 1] = father
      mm[pos, 2] = left_child
      mm[pos, 3] = right_child
      pos=pos+1
    }
  }
  mm
})


list_ML_ordered_tables_by_pop<- lapply(1:length(list_ML_trees), FUN=function(i, list_ML_ordered_tables, list_ML_postorder_tables, list_ML_all_heights_pop1){
  ordered_df<-list_ML_ordered_tables[[i]]
  heights_pop1<-list_ML_all_heights_pop1[[i]]
  ordered_df$population= rep(0,nrow(ordered_df))
  idx=match("population", names(ordered_df))
  for (i in 1:nrow(ordered_df)){
    coal_time_in_pop1= length(which(abs(heights_pop1-ordered_df$coal_time[i])<0.00001))
    if (coal_time_in_pop1 ==1)
      ordered_df[i,idx] = 2
    else
      ordered_df[i,idx] = 1
  }
  
  ordered_df[
    with(ordered_df, order(population, coal_time)),
  ]
  
},list_ML_ordered_tables, list_ML_postorder_tables,   list_ML_all_heights_pop1)


#########################################################################
#### Bayesian inference with stan
#########################################################################
#########################################################################

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

################################################################################################33
#################################################################################################33
# general settings
#####################################################################################################################
Sys.setenv(USE_CXX14 = "1")
Sys.setenv(PKG_CXXFLAGS = StanHeaders:::CxxFlags(as_character = TRUE))
SH <- system.file(ifelse(.Platform$OS.type == "windows", "libs", "lib"), .Platform$r_arch,
                  package = "StanHeaders", mustWork = TRUE)
Sys.setenv(PKG_LIBS = paste0(StanHeaders:::LdFlags(as_character = TRUE),
                             " -L", shQuote(SH), " -lStanHeaders"))
Sys.setenv(PKG_CXXFLAGS = paste0(Sys.getenv("PKG_CXXFLAGS"), " -I",
                                 system.file("include", "src", 
                                             package = "StanHeaders", mustWork = TRUE)))

Sys.setenv(PKG_CXXFLAGS = paste0(Sys.getenv("PKG_CXXFLAGS"), " -I",
                                 system.file("include", "src", 
                                             package = "StanHeaders", mustWork = TRUE)))

path_to_stan_model_one_population = paste0(path_to_save,
                                           "/stan_model_one_pop_GTJC16_asc_16_constant_sites.stan")

K=0.8
num_iter=2000
N=1

require(foreach)
library(doParallel)
library("doFuture")


registerDoFuture()

n.cpus_per_task <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n.cpus_per_task)){
  n.cpus_per_task <- 4
}

ntasks = as.numeric(Sys.getenv("SLURM_NTASKS"))
if (is.na(ntasks)){
  ntasks<- 4
}
print("ntasks")
print(ntasks)
n.cpus <-n.cpus_per_task * ntasks
print("n.cpus")
print(n.cpus)

ncores =ntasks
innerCluster <-
  parallel::makeCluster(ncores, type = "FORK", outfile = "")
on.exit(parallel::stopCluster(innerCluster))
doParallel::registerDoParallel(innerCluster)
#registerDoSEQ()
registerDoParallel(cores=ntasks)

num_chains= 4
n_threads = n.cpus %/% num_chains

Sys.setenv(STAN_NUM_THREADS = n_threads)
#env <- c(callr::rcmd_safe_env(), PROCESSX_NOTIFY_OLD_SIGCHLD = "true")

#list_stan_fit<-foreach(i=start_idx:(start_idx), .export=c("stan", "as.data.frame", "apply", "array", "str_match", "print", "paste0"),
#    .packages = c("stringr"),
#    .combine='c' ) %dopar% {
  i= start_idx
  
      

  file_name = list_IUPACgenotypes_only_filenames[[i]]
  Delta = true_deltas[[i]]
  Torigin = true_torigins[[i]]
  theta= true_thetas[[i]]
  n_var_sites = list_number_var_sites[[i]]
  nMU= true_nMUs[[i]]
  n= true_sample_sizes[[i]]
  res <- str_match(file_name, "true_hap_(.*?)_n=\\s*")
  parameter_text = res[,2]

  j= which(grepl( parameter_text, ML_tree_only_filenames, fixed = TRUE))
  
  print("i=")
  print(i)
  print("j=")
  print(j)

  tree_filename = ML_tree_only_filenames[[j]]
  
  print(paste0("Starting HMC inference for replicate ",i, " IUPAC file: ",  file_name, ", tree file: ",tree_filename ))
  topology <-as.data.frame(list_ML_ordered_tables[[j]][,1:3])
  tipdata <- list_IUPAC_compressed_phyDat[[i]]
  weights <- list_site_weights_compressed_phyDat[[i]]

  tip_association= ML_association[[j]][,3]
 
  max_cols = 10000000
  L=min(max_cols, ncol(tipdata))
  L <- min(max_cols, ncol(tipdata)) 
  if (L==max_cols){
    positions=sample(ncol(tipdata) , size=L, replace=FALSE)
    print(positions)
    tipdata2<-tipdata[, positions]
    weights<-weights[positions]
  }else{
    tipdata2<-tipdata
  }
  
  print("L=", str(L))

  
  topology <- apply(as.matrix.noquote(topology), 
                    2,
                    as.integer)
  print(topology)
  
  tipdata2<-apply(as.matrix(tipdata2), 
                 2,
                 as.integer)
  pop_sizes_proportion = array(c(1), dim=c(1))
  total_sample_size=true_sample_sizes[i]

  number_branches<-2*(total_sample_size-1)
  #check the topology
  coal_times_to_branch_lengths<- matrix(rep(0, number_branches*( total_sample_size-1)), nrow=number_branches, ncol= total_sample_size-1)
  map_internal_node_topology_row= rep(0, 2* total_sample_size-1)
  for(k in 1:(total_sample_size-1)){
    map_internal_node_topology_row[topology[k,1]]=k; 
    coal_times_to_branch_lengths[2*k-1, k]=1.0;
    coal_times_to_branch_lengths[2*k, k]=1.0;
    
  }
  for( k in 1:(total_sample_size-1)){
    if(as.integer(topology[k,2]) >total_sample_size)
    {
      
      coal_times_to_branch_lengths[2*k-1,  map_internal_node_topology_row[as.integer(topology[k,2])]]= -1.0;
    }
    if(as.integer(topology[k,3]) >total_sample_size)
    {
      
      coal_times_to_branch_lengths[2*k,  map_internal_node_topology_row[as.integer(topology[k,3])]]= -1.0;
    }
  }
  print("Checking coal_times_to_branch_lengths equal to 0")
  print(total_sample_size)
  stopifnot(length(which(coal_times_to_branch_lengths!=0))==(3*total_sample_size-4))
  print("Checking coal_times_to_branch_lengths less than  0")
  stopifnot(length(which(coal_times_to_branch_lengths<0))==(total_sample_size-2))
  
  row_sums<-rowSums(coal_times_to_branch_lengths, dims = 1)
  print("Checking row sums equal to  0")
  stopifnot(length(which(row_sums==0))==(total_sample_size-2))
  print("Checking row sums equal equals to   1")
  stopifnot(length(which(row_sums==1))==total_sample_size)
  ("All checking passed!")

  seq_error <- 0.000
  ado_error <- 0.000
  print(n_var_sites)
  number_invariable_sites = 1000000-n_var_sites
  input_stan<-list(K=K,N=N, total_sample_size=total_sample_size,
                   L=L,
                   genotype_tipdata=tipdata2,
                   site_weights=weights,
                   topology = topology,
                   tip_association= tip_association,
                   coal_times_to_branch_lengths = as.matrix(coal_times_to_branch_lengths),
                   seq_amp_error = seq_error,
                   allele_dropout_error = ado_error,
                   number_invariant_sites = number_invariable_sites,
                    n_cores= n_threads)
    
    stopifnot(L==length(weights))
    stopifnot(total_sample_size==length(tip_association))
    print(tip_association)

    check_cmdstan_toolchain()
  
  model <- cmdstan_model(path_to_stan_model_one_population,cpp_options = list(stan_threads = TRUE) )
  
  file_name <- str_remove(file_name, ".phy")
  file_name <- str_remove(file_name, ".txt")
  file_name <- str_remove(file_name, "true_hap_")
  output_path<- paste(
            path_to_save,
            path_all_data_folders,
            "/",
            as.character(i),
            "_",
            file_name,
            "_sp=",
            as.character(L),
            "_true_iter_",
            num_iter,
            ".rds",
            sep = ""
          )
  
 stan_fit <- model$sample(data = input_stan, chains = num_chains, parallel_chains = num_chains,
                                               threads_per_chain = n_threads,
                                               refresh = 100,
                                               iter_warmup = num_iter %/% 2,
                                               iter_sampling = num_iter %/% 2)

stan_fit$save_object(file = output_path)
  
#}




print ("Stan inference finished")
