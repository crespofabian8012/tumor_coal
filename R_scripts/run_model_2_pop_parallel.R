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
    #"phytools",
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
  install.packages("Rcpp", repos="https://rcppcore.github.io/drat", type="source")

install_required_packages = function(list.of.packages) {
  new.packages <-
    list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages) > 0)
    install.packages(new.packages, dependencies = TRUE, repos = "http://cran.us.r-project.org")
  
}


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", dependencies = TRUE, repos = "http://cran.us.r-project.org")

#if (!require("ggtree", quietly = TRUE))
#  BiocManager::install("ggtree")

if (!require("supraHex", quietly = TRUE))
  BiocManager::install("supraHex")

if (!require("Rgraphviz", quietly = TRUE))
  BiocManager::install("Rgraphviz")



install_required_packages(list_packages)

if (!require("rstan", quietly = TRUE))
  install.packages("rstan", dependencies=TRUE, repos = "http://cran.us.r-project.org")

lapply(list_packages,
       library,
       character.only = TRUE)

library(rstan)

getCurrentFileLocation <-  function()
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
path = getCurrentFileLocation()

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
    "-"= 0,
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
foo<-function(tree,...){
  fsize<-36*par()$pin[2]/par()$pin[1]/Ntip(tree)
  plotTree(tree,fsize=fsize,lwd=1,...)
  nodelabels(cex=0.5) 
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
path_to_save = getCurrentFileLocation()


path_all_data_folders="/data"


path_all_IUPACgenotypes<-paste0(path_to_save,path_all_data_folders,"/true_haplotypes_dir")#IUPAC codes


healthycell = "healthycell"
tip_pattern<-"*_C1_*"
healthy_pattern<-"*_C2_*"


list_tip_pattern= list(  tip_pattern, healthy_pattern)
#ML tree topologies


ML_tree_filenames <- list.files(paste0(path_to_save,path_all_data_folders,"/Cellphy_ML_trees/"), pattern=".bestTree", full.names=TRUE, all.files=FALSE)

ML_tree_only_filenames <- list.files(paste0(path_to_save,path_all_data_folders, "/Cellphy_ML_trees/"), pattern=".bestTree", full.names=FALSE, all.files=FALSE)



list_ML_trees<-lapply(ML_tree_filenames, FUN=function(x) {
  castor::read_tree(file=x)
})


list_ML_trees<-lapply(list_ML_trees, FUN=function(x) {
  root(x, outgroup=healthycell)
})

list_ordered_tree_topologies_by_increasing_time<-lapply(1:length(list_ML_trees), FUN = function(i,healthycell){
  
  
  tree=list_ML_trees[[i]]

  tree_without_healthy<- ape::drop.tip(tree, healthycell, rooted=TRUE)
  
  tree_tibble<-dplyr::as_tibble(tree_without_healthy)
  
  tree_tibble<-tree_tibble[!is.na(tree_tibble$branch.length),]
  table<-tree_tibble %>% group_by(parent) %>% summarise(childs = paste0(node, collapse = ","))
  number.tips<-length(tree_without_healthy$tip.label)
  
  list_depths<-node.depth.edgelength(tree_without_healthy)
  
  max_tip_depth= max(list_depths[1:number.tips])
  max_depth= max(list_depths)
  stopifnot(max_tip_depth==max_depth)
  list_heights_internal_nodes<-rep(max_depth,
                                   length(list_depths)-number.tips)-list_depths[(number.tips+1):length(list_depths)]
  table<-as.data.frame(table)
  table<-table %>% tidyr::separate(childs, into=c("left", "right"), ",")
  table$coal_time <-  list_heights_internal_nodes
  table$depth <-  list_depths[(number.tips+1):length(list_depths)]
  table<-table[order(-1*table$depth),]
  return(table)
}, healthycell)


list_idxs_MRCAs_by_pop<-lapply(list_ML_trees, FUN = function(tree,healthycell, list_tip_pattern){
  tree_without_healthy= tree
  #if (!str_detect(healthycell, tip_pattern)){
  tree_without_healthy<- ape::drop.tip(tree, healthycell, rooted=TRUE)
  #}
  
  list_tips<-tree_without_healthy$tip.label
  total_number_tips = length(list_tips)
  
  idxs_mrca_by_pop<-unlist(lapply(list_tip_pattern, function(pattern){
    
    idxs=grepl(pattern, list_tips)
    tips_pop<-list_tips[idxs]
    
    idx_mrca<-ape::getMRCA(tree_without_healthy, tips_pop)
    idx_mrca
  }))
  
  return(unlist(idxs_mrca_by_pop))
  
}, healthycell, list_tip_pattern)


list_order_pop_increasing_time<-lapply(1:length(list_ML_trees), FUN = function(i,healthycell, 
                                                                               list_tip_pattern, 
                                                                               list_idxs_MRCAs_by_pop,
                                                                               list_ordered_tree_topologies_by_increasing_time){
  
  
  tree=list_ML_trees[[i]]
  list_idxs_MRCAs= list_idxs_MRCAs_by_pop[[i]]
  
  table<-list_ordered_tree_topologies_by_increasing_time[[i]]
  pos_MRCAs= unlist(lapply(list_idxs_MRCAs, function(idx_MRCA){
    which(table$parent==idx_MRCA)
    
  }))
  
  return(order(pos_MRCAs))
  
}, healthycell, list_tip_pattern, list_idxs_MRCAs_by_pop, list_ordered_tree_topologies_by_increasing_time)

list_tip_pattern_order<-lapply(list_order_pop_increasing_time, FUN = function(ord, list_tip_pattern){
  
  return(unlist(list_tip_pattern[ord]))
},list_tip_pattern)

list_idxs_MRCAs_by_pop<-lapply(1:length(list_order_pop_increasing_time), FUN = function(i, 
                                                                                        list_tip_pattern_order, 
                                                                                        list_idxs_MRCAs_by_pop,
                                                                                        list_ordered_tree_topologies_by_increasing_time){
  ord=list_order_pop_increasing_time[[i]]
  original_list_MRCAs=list_idxs_MRCAs_by_pop[[i]]
  
  reordered_list<-unlist(original_list_MRCAs[ord])
  
  table_all_tree<-list_ordered_tree_topologies_by_increasing_time[[i]]
  
  if (reordered_list[length(reordered_list)]!=table_all_tree[nrow(table_all_tree), ]$parent){
    
    reordered_list[length(reordered_list)]<-table_all_tree[nrow(table_all_tree), ]$parent
  }
  
  return(reordered_list)
},list_tip_pattern_order, list_idxs_MRCAs_by_pop, list_ordered_tree_topologies_by_increasing_time)

list_IUPACgenotypes_filenames<-list.files(path_all_IUPACgenotypes, all.files=FALSE,
                                          pattern = "*\\.[phy txt]*",
                                          full.names=TRUE)
list_IUPACgenotypes_only_filenames<-list.files(path_all_IUPACgenotypes, all.files=FALSE,
                                               pattern = "*\\.[phy txt]",
                                               full.names=FALSE)

if (length(list_IUPACgenotypes_filenames)==0){
  
  stop("Error: The number of IUPAC files is 0!")
}


list_IUPACgenotypes_data<-lapply(list_IUPACgenotypes_filenames, FUN=function(filename){
  phylotools::read.phylip(filename)
})





list_IUPAC_numeric_genotypes_data<- lapply(1:length(list_IUPACgenotypes_data), FUN=function(i,  healthycell, list_tip_pattern_order){
  list_data= list_IUPACgenotypes_data[[i]]
  ras<-list_data$seq.text
  names<-list_data$seq.name
  
  thealthycell = gsub('\\.', '_', healthycell)
  
  filtered<-vector()
  
  list_tip_pattern= list_tip_pattern_order[[i]]
  
  lapply(list_tip_pattern, function(pattern){
    idxs=grepl(pattern, names)
    idxs[names==thealthycell]=FALSE
    filtered<<-c(filtered,ras[idxs])
  })
  #check filtered has all sequences except the healthy cell
  stopifnot(length(filtered)+1==length(names))
  min_number_loci = min(unlist(lapply(filtered, FUN= function(p) nchar(p))))
  
  list_numeric_seq<-lapply(1:(length(filtered)), FUN=function(j){
    seq<-NLP::as.String(filtered[j])
    from_list_IUPAC_to_numeric_codes(seq)
  })
  matrix(unlist(list_numeric_seq), nrow=(length(filtered)),ncol=length(list_numeric_seq[[1]]),  byrow = T)
},  healthycell, list_tip_pattern_order)

order_of_seq<- lapply(1:length(list_IUPACgenotypes_data), FUN=function(i, healthycell,list_tip_pattern_order){
  
  list_data=list_IUPACgenotypes_data[[i]]
  names<-list_data$seq.name
  
  list_tip_pattern= list_tip_pattern_order[[i]]
  filtered<-vector()
  thealthycell = gsub('\\.', '_', healthycell)
  
  lapply(list_tip_pattern, function(pattern){
    idxs=grepl(pattern, names)
    idxs[names==thealthycell]=FALSE
    filtered<<-c(filtered,names[idxs])
  })

  
  df=data.frame(label =filtered, order_id= 1:length(filtered))
  df
}, healthycell, list_tip_pattern_order)

list_number_var_sites<-unlist(lapply(1:length(list_IUPACgenotypes_filenames), FUN=function(i, list_IUPACgenotypes_filenames){
  res <- str_match(list_IUPACgenotypes_filenames[[i]], "_nS=\\s*(.*?)\\s*_")
  as.numeric(res[,2])
}, list_IUPACgenotypes_filenames))

list_phyDat<- lapply(list_IUPACgenotypes_data, FUN=function(list_data,  healthycell){
  
  
  thealthycell = gsub('\\.', '_', healthycell)
  
  list_data<-list_data[list_data$seq.name !=thealthycell, ]
  ras<-list_data$seq.text
  ras<-list_data$seq.text
  names<-list_data$seq.name
  
  filtered<-ras
  filtered_names=names
  
  list_numeric_seq<-lapply(1:(length(filtered)), FUN=function(j){
    seq<-NLP::as.String(filtered[j])
    list_IUPAC_codes(seq)
  })
  
  res<-phangorn::phyDat(matrix(unlist(list_numeric_seq), nrow=(length(filtered)),ncol=length(list_numeric_seq[[1]]),  byrow = TRUE, 
                               dimnames= 
                                 list(filtered_names,  1:length(list_numeric_seq[[1]]))), 
                        levels=c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K", "N"), compress=TRUE) 
  res
},  healthycell)

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
  compress_dat<-subset(phy_dat, 1:length(list_phyDat[[i]]), compress=TRUE, site.pattern = TRUE)
  weight= attr(phy_dat, "weight")
  weight
})



order_of_seq<- lapply(list_IUPACgenotypes_data, FUN=function(list_data, healthycell, list_tip_pattern_order){
  
  names<-list_data$seq.name
  filtered<-names
  
  
  thealthycell = gsub('\\.', '_', healthycell)
  filtered = filtered[filtered != thealthycell]
  df=data.frame(label =filtered, order_id= 1:length(filtered))
  df
}, healthycell, list_tip_pattern_order)


######################################################################



library(stringr)


list_idx_parent_populations <-lapply(1:length(list_ML_trees), FUN = function(i,healthycell, 
                                                                             list_tip_pattern_order, 
                                                                             list_idxs_MRCAs_by_pop,
                                                                             list_ordered_tree_topologies_by_increasing_time){
  list_tip_pattern<- list_tip_pattern_order[[i]]
  tree<- list_ML_trees[[i]]
  tree_without_healthy<- tree
  
  tree_without_healthy<- ape::drop.tip(tree, healthycell, rooted=TRUE)

  all_tree_table<-list_ordered_tree_topologies_by_increasing_time[[i]]
  
  list_tips<-tree_without_healthy$tip.label
  total_number_tips <- length(list_tips)
  
  idxs_mrca_by_pop<-list_idxs_MRCAs_by_pop[[i]]
  
  indices_parent_populations<- unlist(lapply(idxs_mrca_by_pop,function(current_idx_mrca, idxs_mrca_by_pop)
  {
    
    if (current_idx_mrca==all_tree_table[nrow(all_tree_table),]$parent){ 
      return (0);
    }
    else{
      
      parent_pop_idx=0
      num_tips_smallest_subtree=2*total_number_tips
      lapply(1:length(idxs_mrca_by_pop), function(i, current_idx_mrca, tree_without_healthy) {
        mrca_parent_candidate = idxs_mrca_by_pop[i]
        if (mrca_parent_candidate!= current_idx_mrca){
          descendants=  phytools::getDescendants(tree_without_healthy, mrca_parent_candidate)
          if (current_idx_mrca %in% descendants){
            if(length(descendants) <= num_tips_smallest_subtree){
              parent_pop_idx <<- i
              num_tips_smallest_subtree<<-length(descendants)
            }
          }
        }
      }, current_idx_mrca, tree_without_healthy)
      return (parent_pop_idx);
    }
  }, idxs_mrca_by_pop))
  
  stopifnot(indices_parent_populations[ length(indices_parent_populations) ]==0)

  return(unlist(indices_parent_populations))
}, healthycell, list_tip_pattern_order, list_idxs_MRCAs_by_pop, list_ordered_tree_topologies_by_increasing_time)


ML_association<-lapply(1:length(list_ML_trees), FUN = function(i, healthycell,
                                                               list_tip_pattern_order, 
                                                               order_of_seq,
                                                               list_idx_parent_populations,
                                                               list_ordered_tree_topologies_by_increasing_time){
  ML_tree = list_ML_trees[[i]]
  order_seq= order_of_seq[[i]]
  list_tip_pattern= list_tip_pattern_order[[i]]

    tree_without_healthy<- ape::drop.tip(ML_tree, healthycell, rooted=TRUE)
 
    list_tips<-tree_without_healthy$tip.label
    
  
    tables_by_pop<-lapply(list_tip_pattern, function(pattern){
      list_tips_to_remove<-list_tips[!grepl(pattern, list_tips)]
      sub_tree<-ape::drop.tip(tree_without_healthy, list_tips_to_remove)
      
      idxs=grepl(pattern, list_tips)
      tips_pop<-list_tips[idxs]
      
      idx_mrca<-ape::getMRCA(tree_without_healthy, tips_pop)
      tree_tibble<-tidytree::as_tibble(sub_tree)

      tree_tibble<-tree_tibble[1:length(sub_tree$tip.label),]
      
      table<-as.data.frame(tree_tibble)
     
      table$label= as.character(table$label)
      
      table$label <- gsub('\\.', '_', table$label)
      
      order_seq$label= as.character(order_seq$label)
      
      merged<-   dplyr::inner_join(table, order_seq, by= "label")
 
      merged<-merged[, c(1,2,4, 5)]
      colnames(merged)[which(names(merged) == "order_id")] <- "tip_id"

      merged
      
    })
    final_table = do.call("rbind", tables_by_pop)
    final_table
    
}, healthycell,list_tip_pattern_order, order_of_seq, list_idx_parent_populations, list_ordered_tree_topologies_by_increasing_time)

library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)


list_sample_sizes= lapply(1:length(list_ML_trees), FUN = function(i, healthycell, 
                                                   list_idxs_MRCAs_by_pop,
                                                   list_idx_parent_populations,
                                                   list_tip_pattern_order)
  {
  tree = list_ML_trees[[i]]
  tree_without_healthy= tree
  tree_without_healthy<- ape::drop.tip(tree, healthycell, rooted=TRUE)
  list_tips<-tree_without_healthy$tip.label
  total_sample_size=length(list_tips)
  
  idx_MRCAs_by_pop= list_idxs_MRCAs_by_pop[[i]]
  idx_parent_pop= list_idx_parent_populations[[i]]
  list_tip_pattern<-list_tip_pattern_order[[i]]
  
  sample_sizes<-lapply(list_tip_pattern, function(pattern){
    
    idxs=grepl(pattern, list_tips)
    tips_pop<-list_tips[idxs]
    return(length(tips_pop))
  })


  return(unlist(sample_sizes))

},healthycell, list_idxs_MRCAs_by_pop, list_idx_parent_populations,list_tip_pattern_order)

list_internal_events= lapply(1:length(list_ML_trees), FUN = function(i, healthycell, 
                                                                     list_sample_sizes, 
                                                                     list_idx_parent_populations,
                                                                     list_idxs_MRCAs_by_pop)
{
  tree = list_ML_trees[[i]]
  tree_without_healthy= tree
  tree_without_healthy<- ape::drop.tip(tree, healthycell, rooted=TRUE)
  list_tips<-tree_without_healthy$tip.label
  total_sample_size=length(list_tips)
  idx_MRCAs_by_pop<-list_idxs_MRCAs_by_pop[[i]]
  
  idx_parent_pop= list_idx_parent_populations[[i]]
  sample_sizes= list_sample_sizes[[i]]
  all_tree_table<-list_ordered_tree_topologies_by_increasing_time[[i]]
  
  num_internal_events=unlist(lapply(1:length(idx_MRCAs_by_pop), FUN=function(j){
    
    
    MRCA_pop = idx_MRCAs_by_pop[j]
    
    all_descendants <-phytools::getDescendants(tree_without_healthy, MRCA_pop)
   
    internals<-all_descendants[all_descendants>total_sample_size]
    
    child_pop_MRCAs= idx_MRCAs_by_pop[idx_MRCAs_by_pop %in% all_descendants] 
    
    all_descendants <-c(all_descendants, MRCA_pop)
    all_child_pop_descendants=c()
    
    lapply(child_pop_MRCAs, FUN = function(child_pop_MRCA){
      
      all_child_pop_descendants<<-c(all_child_pop_descendants, phytools::getDescendants(tree_without_healthy, child_pop_MRCA))
      all_child_pop_descendants<<-c(all_child_pop_descendants,child_pop_MRCA )
      
    })
    
    filtered_table <- all_tree_table[ all_tree_table$parent %in% all_descendants, ]
    filtered_table <- filtered_table[ ! filtered_table$parent %in% all_child_pop_descendants, ]
    filtered_table<-filtered_table[order(-1*filtered_table$depth),]
    
    return(nrow(filtered_table))
  }))
  return(unlist(num_internal_events))
  
},healthycell, list_sample_sizes, list_idx_parent_populations,list_idxs_MRCAs_by_pop)



list_idxs_parent_MRCA_nodes_inside_parent_pop_table= lapply(1:length(list_ML_trees), FUN = function(i, healthycell, 
                                                                  list_idxs_MRCAs_by_pop, 
                                                                  list_sample_sizes,
                                                                  list_idx_parent_populations, 
                                                                  list_tip_pattern_order,
                                                                  list_ordered_tree_topologies_by_increasing_time)
{
  tree = list_ML_trees[[i]]
  list_tip_pattern= list_tip_pattern_order[[i]]
  
  tree_without_healthy= tree
  tree_without_healthy<- ape::drop.tip(tree, healthycell, rooted=TRUE)
  
  all_tree_table<-list_ordered_tree_topologies_by_increasing_time[[i]]
  
  list_tips<-tree_without_healthy$tip.label
  total_sample_size=length(list_tips)
  
  idx_MRCAs_by_pop= list_idxs_MRCAs_by_pop[[i]]
  idx_parent_pops= list_idx_parent_populations[[i]]

  sample_sizes= list_sample_sizes[[i]]
  
  child_directions<-vector("list", length(idx_parent_pops))
  
  list_pos=unlist(lapply(1:length(idx_parent_pops), FUN=function(j){
    idx_father_pop= idx_parent_pops[j]
    MRCA_idx_in_parent_pop= idx_MRCAs_by_pop[j]
    if (idx_father_pop==0){
      return(0)
    }
    else{
      MRCA_parent_pop <-idx_MRCAs_by_pop[idx_father_pop]
      father_pop_descendants<- phytools::getDescendants(tree_without_healthy, MRCA_parent_pop)
      father_pop_descendants<-c(father_pop_descendants,MRCA_parent_pop )
      
      idx_all_child_pop=which(idx_parent_pops==idx_father_pop  )
      all_child_pop_descendants=c()
      
      lapply(idx_all_child_pop, FUN = function(idx_child_pop){
        
        all_child_pop_descendants<<-c(all_child_pop_descendants, phytools::getDescendants(tree_without_healthy, idx_MRCAs_by_pop[idx_child_pop]))
        all_child_pop_descendants<<-c(all_child_pop_descendants, idx_MRCAs_by_pop[idx_child_pop])
      
        })
         
      filtered_table <- all_tree_table[ all_tree_table$parent %in% father_pop_descendants, ]
      filtered_table <- filtered_table[ ! filtered_table$parent %in% all_child_pop_descendants, ]
      filtered_table<-filtered_table[order(-1*filtered_table$depth),]
      
      pos_left =which(filtered_table$left==MRCA_idx_in_parent_pop)
      pos_right =which(filtered_table$right==MRCA_idx_in_parent_pop)
      
      if (length(pos_left)==0){
        child_directions[j]<<-1#right
      }else{
        child_directions[j]<<-0
      }
     
      pos=max(pos_left, pos_right)
      
      return(pos)
    }
    
  }))


  return(list_pos)
  
}, healthycell,list_idxs_MRCAs_by_pop, list_sample_sizes, list_idx_parent_populations, list_tip_pattern_order, list_ordered_tree_topologies_by_increasing_time)


list_number_child_pops= lapply(1:length(list_ML_trees), FUN = function(i, list_idx_parent_populations)
{
  idx_parent_pops= list_idx_parent_populations[[i]]
  
  number_child_pops=unlist(lapply(1:length(idx_parent_pops), FUN=function(j){
    length(which(idx_parent_pops==j))
  }))
  return(unlist(number_child_pops))
  
}, list_idx_parent_populations)

list_child_pop_indexes= lapply(1:length(list_ML_trees), FUN = function(i, list_idx_parent_populations, list_number_child_pops)
{
  idx_parent_pops= list_idx_parent_populations[[i]]
  number_idxs_child_pop = list_number_child_pops[[i]]
  
  idx_child_pops=unlist(lapply(1:length(idx_parent_pops), FUN=function(j){
      result=(which(idx_parent_pops==j))
      result
  }))
  return(unlist(idx_child_pops))
  
}, list_idx_parent_populations, list_number_child_pops)

list_trees_without_healthy<-lapply(1:length(list_ML_trees), FUN = function(i,healthycell){
  tree = list_ML_trees[[i]]
  tree_without_healthy<- ape::drop.tip(tree, healthycell, rooted=TRUE)
  return(tree_without_healthy)
  
},healthycell )

list_ML_ordered_tables_by_pop <-lapply(1:length(list_ML_trees), FUN = function(i,healthycell, list_tip_pattern_order, list_idxs_MRCAs_by_pop,
                                                                         list_sample_sizes,
                                                                         list_idx_parent_populations, 
                                                                         list_ordered_tree_topologies_by_increasing_time){

  tree = list_ML_trees[[i]]
  list_tip_pattern=list_tip_pattern_order[[i]]

  tree_without_healthy<- ape::drop.tip(tree, healthycell, rooted=TRUE)
  idx_MRCAs_by_pop= list_idxs_MRCAs_by_pop[[i]]
  idx_parent_pops= list_idx_parent_populations[[i]]

  table_all_tree<-list_ordered_tree_topologies_by_increasing_time[[i]]
  list_tips<-tree_without_healthy$tip.label
  total_sample_size<-length(list_tips)

  tables_by_pop=lapply(1:length(idx_parent_pops), FUN=function(j, idx_parent_pops, idx_MRCAs_by_pop, tree_without_healthy){
     idx_father_pop= idx_parent_pops[j]
     MRCA_idx_in_parent_pop= idx_MRCAs_by_pop[j]

      idx_all_child_pop=which(idx_parent_pops==j)
      all_child_pop_descendants=c()

      lapply(idx_all_child_pop, FUN = function(idx_child_pop){

        all_child_pop_descendants<<-c(all_child_pop_descendants, phytools::getDescendants(tree_without_healthy, idx_MRCAs_by_pop[idx_child_pop]))
        all_child_pop_descendants<<-c(all_child_pop_descendants,idx_MRCAs_by_pop[idx_child_pop])

      })

 

      current_subtree_descendants=c()
      current_subtree_descendants = c(current_subtree_descendants, phytools::getDescendants(tree_without_healthy, idx_MRCAs_by_pop[j]))
      current_subtree_descendants = c(current_subtree_descendants, idx_MRCAs_by_pop[j])
      internal_descendants<-current_subtree_descendants[current_subtree_descendants>total_sample_size]

      filtered_table <-table_all_tree[table_all_tree$parent %in% internal_descendants, ]


      if (!(is.null(all_child_pop_descendants))){
        internal_descendants_child_pop<-all_child_pop_descendants[all_child_pop_descendants>total_sample_size]
        filtered_table <- filtered_table[ ! filtered_table$parent %in% all_child_pop_descendants, ]
      }

      print(filtered_table)
      return(as.data.frame(filtered_table))

  }, idx_parent_pops, idx_MRCAs_by_pop, tree_without_healthy)

  final_table = do.call("rbind", tables_by_pop)
  
 
  
  final_table

}, healthycell, list_tip_pattern_order, list_idxs_MRCAs_by_pop,  list_sample_sizes, list_idx_parent_populations, list_ordered_tree_topologies_by_increasing_time )





#########################################################################

#test
BAD_count<-0
                 
BAD_list<-c()        
test<-lapply(1:length(list_IUPACgenotypes_only_filenames), FUN=function(i) {                         
                         
                         file_name = list_IUPACgenotypes_only_filenames[[i]]
                         
                         res <- str_match(file_name, "true_hap_(.*?).txt")
                         parameter_text = res[,2]
                         
                         idx_j= which(grepl( parameter_text, ML_tree_only_filenames, fixed = TRUE))
                         print(parameter_text)
                         print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
                         print("i=")
                         print(i)
                         if (length(idx_j)==0){
                           j=i
                         }else{
                           j=idx_j
                         }
                         print("j=")
                         print(j)
                         
                         tree_filename = ML_tree_only_filenames[[j]]
                         
                         print(paste0("Starting HMC inference for replicate ",i, " IUPAC file: ",  file_name, ", tree file: ",tree_filename ))
                        
                         file_name = list_IUPACgenotypes_only_filenames[[i]]
                         
                         tree_without_healthy<-list_trees_without_healthy[[i]]
                         
                         tip_pattern_order<-list_tip_pattern_order[[i]]
                         print("tip_pattern_order")
                         print(tip_pattern_order)
                         
                         topology <-as.data.frame(list_ML_ordered_tables_by_pop[[i]][,1:3])
                         print(topology)
                         
                         tipdata <- list_IUPAC_compressed_phyDat[[i]]
                         weights <- list_site_weights_compressed_phyDat[[i]]
                         
                         tip_association= ML_association[[i]][,4]
                         print("tip_associaction")
                         print(tip_association)
                         idxs_MRCA=list_idxs_MRCAs_by_pop[[i]]
                         print("idxs_MRCA")
                         print(idxs_MRCA)
                         indexes_father_populations =list_idx_parent_populations[[i]]
                         print("indexes_father_populations")
                         print(indexes_father_populations)
                         sample_sizes=list_sample_sizes[[i]]
                         print("sample_sizes")
                         print(sample_sizes)
                         total_sample_size=sum(sample_sizes)
                         N=length(sample_sizes)
                         print("N")
                         print(N)
                         number_internal_nodes=list_internal_events[[i]]
                         print("number_internal_nodes")
                         print(number_internal_nodes)
                         child_pop_indexes= as.array(list_child_pop_indexes[[i]])
                         print("child_pop_indexes")
                         print(child_pop_indexes)
                         number_child_pops= list_number_child_pops[[i]]
                         pos_child_pop_MRCA_in_parent_pop=list_idxs_parent_MRCA_nodes_inside_parent_pop_table[[i]]
                         print("pos_child_pop_MRCA_in_parent_pop")
                         print(pos_child_pop_MRCA_in_parent_pop)
                         indexes_parent_populations=unlist(lapply(indexes_father_populations, function(x) {x[x!=0]}))
                         indexes_parent_populations<-c(indexes_parent_populations,0)
                         print("indexes_parent_populations")
                         print(indexes_parent_populations)
                         
                         print("Checking pos_child_pop_MRCA_in_parent_pop")
                         is_BAD<-FALSE
                         for( k in 1:N){
                           idx_father_pop<-indexes_parent_populations[[k]]
                           idx_pop_MRCA<-idxs_MRCA[[k]]
                           
                           if (idx_father_pop!=0){
                             if(pos_child_pop_MRCA_in_parent_pop[k]>number_internal_nodes[idx_father_pop]){
                               is_BAD<-TRUE
                             }
                           }
                           descendants<- phytools::getDescendants(tree_without_healthy, idx_pop_MRCA)
                           descendants<-c(descendants,idx_pop_MRCA )
                           tips_under_MRCA<-descendants[descendants<=total_sample_size]
                           
                           all_child_pop_samples<-0
                           if (length(which(indexes_parent_populations==k)) >0){
                             all_child_pop_samples<-sample_sizes[which(indexes_parent_populations==k)]
                           }
                            all_sample_size_under_MRCA<-sample_sizes[k]+sum(all_child_pop_samples)
                           
                             if(length(tips_under_MRCA)!=all_sample_size_under_MRCA){
                               is_BAD<-TRUE
                             }
                         }
                         if (is_BAD){
                           BAD_count<<-BAD_count+1
                           BAD_list<<-c(BAD_list, i)     
                         }    
                         max_cols <- 500
                         L <- min(max_cols, ncol(tipdata))
                         positions=sample(ncol(tipdata) , size=L, replace=FALSE)
                         
                         tipdata2<-tipdata[, positions]
                         weights<-weights[positions]
                         
                         
                         total_sample_size=sum(sample_sizes)
                         print(total_sample_size)
                         print(total_sample_size)
                         
                         print("Checking idxs_MRCA")
                         stopifnot(length(idxs_MRCA)==N)
                         print("Checking number_internal_nodes")
                         stopifnot(length(number_internal_nodes)==N)
                         print("Checking sample_sizes")
                         stopifnot(length(sample_sizes)==N)
                         print("Checking number_internal_nodes")
                         print(sum(number_internal_nodes))
                         print(total_sample_size-1)
                         stopifnot(sum(number_internal_nodes)==(total_sample_size-1))
                         print("Checking tip_association")
                         stopifnot(length(tip_association)==(total_sample_size))
                         print("Checking tipdata2 col")
                         stopifnot(ncol(tipdata2)==length(weights))
                         print("Checking tipdata2 row")
                         stopifnot(nrow(tipdata2)==total_sample_size)
                         print("Checking topology row")
                         stopifnot(nrow(topology)==(total_sample_size-1))
                         
                         number_branches<-2*(total_sample_size-1)
                         #check the topology
                         coal_times_to_branch_lengths<- matrix(rep(0, number_branches*( total_sample_size-1)), nrow=number_branches, ncol= total_sample_size-1)
                         map_internal_node_topology_row= rep(0, 2* total_sample_size-1)
                         for(i in 1:(total_sample_size-1)){
                           map_internal_node_topology_row[topology[i,1]]=i; 
                           coal_times_to_branch_lengths[2*i-1, i]=1.0;
                           coal_times_to_branch_lengths[2*i, i]=1.0;
                           
                         }
                         for( i in 1:(total_sample_size-1)){
                           if(as.integer(topology[i,2]) >total_sample_size)
                           {
                             
                             coal_times_to_branch_lengths[2*i-1,  map_internal_node_topology_row[as.integer(topology[i,2])]]= -1.0;
                           }
                           if(as.integer(topology[i,3]) >total_sample_size)
                           {
                             
                             coal_times_to_branch_lengths[2*i,  map_internal_node_topology_row[as.integer(topology[i,3])]]= -1.0;
                           }
                         }
                         print("Checking coal_times_to_branch_lengths equal to 0")
                         stopifnot(length(which(coal_times_to_branch_lengths!=0))==(3*total_sample_size-4))
                         print("Checking coal_times_to_branch_lengths less than  0")
                         stopifnot(length(which(coal_times_to_branch_lengths<0))==(total_sample_size-2))
                         
                         row_sums<-rowSums(coal_times_to_branch_lengths, dims = 1)
                         print("Checking row sums equal to  0")
                         stopifnot(length(which(row_sums==0))==(total_sample_size-2))
                         print("Checking row sums equal equals to   1")
                         stopifnot(length(which(row_sums==1))==total_sample_size)
                         ("All checking passed!")
                         input_stan<-list(K=K,N=N, total_sample_size=total_sample_size,
                                          L=L,
                                          ids_MRCA=idxs_MRCA,
                                          sample_sizes=sample_sizes,
                                          number_internal_nodes=number_internal_nodes,
                                          indexes_father_populations=indexes_father_populations,
                                          number_child_pops= number_child_pops,
                                          child_pop_indexes=child_pop_indexes,
                                          pos_parent_child_pop_MRCA_in_parent_pop=pos_child_pop_MRCA_in_parent_pop,
                                          genotype_tipdata=tipdata2,
                                          site_weights=weights,
                                          tip_association=tip_association,
                                          topology=topology,
                                          coal_times_to_branch_lengths=as.matrix(coal_times_to_branch_lengths))
                         
                         output=parameter_text
                         if (length(idx_j)==0){
                           output="result"+ as.character(i)
                           output=parameter_text
                         }
                         
                       print(output)
                       })
all_list<-seq(1,length(list_IUPACgenotypes_only_filenames) , 1)
white_list <- setdiff(all_list, BAD_list)
if (is.na(start_idx)){
  start_idx<- 1
}


if (!(start_idx %in% white_list)){
  
  stop(paste0("The dataset ", start_idx, " has a input tree topology that is not monophyletic!"))
  
}
################################################################################################33
#################################################################################################33
# general settings
#####################################################################################################################
#### Bayesian inference with stan
#########################################################################
#########################################################################
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


if (!require("cmdstanr", quietly = TRUE))
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")),  dependencies = TRUE)

library(cmdstanr)

path_to_stan_model = paste0(path_to_save,
                                           "/stan_model_n_population_physical_time_GTJC6_error_asc_paralllel_refactored.stan")
#stan_model_N_population_JC_genotypes_corrected5.stan es el modelo original(corregido) 

#stan_model_N_population_JC_genotypes_corrected3.stan consider the last population in model time(not multiply by x the last pop)

print(path_to_save)

#expose_stan_functions(path_to_stan_model_one_population)

K=0.8
num_iter=2000

list_idxs_datasets <- seq(start_idx,(start_idx+1),1)


n.cpus_per_task <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n.cpus_per_task)){
  n.cpus_per_task <- 4
}



taskID <- as.numeric(Sys.getenv('SLURM_PROCID'))
if (is.na(taskID)){
  taskID<- 0
}
print("taskID")
print(taskID)

ntasks = as.numeric(Sys.getenv("SLURM_NTASKS"))
if (is.na(ntasks)){
  ntasks<- 4
}
print("ntasks")
print(ntasks)
n.cpus <-n.cpus_per_task * ntasks
print("n.cpus")
print(n.cpus)
list1<-split(list_idxs_datasets, 1:ntasks)[[taskID+1]]
print(list1)



require(foreach)
library(doParallel)
library("doFuture")



registerDoFuture()
innerCluster <-
  parallel::makeCluster(ntasks, type = "FORK", outfile = "")
on.exit(parallel::stopCluster(innerCluster))
doParallel::registerDoParallel(innerCluster)
#registerDoSEQ()
registerDoParallel(cores=ntasks)

getDoParWorkers()
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


list_stan_fit<-foreach(i=list_idxs_datasets[1], .export=c("stan", "as.data.frame", "apply", "array", "str_match", "print", "paste0"),
                       .packages = c("stringr"),
                       .combine='c' ) %dopar% {
                         
                         
                         
                         file_name = list_IUPACgenotypes_only_filenames[[i]]
                         
                         n_var_sites = list_number_var_sites[[i]]
                         
                         res <- str_match(file_name, "true_hap_(.*?).txt")
                         parameter_text = res[,2]
                         
                         idx_j= which(grepl( parameter_text, ML_tree_only_filenames, fixed = TRUE))
                         print(parameter_text)
                         print("i=")
                         print(i)
                         if (length(idx_j)==0){
                           j=i
                         }else{
                           j=idx_j
                         }
                         print("j=")
                         print(j)
                         
                         tree_filename = ML_tree_only_filenames[[j]]
                         
                         print(paste0("Starting HMC inference for replicate ",i, " IUPAC file: ",  file_name, ", tree file: ",tree_filename ))
                         
                         
                         
                         
                         
                         file_name = list_IUPACgenotypes_only_filenames[[i]]
                         
                         
                         topology <-as.data.frame(list_ML_ordered_tables_by_pop[[i]][,1:3])
                         
                         topology <- apply(as.matrix.noquote(topology), 
                                           2,
                                           as.integer)
                         print(topology)
                         
                         tipdata <- list_IUPAC_compressed_phyDat[[i]]
                         weights <- list_site_weights_compressed_phyDat[[i]]
                         
                         tip_association= ML_association[[i]][,4]
                         idxs_MRCA=list_idxs_MRCAs_by_pop[[i]]
                         indexes_father_populations =list_idx_parent_populations[[i]]
                         sample_sizes=list_sample_sizes[[i]]
                         N=length(sample_sizes)
                         number_internal_nodes=list_internal_events[[i]]
                         child_pop_indexes= as.array(list_child_pop_indexes[[i]])
                         number_child_pops= list_number_child_pops[[i]]
                         pos_child_pop_MRCA_in_parent_pop=list_idxs_parent_MRCA_nodes_inside_parent_pop_table[[i]]
                         
                         indexes_parent_populations=unlist(lapply(indexes_father_populations, function(x) {x[x!=0]}))
                         
                         
                         max_cols <- 100000000
                         L <- min(max_cols, ncol(tipdata))
                         positions=sample(ncol(tipdata) , size=L, replace=FALSE)
                         
                         tipdata2<-tipdata[, positions]
                         
                         tipdata2<-apply(as.matrix(tipdata2), 
                                         2,
                                         as.integer)
                         weights<-weights[positions]
                         
                         
                         total_sample_size=sum(sample_sizes)
                         
                         print("Checking idxs_MRCA")
                         stopifnot(length(idxs_MRCA)==N)
                         print("Checking number_internal_nodes")
                         stopifnot(length(number_internal_nodes)==N)
                         print("Checking sample_sizes")
                         stopifnot(length(sample_sizes)==N)
                         print("Checking number_internal_nodes")
                         print(sum(number_internal_nodes))
                         print(total_sample_size-1)
                         stopifnot(sum(number_internal_nodes)==(total_sample_size-1))
                         print("Checking tip_association")
                         stopifnot(length(tip_association)==(total_sample_size))
                         print("Checking tipdata2 col")
                         stopifnot(ncol(tipdata2)==length(weights))
                         print("Checking tipdata2 row")
                         stopifnot(nrow(tipdata2)==total_sample_size)
                         print("Checking topology row")
                         stopifnot(nrow(topology)==(total_sample_size-1))
                         
                         number_branches<-2*(total_sample_size-1)
                         #check the topology
                         coal_times_to_branch_lengths<- matrix(rep(0, number_branches*( total_sample_size-1)), nrow=number_branches, ncol= total_sample_size-1)
                         map_internal_node_topology_row= rep(0, 2* total_sample_size-1)
                         for(i in 1:(total_sample_size-1)){
                           map_internal_node_topology_row[topology[i,1]]=i; 
                           coal_times_to_branch_lengths[2*i-1, i]=1.0;
                           coal_times_to_branch_lengths[2*i, i]=1.0;
                           
                         }
                         for( i in 1:(total_sample_size-1)){
                           if(as.integer(topology[i,2]) >total_sample_size)
                           {
                             
                             coal_times_to_branch_lengths[2*i-1,  map_internal_node_topology_row[as.integer(topology[i,2])]]= -1.0;
                           }
                           if(as.integer(topology[i,3]) >total_sample_size)
                           {
                             
                             coal_times_to_branch_lengths[2*i,  map_internal_node_topology_row[as.integer(topology[i,3])]]= -1.0;
                           }
                         }
                         print("Checking coal_times_to_branch_lengths equal to 0")
                         stopifnot(length(which(coal_times_to_branch_lengths!=0))==(3*total_sample_size-4))
                         print("Checking coal_times_to_branch_lengths less than  0")
                         stopifnot(length(which(coal_times_to_branch_lengths<0))==(total_sample_size-2))
                         
                         row_sums<-rowSums(coal_times_to_branch_lengths, dims = 1)
                         print("Checking row sums equal to  0")
                         stopifnot(length(which(row_sums==0))==(total_sample_size-2))
                         print("Checking row sums equal equals to   1")
                         stopifnot(length(which(row_sums==1))==total_sample_size)
                         print("Checking weights")
                         stopifnot(sum(weights)==n_var_sites)
                         stopifnot(L==length(weights))
                         stopifnot(total_sample_size==length(tip_association))
                         print(tip_association)
                         ("All checking passed!")
                         number_invariable_sites = 1000000-n_var_sites
                         seq_error <- 0.000
                         ado_error <- 0.000
                       
                         
                         input_stan<-list(K=K,N=N,
                                          total_sample_size=total_sample_size,
                                          L=L,
                                          ids_MRCA=idxs_MRCA,
                                          sample_sizes=sample_sizes,
                                          number_internal_nodes=number_internal_nodes,
                                          indexes_father_populations=indexes_father_populations,
                                          number_child_pops= number_child_pops,
                                          child_pop_indexes=child_pop_indexes,
                                          pos_parent_child_pop_MRCA_in_parent_pop=pos_child_pop_MRCA_in_parent_pop,
                                          genotype_tipdata=tipdata2,
                                          site_weights=weights,
                                          tip_association=tip_association,
                                          topology=topology,
                                          coal_times_to_branch_lengths=as.matrix(coal_times_to_branch_lengths),
                                          seq_amp_error=seq_error,
                                          allele_dropout_error=ado_error,
                                          number_invariant_sites=number_invariable_sites,
                                          n_cores= n_threads)
                         
                        
                         
                         check_cmdstan_toolchain()
                         
                         model <- cmdstan_model(path_to_stan_model,force_recompile = TRUE, cpp_options = list(stan_threads = TRUE, STANC2 = TRUE), pedantic = TRUE)
                         
                         model$check_syntax(pedantic = TRUE)
                         
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
                         
                        print("Starting inference")
                   
                         
                         stan_fit <- model$sample(data = input_stan, chains = num_chains, parallel_chains = num_chains,
                                                  threads_per_chain = n_threads,
                                                  refresh = 100,
                                                  iter_warmup = num_iter %/% 2,
                                                  iter_sampling = num_iter %/% 2)
                         
                         stan_fit$save_object(file = output_path)
                       }

print ("Stan inference finished")

output_path<-"/Users/faustofabiancrespofernandez/Downloads/Benchmark_PHD_thesis/codigo/two_pop_test/data/19_G0=100.0000_T0=0.0523_n0=10_m0=0.0093_x0=0.5000_G1=20.0000_T1=0.2239_n1=10_m1=0.0093_x1=0.5000_nMU=1.003228_nS=4647_0002_0001_0001_sp=437_true_iter_2000.rds"
options(max.print = 10000000)
options(cmdstanr_max_rows = 1000000000)

fit_cmdstan <- readRDS(output_path)
fit_cmdstan$summary()
print(fit_cmdstan$summary(), n=1000)
