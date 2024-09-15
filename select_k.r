library(rlist)
library(ape)
library(matrixStats)
library(patchwork)
library(geiger)
library(igraph)
library(dplyr)
library(ggplot2)
library(Seurat)

# source('/home/rx238/rds/hpc-work/aging/code/select_k_functions.r')

build_phylo_tree <- function(x){
  print('Building phylogenetic tree')
  
  # add in dummy program as outgroup to root tree (for depth-first search)
  #W_list = lapply(x$consensus.results, '[[', "W")
  W_list = x
  #W_list = lapply(W_list, t)
  #W_list[["root"]] <- as.matrix(data.frame(root = sample(rowMedians(W_list$R2))))
  W_list[["root"]] <- as.matrix(data.frame(root = rowMedians(W_list$R2)))
  W_list = W_list[c("root", names(x))]
  
  # build phylogenetic tree
  cor_mat = cor(list.cbind(W_list))
  phylo_tree = fastme.bal(1-cor_mat)
  phylo_tree = root(phylo_tree, outgroup = "root", resolve.root = T)
  phylo_tree = drop.tip(phylo_tree, "root")
  # convert negative branches to zero and filter 
  phylo_tree_pruned = phylo_tree
  phylo_tree_pruned$edge.length[phylo_tree_pruned$edge.length<0]=0
  dist_mat = cophenetic.phylo(phylo_tree_pruned)
  
  return(list(phylo_tree = phylo_tree, phylo_tree_pruned = phylo_tree_pruned, distance_matrix = dist_mat))
}

identify_outlier_programs <- function(x, ncells = 50){
  print('Identifying outlier programs')
  test = x
  res = lapply(test, function(y) {colMaxs(y)/apply(y, 2, function(z) {mean(sort(z, decreasing = T)[(1:ncells)+1])})})
  return(res)
}

suggest_dist_thresh <- function(phylo_trees){
  print('Suggesting distance threshold')

  res_list = lapply(phylo_trees, function(x){
    mat = x$distance_matrix
    diag(mat) = NA
    res = list(min = quantile(colMins(mat, na.rm = T), 0.95), max = max(colMins(mat, na.rm = T)))
  })
  common_min_thresh = max(unlist(lapply(res_list, `[[`, "min")))
  common_max_thresh = min(unlist(lapply(res_list, `[[`, "max")))
  res = c(common_min_thresh, common_max_thresh)
  names(res) = c("min", "max")
  return(res)
}

plot_hists = function(phylo_trees, suggested_thresholds = NULL, thresh.use = NULL){
  
  if (is.null(suggested_thresholds)){
    suggested_thresholds = suggest_dist_thresh(phylo_trees)
  }
  if (suggested_thresholds[1]>suggested_thresholds[2]){
    warning("no common threshold found")
  }
  if (is.null(thresh.use)){
    thresh.use = round(max(suggested_thresholds), 1)
  }
  x.max = max(unlist(lapply(phylo_trees, `[[`, "distance_matrix")))
  p = list()
  for (i in names(phylo_trees)){
    print('Plotting histograms')
    x = phylo_trees[[i]]
    df = data.frame(dist = x[["distance_matrix"]][upper.tri(x[["distance_matrix"]])])
    bins = 100
    p[[i]] = ggplot(df, aes(x = dist)) + 
      geom_histogram(aes(y = stat(count)/sum(count)), bins = bins, fill = "grey") + 
      geom_vline(xintercept = thresh.use, linetype = "dashed", color = "red") +
      labs(x = "Pairwise patristic distance", y = "Frequency", title = i) +
      scale_x_continuous(limits = c(0, x.max), expand = c(0,0))
    
    plot_filename = paste0("/home/rx238/rds/hpc-work/aging/NMF_lym_result/hists/k_hist_plot_", i, ".png") 
    ggsave(plot_filename, plot = p[[i]], width = 10, height = 8, dpi = 300)
  }
  return(p)
}

partition_phylo_tree <- function(x, y, dist.thresh = NULL, outlier.thresh = 5){
  print('Partitioning phylogenetic tree')
  if (is.null(dist.thresh)){
    stop("run suggest_dist_thresh")
  }
  tree = x$phylo_tree_pruned
  nodes = 1:tree$Nnode+length(tree$tip.label)
  names(nodes) = nodes
  dist_mat = x$distance_matrix
  res = lapply(nodes, function(x){
    tiplist = tips(tree, x)
    dist = median(dist_mat[tiplist,tiplist][upper.tri(dist_mat[tiplist,tiplist], diag = F)])
  })
  dist = unlist(res)
  ntips<-Ntip(tree)
  cnum <- 0 ## cluster number
  assign <- rep(0,ntips+tree$Nnode) ## cluster assignment
  igraph.tree <- graph.edgelist(tree$edge) ## tree in igraph form
  dfs <- graph.dfs(igraph.tree,root=ntips+1,neimode='out',
                   order=TRUE,dist=TRUE)
  distvec = dist
  ## transverse the tree in depth first order
  for(i in 1:length(dfs$order)){
    node <- dfs$order[i]
    ## skip leaves
    if(node < ntips+1){ next }
    ## If the node's subtree is below the threshold, mark it and
    ## its subtree as members of a new cluster
    if(distvec[node-ntips]<=dist.thresh && assign[node]<=0){
      cnum <- cnum+1
      subtree <- graph.dfs(igraph.tree,node,
                           neimode='out',unreachable=FALSE)$order
      subtree <- subtree[! is.na(subtree)]
      assign[subtree] <- cnum
    }}
  ans <- as.character(assign)
  ans <- ans[1:ntips]
  names(ans) <- tree$tip.label
  
  # identify outlier activity programs
  outliers = unlist(lapply(y, function(z) names(which(z>outlier.thresh))))
  ans[outliers] = "outlier"
  
  # set minimum subtree size to at least 2
  ans = plyr::mapvalues(ans, names(which(table(ans)<=2)), rep("0", length(names(which(table(ans)<=2)))))
  return(ans)
  
}

calculate_K_metric <- function(clusters, K.max = NULL){
  print('Calculating K metric')
  if (is.null(K.max)){
    K.max = max(as.numeric(gsub("R", "", gsub("_.*", "", names(clusters)))))
  }
  df = data.frame(program = names(clusters), clust = clusters)
  df$rank = gsub("_.*", "", df$program)
  # filter subtrees containing outlier programs or fewer than 5 activity programs
  df = data.frame(program = names(clusters), clust = clusters)
  df$rank = gsub("_.*", "", df$program)
  df = subset(df, !clust%in%c("0", "outlier"))
  df = df[df$clust%in%names(which(table(df$clust)>=5)),]
  
  # weight subtrees by the total number of programs in that subtree
  df$clust_weight = as.numeric(plyr::mapvalues(df$clust, from = names(table(df[,c("clust")])), to = table(df[,c("clust")])/(K.max - 1)))
  df$rank_clust = paste(df$rank, df$clust, sep = "_")
  # print(df)
  # only count subtrees once
  duplicated_rows <- which(duplicated(df$rank_clust) | duplicated(df$rank_clust, fromLast = TRUE))
  if(length(duplicated_rows) > 0) {
  print(df[duplicated_rows, ])
  } else {
  print("No duplicated rows found.")
  }
  df = df[!duplicated(df$rank_clust), ]
  # df = df[-which(duplicated(df$rank_clust)),]
  res = df %>%
    group_by(rank) %>%
    summarise(weighted_n_subtrees = sum(clust_weight))
  res = data.frame(res)
  rownames(res) = res$rank
  res = res[paste0("R", 2:K.max),]
}

suggest.k = function(K_metrics){
  non_zero_indices = which(K_metrics$weighted_n_subtrees != 0)

  largest_non_zero_index = max(non_zero_indices)
  ind = largest_non_zero_index
  k.suggest1 = which(round(K_metrics$weighted_n_subtrees)>=round(K_metrics$weighted_n_subtrees[ind]))+1
  k.suggest = intersect(k.suggest1, which(diff(K_metrics$weighted_n_subtrees)<0)+1)[1]
  if (is.na(k.suggest)){
    a = which(diff(K_metrics$weighted_n_subtrees)<=0)+1
    k.suggest = intersect(k.suggest1, which(diff(K_metrics$weighted_n_subtrees)<=0)+1)[1]
  }
  return (k.suggest)
}


color_set = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#1B9E77", "#D95F02", "#7570B3",
              "#E7298A", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#6C9A91", "#DF55BD",
              "#DEE14C", "#D19A83", "#DFBBD6", "#89E055", "#859C59", "#E2A550", "#6AE195", "#8664D5", "#D0E588",
              "#D894D6", "#6CE3DE", "#DDD9AA", "#9DE2BA", "#D1586D", "#7EC3E4", "#A33BE3")


plot_phylo_trees <- function(phylo_trees, phylo_partitions){
  
  # Ensure color_set and alpha function are defined correctly before this function
  for (i in names(phylo_trees)){
    print(paste0('Plotting phylogenetic trees ', i))
    # Open a PDF device for each plot
    pdf(paste0("/home/rx238/rds/hpc-work/aging/NMF_lym_result/trees/phylo_tree_", i, ".pdf"))
    
    tree = phylo_trees[[i]]$phylo_tree_pruned
    clusts = phylo_partitions[[i]]
    phy_sub = drop.tip(tree, names(clusts)[which(clusts %in% c("0", "outlier"))])
    plot(phy_sub, cex = 0.7, show.tip.label = F, type = "phylo", no.margin = F, main = i)
    labels = clusts[!(clusts %in% c("0", "outlier"))]
    labels = factor(labels)
    nclust = sum(!(levels(labels) %in% c("0", "outlier")))
    
    # Assuming color_set is defined and contains enough colors
    # Assuming alpha function is available for use
    tiplabels(pch = 16, cex = rep(1, nclust)[labels], 
              col = color_set[1:nclust][labels])
    tiplabels(pch = 1, lwd = 0.1, cex = rep(1, nclust)[labels], 
              col = rep(alpha("black", 0.3), nclust)[labels])
    
    # Close the PDF device
    dev.off()
  }
}

NMF_results_opt_k <- function(NMF_results_W, NMF_results_H, k, out_score, outlier.thresh = 5) {
  NMF_results_W = NMF_results_W[[paste0("R", k)]]
  NMF_results_H = NMF_results_H[[paste0("R", k)]]
  out_score = out_score[[paste0("R", k)]]
  NMF_results_W = NMF_results_W[,out_score<=outlier.thresh,drop = FALSE]
  NMF_results_W = t(NMF_results_W)
  NMF_results_H = NMF_results_H[,out_score<=outlier.thresh,drop = FALSE]
  NMF_results = list(W = NMF_results_W, H = NMF_results_H)
  return(NMF_results)
}


NMF_results_opt_k1 <- function(NMF_results, k, out_score, outlier.thresh = 5) {
  print('Optimizing K')
  NMF_results = NMF_results[[paste0("R", k)]]
  out_score = out_score[[paste0("R", k)]]
  # NMF_results_list = lapply(NMF_results, function(x) x[,out_score<=outlier.thresh])
  # NMF_results$W = NMF_results$W[out_score<=outlier.thresh,]
  # NMF_results$H = NMF_results$H[,out_score<=outlier.thresh]

  print(out_score)

  if (all(out_score >= outlier.thresh)) {
    # Return NMF_results as an empty list if there are no scores less than the threshold
    return(list())
  }

  values_less_than_5 = out_score[out_score < 5]

# Extract the names of elements where the condition is TRUE
  names_of_elements_less_than_5 <- names(values_less_than_5)

  NMF_results <- sapply(strsplit(names_of_elements_less_than_5, "_"), function(x) x[2])

  return(NMF_results)

}

##read in files

parent_directory <- '/home/rx238/rds/hpc-work/aging/NMF_lym/cNMF'
directories <- list.dirs(path = parent_directory, recursive = FALSE, full.names = FALSE)
print(directories)

new_list <- list()

cell_names <- sub("_cNMF$", "", directories)


for (directory in directories) {
    cellname = sub("_cNMF$", "", directory)
    new_list[[cellname]] <- list()
    filesInfolder <- list.files(paste0(parent_directory, '/', directory, '/tpm'))
    path = paste0(parent_directory, '/', directory, '/tpm/')
    for(file in filesInfolder) {
    k <- gsub(".*k_([0-9]+).*", "\\1", file)
    name <- paste0('R', k)
    print(name)
    print(paste0(path, file))
    temp <- read.table(paste0(path, file), header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE, row.names = 1)
    rownames(temp) <- paste(name, rownames(temp), sep = "_")
    normalized_temp <- apply(temp, 1, function(x) x / sqrt(sum(x^2)))
    #normalized_temp <- L2Dim(temp, , new.dr = NULL, new.key = NULL)
    new_list[[cellname]][[name]] <- as.matrix(normalized_temp)
    }
}

cell_list <- list()

for (directory in directories) {
    cellname = sub("_cNMF$", "", directory) 
    cell_list[[cellname]] <- list()
    filesInfolder <- list.files(paste0(parent_directory, '/', directory, '/usage'))
    path <- paste0(parent_directory, '/', directory, '/usage/')
    print(directory)
    for(file in filesInfolder) {
        k <- gsub(".*k_([0-9]+).*", "\\1", file)
        name <- paste0('R', k)
        print(name)
        temp <- read.table(paste0(path, file), header = TRUE, sep = "\t", na.strings = "NA", stringsAsFactors = FALSE, row.names = 1)
        temp <- as.matrix(temp)
        colnames(temp) <- gsub("X", "", paste0(name, "_", colnames(temp)))
        cell_list[[cellname]][[name]] <- temp
    }
}

saveRDS(new_list, "/home/rx238/rds/hpc-work/aging/NMF_lym_result/K_file_new_list.rds")
saveRDS(cell_list, "/home/rx238/rds/hpc-work/aging/NMF_lym_result/K_file_cell_list.rds")

new_list <- readRDS("/home/rx238/rds/hpc-work/aging/NMF_lym_result/K_file_new_list.rds")
cell_list <- readRDS("/home/rx238/rds/hpc-work/aging/NMF_lym_result/K_file_cell_list.rds")


new_list <- readRDS("/home/rx238/rds/hpc-work/aging/NMF_result_new/K_file_new_list.rds")
cell_list <- readRDS("/home/rx238/rds/hpc-work/aging/NMF_result_new/K_file_cell_list.rds")


phylo_trees <- lapply(new_list, build_phylo_tree)

saveRDS(phylo_trees, "/home/rx238/rds/hpc-work/aging/NMF_lym_result/phylo_trees.rds")

phylo_trees <- readRDS("/home/rx238/rds/hpc-work/aging/NMF_result_new/phylo_trees.rds")

program_outlier_score <- lapply(cell_list, identify_outlier_programs)
suggested_thresholds = suggest_dist_thresh(phylo_trees)
thresh_use = round(max(suggested_thresholds), 2)
thresh_use = 0.05
plot_hists(phylo_trees, suggested_thresholds, thresh.use = thresh_use)
phylo_partitions = mapply(partition_phylo_tree, x = phylo_trees, y = program_outlier_score, dist.thresh = thresh_use, outlier.thresh = 5, SIMPLIFY = F)
plot_phylo_trees(phylo_trees, phylo_partitions)
K_metrics = lapply(phylo_partitions, calculate_K_metric)
k.use = lapply(K_metrics, suggest.k)
NMF_results_atK <- mapply(FUN = NMF_results_opt_k, new_list, cell_list, k.use, program_outlier_score, SIMPLIFY = F)
NMF_results_atK1 <- mapply(FUN = NMF_results_opt_k1, new_list, k.use, program_outlier_score, SIMPLIFY = F)


saveRDS(NMF_results_atK, "/home/rx238/rds/hpc-work/aging/NMF_result_new/NMF_results_atK_all.rds")


saveRDS(NMF_results_atK, "/home/rx238/rds/hpc-work/aging/NMF_lym_result/NMF_results_atK_all.rds")

for (i in 1:length(K_metrics)){
  name = names(K_metrics)[i]
  pdf(paste0("/home/rx238/rds/hpc-work/aging/NMF_lym_result/weights/weighted_sub_", name, ".pdf"))
  plot(2:(dim(K_metrics[[i]])[1]+1), K_metrics[[i]]$weighted_n_subtrees,
       ylab = "Weighted N subtrees", xlab = "NMF Rank K", 
       main = names(K_metrics)[i], pch = 16, type = "b")
  abline(v = k.use[[i]])
  dev.off()
}

for (i in 1:length(program_outlier_score)){
  name = names(program_outlier_score)[i]
  k = k.use[[i]]
  if (is.na(k)){
    next
  }
  k_str = as.character(k.use[[i]])
  print(k)
  program = paste0("R", k_str)
  pdf(paste0("/home/rx238/rds/hpc-work/aging/NMF_lym_result/program_score/outlier_score_", program, '_', name, ".pdf"))
  hist(unlist(program_outlier_score[[i]][[program]]), main = name, xlab = "Outlier Score")
  dev.off()
}

max_length <- max(sapply(NMF_results_atK1, length))
padded_list <- lapply(NMF_results_atK1, function(x) {
  # Ensure x is a vector; if not, create a vector of NAs of max_length
  if (length(x) == 0){
    x <- rep(NA, max_length)
  }  # Adjusted to handle NULLs; choose appropriate type
  
  # Now pad with NAs to ensure length is max_length
  else if (length(x) < max_length) {
    x <- c(x, rep(NA, max_length - length(x)))
  }
  x
})
# Now convert the padded list to a data frame and write to CSV as before
my_df_padded <- data.frame(padded_list)
write.csv(my_df_padded, "/home/rx238/rds/hpc-work/aging/NMF_lym_result/all_p_selected.csv", row.names = FALSE)

k_selected <- data.frame(k.use)
write.csv(k_selected, "/home/rx238/rds/hpc-work/aging/NMF_lym_result/k_selected.csv", row.names = FALSE)


