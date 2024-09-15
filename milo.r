library(miloR)
library(igraph)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Matrix)
library(cowplot)
library(ggsci)
library(enrichR)
library(ggbeeswarm)


milo <- readRDS('/home/rx238/rds/hpc-work/aging/final_data/final_final/adata_all_reannotation_final1.rds')

milo <- milo[, (!milo$ct_level3 %in% c('Hs Nul2_2', 'Hs Uni2_1', 'Mo_TAM', 'Fb3 Uni2_2', 'Myocytes', 'Erythroid', 'Caf', 'Tumour')) & (!milo$Batch == 7) & (!milo$Sample %in% c('ctrl1', 'ctrl2', 'NulControl'))]

milo <- milo[, milo$Parity %in% c('NP', 'Parous')]
milo <- milo[, milo$Genotype %in% c('WT', 'WKBR')]
colData(milo) <- droplevels(colData(milo))


new_parity <- colData(milo)$Parity != 'NP'
colData(milo)$NewParity <- new_parity
milo$Age <- as.numeric(milo$Age)
milo$Genotype <- as.factor(milo$Genotype)
new_genotype <- colData(milo)$Genotype != 'WT'
colData(milo)$NewGenotype <- new_genotype
milo$NewGenotype <- as.factor(milo$NewGenotype)


####Immune

lym <- milo[, milo$ct_level1 == "Lymphoid"]
mye <- milo[, milo$ct_level1 == "Myeloid"]
epi <- milo[, milo$ct_level1 == "Epithelial"]
fib <- milo[, milo$ct_level1 == "Fibroblast"]
str <- milo[, milo$ct_level1 == "Stroma"]


cell_list <- list(lymphoid = lym,
                  myeloid = mye,
                    epithelial = epi,
                    fibroblast = fib,
                    stroma = str)



plotDAbeeswarm1 <- function(da.res, group.by=NULL, alpha=0.1, subset.nhoods=NULL){
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", group.by,")?")
    }
    if (is.numeric(da.res[,group.by])) {
      # stop(group.by, " is a numeric variable. Please bin to use for grouping.")
    }
    da.res <- mutate(da.res, group_by = da.res[,group.by])
  } else {
    da.res <- mutate(da.res, group_by = "g1")
  }

  if (!is.factor(da.res[,"group_by"])) {
    message("Converting group_by to factor...")
    da.res <- mutate(da.res, group_by = factor(group_by, levels=unique(group_by)))
    # anno_vec <- factor(anno_vec, levels=unique(anno_vec))
  }

  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods,]
  }

  # Get position with ggbeeswarm
  beeswarm_pos <- ggplot_build(
    da.res %>%
      mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
      arrange(group_by) %>%
      ggplot(aes(group_by, logFC)) +
      geom_quasirandom()
  )

  pos_x <- beeswarm_pos$data[[1]]$x
  pos_y <- beeswarm_pos$data[[1]]$y

  n_groups <- unique(da.res$group_by) %>% length()

  da.res %>%
    mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
    mutate(logFC_color = ifelse(is_signif==1, logFC, 0)) %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
    mutate(pos_x = pos_x, pos_y=pos_y) %>%
    ggplot(aes(pos_x, pos_y)) +
    geom_point(color="darkgrey") +
    #scale_color_gradient2() +
    #guides(color="none") +
    xlab(group.by) + ylab("Log Fold Change") +
    scale_x_continuous(
      breaks = seq(1,n_groups),
      labels = setNames(levels(da.res$group_by), seq(1,n_groups))
      ) +
    #eom_point() +
    coord_flip() +
    theme_bw(base_size=22) +
    theme(strip.text.y =  element_text(angle=0))

}

set.seed(123)
milo_func <- function(adata, celltype, k_value, comparison){
    
    set.seed(123)
    print('count')
    milo <- Milo(adata)
    names(assays(milo))=c("counts","logcounts", "raw")
    milo <- buildGraph(milo, k = k_value, d = 20, reduced.dim = "X_scVI", transposed = TRUE)
    milo <- makeNhoods(milo, prop = 0.1, k = k_value, d=20, refined = TRUE, reduced_dims = "X_scVI", refinement_scheme = "graph")
    plotNhoodSizeHist(milo)
    ggsave(paste0("/home/rx238/rds/hpc-work/aging/report_graph/milo/",celltype,'/',comparison,'/', comparison, '_', celltype, ' hist_new.png'),height=6, width=5.15)
    milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="Sample")

    if (comparison == "Parity"){
        print('parity')
        # adata <- adata[,(adata$Age >90) & (adata$Genotype == "WT")]
        # colData(adata) <- droplevels(colData(adata))
        # milo <- Milo(adata)
        # names(assays(milo))=c("counts","logcounts", "raw")
        # milo <- buildGraph(milo, k = k_value, d = 20, reduced.dim = "X_scVI", transposed = TRUE)
        # milo <- makeNhoods(milo, prop = 0.1, k = k_value, d=20, refined = TRUE, reduced_dims = "X_scVI", refinement_scheme = "graph")
        # plotNhoodSizeHist(milo)
        # ggsave(paste0("/home/rx238/rds/hpc-work/aging/milo/",celltype,'/',comparison,'/', comparison, '_', celltype, ' hist_new.png'),height=6, width=5.15)
        # milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="Sample")
        design_df <- data.frame(colData(milo)) %>%
        filter((Age > 80) & (Genotype == "WT")) %>%
        select(Sample, NewParity)
        design_df <- distinct(design_df)
        design_df <- droplevels(design_df)
        rownames(design_df) = design_df$Sample
        print('design3')
        milo.res <- testNhoods(milo, design = ~ NewParity, design.df = design_df, reduced.dim='X_scVI', fdr.weighting="graph-overlap")
    }
    else if (comparison == "Age"){
        print('age')
        design_df <- data.frame(colData(milo)) %>%
        filter((Genotype == 'WT') & (Parity == "NP")) %>%
        select(Sample, Age)
        design_df <- distinct(design_df)
        design_df <- droplevels(design_df)
        rownames(design_df) = design_df$Sample
        milo.res <- testNhoods(milo, design = ~ Age, design.df = design_df, reduced.dim='X_scVI', fdr.weighting="graph-overlap")
    }
    
    else if (comparison == "Genotype"){
        print('genotype')
        design_df <- data.frame(colData(milo)) %>%
        filter((Parity == "NP") & (Age < 60)) %>%
        select(Sample, Age, NewGenotype)
        design_df <- distinct(design_df)
        design_df <- droplevels(design_df)
        rownames(design_df) = design_df$Sample
        milo.res <- testNhoods(milo, design = ~ Age + NewGenotype, design.df = design_df, reduced.dim='X_scVI', fdr.weighting="graph-overlap")
    }
    
    ggplot(milo.res, aes(PValue)) + geom_histogram(bins=50)
    ggsave(paste0("/home/rx238/rds/hpc-work/aging/report_graph/milo/", celltype,'/',comparison,'/', comparison, '_', celltype, ' p_distribution_new.png'),height=6, width=5.15)
    
    ggplot(milo.res, aes(logFC, -log10(SpatialFDR))) + 
    geom_point() +
    geom_hline(yintercept = 1)
    ggsave(paste0("/home/rx238/rds/hpc-work/aging/report_graph/milo/", celltype,'/',comparison,'/', comparison, '_', celltype, ' volcano_new.png'),height=6, width=5.15)
    
    milo.res$Diff <- sign(milo.res$logFC)
    milo.res$Diff[milo.res$SpatialFDR > 0.1] <- 0
    
    milo <- buildNhoodGraph(milo, overlap=5)
    milo.res <- annotateNhoods(milo, milo.res, 'ct_level3')
    milo.res$CellType <- ifelse(milo.res$ct_level3_fraction < 0.5, "Mixed", milo.res$ct_level3)

    milo.res_new = milo.res[which(milo.res$CellType != "Mixed"),]
    milo_sorted <- milo.res_new[order(milo.res_new$CellType), ]

    print(unique(milo_sorted$CellType))
    
    plotReducedDim(milo, dimred = "UMAP", text_by = "ct_level3", colour_by=comparison, point_size=1) + guides(fill="none") 
    ggsave(paste0("/home/rx238/rds/hpc-work/aging/report_graph/milo/",celltype,'/',comparison,'/', comparison, '_', celltype, ' umap_new.png'),height=6, width=5.15)
    
    plotNhoodGraphDA(milo, layout="UMAP", milo_res=milo_sorted, alpha=0.1) 
    ggsave(paste0("/home/rx238/rds/hpc-work/aging/report_graph/milo/", celltype,'/',comparison,'/', comparison, '_', celltype, ' umap_DA_new.png'),height=6, width=5.15)
    
    if (!any(milo_sorted$Diff == 1)){
        plotDAbeeswarm1(milo_sorted, group.by = "CellType")
        ggsave(paste0("/home/rx238/rds/hpc-work/aging/report_graph/milo/", celltype,'/',comparison,'/', comparison, '_', celltype, ' bee_DA_new.png'),height=6, width=9)
    }
    
    else{
        plotDAbeeswarm(milo_sorted, group.by = "CellType")
        ggsave(paste0("/home/rx238/rds/hpc-work/aging/report_graph/milo/", celltype,'/',comparison,'/', comparison, '_', celltype, ' bee_DA_new.png'),height=6, width=9)
    }
    
    list(milo, milo.res)
}

result_list <- list()

cell_list <- result_list

for (cell in names(cell_list)){
    set.seed(123)
    print(cell)
    # print("Parity")
    # result_list[[cell]][['Parity']] <- milo_func(cell_list[[cell]], cell, 50, "Parity")
    print("Age")
    result_list[[cell]][['Age']] <- milo_func(cell_list[[cell]], cell, 50, "Age")
    # print("Genotype")
    # result_list[[cell]][['Genotype']] <- milo_func(cell_list[[cell]], cell, 50, "Genotype")
}

saveRDS(result_list, '/home/rx238/rds/hpc-work/aging/report_graph/milo/milo_list_result.rds')


milo <- result_list[['lymphoid']][['Age']][[1]]
milo.res <- result_list[['lymphoid']][['Age']][[2]]

set.seed(123)
milo.res <- groupNhoods(milo, milo.res, max.lfc.delta = 0.01, overlap = 5)

plotDAbeeswarm(milo.res , group.by = "NhoodGroup")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/age_t_reg/nhood_groups_new_bee.png')
plotNhoodGroups(milo, milo.res, layout="UMAP")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/age_t_reg/nhood_groups_new1.png')

milo.res_reg = milo.res[which(milo.res$CellType%in% c('CD4_Treg')),]

keep.rows <- rowSums(logcounts(milo)) != 0
milo <- milo[keep.rows,]
    
dec <- modelGeneVar(milo)
hvgs <- getTopHVGs(dec, n = 2000)

milo.markers <- findNhoodGroupMarkers(milo, milo.res,
                                    subset.row = hvgs,
                                    subset.nhoods = (milo.res$CellType %in% c('CD4_Treg') & (milo.res$NhoodGroup %in% c('14', '5'))),
                                    aggregate.samples = TRUE,
                                    sample_col = 'Sample',
                                    gene.offset = FALSE)

write.csv(milo.markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/age_t_reg/CD4_Treg_markers_new_1.csv')

markers_all <- milo.markers[which(milo.markers$adj.P.Val_14 < 0.05 & milo.markers$logFC_14 > 0),]
sorted_markers <- markers_all[order(markers_all$logFC_14, decreasing = TRUE), ]
top_30_genes <- head(sorted_markers$GeneID, 30)
markers <- top_30_genes

milo <- calcNhoodExpression(milo, subset.row=markers)

plotNhoodExpressionGroups(milo, milo.res, features=markers,
                          subset.nhoods = (milo.res$CellType %in% c('CD4_Treg')) & (milo.res$NhoodGroup %in% c('14', '5')),
                          scale=TRUE,
                          grid.space = "fixed")


ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/age_t_reg/CD4_TregG14vsG5_G14_new.png')

##### NKT




milo <- result_list[['lymphoid']][['Parity']][[1]]
milo.res <- result_list[['lymphoid']][['Parity']][[2]]

set.seed(123)
milo.res <- groupNhoods(milo, milo.res, max.lfc.delta = 2)

plotDAbeeswarm(milo.res , group.by = "NhoodGroup")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/parity_NKT/nhood_groups_new_bee.png')
plotNhoodGroups(milo, milo.res, layout="UMAP")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/parity_NKT/nhood_groups_new1.png')

milo.res_reg = milo.res[which(milo.res$CellType%in% c('CD8_naive')),]

keep.rows <- rowSums(logcounts(milo)) != 0
milo <- milo[keep.rows,]
    
dec <- modelGeneVar(milo)
hvgs <- getTopHVGs(dec, n = 2000)

milo.markers <- findNhoodGroupMarkers(milo, milo.res,
                                    subset.row = hvgs,
                                    subset.nhoods = (milo.res$CellType %in% c('CD8_naive') & (milo.res$NhoodGroup %in% c('16', '3'))),
                                    aggregate.samples = TRUE,
                                    sample_col = 'Sample',
                                    gene.offset = FALSE)

write.csv(milo.markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/parity_cd8_naive/CD8_naive_markers_new.csv')

###NAIVE CD8 PARITY AND NKT

markers_all <- milo.markers[which(milo.markers$adj.P.Val_16 < 0.01 & milo.markers$logFC_16 > 0),]
sorted_markers <- markers_all[order(markers_all$adj.P.Val_16), ]
top_30_genes <- head(sorted_markers$GeneID, 30)
markers <- top_30_genes

milo <- calcNhoodExpression(milo, subset.row=markers)

plotNhoodExpressionGroups(milo, milo.res, features=markers,
                          subset.nhoods = (milo.res$CellType %in% c('CD8_naive')) & (milo.res$NhoodGroup %in% c('16', '3')),
                          scale=TRUE,
                          grid.space = "fixed")


ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/parity_cd8_naive/CD8_naiveG16vsG3_G16_new.png')


###Epi Age Lp

result_list <- readRDS('/home/rx238/rds/hpc-work/aging/report_graph/milo/milo_list_result.rds')

milo <- result_list[['epithelial']][['Age']][[1]]
milo.res <- result_list[['epithelial']][['Age']][[2]]

set.seed(123)
milo.res <- groupNhoods(milo, milo.res, max.lfc.delta = 0.05)

plotDAbeeswarm(milo.res , group.by = "NhoodGroup")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_Lp/nhood_groups_new_bee.png')
plotNhoodGroups(milo, milo.res, layout="UMAP")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_Lp/nhood_groups_new1.png')

milo.res1 <- groupNhoods(milo, milo.res, max.lfc.delta = 0.05, subset.nhoods = (milo.res$CellType %in% c('Lp')))
plotNhoodGroups(milo, milo.res1, layout="UMAP")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_Lp/nhood_groups_checklp.png')


milo.res_reg = milo.res[which(milo.res$CellType%in% c('Lp')),]

keep.rows <- rowSums(logcounts(milo)) != 0
milo <- milo[keep.rows,]
    
dec <- modelGeneVar(milo)
hvgs <- getTopHVGs(dec, n = 2000)

milo.markers <- findNhoodGroupMarkers(milo, milo.res,
                                    subset.row = hvgs,
                                    subset.nhoods = (milo.res$CellType %in% c('Lp') & (milo.res$NhoodGroup %in% c('11', '1'))),
                                    aggregate.samples = TRUE,
                                    sample_col = 'Sample',
                                    gene.offset = FALSE)

write.csv(milo.markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_Lp/Lp_markers_new_1.csv')

markers_all <- milo.markers[which(milo.markers$adj.P.Val_1 < 0.01 & milo.markers$logFC_1 > 0),]
sorted_markers <- markers_all[order(markers_all$adj.P.Val_1), ]
top_30_genes <- head(sorted_markers$GeneID, 30)
markers <- top_30_genes

milo <- calcNhoodExpression(milo, subset.row=markers)

plotNhoodExpressionGroups(milo, milo.res, features=c('Areg', 'Mki67', 'Wnt4', 'Ly6a', 'Pdgfb'),
                          subset.nhoods = (milo.res$CellType %in% c('Lp')) & (milo.res$NhoodGroup %in% c('11', '1')),
                          scale=TRUE,
                          grid.space = "fixed")


ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_Lp/LpG11vsG1_G1_new.png')

sorted_markers <- as.data.frame(sorted_markers)
write.csv(sorted_markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_Lp/Lp_markers_G1_new_1.csv')



##Epi Parity Lp

milo <- result_list[['epithelial']][['Parity']][[1]]
milo.res <- result_list[['epithelial']][['Parity']][[2]]

set.seed(123)
milo.res <- groupNhoods(milo, milo.res, max.lfc.delta = 2.5, overlap = 8)

plotDAbeeswarm(milo.res , group.by = "NhoodGroup")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Parity_Lp/nhood_groups_new_bee.png')
plotNhoodGroups(milo, milo.res, layout="UMAP")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Parity_Lp/nhood_groups_new1.png')

# milo.res1 <- groupNhoods(milo, milo.res, max.lfc.delta = 0.05, subset.nhoods = (milo.res$CellType %in% c('Lp')))
# plotNhoodGroups(milo, milo.res1, layout="UMAP")
# ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_Lp/nhood_groups_checklp.png')


milo.res_reg = milo.res[which(milo.res$CellType%in% c('Lp')),]

keep.rows <- rowSums(logcounts(milo)) != 0
milo <- milo[keep.rows,]
    
dec <- modelGeneVar(milo)
hvgs <- getTopHVGs(dec, n = 2000)

milo.markers <- findNhoodGroupMarkers(milo, milo.res,
                                    subset.row = hvgs,
                                    subset.nhoods = (milo.res$CellType %in% c('Lp') & (milo.res$NhoodGroup %in% c('7', '13'))),
                                    aggregate.samples = TRUE,
                                    sample_col = 'Sample',
                                    gene.offset = FALSE)

write.csv(milo.markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Parity_Lp/Lp_markers_new_1.csv')

markers_all <- milo.markers[which(milo.markers$adj.P.Val_7 < 0.01 & milo.markers$logFC_7 > 0),]
sorted_markers <- markers_all[order(markers_all$adj.P.Val_7), ]
top_30_genes <- head(sorted_markers$GeneID, 30)
markers <- top_30_genes

milo <- calcNhoodExpression(milo, subset.row=markers)

plotNhoodExpressionGroups(milo, milo.res, features=markers,
                          subset.nhoods = (milo.res$CellType %in% c('Lp')) & (milo.res$NhoodGroup %in% c('7', '13')),
                          scale=TRUE,
                          grid.space = "fixed")


ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Parity_Lp/LpG7vsG13_G7_new.png')

sorted_markers <- as.data.frame(sorted_markers)
write.csv(sorted_markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Parity_Lp/Lp_markers_G7_new_1.csv')



##B cell Age

milo <- result_list[['lymphoid']][['Age']][[1]]
milo.res <- result_list[['lymphoid']][['Age']][[2]]

set.seed(123)
milo.res <- groupNhoods(milo, milo.res, max.lfc.delta = 0.01, overlap = 8)

plotDAbeeswarm(milo.res , group.by = "NhoodGroup")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_B/nhood_groups_new_bee.png')
plotNhoodGroups(milo, milo.res, layout="UMAP")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_B/nhood_groups_new1.png')

# milo.res1 <- groupNhoods(milo, milo.res, max.lfc.delta = 0.05, subset.nhoods = (milo.res$CellType %in% c('Lp')))
# plotNhoodGroups(milo, milo.res1, layout="UMAP")
# ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_Lp/nhood_groups_checklp.png')


milo.res_reg = milo.res[which(milo.res$CellType%in% c('B_cell')),]

keep.rows <- rowSums(logcounts(milo)) != 0
milo <- milo[keep.rows,]
    
dec <- modelGeneVar(milo)
hvgs <- getTopHVGs(dec, n = 2000)

milo.markers <- findNhoodGroupMarkers(milo, milo.res,
                                    subset.row = hvgs,
                                    subset.nhoods = (milo.res$CellType %in% c('B_cell') & (milo.res$NhoodGroup %in% c('26', '20', '4'))),
                                    aggregate.samples = TRUE,
                                    sample_col = 'Sample',
                                    gene.offset = FALSE)

write.csv(milo.markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_B/B_markers_new_1.csv')

markers_all <- milo.markers[which(milo.markers$adj.P.Val_20 < 0.01 & milo.markers$logFC_20 > 0),]
sorted_markers <- markers_all[order(markers_all$adj.P.Val_20), ]
top_30_genes <- head(sorted_markers$GeneID, 30)
markers <- top_30_genes

milo <- calcNhoodExpression(milo, subset.row=markers)

plotNhoodExpressionGroups(milo, milo.res, features=markers,
                          subset.nhoods = (milo.res$CellType %in% c('B_cell')) & (milo.res$NhoodGroup %in% c('20', '26', '4')),
                          scale=TRUE,
                          grid.space = "fixed")


ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_B/B_G20vsG26&G4_G20_new.png', width = 10)

sorted_markers <- as.data.frame(sorted_markers)
write.csv(sorted_markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_B/B_markers_G20_new_1.csv')


##G26

markers_all <- milo.markers[which(milo.markers$adj.P.Val_26 < 0.01 & milo.markers$logFC_26 > 0),]
sorted_markers <- markers_all[order(markers_all$adj.P.Val_26), ]
top_30_genes <- head(sorted_markers$GeneID, 30)
markers <- top_30_genes

milo <- calcNhoodExpression(milo, subset.row=markers)

plotNhoodExpressionGroups(milo, milo.res, features=markers,
                          subset.nhoods = (milo.res$CellType %in% c('B_cell')) & (milo.res$NhoodGroup %in% c('20', '26', '4')),
                          scale=TRUE,
                          grid.space = "fixed")


ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_B/B_G26vsG26&G4_G26_new.png', width = 10)

sorted_markers <- as.data.frame(sorted_markers)
write.csv(sorted_markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_B/B_markers_G26_new_1.csv')


##B cell Parity

milo <- result_list[['lymphoid']][['Parity']][[1]]
milo.res <- result_list[['lymphoid']][['Parity']][[2]]

set.seed(123)
milo.res <- groupNhoods(milo, milo.res, max.lfc.delta = 1)

plotDAbeeswarm(milo.res , group.by = "NhoodGroup")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Parity_B/nhood_groups_new_bee.png')
plotNhoodGroups(milo, milo.res, layout="UMAP")
ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Parity_B/nhood_groups_new1.png')

# milo.res1 <- groupNhoods(milo, milo.res, max.lfc.delta = 0.05, subset.nhoods = (milo.res$CellType %in% c('Lp')))
# plotNhoodGroups(milo, milo.res1, layout="UMAP")
# ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_Lp/nhood_groups_checklp.png')


milo.res_reg = milo.res[which(milo.res$CellType%in% c('B_cell')),]

keep.rows <- rowSums(logcounts(milo)) != 0
milo <- milo[keep.rows,]
    
dec <- modelGeneVar(milo)
hvgs <- getTopHVGs(dec, n = 2000)

milo.markers <- findNhoodGroupMarkers(milo, milo.res,
                                    subset.row = hvgs,
                                    subset.nhoods = (milo.res$CellType %in% c('B_cell') & (milo.res$NhoodGroup %in% c('14', '4'))),
                                    aggregate.samples = TRUE,
                                    sample_col = 'Sample',
                                    gene.offset = FALSE)

write.csv(milo.markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Parity_B/B_markers_new_1.csv')

markers_all <- milo.markers[which(milo.markers$adj.P.Val_14 < 0.01 & milo.markers$logFC_14 > 0),]
sorted_markers <- markers_all[order(markers_all$adj.P.Val_14), ]
top_30_genes <- head(sorted_markers$GeneID, 30)
markers <- top_30_genes

milo <- calcNhoodExpression(milo, subset.row=markers)

plotNhoodExpressionGroups(milo, milo.res, features=markers,
                          subset.nhoods = (milo.res$CellType %in% c('B_cell')) & (milo.res$NhoodGroup %in% c('14','4')),
                          scale=TRUE,
                          grid.space = "fixed")


ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Parity_B/B_G14vsG4_G14_new.png', width = 10)

sorted_markers <- as.data.frame(sorted_markers)
write.csv(sorted_markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Parity_B/B_markers_G14_new_1.csv')


##G26

markers_all <- milo.markers[which(milo.markers$adj.P.Val_26 < 0.01 & milo.markers$logFC_26 > 0),]
sorted_markers <- markers_all[order(markers_all$adj.P.Val_26), ]
top_30_genes <- head(sorted_markers$GeneID, 30)
markers <- top_30_genes

milo <- calcNhoodExpression(milo, subset.row=markers)

plotNhoodExpressionGroups(milo, milo.res, features=markers,
                          subset.nhoods = (milo.res$CellType %in% c('B_cell')) & (milo.res$NhoodGroup %in% c('20', '26', '4')),
                          scale=TRUE,
                          grid.space = "fixed")


ggsave('/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_B/B_G26vsG26&G4_G26_new.png', width = 10)

sorted_markers <- as.data.frame(sorted_markers)
write.csv(sorted_markers, file = '/home/rx238/rds/hpc-work/aging/report_graph/milo/neighbourhood/Age_B/B_markers_G26_new_1.csv')




### Milo obj

object <- readRDS('/home/rx238/rds/hpc-work/aging/report_graph/milo/milo_list_result.rds')

object <- result_list

parity_res_total <- c()

for (cell in names(object)){
  parity_res <- object[[cell]][['Parity']][[2]]
  parity_res$ct <- cell
  parity_res_total <- rbind(parity_res_total, parity_res)
}

age_res_total <- c()  

for (cell in names(object)){
  age_res <- object[[cell]][['Age']][[2]]
  age_res$ct <- cell
  age_res_total <- rbind(age_res_total, age_res)
}

genotype_res_total <- c()

for (cell in names(object)){
  genotype_res <- object[[cell]][['Genotype']][[2]]
  genotype_res$ct <- cell
  genotype_res_total <- rbind(genotype_res_total, genotype_res)
}

graph_res <- function(milo.res, test){
  
  
  milo.res$Diff <- sign(milo.res$logFC)
  milo.res$Diff[milo.res$SpatialFDR > 0.1] <- 0
    


  milo.res_new = milo.res[which(milo.res$CellType != "Mixed"),]
  milo_sorted <- milo.res_new[order(milo.res_new$ct, decreasing = TRUE), ]

  print(unique(milo_sorted$CellType))
    
  
  if (!any(milo_sorted$Diff == 1)){
      plotDAbeeswarm1(milo_sorted, group.by = "CellType")
      ggsave(paste0("/home/rx238/rds/hpc-work/aging/report_graph/milo/sum/", test, '_bee_DA_new.png'), height=15, width=8)
    }
    
  else{
      plotDAbeeswarm(milo_sorted, group.by = "CellType")
      ggsave(paste0("/home/rx238/rds/hpc-work/aging/report_graph/milo/sum/", test, '_bee_DA_new.png'), height=15, width=8)
    }
  
}

graph_res(parity_res_total, 'parity_90')

graph_res(age_res_total, 'age')

graph_res(genotype_res_total, 'genotype_60')