library(ggrepel)
library(ggplot2)
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
library(edgeR)
library(EnhancedVolcano)


final <- readRDS('/home/rx238/rds/hpc-work/aging/final_data/final_final/adata_all_reannotation_final1.rds')



final <- final[, !final$ct_level3 %in% c("Myocytes", "Erythroid", "Tumour", "Caf", "Fb_proliferating", "Mo_TAM", "Fb3 Uni2_2", 'Hs Uni2_1', 'Hs Nul2_2', 'Bsl_proliferating')]

final <- final[,(final$Batch != 7) & (!final$Sample %in% c('ctrl2', 'ctrl1', 'NulControl'))]


final$tage <- (final$Age - min(final$Age)) / (max(final$Age) - min(final$Age))
final$parity_bool <- ifelse(final$Parity == 'NP', 0, 1)
final$tumour_bool <- ifelse(final$Genotype == 'WT', 0, 1)
final$ct_level3[final$ct_level3 == 'Avd'] <- 'Lp'


names(assays(final)) <- c('X', 'counts', 'raw')

### subset

age <- final[, (final$Genotype == "WT") & (final$Parity == 'NP')]

wkbr <- final[, (((final$Genotype == 'WKBR') & (final$Parity == 'NP'))| ((final$Genotype == "WT") & (final$Parity == 'NP')))]

parity <- final[, (final$Parity %in% c('NP', 'Parous')) & (final$Age > 80) & (final$Genotype == 'WT')]


#### use old annotation






###Function

set_up_DGE(age, 'tage')
set_up_DGE(parity, 'parity_bool')
set_up_DGE(wkbr, 'tumour_bool')

#Age_(Parity)


###Parity (Age)



set_up_DGE <- function(data, test_var){
    dge_pwd <- '/home/rx238/rds/hpc-work/aging/report_graph/dge'
    dir_path <- paste0(dge_pwd, '/', test_var)
    celltypes <- unique(data$ct_level3)
    result_list <- list()
    for (cell in celltypes){
        cell_name <- gsub("/", ";", cell)
        print(cell)
        dir_path <- paste0(dge_pwd, '/', test_var, '/', cell_name)
        if (!file.exists(dir_path)){
        dir.create(dir_path, recursive = TRUE)
        }
        data1 <- data[, data$ct_level3 == cell]
        DEGs <- DGE_func(data1, cell, test_var, dir_path)
        result_list[[cell]] <- DEGs
        if (length(DEGs) != 0){
            plot_df(DEGs, cell, test_var, FCcutoff = 1, maxoverlapsConnectors = 20, dir_path)
        }
        
    }
    return(result_list)
}


DGE_func <- function(data, label, test_var, dge_pwd){

    label_name <- gsub("/", ";", label)

    print(paste0('Running DGE for: ', label))

    summed <- aggregateAcrossCells(data, id = colData(data)[, c("ct_level3", "Sample")])
    summed.filt <- summed[,summed$ncells >= 5]
   
    current <- summed[,label==summed$ct_level3]

    y <- DGEList(counts(current), samples=colData(current))

    if (dim(current)[2] == 0){
    print('Too few cells of this type')
    return()}

    discarded <- current$ncells < 5
    y <- y[,!discarded]

    y$samples <- droplevels(y$samples)


    print('creating design matrix')

    factors <- list(
        tumour_bool = factor(y$samples$tumour_bool),
        parity_bool = factor(y$samples$parity_bool),
        tage = y$samples$tage
    )

    factors <- factors[setdiff(names(factors), test_var)]
    if (test_var == 'parity_bool'){
        factors <- factors[setdiff(names(factors), 'tage')]
    }

    formula_terms <- c()
    for (factor_name in names(factors)) {
        factor_values <- factors[[factor_name]]
        if (length(unique(factor_values)) > 1) {
            formula_terms <- c(formula_terms, factor_name)
        } else {
            print(paste0("Excluding unique factor:", factor_name))
        }
    }

    print(test_var)
    if (length(unique(y$samples[[test_var]])) > 1) {
        formula_terms <- c(formula_terms, test_var)
    } else {
        print(paste("test variable is unique end:", test_var))
        return()
    }

    formula <- as.formula(paste("~", paste(formula_terms, collapse = " + ")))
    print(formula)
    design <- model.matrix(formula, y$samples)
    
    print("Design matrix:")
    print(design)

    print('filtering by expression')
    keep <- filterByExpr(y, design=design)
  
    y <- y[keep,]
    y <- calcNormFactors(y)
  
    y <- estimateDisp(y, design)
    if (is.na(y$common.dispersion)){
        print('Common dispersion is NA return')
        return()
    }
    fit <- glmQLFit(y, design, robust=TRUE)
    res <- glmQLFTest(fit, coef=ncol(design))

    res$table$FDR <- p.adjust(res$table$PValue, method="BH")

    DEGs <- res$table[order(res$table$FDR),]
    DEGs_up <- DEGs[DEGs$logFC > 0,]
    DEGs_down <- DEGs[DEGs$logFC < 0,]
    print(paste0('Number of DEGs: ', sum(DEGs$FDR < 0.1))) 
  
    dir.create(dge_pwd, showWarnings = FALSE, recursive = TRUE)

    write.csv(DEGs, file = paste0(dge_pwd, '/dge-', label_name, '_',test_var,'-all.csv'))
    write.csv(DEGs_up, file = paste0(dge_pwd, '/dge-', label_name, '_',test_var,'-up.csv'))
    write.csv(DEGs_down, file = paste0(dge_pwd, '/dge-', label_name, '_',test_var,'-down.csv'))
    return(DEGs)
}





df <- read.csv('/home/rx238/rds/hpc-work/aging/report_graph/dge/parity_bool_80/Lp/dge-Lp_parity_bool-all.csv')

plot_df(df, 'Lp', 'parity_bool', FCcutoff = 1, maxoverlapsConnectors = 20, '/home/rx238/rds/hpc-work/aging/report_graph/dge/parity_bool_80/Lp')


plot_df <- function(DEGs, label, test_var, FCcutoff = 1, maxoverlapsConnectors = 20, dge_pwd){

    rownames(DEGs) <- DEGs$X

    label_name <- gsub("/", ";", label)
    x_buffer <- 1.5

    DEGs <- as.data.frame(DEGs)
    if (test_var == 'donor_age'){
        FCcutoff <- 0.01
        x_buffer <- 0.1
    }
    
    ttl <- paste0(label, ' ', test_var)

    DEGs <- DEGs[order(DEGs$FDR),]

    top_20_genes <- head(DEGs, 20)

    top_20_cd1d <- c(top_20_genes$X, 'Cd1d1')

    Volcano <- EnhancedVolcano(DEGs,
                           lab = rownames(DEGs),
                           x = 'logFC',
                           y = 'FDR',
                           xlim = c(min(DEGs[['logFC']], na.rm = TRUE) - x_buffer,
                                    max(DEGs[['logFC']], na.rm = TRUE) + x_buffer),
                           title = ttl,
                           subtitle = '',
                           subtitleLabSize = 3,
                           legendPosition = "bottom",
                           pointSize = 1.0,
                           labSize = 3,
                           labFace = 'bold',
                           labCol = 'red2',
                           FCcutoff = FCcutoff,
                           pCutoff = 0.05,
                           col = c("grey30", "forestgreen", "royalblue", "purple"),
                           drawConnectors = TRUE,
                           typeConnectors = 'open',
                           maxoverlapsConnectors = maxoverlapsConnectors,
                           selectLab = top_20_cd1d,
                           )
   
# Print the plot
    pdf(paste0(dge_pwd, '/dge-', label_name, '_',test_var,'-volcano_new.pdf'))
    print(Volcano)
    dev.off() 
}


