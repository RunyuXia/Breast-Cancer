library(nlme)
library(lme4)
library(SingleCellExperiment)
library(ggplot2)
library(dplyr)
library(broom)
library(tidyverse)
library(broom.mixed)
library(lmerTest)
library(ggpubr)
library(rstatix)


immune_receptors_inhi <- c('Pdcd1', 'Ctla4','Lag3', 'Tigit', 'Cd96', 'Havcr2', 'Cd244a', 'Cd160', 'Klra9')
immune_receptors_act <- c('Cd226', 'Cd28')
immune_ligands <- c('Cd274', 'Cd80', 'Cd86', 'Clec4g', 'Fgl1', 'Nectin2', 'Pvr', 'Tnfrsf14', 'Ceacam1', 'Lgals9', 'Hmgb1', 'Ptdss1', 'Ptdss2', 'Pdcd1lg2', 'Vtcn1')


all_markers <- c(immune_receptors_inhi, immune_receptors_act, immune_ligands)
receptor_markers <- c(immune_receptors_inhi, immune_receptors_act)
ligand_markers <- immune_ligands


milo <- readRDS('/home/rx238/rds/hpc-work/aging/final_data/final_final/adata_all_reannotation_final1.rds')


milo <- milo[, (!milo$ct_level3 %in% c('Hs Nul2_2', 'Hs Uni2_1', 'Mo_TAM', 'Fb3 Uni2_2', 'Myocytes', 'Erythroid', 'Caf', 'Tumour', 'Bsl_proliferating', 'Proliferating Fibroblast', 'ECAT', 'Fb_proliferating', 'cycling_Treg')) & (!milo$Batch == 7) & (!milo$Sample %in% c('ctrl1', 'ctrl2', 'NulControl'))]
milo <- milo[, milo$Parity %in% c('NP', 'Parous')]
milo <- milo[, milo$Genotype %in% c('WT', 'WKBR')]
milo$Index <- colnames(milo)

# old_meta <- read_csv('/home/rx238/rds/hpc-work/aging/report_graph/markers/old_immune_meta.csv')
# old_meta <- old_meta[old_meta$ct_level1 == 'Lymphoid',]
# match_indices <- match(milo$Index, old_meta$Index)
# milo$ct_level3 <- as.character(milo$ct_level3)
# old_meta$ct_level3 <- as.character(old_meta$ct_level3)
# milo$ct_level3_new <- ifelse(!is.na(match_indices), old_meta$ct_level3[match_indices], milo$ct_level3)
# milo$ct_level3 <- milo$ct_level3_new


colData(milo) <- droplevels(colData(milo))
col_data <- as.data.frame(colData(milo))
col_data$ct_level2 <- ''

col_data <- col_data %>%
  mutate(ct_level2 = if_else(ct_level3 %in% c('cDC1', 'cDC2'), "cDC", ct_level2))

col_data <- col_data %>%
  mutate(ct_level2 = if_else(ct_level3 %in% c('CD8_Tem','CD8_Trm', 'CD8_NKT', 'CD8_naive'), "CD8", ct_level2))

col_data <- col_data %>%
  mutate(ct_level2 = if_else(ct_level3 %in% c('CD4_Th1/17','CD4_Th2', 'CD4_Treg', 'CD4_Th1', 'CD4_naive'), "CD4", ct_level2))

col_data <- col_data %>%
  mutate(ct_level2 = if_else(ct_level3 %in% c('Avd', 'Lp'), "Lp_Avd", ct_level2))

col_data <- col_data %>%
  mutate(ct_level2 = if_else(ct_level3 %in% c('Mo2', 'Mo1'), "Macrophage", ct_level2))

col_data <- col_data %>%
  mutate(ct_level2 = if_else(ct_level3 %in% c('LEC1', 'LEC2'), "LEC", ct_level2))

col_data <- col_data %>%
  mutate(ct_level2 = if_else(ct_level3 %in% c('VM1', 'VM2'), "VM", ct_level2))


Not_immune <- col_data[col_data$ct_level1 %in% c('Epithelial', 'Fibroblast', 'Stroma'), ]

ligand_cts <- as.character(unique(Not_immune$ct_level3))

ligand_cts <- c(ligand_cts, 'Lp_Avd', 'LEC', 'VM')

Immune_lym <- col_data[col_data$ct_level1 %in% c('Lymphoid'),]

receptor_cts <- as.character(unique(Immune_lym$ct_level3))

receptor_cts <- c(receptor_cts, 'CD4', 'CD8')

Immune_mye <- col_data[col_data$ct_level1 %in% c('Myeloid'),]

receptor_ligand_cts <- as.character(unique(Immune_mye$ct_level3))

receptor_ligand_cts <- c(receptor_ligand_cts, 'Macrophage', 'cDC')

cts_dict <- list()

celltypes <- unique(col_data$ct_level3)
for (cell in celltypes){
    cts_dict[[cell]] <- cell
}


level2_cts <- unique(col_data$ct_level2)
for (cell in level2_cts){
    if (cell == ''){
        next}
    subset <- col_data[col_data$ct_level2 == cell, ]
    cts_dict[[cell]] <- as.character(unique(subset$ct_level3))
}



extract_data <- function(data, cts, gene_list){
    print('start')
    data <- data[, data$ct_level3 %in% cts]
    print('getsample')
    expression_matrix <- assay(data, "lognorm")  # or "logcounts" depending on your data
    expression_subset <- expression_matrix[rownames(expression_matrix) %in% gene_list, ]
    expression_subset_matrix <- as.matrix(expression_subset)
    expression_df <- as.data.frame(t(expression_subset_matrix))
    rownames(expression_df) <- colnames(data)
    colnames(expression_df) <- rownames(expression_subset)
    meta <- colData(data)
    meta <- meta[, c('Parity', 'Genotype', 'Sample', 'Age', 'Batch')]
    expression_df <- cbind(expression_df, meta)
    mean_expression_df <- expression_df %>% 
    group_by(Sample) %>% 
    summarise(across(all_of(gene_list), mean),
        across(everything(), ~ if(is.numeric(.)) mean(.) else first(.))) 
    print(mean_expression_df)
    mean_expression_df$Genotype <- as.factor(mean_expression_df$Genotype)
    mean_expression_df$Parity <- as.factor(mean_expression_df$Parity)
    mean_expression_df$Age <- as.numeric(mean_expression_df$Age)
    return(mean_expression_df)
}


models <- function(df, cell, gene_list, test){
  cell_name <- gsub("/", ";", cell)
  new_dir2 <- paste0('/home/rx238/rds/hpc-work/aging/report_graph/markers/new_anno/', cell_name, '/',test)
  if (!dir.exists(new_dir2)){dir.create(new_dir2, recursive = TRUE)}
  df1 <-  df %>%
    pivot_longer(cols = all_of(gene_list), 
    names_to = "Gene", 
    values_to = "Expression")
  if (test == 'age'){
    print('age')
    df1 <- df1[(df1$Genotype == 'WT'), ]
    df1 <- df1[(df1$Parity == 'NP'),]
    df1 <- droplevels(df1)
    result <- perform_cor_test(df1, cell, test)
    # result <- perform_t_test_age(df1, cell, test)
    } else if (test == 'parity'){
      df1 <- df1[(df1$Genotype == 'WT'), ]
      df1 <- df1[(df1$Age > 80),]
      df1 <- droplevels(df1)
      result <- perform_t_test(df1, cell, test)
    } else if (test == 'risk') {
      df1 <- df1[(df1$Parity == 'NP') & (df1$Age < 80),]
      df1 <- droplevels(df1)
      result <- perform_t_test_wkbr(df1, cell, test)
      # result <- perform_lm_test(df1, cell, test)
    }
  if (length(result) > 0) {
    result$cell <- cell
    write.csv(result, file = paste0(new_dir2, '/', test, 'after_result.csv'))
    return(result)
  } else {
    message("No models were successfully fitted.")
    return(NULL)
  }}


perform_cor_test <- function(data, cell, test) {
    print('cor_test')
   cell_name <- gsub("/", ";", cell)

   new_dir2 <- paste0('/home/rx238/rds/hpc-work/aging/report_graph/markers/new_anno/', cell_name, '/',test)
   if (!dir.exists(new_dir2)){dir.create(new_dir2, recursive = TRUE)}

    print('stats')
   stat.test <- data %>%
   group_by(Gene) %>%
   summarise(
    cor = cor.test(Age, Expression, method = "pearson")$estimate,
    p.value = cor.test(Age, Expression, method = "pearson")$p.value) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"),
    label = paste("R:", round(cor, 2), "P:", format(p.adj, digits = 3)),
    Label = paste0(Gene, ' ', label))
    
    labeller_custom <- function(string){
    label <- stat.test[stat.test$Gene == string,]$Label
    return(label)}
    
    print('plot')
    corr_plt <- ggscatter(data, x = 'Age', y = 'Expression', add = 'reg.line',conf.int = TRUE,
    add.params = list(color = "blue", fill = "lightgray")) +
    facet_wrap(~ Gene, ncol = 3,scales = "free", labeller = labeller(Gene = labeller_custom)) 

    nrow = nrow(stat.test)/3
    
    ggsave(corr_plt, filename = paste0(new_dir2, '/', cell_name,'_', test, 'corrafter.png'), width = 10, height = nrow*2.5)

  return(stat.test)
}


perform_t_test <- function(data, cell, test){
    print('t_test')
   cell_name <- gsub("/", ";", cell)
   new_dir2 <- paste0('/home/rx238/rds/hpc-work/aging/report_graph/markers/new_anno/', cell_name, '/',test)
   if (!dir.exists(new_dir2)){dir.create(new_dir2, recursive = TRUE)}

    print('stats')
   stat.test <- data %>%
   group_by(Gene) %>%
   t_test(Expression ~ Parity) %>%
   adjust_pvalue(method = 'BH') %>%
   add_significance()
    
    stat.test <- stat.test %>% add_x_position(x = "Parity")
    stat.test$p.adj_r <- round(stat.test$p.adj, 3)
    max_expressions <- data %>% 
     group_by(Gene) %>% 
     summarise(y.position = max(Expression, na.rm = TRUE) * 1.1)
    
    stat.test <- left_join(stat.test, max_expressions, by = "Gene")
    print('plot')
    bxp <- ggboxplot(data, x = "Parity", y = "Expression", fill = "Parity", add = "jitter", palette = c("#B3DFD6", "#F4D166"))+
    facet_wrap(~ Gene, ncol = 3, scales = "free")+
    scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.20))) +
    stat_pvalue_manual(stat.test, label = "{p.adj_r} {p.adj.signif}") 

    nrow = nrow(stat.test)/3
    
    ggsave(bxp, filename = paste0(new_dir2, '/', cell_name,'_', test, '80after.png'), width = 10, height = nrow*2.5)
    
    return(stat.test)
}


perform_t_test_wkbr <- function(data, cell, test){
    print('t_test')
   cell_name <- gsub("/", ";", cell)
   new_dir2 <- paste0('/home/rx238/rds/hpc-work/aging/report_graph/markers/new_anno/', cell_name, '/',test)
   if (!dir.exists(new_dir2)){dir.create(new_dir2, recursive = TRUE)}
    
    data <- data %>%
    mutate(Genotype = factor(Genotype, levels = c("WT", "WKBR")))
    print('stats')
        
    stat.test <- tryCatch({
        data %>%
        group_by(Gene) %>%
        t_test(Expression ~ Genotype) %>%
        adjust_pvalue(method = 'BH') %>%
        add_significance()
    }, error = function(e) {
        message("Error in t_test: ", e)
        return(NULL)  # Return NULL or an empty tibble/dataframe
    })
    
    # Check if stat.test is NULL (due to error in t_test)
    if (is.null(stat.test)) {
        return(NULL)  # Or return a default value if you prefer
    }

    stat.test <- stat.test %>% add_x_position(x = "Genotype")
    stat.test$p.adj_r <- round(stat.test$p.adj, 3)
    max_expressions <- data %>% 
     group_by(Gene) %>% 
     summarise(y.position = max(Expression, na.rm = TRUE) * 1.1)
    
    stat.test <- left_join(stat.test, max_expressions, by = "Gene")
    print('plot')
    bxp <- ggboxplot(data, x = "Genotype", y = "Expression", fill = "Genotype", add = "jitter", palette = c('#D5D5D5', '#E9BFDD'))+
    facet_wrap(~ Gene, ncol = 3, scales = "free")+
    scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.20))) +
    stat_pvalue_manual(stat.test, label = "{p.adj_r} {p.adj.signif}") 

    nrow = nrow(stat.test)/3
    
    ggsave(bxp, filename = paste0(new_dir2, '/', cell_name,'_', test, 'risk_80_ttest_after.png'), width = 10, height = nrow*2.5)
    return(stat.test)
}


perform_t_test_age <- function(data, cell, test){
    print('t_test')
   cell_name <- gsub("/", ";", cell)
   new_dir2 <- paste0('/home/rx238/rds/hpc-work/aging/report_graph/markers/new_anno/', cell_name, '/',test)
   if (!dir.exists(new_dir2)){dir.create(new_dir2, recursive = TRUE)}

    data$age_bin <- cut(data$Age, breaks = c(0, 50, Inf), labels = c('young', 'old'))
    
    data <- data %>%
    mutate(age_bin = factor(age_bin, levels = c("young", "old")))
    print('stats')

        
    stat.test <- tryCatch({
        data %>%
            group_by(Gene) %>%
            t_test(Expression ~ age_bin) %>%
            adjust_pvalue(method = 'BH') %>%
            add_significance()
    }, error = function(e) {
        message("Error in t_test: ", e)
        return(NULL)  # Return NULL or an empty tibble/dataframe
    })
    
    # Check if stat.test is NULL (due to error in t_test)
    if (is.null(stat.test)) {
        return(NULL)  # Or return a default value if you prefer
    }

    stat.test <- stat.test %>% add_x_position(x = "age_bin")
    stat.test$p.adj_r <- round(stat.test$p.adj, 3)
    max_expressions <- data %>% 
     group_by(Gene) %>% 
     summarise(y.position = max(Expression, na.rm = TRUE) * 1.1)
    
    stat.test <- left_join(stat.test, max_expressions, by = "Gene")
    print('plot')
    bxp <- ggboxplot(data, x = "age_bin", y = "Expression", fill = "age_bin", add = "jitter", palette = c('#D5D5D5', '#E9BFDD'))+
    facet_wrap(~ Gene, ncol = 3, scales = "free")+
    scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.20))) +
    stat_pvalue_manual(stat.test, label = "{p.adj_r} {p.adj.signif}") 

    nrow = nrow(stat.test)/3
    
    ggsave(bxp, filename = paste0(new_dir2, '/', cell_name,'_', test, 'age_bin_ttest.png'), width = 10, height = nrow*2.5)
    return(stat.test)
}


perform_lm_test <- function(data, cell, test) {
    print('lm_test')
    cell_name <- gsub("/", ";", cell)
    new_dir2 <- paste0('/home/rx238/rds/hpc-work/aging/report_graph/markers/new_anno/', cell_name, '/',test)
    if (!dir.exists(new_dir2)){dir.create(new_dir2, recursive = TRUE)}

    data <- data %>%
    mutate(Genotype = factor(Genotype, levels = c("WT", "WKBR")))
    print('stats')
    stat.test <- data %>%
    group_by(Gene) %>%
    do(model = lm(Expression ~ Age * Genotype, data = .)) %>%
    mutate(tidied = list(tidy(model))) %>%
    unnest(tidied) %>%
    filter(term == "GenotypeWKBR") %>%  # Now we're interested in the WKBR effect
    select(Gene, term, estimate, std.error, statistic, p.value) %>%
    adjust_pvalue(p.col = "p.value", method = "BH") %>%
    add_significance()

    labeller_custom <- function(string) {
    adj <- stat.test[stat.test$Gene == string,]$p.value.adj
    label <- paste0(string, ' ', "P:", format(adj, digits = 3))
    return(label)}
    
    print('plot')
    plt <- ggplot(data, aes(x = Age, y = Expression, color = Genotype)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ Gene, ncol = 3, scales = "free", labeller = labeller(Gene = labeller_custom)) +
    scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.20)))

    nrow = nrow(stat.test)/3

    ggsave(plt, filename = paste0(new_dir2, '/', cell_name, '_', test, 'regression.png'), width = 10, height = nrow*2.5)

    stat.test <- stat.test %>%
        mutate(
            p.adj_r = round(p.value.adj, 3),
            group1 = "WT",
            group2 = "WKBR",
        )
    
    max_expressions <- data %>% 
     group_by(Gene) %>% 
     summarise(y.position = max(Expression, na.rm = TRUE) * 1.1)
    
    stat.test <- left_join(stat.test, max_expressions, by = "Gene")


    plt <- ggboxplot(data, x = "Genotype", y = "Expression", fill = "Genotype", add = "jitter", palette = c('#D5D5D5', '#E9BFDD'))+
    facet_wrap(~ Gene, ncol = 3, scales = "free")+
    scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.20))) +
    stat_pvalue_manual(stat.test, label = "{p.adj_r} {p.value.adj.signif}") 

    ggsave(plt, filename = paste0(new_dir2, '/', cell_name, '_', test, 'boxplot_lm_80.png'), width = 10, height = nrow*2.5)

    return(stat.test)
}
    
process_data <- function(data, cts_dict){
  age_result <- list()
  parity_result <- list()
  risk_result <- list()
  for (key in names(cts_dict)){
    cts <- cts_dict[[key]]

    if (all(cts %in% receptor_cts)){
      marker <- receptor_markers
    } else if (all(cts %in% ligand_cts)){
      marker <- ligand_markers
    } else if (all(cts %in% receptor_ligand_cts)){
      marker <- all_markers
    } else {
      marker <- cts
    }
    marker <- intersect(marker, rownames(assay(data, "lognorm")))
    df <- extract_data(data, cts, marker)
    print(paste('process_age', key))
    age_result1  <- models(df, key, marker, 'age')
    print(paste('process_parity', key))
    parity_result1 <- models(df, key, marker, 'parity')
    print(paste('process_risk', key))
    risk_result1 <- models(df, key, marker, 'risk')
    age_result[[key]] <- age_result1
    parity_result[[key]] <- parity_result1
    risk_result[[key]] <- risk_result1
  }
  return(list(age_result, parity_result, risk_result))
}

result <- process_data(milo, cts_dict)

age_result_all <- result[[1]]
parity_result_all <- result[[2]]
risk_result_all <- result[[3]]

age_df_list <- list()
parity_df_list <- list()
risk_df_list <- list()


for (i in 1:length(age_result_all)){
    ind_df <- age_result_all[[i]]
    age_df_list <- c(age_df_list, list(ind_df))
}

for (i in 1:length(parity_result_all)){
    ind_df <- parity_result_all[[i]]
    parity_df_list <- c(parity_df_list, list(ind_df))
}

for (i in 1:length(risk_result_all)){
    ind_df <- risk_result_all[[i]]
    risk_df_list <- c(risk_df_list, list(ind_df))
}

age_result_df <- do.call(rbind, age_df_list)
parity_result_df <- do.call(rbind, parity_df_list)
risk_result_df <- do.call(rbind, risk_df_list)

print('writing')
write.csv(age_result_df, file = '/home/rx238/rds/hpc-work/aging/report_graph/markers/age_result_new_anno_corr.csv')
write.csv(parity_result_df, file = '/home/rx238/rds/hpc-work/aging/report_graph/markers/parity_result_new_anno.csv')
write.csv(risk_result_df, file = '/home/rx238/rds/hpc-work/aging/report_graph/markers/risk_result_new_80_anno_lm.csv')


cts <- cts_dict[['CD8']]

genes <- c('Pdcd1', 'Lag3', 'Ctla4', 'Cd28')



df <- extract_data(milo, cts, genes)

selected_df <- df %>% filter(Genotype == 'WT')

selected_df <- selected_df %>% filter(Age < 20)

selected_df <- selected_df[, c('Sample', 'Cd274', 'Age', 'Parity', 'Genotype', 'Batch')]

library(ggplot2)

pdf('/home/rx238/rds/hpc-work/aging/report_graph/markers/Cd274_parity.pdf')
ggplot(selected_df, aes(x = Parity, y = Cd274)) +
  geom_point() +
  labs(x = "Parity", y = "Expression of Cd274") +
  theme_minimal()
dev.off()



remake_graph <- function(data, ct, gene_list, test, plot_gene_list){

    cts <- cts_dict[[ct]]

    df <- extract_data(data, cts, gene_list)
    models(df, ct, gene_list, test, plot_gene_list)

}


models <- function(df, cell, gene_list, test, plot_gene_list){
  cell_name <- gsub("/", ";", cell)
  new_dir2 <- paste0('/home/rx238/rds/hpc-work/aging/report_graph/markers/new_anno/', cell_name, '/',test)
  if (!dir.exists(new_dir2)){dir.create(new_dir2, recursive = TRUE)}
  df1 <-  df %>%
    pivot_longer(cols = all_of(gene_list), 
    names_to = "Gene", 
    values_to = "Expression")
  if (test == 'age'){
    print('age')
    df1 <- df1[(df1$Genotype == 'WT'), ]
    df1 <- df1[(df1$Parity == 'NP'),]
    df1 <- droplevels(df1)
    result <- perform_cor_test(df1, cell, test, plot_gene_list)
    # result <- perform_t_test_age(df1, cell, test)
    } else if (test == 'parity'){
      df1 <- df1[(df1$Genotype == 'WT'), ]
      df1 <- df1[(df1$Age > 80),]
      df1 <- droplevels(df1)
      result <- perform_t_test(df1, cell, test, plot_gene_list)
    } else if (test == 'risk') {
      df1 <- df1[(df1$Parity == 'NP') & (df1$Age < 80),]
      df1 <- droplevels(df1)
      result <- perform_t_test_wkbr(df1, cell, test, plot_gene_list)
      # result <- perform_lm_test(df1, cell, test)
    }
  if (length(result) > 0) {
    result$cell <- cell
    # write.csv(result, file = paste0(new_dir2, '/', test, '_multi_risk_80_paper.csv'))
    return(result)
  } else {
    message("No models were successfully fitted.")
    return(NULL)
  }}


perform_cor_test_new <- function(data, cell, test, gene_list) {
    print('cor_test')
   cell_name <- gsub("/", ";", cell)

   new_dir2 <- paste0('/home/rx238/rds/hpc-work/aging/report_graph/marker_paper/', cell_name, '/',test)
   if (!dir.exists(new_dir2)){dir.create(new_dir2, recursive = TRUE)}

    print('stats')
   stat.test <- data %>%
   group_by(Gene) %>%
   summarise(
    cor = cor.test(Age, Expression, method = "pearson")$estimate,
    p.value = cor.test(Age, Expression, method = "pearson")$p.value) %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"),
    label = paste("R:", round(cor, 2), "P:", format(p.adj, digits = 3)),
    Label = paste0(Gene, ' ', label))
    
    labeller_custom <- function(string){
    label <- stat.test[stat.test$Gene == string,]$Label
    return(label)}

    data <- data %>% filter(Gene %in% gene_list)
    
    print('plot')
    corr_plt <- ggscatter(data, x = 'Age', y = 'Expression', add = 'reg.line',conf.int = TRUE,
    add.params = list(color = "blue", fill = "lightgray")) +
    facet_wrap(~ Gene, ncol = 3,scales = "free", labeller = labeller(Gene = labeller_custom)) 

    nrow = nrow(stat.test)/3
    
    ggsave(corr_plt, filename = paste0(new_dir2, '/', cell_name,'_', test, 'corr.png'), width = 10, height = nrow*2.5)

  return(stat.test)
}

perform_t_test <- function(data, cell, test, gene_list){
    print('t_test')
   cell_name <- gsub("/", ";", cell)
   new_dir2 <- paste0('/home/rx238/rds/hpc-work/aging/report_graph/marker_paper/', cell_name, '/',test)
   if (!dir.exists(new_dir2)){dir.create(new_dir2, recursive = TRUE)}

    print('stats')
   stat.test <- data %>%
   group_by(Gene) %>%
   t_test(Expression ~ Parity) %>%
   adjust_pvalue(method = 'BH') %>%
   add_significance()
    
    stat.test <- stat.test %>% add_x_position(x = "Parity")
    stat.test$p.adj_r <- round(stat.test$p.adj, 3)
    max_expressions <- data %>% 
     group_by(Gene) %>% 
     summarise(y.position = max(Expression, na.rm = TRUE) * 1.1)
    
    stat.test <- left_join(stat.test, max_expressions, by = "Gene")
    print('plot')

    data <- data %>% filter(Gene %in% gene_list)
    stat.test <- stat.test %>% filter(Gene %in% gene_list)

    bxp <- ggboxplot(data, x = "Parity", y = "Expression", fill = "Parity", add = "jitter", palette = c("#B3DFD6", "#F4D166"))+
    facet_wrap(~ Gene, ncol = 3, scales = "free")+
    scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.20))) +
    stat_pvalue_manual(stat.test, label = "{p.adj_r} {p.adj.signif}") 

    nrow = nrow(stat.test)/3
    
    ggsave(bxp, filename = paste0(new_dir2, '/', cell_name,'_', test, '80.png'), width = 10, height = nrow*2.5)
    
    return(stat.test)
}


perform_t_test_wkbr <- function(data, cell, test, gene_list){
    print('t_test')
   cell_name <- gsub("/", ";", cell)
   new_dir2 <- paste0('/home/rx238/rds/hpc-work/aging/report_graph/marker_paper/', cell_name, '/',test)
   if (!dir.exists(new_dir2)){dir.create(new_dir2, recursive = TRUE)}
    
    data <- data %>%
    mutate(Genotype = factor(Genotype, levels = c("WT", "WKBR")))
    print('stats')
        
    stat.test <- tryCatch({
        data %>%
        group_by(Gene) %>%
        t_test(Expression ~ Genotype) %>%
        adjust_pvalue(method = 'BH') %>%
        add_significance()
    }, error = function(e) {
        message("Error in t_test: ", e)
        return(NULL)  # Return NULL or an empty tibble/dataframe
    })
    
    # Check if stat.test is NULL (due to error in t_test)
    if (is.null(stat.test)) {
        return(NULL)  # Or return a default value if you prefer
    }

    stat.test <- stat.test %>% add_x_position(x = "Genotype")
    stat.test$p.adj_r <- round(stat.test$p.adj, 3)
    max_expressions <- data %>% 
     group_by(Gene) %>% 
     summarise(y.position = max(Expression, na.rm = TRUE) * 1.1)
    
    data <- data %>% filter(Gene %in% gene_list)
    stat.test <- stat.test %>% filter(Gene %in% gene_list)
    
    stat.test <- left_join(stat.test, max_expressions, by = "Gene")
    print('plot')
    bxp <- ggboxplot(data, x = "Genotype", y = "Expression", fill = "Genotype", add = "jitter", palette = c('#D5D5D5', '#E9BFDD'))+
    facet_wrap(~ Gene, ncol = 3, scales = "free")+
    scale_y_continuous(limits = c(0, NA),expand = expansion(mult = c(0, 0.20))) +
    stat_pvalue_manual(stat.test, label = "{p.adj_r} {p.adj.signif}") 

    nrow = nrow(stat.test)/3
    
    ggsave(bxp, filename = paste0(new_dir2, '/', cell_name,'_', test, 'risk_80_ttest.png'), width = 10, height = nrow*2.5, dpi = 300)
    return(stat.test)
}