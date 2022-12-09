#' @title Candidate genes specificity
#' @description Calculate and sequence the percentage of candidate gene expression in each cell type
#' @details Input a Seurat format single-cell expression matrix, cell type annotation column names, candidate gene list, and correlation lists of genes and phenotypes
#' @param exp A Seurat format single-cell expression matrix
#' @param candidate Candidate gene and trait-relavent correlation list with colnames 'gene' and 'cor'
#' @return A list with specificity percent and the rank of specificity score and correlation for candidate genes
#' @export
#'
#' @import dplyr
#' @import tidyverse
#' @author Yijun Zhou
Candidate_genes_specificity <- function(exp, candidate) {
  exp1 <- as.data.frame(exp@assays[["RNA"]]@data)
  exp1$gene <- rownames(exp1)
  exp1 <- subset(exp1, gene %in% candidate$gene)
  #############################
  ### 1.Only keep genes with a unique name and tidy data.
  ############################
  exp1 <- exp1 %>%
    add_count(gene) %>%
    filter(n == 1) %>%
    select(-n) %>%
    gather(key = column, value = expr, -gene) %>%
    as_tibble()
  #############################
  ### 2.Add cell type names.
  ############################
  cell_types <- as.data.frame(exp@meta.data) %>% select(Cell_type)
  cell_types$column <- rownames(cell_types)

  exp1 <- inner_join(exp1, cell_types, by = c("column")) %>% ungroup()

  #############################
  ### 3.Remove not expressed genes and calculate gene expression in each cell type
  ############################
  not_expressed <- exp1 %>%
    group_by(gene) %>%
    summarise(total_sum = sum(expr)) %>%
    filter(total_sum == 0) %>%
    select(gene) %>%
    unique()

  exp1 <- filter(exp1, !gene %in% not_expressed$gene)
  exp1 <- exp1 %>%
    group_by(Cell_type, gene) %>%
    summarise(expr_sum_mean = mean(expr))

  #############################
  ### 4.Each cell type is scaled to the same total number of molecules.
  ############################
  exp1 <- exp1 %>%
    group_by(Cell_type) %>%
    mutate(expr_sum_mean = expr_sum_mean * 1e6 / sum(expr_sum_mean))

  ##############################
  ### 5.Specificity calculation
  #####################################
  # The specifitiy is defined as the proportion of total exp1ression performed by the cell type of interest (x/sum(x)).
  exp1 <- exp1 %>%
    group_by(gene) %>%
    mutate(specificity = expr_sum_mean / sum(expr_sum_mean)) %>%
    ungroup() # ungroup():取消分组

  exp2 <- inner_join(exp1, candidate, by = c("gene"))

  ##############################
  ### 6.Percent rank calculation of specificity score and gene correlation score in each cell type
  #####################################
  exp2 <- exp2 %>%
    mutate(specificity_rank = rank(-desc(specificity)) * 10 / length(specificity))

  exp2 <- exp2 %>%
    group_by(Cell_type) %>%
    mutate(correlation_rank = rank(-desc(cor)) * 10 / length(cor)) %>%
    ungroup()

  return(exp2)
}
