#' @title Single Cell Drug Hunter
#' @description Calculate single cell Drug Gene Score(scDGS) for candidate genes' drug-gene pairs based on DGIdb 2021.5 database
#' @details Input a Seurat format single-cell expression matrix, cell type annotation column names, candidate gene list, correlation lists of genes and phenotypes, significance pvalue for genes
#' @param seurat_obj A Seurat format single-cell expression matrix with 'Cell_type' colname in meta.data
#' @param gene_cor Correlation lists of genes and phenotypes with colnames 'gene' and 'cor'
#' @param gene_pvalue P-value lists of genes and phenotypes with colnames 'gene' and 'pvalue'
#' @param candidate Candidate gene list
#' @return A list with scDGS and the rank of specificity score, correlation, drug-gene interaction and pvalue for candidate drug-gene pairs
#' @export
#'
#' @import dplyr
#' @import tidyverse
#' @author Yijun Zhou

scDGS <- function(seurat_obj, gene_cor, gene_pvalue,candidate){

  #############################
  ###1.Load DGIdb database
  ############################
  DGIdb_interaction <- readr::read_tsv(system.file('extdata','DGIdb_interactions.tsv',package = 'scDrugHunter'))
  gene_drug_pairs <- filter(DGIdb_interaction, !is.na(interaction_group_score)) %>% select(gene_name,drug_name,interaction_group_score)
  colnames(gene_drug_pairs) <- c('gene','drug_name','interaction_group_score')

  #############################
  ###2.Calcualate candidate genes specificity.
  ############################
  candidate_cor <- subset(gene_cor, gene %in% candidate)
  cts <- Candidate_genes_specificity(exp = seurat_obj,candidate = candidate_cor)
  gc()

  #############################
  ###3.Link candaidate gene and DGIdb database.
  ############################
  pairs <- inner_join(cts,gene_drug_pairs,by = c('gene'))
  pairs$pvalue <- gene_pvalue$pvalue[match(pairs$gene,gene_pvalue$gene)]

  #############################
  ###4.Calcualate drug-gene interaction rank, genes-trait significance rank and scDGS.
  ############################
  pairs <- pairs %>% group_by(Cell_type) %>%
    mutate(pair_name = paste(gene,drug_name,sep = '_')) %>%
    distinct() %>%
    mutate(interaction_rank=rank(desc(interaction_group_score))*10/length(interaction_group_score)) %>%
    mutate(pvalue_rank=rank(desc(pvalue))*10/length(pvalue)) %>%
    mutate(scDGS = 0.5*(specificity_rank*pvalue_rank+specificity_rank*interaction_rank+correlation_rank*interaction_rank+correlation_rank*pvalue_rank)) %>%
    select(gene, drug_name, pair_name, specificity_rank, correlation_rank, interaction_rank, pvalue_rank,scDGS) %>%
    ungroup()
  pairs <- pairs[order(-pairs$scDGS),]

  return(pairs)

}

