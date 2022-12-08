#' @title Single Cell Drug Hunter
#' @description Calculate single cell Drug Gene Score(scDGS) for candidate genes' drug-gene pairs based on DGIdb 2021.5 database
#' @details Input a Seurat format single-cell expression matrix, cell type annotation column names, candidate gene list, correlation lists of genes and phenotypes, significance pvalue for genes
#' @param seurat_obj A Seurat format single-cell expression matrix
#' @param gene_cor Correlation lists of genes and phenotypes with colnames 'gene' and 'cor'
#' @param gene_pvalue P-value lists of genes and phenotypes with colnames 'gene' and 'pvalue'
#' @param candidate Candidate gene list
#' @param annotation Cell type annotation column names in seurat meta.data part
#' @return A list with scDGS and the rank of specificity score, correlation, drug-gene interaction and pvalue for candidate drug-gene pairs
#' @export
#'
#' @import dplyr
#' @import tidyverse
#' @author Yijun Zhou

scDGS <- function(seurat_obj, gene_cor, gene_pvalue,candidate,annotation){

  #############################
  ###1.Load DGIdb database
  ############################
  DGIdb_interaction <- readr::read_table(system.file('extdata','DGIdb_interaction.tsv',package = 'scDrugHunter'))
  gene_drug_pairs <- filter(DGIdb_interaction, !is.na(interaction_group_score)) %>% select(gene_name,drug_name,interaction_group_score)
  colnames(gene_drug_pairs) <- c('gene','drug_name','interaction_group_score')

  #############################
  ###2.Calcualate candidate genes specificity.
  ############################
  candidate_cor <- subset(gene_cor, gene %in% candidate)
  cts <- Candidate_genes_specificity(exp = seurat_obj,candidate = candidate_cor, annotation = annotation,correlation = 'cor')

  #############################
  ###3.Link candaidate gene and DGIdb database.
  ############################
  pairs <- inner_join(cts,gene_drug_pairs,by = c('gene'))
  pairs$pvalue <- gene_pvalue$pvalue[match(pairs$gene,gene_pvalue$gene_name)]

  #############################
  ###4.Calcualate drug-gene interaction rank, genes-trait significance rank and scDGS.
  ############################
  pairs <- pairs %>% group_by(!!ensym(annotation)) %>%
    mutate(pair_name = paste(gene,drug_name,sep = '_')) %>%
    distinct() %>%
    mutate(interaction_rank=rank(desc(interaction_group_score))*10/length(interaction_group_score)) %>%
    mutate(pvalue_rank=rank(desc(pvalue))*10/length(pvalue)) %>%
    mutate(scDGS = 0.5*(specificity_rank*pvalue_rank+specificity_rank*interaction_rank+correlation_rank*interaction_rank+correlation_rank*pvalue_rank)) %>%
    select(gene, drug_name, pair_name, specificity_rank, correlation_rank, interaction_rank, pvalue_rank,scDGS) %>%
    ungroup()
  pairs <- pairs[order(-pairs$radar_score),]

  return(pairs)

}
