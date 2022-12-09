#' @title scDH_radar_plot
#' @description Plot the top n drug-gene pairs' rank in specific cell type
#' @details Input scDGS format results
#' @param DG_pairs A list of scDH results
#' @param annotation Cell type annotation column names in scDH results
#' @param cell_type Candidate cell type
#' @param top Number of top pairs to be ploted in radar chart
#' @return A radar chart for top n drug-gene pairs in candidate cell type
#' @export
#'
#' @import dplyr
#' @import tidyverse
#' @import ggradar
#' @author Yijun Zhou

scDH_radar_plot <- function(DG_pairs,  cell_type, top,unique_g) {
  if(unique_g){
  sub_pairs <- DG_pairs %>%
    subset(Cell_type == cell_type) %>%
    group_by(gene) %>%
    do(head(., n = 1)) %>%
    ungroup() %>%
    arrange(-scDGS)  %>%
    select(pair_name, scDGS) %>%
    head(n = top)

  }else{
  sub_pairs <- DG_pairs %>%
    subset(Cell_type == cell_type) %>%
    select(pair_name, scDGS) %>%
    head(n = top)
  }
  sub_pairs[, -1] <- sapply(sub_pairs[, -1], as.numeric)
  rank_pairs <- sub_pairs %>%
    mutate(scDGS_rank = rank(scDGS)) %>%
    select(pair_name, scDGS_rank)
  top_pairs <- rank_pairs$pair_name
  rank_pairs <- as.data.frame(t(rank_pairs[, 2]))
  colnames(rank_pairs) <- top_pairs
  v1_name <- glue::glue("{cell_type}_top{top}_scDGS")
  rank_pairs <- cbind(v1_name, rank_pairs)
  p <- ggradar(rank_pairs,
    grid.min = 0,
    grid.mid = 0.5 * top,
    grid.max = top
  )

  return(p)
}
