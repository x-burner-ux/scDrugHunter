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

scDH_radar_plot <- function(DG_pairs, annotation, cell_type, top) {
  sub_pairs <- DG_pairs %>%
    subset(!!ensym(annotation) == cell_type) %>%
    select(pair_name, scDGS) %>%
    head(n = top)
  sub_pairs[, -1] <- sapply(sub_pairs[, -1], as.numeric)
  rank_pairs <- sub_pairs %>%
    mutate(scDGS_rank = rank(scDGS)) %>%
    select(pair_name, scDGS_rank)
  top_pairs <- rank_pairs$pair_name
  top_pairs <- as.data.frame(t(top_pairs[, 2]))
  colnames(top_pairs) <- top_pairs
  v1_name <- glue::glue("{cell_type}_top{top}_scDGS")
  top_pairs <- cbind(v1_name, top_pairs)
  p <- ggradar(top_pairs,
    grid.min = 0,
    grid.mid = 0.5 * top,
    grid.max = top
  )

  return(p)
}
