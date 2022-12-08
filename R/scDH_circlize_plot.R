#' @title scDH_circlize_plot
#' @description Plot the top drug-gene pairs' scDGS network in specific cell type
#' @details Input scDH format results
#' @param DG_pairs A list of scDH results
#' @param annotation Cell type annotation column names in scDH results
#' @param cell_type Candidate cell type
#' @param top_g Number of top genes to be ploted in circlize
#' @param top_d Number of top drugs for gene specific drug_gene pairs to be ploted in circlize
#' @return A circlize chart for top drug-gene pairs in candidate cell type
#' @export
#'
#' @import dplyr
#' @import tidyverse
#' @import circlize
#' @author Yijun Zhou

scDGS_circlize_plot <- function(DG_pairs, annotation, cell_type, top_g, top_d) {
  sub_pairs <- DG_pairs %>% subset(!!ensym(annotation) == cell_type)
  top_genes <- unique(sub_pairs$gene) %>% head(n = top_g)
  top_pairs <- subset(sub_pairs, gene %in% top_genes) %>%
    group_by(gene) %>%
    do(head(., n = top_d)) %>%
    ungroup() %>%
    arrange(-scDGS) %>%
    select(gene, drug_name, radar_score)
  top_pairs <- top_pairs %>% mutate(scDGS = (scDGS - 1.01 * min(scDGS) + 0.01 * max(scDGS)) / (1.01 * (max(scDGS) - min(scDGS))))


  circos.par("track.height" = 0.1)
  df_col <- c(rand_color(top_g), rep(c("grey"), each = length(unique(top_pairs$drug_name))))
  names(df_col) <- c(unique(top_pairs$gene), unique(top_pairs$drug_name))


  circos.par(start.degree = 90, clock.wise = FALSE)
  chordDiagram(top_pairs,
    grid.col = df_col,
    annotationTrack = "grid",
    annotationTrackHeight = c(0.03, 0.01)
  )


  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
      facing = "clockwise", niceFacing = TRUE, adj = c(-0.4, 0)
    )
  }, bg.border = NA)
}
