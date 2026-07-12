# ==============================================================================
# collinearity_analysis.R
#
# Within-compartment transcriptional collinearity from the deconvolution
# reference (phi_cellState).  Motivation: niche pressure can only redistribute
# reads between states whose reference profiles are hard to tell apart, so the
# states that absorb mass should be the collinear attractors of their
# compartment (high mean correlation to their compartment-mates).  This measures
# that directly, on the same phi the Gibbs sampler sees.
#
# Note on the correlation space:
#   transform = "log" (default) correlates log10 profiles -- the conventional
#     expression-similarity readout, not dominated by a few ultra-high genes.
#   transform = "none" correlates raw phi rows -- linear collinearity in the
#     exact space of the multinomial likelihood, i.e. the identifiability-
#     relevant version.  Worth running both; the argument is strongest if the
#     attractors rank top under either.
#
# Usage:
#   source(here("R/collinearity_analysis.R"))
#   phi <- myPrism@phi_cellState@phi            # states x genes
#   collin <- state_collinearity(phi, state_to_compartment,
#                                compartment_levels = names(compartment_to_state))
#   collin$summary                              # per-state within/between means
#
#   plot_collinearity_heatmap(collin, compartment_colors = compartment_colors)
#   plot_collinearity_centrality(
#     collin, compartment_colors = compartment_colors,
#     highlight = c("Oligo_01","Oligo_03","Excit_09","Excit_10","Excit_13","EndoMural_02"))
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a


# ------------------------------------------------------------------------------
# Compute state x state collinearity + per-state within/between summaries
# ------------------------------------------------------------------------------

#' @param phi states x genes reference matrix (e.g. prism@phi_cellState@phi).
#' @param state_to_compartment named chr vector: names = states, values = compartment.
#' @param compartment_levels compartment order (defaults to first appearance).
#' @param method "pearson" (default) or "spearman".
#' @param transform "log" (default, log10 profiles) or "none" (raw/linear).
#' @param genes optional gene subset to restrict the correlation to (e.g. HVGs
#'   or markers) -- sharpens contrasts; defaults to all genes in phi.
#' @param order_within how to order states inside each compartment block for the
#'   heatmap: "hclust" (default), "collinearity" (desc within_mean), or "input".
#' @return list(cor, comp, summary, state_order, blocks, params).
state_collinearity <- function(phi, state_to_compartment,
                               compartment_levels = NULL,
                               method    = c("pearson", "spearman"),
                               transform = c("log", "none"),
                               genes = NULL, pseudocount = NULL,
                               order_within = c("hclust", "collinearity", "input")) {
  
  method       <- match.arg(method)
  transform    <- match.arg(transform)
  order_within <- match.arg(order_within)
  
  phi <- as.matrix(phi)                                   # states x genes
  if (!is.null(genes)) {
    keep <- intersect(genes, colnames(phi))
    if (length(keep) < 2) stop("state_collinearity: <2 of `genes` found in phi columns.")
    phi <- phi[, keep, drop = FALSE]
  }
  
  M <- phi
  if (transform == "log") {
    if (is.null(pseudocount))
      pseudocount <- if (any(M == 0)) min(M[M > 0]) / 2 else 0
    M <- log10(M + pseudocount)
  }
  
  C <- cor(t(M), method = method)                          # states x states
  states <- rownames(C)
  
  comp <- state_to_compartment[states]
  if (anyNA(comp))
    stop("state_collinearity: states without a compartment: ",
         paste(states[is.na(comp)], collapse = ", "))
  if (is.null(compartment_levels)) compartment_levels <- unique(comp)
  
  # per-state within- and between-compartment mean correlation (self excluded)
  within_mean <- vapply(states, function(s) {
    mates <- setdiff(states[comp == comp[[s]]], s)
    if (!length(mates)) NA_real_ else mean(C[s, mates])
  }, numeric(1))
  between_mean <- vapply(states, function(s) {
    others <- states[comp != comp[[s]]]
    if (!length(others)) NA_real_ else mean(C[s, others])
  }, numeric(1))
  
  summary <- data.frame(
    state        = states,
    compartment  = factor(unname(comp), levels = compartment_levels),
    within_mean  = unname(within_mean),
    between_mean = unname(between_mean),
    contrast     = unname(within_mean - between_mean),
    stringsAsFactors = FALSE
  )
  
  # ---- ordering for the heatmap: by compartment, then within-block ----
  order_block <- function(s_in) {
    if (length(s_in) <= 2 || order_within == "input") return(s_in)
    if (order_within == "collinearity")
      return(s_in[order(-within_mean[s_in])])
    # hclust on 1 - correlation within the block
    d  <- as.dist(1 - C[s_in, s_in])
    hc <- hclust(d, method = "average")
    s_in[hc$order]
  }
  state_order <- unlist(lapply(compartment_levels,
                               function(cl) order_block(states[comp == cl])),
                        use.names = FALSE)
  
  # block index ranges (for rectangles / annotation on the ordered axis)
  ord_comp <- comp[state_order]
  runs <- rle(as.character(ord_comp))
  ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
  blocks <- data.frame(compartment = runs$values, start = starts, end = ends,
                       stringsAsFactors = FALSE)
  
  list(cor = C, comp = comp, summary = summary,
       state_order = state_order, blocks = blocks,
       compartment_levels = compartment_levels,
       params = list(method = method, transform = transform,
                     n_genes = ncol(M)))
}


# ------------------------------------------------------------------------------
# Heatmap: state x state correlation, blocked by compartment
# ------------------------------------------------------------------------------

plot_collinearity_heatmap <- function(collin, compartment_colors = NULL,
                                      title = "Transcriptional collinearity of cell states",
                                      subtitle = NULL, limits = c(-1, 1)) {
  C   <- collin$cor
  ord <- collin$state_order
  
  if (is.null(subtitle))
    subtitle <- sprintf("%s correlation of %s reference profiles over %d genes",
                        tools::toTitleCase(collin$params$method),
                        collin$params$transform, collin$params$n_genes)
  
  long <- as.data.frame(as.table(C), stringsAsFactors = FALSE)
  colnames(long) <- c("state_i", "state_j", "corr")
  long$state_i <- factor(long$state_i, levels = ord)
  long$state_j <- factor(long$state_j, levels = rev(ord))   # matrix-style y
  
  # compartment color for each ordered state (axis text)
  comp_ord <- collin$comp[ord]
  axis_cols_x <- if (!is.null(compartment_colors)) unname(compartment_colors[comp_ord]) else "black"
  axis_cols_y <- if (!is.null(compartment_colors)) unname(compartment_colors[rev(comp_ord)]) else "black"
  
  p <- ggplot(long, aes(state_i, state_j, fill = corr)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = limits, name = "Correlation") +
    # block-diagonal rectangles per compartment
    annotate("rect",
             xmin = collin$blocks$start - 0.5, xmax = collin$blocks$end + 0.5,
             ymin = nrow(C) - collin$blocks$end + 0.5,
             ymax = nrow(C) - collin$blocks$start + 1.5,
             fill = NA, color = "grey20", linewidth = 0.5) +
    coord_equal() +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 color = axis_cols_x, size = 6),
      axis.text.y = element_text(color = axis_cols_y, size = 6),
      panel.grid  = element_blank()
    )
  p
}


# ------------------------------------------------------------------------------
# Centrality: per-state mean within-compartment collinearity (the key test)
# ------------------------------------------------------------------------------

#' Ranked within-compartment collinearity per state, faceted by compartment.
#' Mass-absorbing states passed to `highlight` are marked; the hypothesis is
#' that they sit at the top of their compartment.
plot_collinearity_centrality <- function(collin, highlight = character(0),
                                         compartment_colors = NULL,
                                         title = "Within-compartment collinearity centrality",
                                         subtitle = NULL) {
  df <- collin$summary
  if (is.null(subtitle))
    subtitle <- "Mean correlation to compartment-mates (self excluded); marked = mass-absorbing states"
  
  df$is_hi <- df$state %in% highlight
  # order states within each compartment by within_mean (desc)
  df <- df[order(df$compartment, -df$within_mean), ]
  df$state <- factor(df$state, levels = rev(unique(df$state)))
  
  p <- ggplot(df, aes(within_mean, state)) +
    geom_segment(aes(x = 0, xend = within_mean, yend = state,
                     color = compartment), linewidth = 0.6, alpha = 0.5) +
    geom_point(aes(color = compartment, size = is_hi, shape = is_hi)) +
    facet_grid(compartment ~ ., scales = "free_y", space = "free_y") +
    scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 3.2), guide = "none") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    labs(title = title, subtitle = subtitle,
         x = "Mean within-compartment correlation", y = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "none",
          strip.text.y = element_text(angle = 0, face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())
  
  if (!is.null(compartment_colors))
    p <- p + scale_color_manual(values = compartment_colors)
  
  # bold the highlighted state labels
  hi_states <- rev(levels(df$state))
  lab_face  <- ifelse(hi_states %in% highlight, "bold", "plain")
  p + theme(axis.text.y = element_text(face = rev(lab_face)))
}