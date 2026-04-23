plot_pca <- function(
    pca_result,
    meta,
    colour_var,
    shape_var  = NULL,
    title      = "",
    save_dir   = NULL,
    filename   = "pca",
    width      = 14,
    height     = 5,
    show       = TRUE
) {
    var_expl <- summary(pca_result)$importance[2, 1:2] * 100

    df      <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2],
                          group = meta[[colour_var]])
    is_cont <- is.numeric(meta[[colour_var]])
    if (!is_cont)            df$group <- factor(df$group)
    if (!is.null(shape_var)) df$shape <- factor(meta[[shape_var]])

    colour_scale <- if (colour_var == COL_PLATE) {
        scale_colour_manual(values = LAB_COLORS, name = colour_var)
    } else if (colour_var == COL_CELL_TYPE) {
        scale_colour_manual(values = CONDITION_COLORS, name = colour_var)
    } else if (is_cont) {
        scale_colour_distiller(palette = "RdYlBu", name = colour_var)
    } else {
        n   <- nlevels(factor(df$group))
        pal <- colorRampPalette(brewer.pal(min(n, 9), "Set1"))(n)
        scale_colour_manual(values = pal, name = colour_var)
    }

    p <- ggplot(df, aes(x = PC1, y = PC2, colour = group)) +
        { if (!is.null(shape_var)) geom_point(aes(shape = shape), alpha = 0.7, size = 1.8)
          else                     geom_point(alpha = 0.7, size = 1.8) } +
        colour_scale +
        { if (!is.null(shape_var))
            scale_shape_manual(values = seq_len(nlevels(factor(meta[[shape_var]]))),
                               name   = shape_var)
        } +
        labs(title = title,
             x     = sprintf("PC1 (%.1f%%)", var_expl[1]),
             y     = sprintf("PC2 (%.1f%%)", var_expl[2])) +
        theme(plot.title = element_text(face = "bold"))

    if (!is.null(save_dir)) {
        save_figure(p, file.path(save_dir, filename), width = width, height = height)
    }

    if (show) {
        options(repr.plot.width = width, repr.plot.height = height)
        print(p)
    }

    invisible(p)
}