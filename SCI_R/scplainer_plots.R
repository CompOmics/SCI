#' Plot QC scatter panel for a QFeatures object
#'
#' Produces a 2x3 panel of n_precursors vs medianIntensityPeptide scatter plots,
#' each coloured by a different QC variable.
#'
#' @param qf        A QFeatures object with colData containing the expected columns.
#' @param conditions Optional vector of source.name values to filter the data by.
#' @param save_dir  Directory to save the figure. NULL skips saving.
#' @param filename  Output filename (without extension).
#' @param width     Figure width in inches.
#' @param height    Figure height in inches.
#'
#' @return A patchwork figure object (invisible).
plot_qc_scatter_panel <- function(
    qf,
    conditions = NULL,
    save_dir   = NULL,
    filename   = 'qc_scatter_panel',
    width      = 18,
    height     = 10,
    show = FALSE
) {
    cd <- data.frame(colData(qf))

    if (!is.null(conditions)) {
        cd <- cd[cd$source.name %in% conditions, ]
    }

    # enforce factor ordering
    cd$source.name <- factor(cd$source.name, levels = CONDITION_ORDER)
    cd$plate_idx   <- factor(cd$plate_idx,   levels = LAB_ORDER)

    base_aes <- aes(x = n_precursors, y = medianIntensityPeptide)

    make_plot <- function(color_var, title, continuous = FALSE,
                          color_scale = NULL, legend_title = color_var) {
        p <- ggplot(cd) +
            base_aes +
            geom_point(aes(color = .data[[color_var]]), size = 2) +
            labs(title = title, color = legend_title) +
            theme(legend.position = 'right')

        if (continuous) {
            p <- p + scale_color_viridis_c()
        } else if (!is.null(color_scale)) {
            p <- p + color_scale
        }
        p
    }

    panels <- list(
        make_plot('medianCV',    'Median CV',       continuous  = TRUE,
                  legend_title  = 'Median CV'),
        make_plot('source.name', 'Cell type',
                  color_scale   = scale_color_manual(values = CONDITION_COLORS,
                                                     drop   = FALSE),
                  legend_title  = 'Condition'),
        make_plot('comment.X.position.', 'Column position',
                  legend_title  = 'Column'),
        make_plot('plate_idx',   'Lab',
                  color_scale   = scale_color_manual(values = LAB_COLORS,
                                                     drop   = FALSE),
                  legend_title  = 'Lab'),
        make_plot('Reinjection', 'Reinjection',
                  legend_title  = 'Reinjection'),
        make_plot('instrument_type', 'Instrument',
                  legend_title  = 'Instrument')
    )

    fig <- patchwork::wrap_plots(panels, ncol = 3, nrow = 2)

    if (!is.null(save_dir)) {
        save_figure(fig, file.path(save_dir, filename), width = width, height = height)
    }

    if (show) {
        options(repr.plot.width = width, repr.plot.height = height)
        print(fig)
    }
    invisible(fig)
}

plot_scp_model_filter <- function(
    sce_list,
    np_threshold = 20,
    assay_name   = "precursors_log",
    save_dir     = NULL,
    filename     = 'scp_model_filter',
    width        = 14,
    height       = 8,
    show         = TRUE
) {
    # apply threshold to all
    sce_list <- lapply(sce_list, function(obj) {
        scpModelFilterThreshold(obj, assay_name) <- np_threshold
        obj
    })

    # plot
    panels <- mapply(
        function(obj, title) scpModelFilterPlot(obj) + ggtitle(title),
        sce_list, names(sce_list),
        SIMPLIFY = FALSE
    )

    ncols <- min(2, length(panels))
    nrows <- ceiling(length(panels) / ncols)

    fig <- patchwork::wrap_plots(panels, ncol = ncols, nrow = nrows) +
        patchwork::plot_layout(guides = 'keep') &
        theme(plot.margin = margin(10, 15, 10, 15))

    if (!is.null(save_dir)) {
        save_figure(fig, file.path(save_dir, filename), width = width, height = height)
    }

    if (show) {
        options(repr.plot.width = width, repr.plot.height = height)
        print(fig)
    }

    invisible(fig)
}

plot_scp_variance <- function(
    va_list,
    colour_maps,
    save_dir    = NULL,
    filename    = 'scp_variance',
    width       = 14,
    height      = 8,
    show        = TRUE
) {
    panels <- mapply(
        function(va, colour_map, title) {
            scpVariancePlot(va) +
                scale_fill_manual(values = colour_map) +
                ggtitle(title)
        },
        va_list, colour_maps, names(va_list),
        SIMPLIFY = FALSE
    )

    ncols <- min(2, length(panels))
    nrows <- ceiling(length(panels) / ncols)

    fig <- patchwork::wrap_plots(panels, ncol = ncols, nrow = nrows) +
        patchwork::plot_layout(guides = 'keep') &
        theme(plot.margin = margin(10, 15, 10, 15))

    if (!is.null(save_dir)) {
        save_figure(fig, file.path(save_dir, filename), width = width, height = height)
    }

    if (show) {
        options(repr.plot.width = width, repr.plot.height = height)
        print(fig)
    }

    invisible(fig)
}

plot_variance_violin <- function(
    vaRes,
    title      = "Variance explained per precursor",
    colour_map = NULL,
    save_dir   = NULL,
    filename   = 'variance_violin',
    width      = 8,
    height     = 5,
    show       = TRUE
) {
    var_long <- lapply(names(vaRes), function(v) {
        data.frame(variable = v, pct = vaRes[[v]]$percentExplainedVar)
    }) %>% bind_rows()

    # order by median
    var_order <- var_long %>%
        group_by(variable) %>%
        summarise(med = median(pct, na.rm = TRUE), .groups = 'drop') %>%
        arrange(desc(med)) %>%
        pull(variable)

    var_long$variable <- factor(var_long$variable, levels = var_order)

    if (is.null(colour_map)) {
        n_vars     <- length(var_order)
        colour_map <- setNames(
            colorRampPalette(brewer.pal(min(n_vars, 9), 'Set1'))(n_vars),
            var_order
        )
    }

    fig <- ggplot(var_long, aes(x = variable, y = pct / 100, fill = variable)) +
        geom_violin(scale = 'width', alpha = 0.8, show.legend = FALSE) +
        geom_boxplot(width = 0.08, outlier.shape = NA,
                     fill = 'white', show.legend = FALSE) +
        scale_fill_manual(values = colour_map) +
        scale_y_continuous(
            labels = function(x) paste0(round(x * 100), '%'),
            limits = c(0, 1)
        ) +
        labs(
            title    = title,
            subtitle = sprintf('n = %d precursors', nrow(vaRes[[1]])),
            x        = NULL,
            y        = 'Fraction of variance explained'
        ) +
        theme(
            axis.text.x = element_text(angle = 35, hjust = 1),
            plot.title  = element_text(face = 'bold')
        )

    if (!is.null(save_dir)) {
        save_figure(fig, file.path(save_dir, filename), width = width, height = height)
    }

    if (show) {
        options(repr.plot.width = width, repr.plot.height = height)
        print(fig)
    }

    invisible(fig)
}


plot_gmf_embedding <- function(
    sce,
    title       = NULL,
    colour_dims = c("plate_idx", "source.name", "comment.X.position.",
                    "n_precursors", "medianIntensityPeptide", "Reinjection", "Diameter"),
    dimred      = "GMF",
    save_dir    = NULL,
    filename    = "gmf_embedding",
    width       = 30,
    height      = 6,
    show        = TRUE
) {
    if ("Diameter" %in% colnames(colData(sce))) {
        colData(sce)$Diameter[colData(sce)$Diameter == -1] <- NA
    }
    panels <- lapply(colour_dims, function(col_by) {
        p <- tryCatch(
            plotReducedDim(sce, dimred = dimred, colour_by = col_by) +
                labs(title = col_by) +
                theme(legend.position = 'right',
                      plot.title      = element_text(size = 9)),
            error = function(e) {
                message(sprintf("Skipping %s: %s", col_by, e$message))
                NULL
            }
        )

        if (!is.null(p)) {
            if (col_by == "source.name") {
                p <- p + scale_colour_manual(values = CONDITION_COLORS, na.value = "grey80")
            } else if (col_by == "plate_idx") {
                p <- p + scale_colour_manual(values = LAB_COLORS, na.value = "grey80")
            }
        }
        p
    })

    panels <- Filter(Negate(is.null), panels)

    fig <- patchwork::wrap_plots(panels, nrow = 1)

    if (!is.null(title)) {
        fig <- fig + patchwork::plot_annotation(
            title = title,
            theme = theme(plot.title = element_text(face = 'bold', size = 11))
        )
    }

    if (!is.null(save_dir)) {
        save_figure(fig, file.path(save_dir, filename), width = width, height = height)
    }

    if (show) {
        options(repr.plot.width = width, repr.plot.height = height)
        print(fig)
    }

    invisible(fig)
}


plot_silhouette_scores <- function(
    stage_list,
    col_bio   = COL_BIO,
    col_batch = COL_BATCH,
    save_dir  = NULL,
    filename  = "silhouette_scores",
    width     = 9,
    height    = 5,
    show      = TRUE
) {
    sil_results <- lapply(names(stage_list), function(nm) {
        cat(sprintf("  Computing silhouettes: %s...\n", nm))
        flush.console()
        res       <- compute_silhouettes(stage_list[[nm]])
        res$stage <- nm
        res
    }) %>% bind_rows()

    sil_results$stage <- factor(sil_results$stage, levels = names(stage_list))

    cat("\nSilhouette scores per correction stage:\n")
    print(sil_results[, c("stage", "bio_silhouette", "batch_silhouette")])

    sil_long <- rbind(
        data.frame(stage      = sil_results$stage,
                   metric     = "Biological separation (source.name)",
                   silhouette = sil_results$bio_silhouette),
        data.frame(stage      = sil_results$stage,
                   metric     = "Batch structure (plate_idx)",
                   silhouette = sil_results$batch_silhouette)
    )
    sil_long$stage <- factor(sil_long$stage, levels = names(stage_list))

    fig <- ggplot(sil_long, aes(x = stage, y = silhouette, fill = metric)) +
        geom_col(position = position_dodge(width = 0.7), width = 0.6) +
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
        scale_fill_manual(values = c(
            "Biological separation (source.name)" = col_bio,
            "Batch structure (plate_idx)"         = col_batch
        )) +
        labs(
            title    = "Silhouette scores across correction stages",
            subtitle = "First 4 GMF components; higher bio + lower batch = better correction",
            x        = NULL,
            y        = "Mean silhouette width",
            fill     = NULL
        ) +
        theme(
            axis.text.x     = element_text(angle = 25, hjust = 1),
            legend.position = "top"
        )

    if (!is.null(save_dir)) {
        save_figure(fig, file.path(save_dir, filename), width = width, height = height)
    }

    if (show) {
        options(repr.plot.width = width, repr.plot.height = height)
        print(fig)
    }

    invisible(fig)
}