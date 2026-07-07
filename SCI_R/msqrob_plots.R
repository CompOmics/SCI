plot_top_proteins <- function(
    res,
    qf,
    experiment_name = "proteins_medianPolish",
    assay_name      = "proteins_medianPolish",
    n_show          = 6,
    colour_var      = COL_CELL_TYPE,
    jitter_var      = COL_PLATE,
    share_y         = FALSE,
    title           = NULL,
    save_dir        = NULL,
    filename        = "top_proteins",
    width           = 14,
    height          = 6,
    show            = TRUE
) {
    if (is(qf, "QFeatures")) {
        sce <- qf[[experiment_name]]
        cd  <- as.data.frame(colData(qf))
    } else {
        sce <- qf
        cd  <- as.data.frame(colData(qf))
    }

    top_prots <- res %>%
        filter(!is.na(adjPval)) %>%
        arrange(pval) %>%
        slice_head(n = n_show) %>%
        pull(protein)

    fig <- assay(sce, assay_name)[top_prots, , drop = FALSE] %>%
        as.data.frame() %>%
        rownames_to_column("protein") %>%
        pivot_longer(-protein, names_to = "cell", values_to = "intensity") %>%
        left_join(
            cd %>% rownames_to_column("cell") %>%
                select(cell, all_of(c(colour_var, jitter_var))),
            by = "cell"
        ) %>%
        filter(!is.na(intensity)) %>%
        mutate(protein = factor(protein, levels = top_prots)) %>%
        ggplot(aes(x = .data[[colour_var]], y = intensity, fill = .data[[colour_var]])) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
        geom_jitter(aes(colour = .data[[jitter_var]]), width = 0.15, size = 0.8, alpha = 0.5) +
        facet_wrap(~ protein, nrow = 1,
                   scales = if (share_y) "fixed" else "free_y") +
        scale_fill_manual(values  = CONDITION_COLORS, name = colour_var) +
        scale_colour_manual(values = LAB_COLORS,      name = jitter_var) +
        theme(plot.title  = element_text(face = "bold"),
              axis.text.x = element_text(angle = 35, hjust = 1),
              strip.text  = element_text(face = "bold.italic", size = 9)) +
        labs(
            title = title %||% sprintf("Top %d most significant proteins", n_show),
            x     = NULL,
            y     = "log2 Intensity"
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

plot_volcano <- function(
    res,
    title         = "",
    fdr_threshold = .001,
    fc_threshold  = .5,
    num_prot      = 10,
    x_lim         = NULL,
    y_lim         = NULL,
    col_sig       = "#2166ac",
    col_ns        = "grey60",
    save_dir      = NULL,
    filename      = "volcano",
    width         = 7,
    height        = 6,
    show          = FALSE
) {
    res <- res %>%
        filter(is.finite(logFC), is.finite(pval), pval > 0) %>%
        mutate(
            neg_log10_p = -log10(pval),
            significant = adjPval < fdr_threshold & abs(logFC) > fc_threshold,
            direction   = case_when(
                significant & logFC >  fc_threshold ~ "up",
                significant & logFC < -fc_threshold ~ "down",
                TRUE                                ~ "ns"
            )
        )

    n_up   <- sum(res$direction == "up",   na.rm = TRUE)
    n_down <- sum(res$direction == "down", na.rm = TRUE)

    res_ns  <- res[res$direction == "ns", ]
    res_sig <- res[res$direction != "ns", ]

    p <- ggplot(res, aes(x = logFC, y = neg_log10_p)) +
        geom_point(data = res_ns,  colour = col_ns,  size = 0.7, alpha = 0.4) +
        geom_point(data = res_sig, colour = col_sig, size = 1.2, alpha = 0.85) +
        geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                   linetype = "dashed", colour = "grey40", linewidth = 0.4) +
        geom_hline(yintercept = -log10(fdr_threshold),
                   linetype = "dashed", colour = "grey40", linewidth = 0.4) +
        annotate("text", x = -Inf, y = Inf,
                 label   = sprintf("down: %d", n_down),
                 hjust   = -0.1, vjust = 1.5,
                 size    = 3, colour = col_sig) +
        annotate("text", x =  Inf, y = Inf,
                 label   = sprintf("up: %d", n_up),
                 hjust   = 1.1, vjust = 1.5,
                 size    = 3, colour = col_sig) +
        labs(
            title = title,
            x     = expression(log[2]~FC),
            y     = expression(-log[10]~p)
        )

    # label top proteins by significance
    if (num_prot > 0) {
        top <- res %>%
            filter(significant) %>%
            arrange(adjPval) %>%
            head(num_prot)

        if (nrow(top) > 0) {
            p <- p + ggrepel::geom_label_repel(
                data        = top,
                aes(label   = protein),
                size        = 2.5,
                max.overlaps= 20,
                box.padding = 0.3
            )
        }
    }

    if (!is.null(x_lim)) p <- p + xlim(x_lim)
    if (!is.null(y_lim)) p <- p + ylim(y_lim)

    if (!is.null(save_dir)) {
        save_figure(p, file.path(save_dir, filename), width = width, height = height)
    }

    if (show) {
        options(repr.plot.width = width, repr.plot.height = height)
        print(p)
    }

    invisible(p)
}

plot_volcano_grid <- function(
    model_configs,
    fdr_threshold = .001,
    fc_threshold  = .5,
    save_dir      = NULL,
    filename      = "volcano_grid",
    width         = 14,
    height        = 12,
    show          = TRUE
) {
    all_res <- bind_rows(lapply(model_configs, function(cfg) {
        get(paste0("res_", cfg$name)) %>% mutate(model = cfg$name)
    }))

    # compute shared axis limits, guard against Inf and NaN
    logfc_vals <- all_res$logFC[is.finite(all_res$logFC)]
    pval_vals  <- -log10(all_res$pval[is.finite(all_res$pval) & all_res$pval > 0])

    x_lim <- c(floor(min(logfc_vals, na.rm = TRUE)),
                ceiling(max(logfc_vals, na.rm = TRUE)))
    y_lim <- c(0, ceiling(max(pval_vals, na.rm = TRUE) * 1.05))

    cat(sprintf("Shared x_lim: [%.1f, %.1f]  y_lim: [%.1f, %.1f]\n",
                x_lim[1], x_lim[2], y_lim[1], y_lim[2]))

    panels <- lapply(model_configs, function(cfg) {
        plot_volcano(
            get(paste0("res_", cfg$name)),
            title   = cfg$label,
            num_prot = 0,
            x_lim   = x_lim,
            y_lim   = y_lim
        )
    })

    fig <- wrap_plots(panels, ncol = 2) +
        plot_annotation(
            title    = "DEA — K562 vs HEK293T across model complexity",
            subtitle = sprintf("FDR < %.3f | |log2FC| > %.1f | msqrob2 ridge + robust",
                               fdr_threshold, fc_threshold),
            theme    = theme(
                plot.title    = element_text(face = "bold", size = 14),
                plot.subtitle = element_text(size = 11, colour = "grey40")
            )
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

plot_upset_sig <- function(
    model_configs,
    fdr_threshold = FDR_THRESHOLD,
    fc_threshold  = FC_THRESHOLD,
    direction     = "both",
    bar_color     = "#2166ac",
    save_dir      = NULL,
    filename      = "upset_sig",
    width         = 10,
    height        = 5,
    show          = TRUE
) {
    get_sig <- function(res, fdr, fc, direction) {
        res <- res %>% filter(!is.na(adjPval), adjPval < fdr)
        if (direction == "up")   res <- res %>% filter(logFC >  fc)
        if (direction == "down") res <- res %>% filter(logFC < -fc)
        if (direction == "both") res <- res %>% filter(abs(logFC) > fc)
        res$protein
    }

    sig_lists        <- lapply(model_configs, function(cfg) {
        get_sig(get(paste0("res_", cfg$name)), fdr_threshold, fc_threshold, direction)
    })
    names(sig_lists) <- sapply(model_configs, `[[`, "label")

    cat("Significant proteins per model:\n")
    print(sapply(sig_lists, length))

    fig <- UpSetR::upset(
        UpSetR::fromList(sig_lists),
        order.by        = "freq",
        point.size      = 3,
        line.size       = 1,
        text.scale      = 1.3,
        main.bar.color  = bar_color,
        sets.bar.color  = bar_color,
        mainbar.y.label = "Intersection size",
        sets.x.label    = "Proteins per model"
    )

    if (!is.null(save_dir)) {
        # UpSetR plots are base graphics so we need to wrap for save_figure
        for (fmt in c("png", "svg", "pdf")) {
            out_dir <- file.path(save_dir, fmt)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
            out     <- file.path(out_dir, paste0(filename, ".", fmt))
            dev     <- switch(fmt,
                png  = function(f) png(f,  width = width, height = height, units = "in", res = 150),
                svg  = function(f) svglite::svglite(f, width = width, height = height),
                pdf  = function(f) cairo_pdf(f, width = width, height = height)
            )
            dev(out)
            print(fig)
            dev.off()
        }
    }

    if (show) {
        options(repr.plot.width = width, repr.plot.height = height)
        print(fig)
    }

    invisible(fig)
}

plot_model_comparison_grid <- function(
    model_configs,
    save_dir = NULL,
    filename = "model_comparison_grid",
    width    = 14,
    height   = 5,
    show     = TRUE
) {
    # consecutive pairs: norm→plate, plate→batch, batch→all
    pairs <- lapply(seq_len(length(model_configs) - 1), function(i) {
        list(
            simple = model_configs[[i]],
            full   = model_configs[[i + 1]]
        )
    })

    fc_panels <- lapply(pairs, function(pair) {
        res_simple <- get(paste0("res_", pair$simple$name))
        res_full   <- get(paste0("res_", pair$full$name))

        comp <- data.frame(
            logFC_simple = res_simple$logFC,
            logFC_full   = res_full$logFC[match(res_simple$protein, res_full$protein)]
        ) %>% filter(!is.na(logFC_simple), !is.na(logFC_full))

        r <- cor(comp$logFC_simple, comp$logFC_full, use = "complete.obs")

        ggplot(comp, aes(x = logFC_simple, y = logFC_full)) +
            geom_point(alpha = 0.3, size = 0.6) +
            geom_abline(slope = 1, intercept = 0,
                        colour = "firebrick", linetype = "dashed") +
            labs(title    = sprintf("%s → %s", pair$simple$label, pair$full$label),
                 subtitle = sprintf("r = %.4f", r),
                 x        = sprintf("logFC (%s)", pair$simple$name),
                 y        = sprintf("logFC (%s)", pair$full$name))
    })

    se_panels <- lapply(pairs, function(pair) {
        res_simple <- get(paste0("res_", pair$simple$name))
        res_full   <- get(paste0("res_", pair$full$name))

        comp <- data.frame(
            se_simple = res_simple$se,
            se_full   = res_full$se[match(res_simple$protein, res_full$protein)]
        ) %>% filter(!is.na(se_simple), !is.na(se_full))

        se_red <- 100 * (1 - mean(comp$se_full,   na.rm = TRUE) /
                              mean(comp$se_simple, na.rm = TRUE))

        ggplot(comp, aes(x = se_simple, y = se_full)) +
            geom_point(alpha = 0.3, size = 0.6) +
            geom_abline(slope = 1, intercept = 0,
                        colour = "firebrick", linetype = "dashed") +
            labs(title    = sprintf("%s → %s", pair$simple$label, pair$full$label),
                 subtitle = sprintf("SE reduction: %.1f%%", se_red),
                 x        = sprintf("SE (%s)", pair$simple$name),
                 y        = sprintf("SE (%s)", pair$full$name))
    })

    panels <- unlist(
        lapply(seq_along(pairs), function(i) list(fc_panels[[i]], se_panels[[i]])),
        recursive = FALSE
    )

    fig <- wrap_plots(panels, nrow = 2, byrow = FALSE) +
        plot_annotation(
            title    = "Incremental model comparison",
            subtitle = "Top row: logFC agreement  |  Bottom row: SE reduction",
            theme    = theme(
                plot.title    = element_text(face = "bold", size = 13),
                plot.subtitle = element_text(size = 10, colour = "grey40")
            )
        )

    n_pairs <- length(pairs)

    if (!is.null(save_dir)) {
        save_figure(fig, file.path(save_dir, filename),
                    width  = width,
                    height = height)
    }

    if (show) {
        options(repr.plot.width  = width,
                repr.plot.height = height)
        print(fig)
    }

    invisible(fig)
}

# ── 1. Per-plate volcanos colored by integrated model significance ─────────────
plot_volcano_grid_per_plate <- function(
    results_per_plate,
    res_integrated,
    fdr_threshold   = .001,
    fc_threshold    = 0.5,
    ncol            = 3,
    point_size      = 1.8,
    point_alpha     = 0.6,
    title_size      = 22,
    annot_size      = 6.5,
    axis_title_size = 18,
    axis_text_size  = 13,
    save_dir        = NULL,
    filename        = "volcano_grid_per_plate",
    width           = 18,
    height          = 12,
    show            = TRUE
) {
    # proteins significant in integrated model
    sig_integrated <- res_integrated %>%
        filter(!is.na(adjPval), adjPval < fdr_threshold, abs(logFC) > fc_threshold) %>%
        pull(protein)

    all_res    <- bind_rows(results_per_plate)
    logfc_vals <- all_res$logFC[is.finite(all_res$logFC)]
    pval_vals  <- -log10(all_res$pval[is.finite(all_res$pval) & all_res$pval > 0])

    x_lim <- c(floor(min(logfc_vals,   na.rm = TRUE)),
                ceiling(max(logfc_vals, na.rm = TRUE)))
    y_lim <- c(0, ceiling(max(pval_vals, na.rm = TRUE) * 1.05))

    panels <- imap(results_per_plate, function(res, plate) {
        res <- res %>%
            filter(is.finite(logFC), is.finite(pval), pval > 0) %>%
            mutate(
                neg_log10_p = -log10(pval),
                sig_plate   = !is.na(adjPval) & adjPval < fdr_threshold &
                              abs(logFC) > fc_threshold,
                sig_integrated = protein %in% sig_integrated,
                colour_group = case_when(
                    sig_plate & sig_integrated  ~ "both",
                    sig_plate & !sig_integrated ~ "plate only",
                    !sig_plate & sig_integrated ~ "integrated only",
                    TRUE                        ~ "ns"
                )
            )

        n_both  <- sum(res$colour_group == "both")
        n_plate <- sum(res$colour_group == "plate only")

        res_ns  <- res[res$colour_group == "ns", ]
        res_col <- res[res$colour_group != "ns", ]

        ggplot(res, aes(x = logFC, y = neg_log10_p)) +
            geom_point(data = res_ns,  colour = "grey85",
                       size = point_size * 0.75, alpha = point_alpha * 0.5) +
            geom_point(data = res_col, aes(colour = colour_group),
                       size = point_size, alpha = min(point_alpha + 0.25, 1)) +
            scale_colour_manual(
                values = c(
                    "both"            = "#d73027",
                    "plate only"      = "#f4a442",
                    "integrated only" = "#4393c3"
                ),
                name   = NULL,
                guide  = guide_legend(
                    override.aes = list(size = 4, alpha = 1)
                )
            ) +
            geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                       linetype = "dashed", colour = "grey40", linewidth = 0.3) +
            geom_hline(yintercept = -log10(fdr_threshold),
                       linetype = "dashed", colour = "grey40", linewidth = 0.3) +
            annotate("text", x = Inf, y = Inf,
                     label = sprintf("both: %d\nplate: %d", n_both, n_plate),
                     hjust = 1.05, vjust = 1.35, size = annot_size,
                     lineheight = 0.9, colour = "grey30") +
            xlim(x_lim) + ylim(y_lim) +
            labs(title = plate,
                 x     = expression(log[2]~FC),
                 y     = expression(-log[10]~p)) +
            theme(plot.title      = element_text(face = "bold", size = title_size,
                                                 hjust = 0.5),
                  axis.title      = element_text(size = axis_title_size),
                  axis.text       = element_text(size = axis_text_size),
                  legend.position = "none")
    })

    # Shared x/y axes: collect the duplicate axis text and titles onto the outer
    # edge only, so inner panels show no repeated ticks or labels.
    # Requires patchwork >= 1.2.0 (the `axes` / `axis_titles` arguments).
    fig <- wrap_plots(panels, ncol = ncol) +
        plot_layout(axes = "collect", axis_titles = "collect") +
        plot_annotation(
            title    = "Per-plate DEA — colored by integrated model",
            subtitle = sprintf("FDR < %.3f | |log2FC| > %.1f", fdr_threshold, fc_threshold),
            theme    = theme(
                plot.title    = element_text(face = "bold", size = title_size),
                plot.subtitle = element_text(size = axis_title_size, colour = "grey40")
            )
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


# ── Per-plate volcanos colored by cross-lab overlap of hits ────────────────────
# The category of each protein is assigned from how many plates it is
# significant in:
#   all plates  -> "Hit in all labs"
#   >1 & < all  -> "Hit in few labs"
#   exactly 1   -> "Hit in one lab"
# A point is coloured by that category only in the plates where it is itself
# significant; elsewhere it is "Other" (grey). Works on any per-plate result
# list (filtered or unfiltered). Significance uses adjPval and logFC; the y-axis
# is -log10(adjPval).
plot_volcano_grid_per_plate_overlap <- function(
    results_per_plate,
    fdr_threshold   = .001,
    fc_threshold    = 0.5,
    ncol            = 3,
    point_size      = 1.8,
    hit_alpha       = 0.9,
    other_alpha     = 0.3,
    title_size      = 22,
    axis_title_size = 18,
    axis_text_size  = 13,
    save_dir        = NULL,
    filename        = "volcano_grid_per_plate_overlap",
    width           = 18,
    height          = 12,
    show            = TRUE
) {
    n_plates <- length(results_per_plate)

    all_res <- bind_rows(results_per_plate) %>%
        mutate(sig = !is.na(adjPval) & adjPval < fdr_threshold &
                     abs(logFC) > fc_threshold)

    # number of plates in which each protein is significant
    # (dplyr::count qualified to avoid masking by ProtGenerics/QFeatures)
    hit_counts <- all_res %>%
        filter(sig) %>%
        dplyr::count(protein, name = "n_hits")

    hit_in_all <- hit_counts$protein[hit_counts$n_hits == n_plates]
    hit_in_one <- hit_counts$protein[hit_counts$n_hits == 1]
    hit_in_few <- hit_counts$protein[hit_counts$n_hits > 1 &
                                     hit_counts$n_hits <  n_plates]

    cat_levels <- c("Hit in all labs", "Hit in few labs", "Hit in one lab", "Other")
    cat_colors <- c("Hit in all labs" = "#8A2BE2",   # purple
                    "Hit in few labs" = "#E84855",   # red
                    "Hit in one lab"  = "#FFA500",    # orange
                    "Other"           = "#AAAAAA")    # grey
    draw_order <- c("Other", "Hit in one lab", "Hit in few labs", "Hit in all labs")

    classify <- function(protein, sig) {
        cls <- rep("Other", length(protein))
        cls[sig & protein %in% hit_in_few] <- "Hit in few labs"
        cls[sig & protein %in% hit_in_one] <- "Hit in one lab"
        cls[sig & protein %in% hit_in_all] <- "Hit in all labs"
        cls
    }

    logfc_vals <- all_res$logFC[is.finite(all_res$logFC)]
    pval_vals  <- -log10(pmax(all_res$adjPval[is.finite(all_res$adjPval)], 1e-300))
    x_lim <- c(floor(min(logfc_vals, na.rm = TRUE)),
               ceiling(max(logfc_vals, na.rm = TRUE)))
    y_lim <- c(0, ceiling(max(pval_vals, na.rm = TRUE) * 1.05))

    panels <- imap(results_per_plate, function(res, plate) {
        res <- res %>%
            filter(is.finite(logFC), !is.na(adjPval), adjPval > 0) %>%
            mutate(
                neg_log10_p = -log10(pmax(adjPval, 1e-300)),
                sig         = adjPval < fdr_threshold & abs(logFC) > fc_threshold,
                category    = factor(classify(protein, sig), levels = cat_levels)
            ) %>%
            arrange(match(as.character(category), draw_order))

        res_other <- res[res$category == "Other", ]
        res_hit   <- res[res$category != "Other", ]

        ggplot(res, aes(x = logFC, y = neg_log10_p)) +
            geom_point(data = res_other, colour = cat_colors[["Other"]],
                       size = point_size * 0.7, alpha = other_alpha, stroke = 0) +
            geom_point(data = res_hit, aes(colour = category),
                       size = point_size, alpha = hit_alpha, stroke = 0) +
            scale_colour_manual(values = cat_colors, limits = cat_levels,
                                drop = FALSE, name = NULL,
                                guide = guide_legend(
                                    override.aes = list(size = 4, alpha = 1))) +
            geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                       linetype = "dashed", colour = "grey40", linewidth = 0.3) +
            geom_hline(yintercept = -log10(fdr_threshold),
                       linetype = "dashed", colour = "grey40", linewidth = 0.3) +
            xlim(x_lim) + ylim(y_lim) +
            labs(title = plate,
                 x     = expression(log[2]~FC),
                 y     = expression(-log[10]~adjP)) +
            theme(plot.title      = element_text(face = "bold", size = title_size,
                                                 hjust = 0.5),
                  axis.title      = element_text(size = axis_title_size),
                  axis.text       = element_text(size = axis_text_size),
                  legend.position = "none")
    })

    fig <- wrap_plots(panels, ncol = ncol) +
        plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") +
        plot_annotation(
            title    = "Per-plate DEA — colored by cross-lab overlap",
            subtitle = sprintf("FDR < %.3f | |log2FC| > %.1f", fdr_threshold, fc_threshold),
            theme    = theme(
                plot.title    = element_text(face = "bold", size = title_size),
                plot.subtitle = element_text(size = axis_title_size, colour = "grey40")
            )
        ) &
        theme(legend.position = "bottom",
              legend.text     = element_text(size = axis_title_size))

    if (!is.null(save_dir)) {
        save_figure(fig, file.path(save_dir, filename), width = width, height = height)
    }
    if (show) {
        options(repr.plot.width = width, repr.plot.height = height)
        print(fig)
    }
    invisible(fig)
}


plot_upset_per_plate <- function(
    results_per_plate,
    res_integrated,
    fdr_threshold = FDR_THRESHOLD,
    fc_threshold  = FC_THRESHOLD,
    bar_color     = "#2166ac",
    save_dir      = NULL,
    filename      = "upset_per_plate",
    width         = 16,
    height        = 7,
    show          = TRUE
) {
    get_sig <- function(res, direction) {
        res <- res %>% filter(!is.na(adjPval), adjPval < fdr_threshold)
        if (direction == "up")   res <- res %>% filter(logFC >  fc_threshold)
        if (direction == "down") res <- res %>% filter(logFC < -fc_threshold)
        if (direction == "both") res <- res %>% filter(abs(logFC) > fc_threshold)
        res$protein
    }

    dir_label <- c(up = "Up", down = "Down", both = "Differentially")

    for (direction in c("up", "down", "both")) {
        sig_lists <- c(
            list("Integrated" = get_sig(res_integrated, direction)),
            lapply(results_per_plate, get_sig, direction = direction)
        )

        cat(sprintf("\n%sregulated proteins per model:\n", dir_label[[direction]]))
        print(sapply(sig_lists, length))

        fig <- UpSetR::upset(
            UpSetR::fromList(sig_lists),
            nsets           = length(sig_lists),
            order.by        = "freq",
            point.size      = 3,
            line.size       = 1,
            text.scale      = c(1.3, 1.2, 1.2, 1.1, 1.2, 1.2),
            main.bar.color  = bar_color,
            sets.bar.color  = bar_color,
            mainbar.y.label = "Intersection size",
            sets.x.label    = sprintf("%sregulated proteins", dir_label[[direction]])
        )

        fname <- paste0(filename, "_", direction)

        if (!is.null(save_dir)) {
            for (fmt in c("png", "svg", "pdf")) {
                out_dir <- file.path(save_dir, fmt)
                dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                out <- file.path(out_dir, paste0(fname, ".", fmt))
                dev <- switch(fmt,
                    png = function(f) png(f, width = width, height = height,
                                          units = "in", res = 150),
                    svg = function(f) svglite::svglite(f, width = width, height = height),
                    pdf = function(f) cairo_pdf(f, width = width, height = height)
                )
                dev(out)
                print(fig)
                dev.off()
            }
        }

        if (show) {
            options(repr.plot.width = width, repr.plot.height = height)
            print(fig)
            flush.console()
        }
    }
}


# ── 3. Overlap heatmap with clustering ────────────────────────────────────────
plot_overlap_heatmap <- function(
    results_per_plate,
    res_integrated,
    fdr_threshold = FDR_THRESHOLD,
    fc_threshold  = FC_THRESHOLD,
    metric        = c("jaccard", "overlap", "logfc_cor"),
    save_dir      = NULL,
    filename      = "overlap_heatmap",
    width         = 8,
    height        = 7,
    show          = TRUE
) {
    metric <- match.arg(metric)

    get_sig <- function(res) {
        res %>%
            filter(!is.na(adjPval), adjPval < fdr_threshold, abs(logFC) > fc_threshold) %>%
            pull(protein)
    }

    all_results <- c(list("Integrated" = res_integrated), results_per_plate)
    models      <- names(all_results)
    n           <- length(models)

    # ── metric: logFC correlation ─────────────────────────────────────────────
    if (metric == "logfc_cor") {

        # union of significant proteins across all models
        sig_any      <- Reduce(union, lapply(all_results, get_sig))
        cat(sprintf("logFC correlation computed on %d proteins (union of significant)\n",
                    length(sig_any)))

        # build protein x model logFC matrix, restricted to sig_any
        logfc_mat <- sapply(all_results, function(res) {
            v <- setNames(res$logFC, res$protein)
            v[sig_any]
        })
        rownames(logfc_mat) <- sig_any

        sim_mat      <- cor(logfc_mat, use = "pairwise.complete.obs", method = "spearman")
        legend_label <- "Spearman\ncorrelation\n(logFC)"
        limits       <- c(-1, 1)

    # ── metric: jaccard or overlap coefficient ────────────────────────────────
    } else {

        sig_lists <- lapply(all_results, get_sig)

        sim_mat <- matrix(0, nrow = n, ncol = n, dimnames = list(models, models))
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                a <- sig_lists[[i]]
                b <- sig_lists[[j]]
                inter <- length(intersect(a, b))
                if (metric == "jaccard") {
                    sim_mat[i, j] <- inter / length(union(a, b))
                } else if (metric == "overlap") {
                    denom <- min(length(a), length(b))
                    sim_mat[i, j] <- if (denom == 0) 0 else inter / denom
                }
            }
        }
        legend_label <- switch(metric,
            jaccard = "Jaccard\nsimilarity",
            overlap = "Overlap\ncoefficient"
        )
        limits <- c(0, 1)
    }

    # ── cluster and reshape ───────────────────────────────────────────────────
    hc     <- hclust(as.dist(1 - sim_mat))
    ord    <- hc$order
    sm_ord <- sim_mat[ord, ord]

    df_long <- as.data.frame(sm_ord) %>%
        rownames_to_column("model_x") %>%
        pivot_longer(-model_x, names_to = "model_y", values_to = "sim") %>%
        mutate(
            model_x = factor(model_x, levels = rownames(sm_ord)),
            model_y = factor(model_y, levels = rownames(sm_ord))
        )

    # annotation: set sizes (not shown for logfc_cor)
    subtitle <- if (metric == "logfc_cor") {
        sprintf("Spearman correlation of logFC restricted to %d proteins significant in ≥1 model",
                length(sig_any))
    } else {
        "logFC Spearman correlation across all tested proteins"
    }

    # ── plot ──────────────────────────────────────────────────────────────────
    title <- switch(metric,
        jaccard  = "Overlap of significant proteins (Jaccard)",
        overlap  = "Overlap of significant proteins (Overlap coefficient)",
        logfc_cor = "logFC correlation across plates (Spearman)"
    )

    fig <- ggplot(df_long, aes(x = model_x, y = model_y, fill = sim)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        geom_text(aes(label = sprintf("%.2f", sim)), size = 2.8) +
        scale_fill_distiller(palette = "RdYlBu", direction = 1,
                             name    = legend_label,
                             limits  = limits) +
        labs(title    = title,
             subtitle = subtitle,
             x = NULL, y = NULL) +
        theme(
            axis.text.x  = element_text(angle = 45, hjust = 1),
            plot.title   = element_text(face = "bold"),
            plot.subtitle = element_text(size = 7, colour = "grey40")
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