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

    p <- ggplot(res, aes(x = logFC, y = neg_log10_p)) +
        geom_point(aes(colour = direction), size = 0.8, alpha = 0.6) +
        scale_colour_manual(
            values = c(up = col_sig, down = col_sig, ns = col_ns),
            guide  = "none"
        ) +
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

    fig <- upset(
        fromList(sig_lists),
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
    fdr_threshold = .001,
    fc_threshold  = 0.5,
    ncol          = 3,
    save_dir      = NULL,
    filename      = "volcano_grid_per_plate",
    width         = 18,
    height        = 12,
    show          = TRUE
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

        ggplot(res, aes(x = logFC, y = neg_log10_p, colour = colour_group)) +
            geom_point(size = 0.7, alpha = 0.6) +
            scale_colour_manual(
                values = c(
                    "both"             = "#d73027",   # red — significant in both
                    "plate only"       = "#f4a442",   # orange — plate only
                    "integrated only"  = "#4393c3",   # blue — integrated only
                    "ns"               = "grey85"     # light grey — not significant
                ),
                name   = NULL,
                guide  = guide_legend(
                    override.aes = list(size = 3, alpha = 1)  # make legend points larger
                )
            ) +
            geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                       linetype = "dashed", colour = "grey40", linewidth = 0.3) +
            geom_hline(yintercept = -log10(fdr_threshold),
                       linetype = "dashed", colour = "grey40", linewidth = 0.3) +
            annotate("text", x = Inf, y = Inf,
                     label = sprintf("both: %d\nplate: %d", n_both, n_plate),
                     hjust = 1.1, vjust = 1.5, size = 2.8, colour = "grey30") +
            xlim(x_lim) + ylim(y_lim) +
            labs(title = plate,
                 x     = expression(log[2]~FC),
                 y     = expression(-log[10]~p)) +
            theme(plot.title    = element_text(face = "bold", size = 9),
                  legend.position = "none")
    })

    fig <- wrap_plots(panels, ncol = ncol) +
        plot_annotation(
            title    = "Per-plate DEA — colored by integrated model",
            subtitle = sprintf("FDR < %.3f | |log2FC| > %.1f", fdr_threshold, fc_threshold),
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
        res$protein
    }

    for (direction in c("up", "down")) {
        sig_lists <- c(
            list("Integrated" = get_sig(res_integrated, direction)),
            lapply(results_per_plate, get_sig, direction = direction)
        )

        cat(sprintf("\n%sregulated proteins per model:\n",
                    ifelse(direction == "up", "Up", "Down")))
        print(sapply(sig_lists, length))

        fig <- upset(
            fromList(sig_lists),
            nsets           = length(sig_lists),
            order.by        = "freq",
            point.size      = 3,
            line.size       = 1,
            text.scale      = c(1.3, 1.2, 1.2, 1.1, 1.2, 1.2),
            main.bar.color  = bar_color,
            sets.bar.color  = bar_color,
            mainbar.y.label = "Intersection size",
            sets.x.label    = sprintf("%sregulated proteins",
                                      ifelse(direction == "up", "Up", "Down"))
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
    save_dir      = NULL,
    filename      = "overlap_heatmap",
    width         = 8,
    height        = 7,
    show          = TRUE
) {
    get_sig <- function(res) {
        res %>%
            filter(!is.na(adjPval), adjPval < fdr_threshold, abs(logFC) > fc_threshold) %>%
            pull(protein)
    }

    sig_lists <- c(
        list("Integrated" = get_sig(res_integrated)),
        lapply(results_per_plate, get_sig)
    )

    # compute pairwise Jaccard similarity
    models  <- names(sig_lists)
    n       <- length(models)
    overlap <- matrix(0, nrow = n, ncol = n, dimnames = list(models, models))

    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            a <- sig_lists[[i]]
            b <- sig_lists[[j]]
            overlap[i, j] <- length(intersect(a, b)) /
                              length(union(a, b))
        }
    }

    # cluster
    hc       <- hclust(as.dist(1 - overlap))
    ord      <- hc$order
    ov_ord   <- overlap[ord, ord]

    df_long <- as.data.frame(ov_ord) %>%
        rownames_to_column("model_x") %>%
        pivot_longer(-model_x, names_to = "model_y", values_to = "jaccard") %>%
        mutate(
            model_x = factor(model_x, levels = rownames(ov_ord)),
            model_y = factor(model_y, levels = rownames(ov_ord))
        )

    fig <- ggplot(df_long, aes(x = model_x, y = model_y, fill = jaccard)) +
        geom_tile(colour = "white", linewidth = 0.4) +
        geom_text(aes(label = sprintf("%.2f", jaccard)), size = 2.8) +
        scale_fill_distiller(palette = "RdYlBu", direction = 1,
                             name = "Jaccard\nsimilarity",
                             limits = c(0, 1)) +
        labs(
            title = "Overlap of significant proteins across plates",
            x     = NULL,
            y     = NULL
        ) +
        theme(
            axis.text.x  = element_text(angle = 45, hjust = 1),
            plot.title   = element_text(face = "bold")
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