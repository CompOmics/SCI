# plotting_utils.R
library(ggplot2)

# ── Colors ────────────────────────────────────────────────────────────────────
CONDITION_COLORS <- c(
    'K562'                 = '#9B72CF',
    'HEK293T'              = '#F4A7B9',
    'QC - Hela'            = '#F9D776',
    'Site specific QC'     = '#2E8B57',
    'Cellenone blank'      = '#E8735A',
    'Site specific blanks' = '#6DBF8A',
    'Processing control'   = '#5B8DB8',
    'Sorting control'      = '#A8CCE8',
    'QC - HYE A'           = '#D87A21',
    'QC - HYE B'           = '#AF6D09'
)

LAB_COLORS <- c(
    'Astral_1'  = '#6aaed6',
    'Astral_2'  = '#f0b84a',
    'Astral_3'  = '#5bbfaa',
    'Astral_4'  = '#f0956e',
    'Astral_5'  = '#dda0cc',
    'Astral_6'  = '#d4b48c',
    'timsTOF_1' = '#fcd4ee',
    'timsTOF_2' = '#b8b8b8',
    'timsTOF_3' = '#f0e87a'
)

# ── Orderings ─────────────────────────────────────────────────────────────────
CONDITION_ORDER <- c(
    'K562', 'HEK293T', 'QC - Hela', 'Site specific QC',
    'Sorting control', 'Processing control',
    'Cellenone blank', 'Site specific blanks'
)

LAB_ORDER <- c(
    'Astral_1', 'Astral_2', 'Astral_3', 'Astral_4',
    'Astral_5', 'Astral_6', 'timsTOF_1', 'timsTOF_2', 'timsTOF_3'
)

# ── Lab name mapping ──────────────────────────────────────────────────────────
LAB_NAME_MAPPING <- c(
    'DTU_Astral_1'   = 'Astral_1',
    'LSMBO_Astral_1' = 'Astral_2',
    'KI_Astral_1'    = 'Astral_3',
    'Wyss_Astral_1'  = 'Astral_4',
    'IMP_Astral_1'   = 'Astral_5',
    'KI_Astral_2'    = 'Astral_6',
    'IMP_timsTOF_1'  = 'timsTOF_1',
    'Broad_timsTOF_1'= 'timsTOF_2',
    'BYU_timsTOF_1'  = 'timsTOF_3'
)

# ── ggplot2 theme ─────────────────────────────────────────────────────────────
theme_sci <- function(base_size = 11) {
    theme_classic(base_size = base_size) %+replace%
    theme(
        panel.grid        = element_blank(),
        strip.background  = element_blank(),
        strip.text        = element_text(face = 'bold'),
        plot.title        = element_text(face = 'bold', hjust = 0),
        legend.key        = element_blank(),
        legend.background = element_blank(),
        axis.line         = element_line(colour = 'black', linewidth = 0.4),
        axis.ticks        = element_line(colour = 'black', linewidth = 0.4),
    )
}

theme_set(theme_sci())

# ── Save utility ──────────────────────────────────────────────────────────────
save_figure <- function(plot = last_plot(), path, width = 7, height = 5) {
    base <- tools::file_path_sans_ext(basename(path))
    root <- dirname(path)

    for (fmt in c('png', 'svg', 'pdf')) {
        out_dir <- file.path(root, fmt)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        out <- file.path(out_dir, paste0(base, '.', fmt))
        dev <- switch(fmt,
            pdf = cairo_pdf,
            svg = svglite::svglite,
            fmt
        )
        ggplot2::ggsave(
            filename = out,
            plot     = plot,
            width    = width,
            height   = height,
            bg       = 'transparent',
            device   = dev
        )
    }
}

# ── Save legend ───────────────────────────────────────────────────────────────
save_legend <- function(color_mapping, save_dir, title = 'Legend',
                        filename = 'legend', ncol = 1) {
    df <- data.frame(
        label = factor(names(color_mapping), levels = names(color_mapping)),
        x     = 1
    )
    p <- ggplot(df, aes(x = x, fill = label)) +
        geom_point(aes(colour = label), shape = 15, size = 4, alpha = 0) +
        scale_colour_manual(values = color_mapping, name = title) +
        guides(colour = guide_legend(ncol = ncol, override.aes = list(alpha = 1, size = 4))) +
        theme_void() +
        theme(legend.position = 'center',
              legend.background = element_blank())

    legend_only <- cowplot::get_legend(p)
    h           <- ceiling(length(color_mapping) / ncol) * 0.35 + 0.3

    out <- cowplot::plot_grid(legend_only)
    save_figure(out, file.path(save_dir, filename), width = 3, height = h)
}