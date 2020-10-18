library(tidyverse)
library(corrplot)
library(skimr)
library(RColorBrewer)

root_dir <- "~/Desktop/github_repo/blumberg_et_al"
hl_published_dir <-
    file.path(root_dir, "data/published_half_life")
hl_dir <- file.path(root_dir, "output/half_life/table")
id_dir <- file.path(root_dir, "data/id_mapper")

hl_figure_dir <- file.path(root_dir, "output/half_life/figure")
elg_dir <- file.path(root_dir, "data/elongation_rate/veloso_et_al")

genetype <- read_csv(file.path(id_dir, "genetype.csv"))
tx2gene <- read_csv(file.path(id_dir, "tx2gene.csv"))

# read in half-life dfs
hl_df <- read_csv(file.path(hl_dir, "half_life_concat_libraries.csv"))
hl_published_df <-
    read_csv(file.path(hl_published_dir, "published_half_lives.csv")) %>%
    rename("ensembl_gene_id" = X1)
hl_df <- hl_df %>%
    full_join(hl_published_df, by = "ensembl_gene_id") %>%
    left_join(genetype, by = "ensembl_gene_id")

# define handy functions
# corrplot
corrplot_partial <-
    partial(
        corrplot,
        method = "ellipse",
        addrect = 2,
        tl.srt = 45,
        cl.lim = c(0, 1),
        addCoef.col = "white",
        tl.col = "black",
        diag = FALSE
    )
# save output figure
save_fig <- function(file_path, fig_fun) {
    png(file_path, width = 1000, height = 1000)
    fig_fun
    dev.off()
}

# correlation
cor_spearman <-
    partial(cor, use = "pairwise.complete.obs", method = "spearman")

# draw heatmap for half-life comparisons
hl_corr <-
    cor_spearman(hl_df %>% select(-c(ensembl_gene_id, biotype)))
save_fig(
    file.path(hl_figure_dir, "half_life_all_comparison.png"),
    corrplot_partial(
        hl_corr,
        tl.cex = 1.2,
        cl.cex = 1.5,
        order = "hclust"
    )
)

hl_K562_corr <-
    cor_spearman(hl_df %>%
                     filter(biotype == "protein_coding") %>%
                     select(contains("K562")) %>%
                     # select(-contains("_IE"), -contains("2014NG")) %>% na.omit())
                     select(-contains("_IE"), -contains("Amit")) %>% na.omit())

nice_names <-
    c(
        "ENCODE polyA RNA-seq \n   + Core et al PRO-seq",
        #"ENCODE total RNA-seq \n   + Core et al PRO-seq",
        # "Blumberg et al (PRO-seq + RNA-seq)",
        # "ENCODE polyA RNA-seq \n   exon + intron",
        # "ENCODE total RNA-seq \n   exon + intron",
        # "total RNA-seq \n   exon + intron (this paper)",
        "Schofield et al (TimeLapse-seq)",
        "Mele et al (Actinomycin D)",
        "Wu et al (SLAM-seq)",
        "Wachutka et al (TT-seq)"
    )

names(nice_names) <- colnames(hl_K562_corr)
colnames(hl_K562_corr) <- nice_names
rownames(hl_K562_corr) <- nice_names

save_fig(
    file.path(
        hl_figure_dir,
        "half_life_K562_same_set_protein_coding_comparison.png"
    ),
    corrplot_partial(
        hl_K562_corr,
        tl.cex = 1.3,
        cl.cex = 1.8,
        number.cex = 1.5,
        order = NULL
    )
)

# sel_cols <-
#     c(
#         "K562_Amit_PR",
#         "K562_2014NG_polyA_PR",
#         "K562_2014NG_total_PR",
#         "K562_Amit_IE",
#         "K562_2014NG_polyA_IE",
#         "K562_2014NG_total_IE",
#         "Schofield_et_al_K562",
#         "Mele_et_al_K562"
#     )
# comb_cols <- t(combn(sel_cols, 2))
#
# create_scatter <- function(df, slope, corr) {
#     p <-
#         ggplot(df, aes_string(x = colnames(df)[1], y = colnames(df)[2])) +
#         geom_point(alpha = 0.1, color = "steelblue3") +
#         geom_smooth(method = "lm",
#                     se = TRUE,
#                     color = "red") +
#         coord_fixed(ratio = 1,
#                     xlim = c(-6, 4),
#                     ylim = c(-4, 6)) +
#         annotate(
#             geom = "text",
#             x = -4,
#             y = 5,
#             label = paste("ρ = ", round(corr, 2), "\n", "slope = ", round(slope, 2)),
#             size = 5
#         ) +
#         theme_classic(base_size = 15) +
#         theme(plot.title = element_text(hjust = 0.5))
#     p
# }
#
# for (i in seq_len(nrow(comb_cols))) {
#     df <- hl_df[comb_cols[i,]] %>% na.omit() %>% log2()
#     print(colnames(df))
#     spearman_r <- cor(df, method = "spearman")[1, 2]
#     slope <- unname(coef(lm(df[[1]] ~ df[[2]], data = df))[2])
#
#     ggsave(
#         filename = str_c(comb_cols[i, 1], "_vs_", comb_cols[i, 2], "_half_life.png"),
#         plot = create_scatter(df, slope, spearman_r),
#         path = hl_figure_dir,
#         width = 5,
#         height = 5
#     )
# }

# Comparison for Amit's data
# amit_PR <- c("K562_Amit_PR", "Schofield_et_al_K562")
# amit_IE <- c("K562_Amit_IE", "Schofield_et_al_K562")

NG_PR <- c("K562_2014NG_PR", "Schofield_et_al_K562")
NG_IE <- c("K562_2014NG_IE", "Schofield_et_al_K562")

scatter_contour_plot <- function(amit_cols, xlab, label, c1 = -6, c2 = 4, halflife_df = hl_df) {
    halflife_df <- halflife_df %>%
        select(amit_cols) %>% na.omit() %>% log2()

    print(paste("The number of observations:", NROW(halflife_df)))

    halflife_df %>%
        ggplot(mapping = aes_string(x = amit_cols[1], y = amit_cols[2])) +
        geom_point(alpha = 0.1, color = "red") +
        geom_density2d(aes(colour = ..level..), show.legend = FALSE) +
        coord_fixed(ratio = 1,
                    xlim = c(c1, c2),
                    ylim = c(-4, 6)) +
        annotate(
            geom = "text",
            x = c1 + 2,
            y = 5,
            label = label,
            size = 8
        ) +
        labs(x = bquote(paste(log[2], "(T"["1/2"] ^ .(xlab), ")", sep = "")),
             y = expression(paste(log[2], "(T"["1/2"] ^ "TLS", ", hours)", sep = ""))) +
        theme_classic(base_size = 20) +
        theme(plot.title = element_text(hjust = 0.5))
}

# ggsave(file.path(hl_figure_dir,
#                  "K562_Amit_total_vs_Schofield_et_al_K562_PR_paper_figure.png"
#     ), plot = scatter_contour_plot(tidyselect::all_of(amit_PR), "PR", "ρ = 0.59\nn=5068 "), width = 6, height = 5)
#
# ggsave(file.path(hl_figure_dir,
#                  "K562_Amit_total_vs_Schofield_et_al_K562_IE_paper_figure.png"
# ), plot = scatter_contour_plot(tidyselect::all_of(amit_IE), "IE", "ρ = 0.22\nn=3090", c1 = -4, c2 = 6), width = 5, height = 5)

ggsave(
    file.path(hl_figure_dir, "K562_2014NG_vs_Schofield_et_al_K562_PR_paper_figure.png"),
    plot = scatter_contour_plot(tidyselect::all_of(NG_PR), "PR", "n = 4351\nρ = 0.71"),
    width = 6, height = 5)

ggsave(
    file.path(hl_figure_dir, "K562_2014NG_vs_Schofield_et_al_K562_IE_paper_figure.png"),
    plot = scatter_contour_plot(tidyselect::all_of(NG_IE), "IE", "n = 3629\nρ = 0.47",
                                c1 = -4, c2 = 6), width = 5, height = 5)

# Considering Elongation rate
elg_rate <- read_csv(file.path(elg_dir, "K562_elongation_rate.csv"))
colnames(elg_rate) <- c("ensembl_gene_id", "elongation_rate", "expression")
elg_rate <- elg_rate[elg_rate$elongation_rate > 0, ]

elg_rate_alex <- read_csv(file.path(root_dir,
                                    "data/elongation_rate/alex/elongation_rates.csv"))
elg_rate_alex <- elg_rate_alex %>%
    left_join(tx2gene, by = "ensembl_transcript_id") %>%
    select(ensembl_gene_id, everything(), -ensembl_transcript_id)

hl_elg_df <- hl_df %>%
    # select(c("ensembl_gene_id", tidyselect::all_of(amit_PR))) %>%
    select(c("ensembl_gene_id", tidyselect::all_of(NG_PR))) %>%
    left_join(elg_rate, by = "ensembl_gene_id") %>%
    left_join(elg_rate_alex, by = "ensembl_gene_id") %>%
    # mutate(K562_Amit_total_PR_elg = K562_Amit_PR * 1000 / elongation_rate,
    #        K562_Amit_total_PR_60to120min = K562_Amit_PR * 1000 / Rate_60to120min)
    mutate(K562_2014NG_PR_elg = K562_2014NG_PR * 1000 / elongation_rate,
       K562_2014NG_PR_60to120min = K562_2014NG_PR * 1000 / Rate_60to120min) %>%
    select(-expression, -contains("rate"))

# amit_PR <- c("K562_Amit_PR", "Schofield_et_al_K562")
# amit_PR_elg <- c("K562_Amit_total_PR_elg", "Schofield_et_al_K562")
# amit_PR_elg_60to120min <- c("K562_Amit_total_PR_60to120min", "Schofield_et_al_K562")

NG_PR <- c("K562_2014NG_PR", "Schofield_et_al_K562")
NG_PR_elg <- c("K562_2014NG_PR_elg", "Schofield_et_al_K562")
NG_PR_elg_60to120min <- c("K562_2014NG_PR_60to120min", "Schofield_et_al_K562")

hl_elg_corr <-  cor(log2(hl_elg_df[2:ncol(hl_elg_df)]), use = "pairwise.complete.obs",
                    method = "spearman")

# ggsave(
#     file.path(hl_figure_dir,
#               "K562_Amit_total_vs_Schofield_et_al_K562_PR_without_elongation_rate_paper_figure.png"),
#     plot = scatter_contour_plot(amit_PR, "PR", paste0("ρ = ", as.character(round(hl_elg_corr["K562_Amit_PR", "Schofield_et_al_K562"], 2))), halflife_df = hl_elg_df), width = 5, height = 5)
#
# ggsave(file.path(hl_figure_dir,
#                  "K562_Amit_total_vs_Schofield_et_al_K562_PR_with_elongation_rate_paper_figure.png"),
#        plot = scatter_contour_plot(amit_PR_elg, "PR", paste0("ρ = ", as.character(round(hl_elg_corr["K562_Amit_total_PR_elg", "Schofield_et_al_K562"], 2)), "\nn=1573"), halflife_df = hl_elg_df), width = 5, height = 5)
#
# ggsave(file.path(hl_figure_dir,
#                  "K562_Amit_total_vs_Schofield_et_al_K562_PR_with_elongation_rate_60to120min_paper_figure.png"),
#        plot = scatter_contour_plot(amit_PR_elg_60to120min, "PR", paste0("ρ = ", as.character(round(hl_elg_corr["K562_Amit_total_PR_60to120min", "Schofield_et_al_K562"], 2)), "\nn=376"), halflife_df = hl_elg_df), width = 5, height = 5)

hl_elg_narm_df <- hl_elg_df %>% select(K562_2014NG_PR, Schofield_et_al_K562, K562_2014NG_PR_elg) %>% na.omit()

ggsave(
    file.path(
        hl_figure_dir,
        "K562_2014NG_polyA_vs_Schofield_et_al_K562_PR_without_elongation_rate_paper_figure.png"),
    plot = scatter_contour_plot(NG_PR,
                                "PR", paste0("n = 1286\n",
                                             "ρ = ",
                                             as.character(round(hl_elg_corr["K562_2014NG_PR",
                                                                            "Schofield_et_al_K562"], 2))
                                             ),
                                halflife_df = hl_elg_narm_df), width = 5, height = 5)

ggsave(
    file.path(hl_figure_dir,
            "K562_2014NG_polyA_vs_Schofield_et_al_K562_PR_with_elongation_rate_paper_figure.png"),
       plot = scatter_contour_plot(NG_PR_elg, "PR", paste0("n = 1286\n", "ρ = ", as.character(round(hl_elg_corr["K562_2014NG_PR_elg", "Schofield_et_al_K562"], 2))), halflife_df = hl_elg_narm_df), width = 5, height = 5)

ggsave(
    file.path(hl_figure_dir,
            "K562_2014NG_polyA_vs_Schofield_et_al_K562_PR_with_elongation_rate_60to120min_paper_figure.png"),
       plot = scatter_contour_plot(NG_PR_elg_60to120min, "PR", paste0("n = 282\n", "ρ = ", as.character(round(hl_elg_corr["K562_2014NG_PR_60to120min", "Schofield_et_al_K562"], 2))), halflife_df = hl_elg_df), width = 5, height = 5)

# cv and quantile
cv <- function(vec) {
    sd(vec, na.rm = TRUE) / mean(vec, na.rm = TRUE)
}

hl_el_cv <- hl_df %>%
    # select(c("ensembl_gene_id", tidyselect::all_of(amit_PR))) %>%
    select(c("ensembl_gene_id", tidyselect::all_of(NG_PR))) %>%
    left_join(elg_rate, by = "ensembl_gene_id") %>%
    left_join(elg_rate_alex, by = "ensembl_gene_id")

hl_el_cv %>%
    # select(K562_Amit_PR, elongation_rate) %>%
    select(K562_2014NG_PR, elongation_rate) %>%
    na.omit() %>%
    # map(cv)
    map(quantile, probs = c(0, 0.05, 0.1, 0.5, 0.9, 0.95, 1))

hl_el_cv %>%
    select(K562_2014NG_PR, Rate_60to120min) %>%
    na.omit() %>%
    # map(cv)
    map(quantile, probs = c(0, 0.05, 0.1, 0.5, 0.9, 0.95, 1))

# top 50% genes
# amit_hl <- read_csv(file.path(root_dir, "output/half_life/table/half_life_K562_Amit_total_RNA.csv"))

NG_hl <- read_csv(file.path(root_dir, "output/half_life/table/half_life_K562_2014NG_polyA_RNA.csv"))

NG_hl %>%
    filter(TPM_PROseq > median(TPM_PROseq)) %>%
    select(ensembl_gene_id, half_life) %>%
    left_join(hl_published_df, by = "ensembl_gene_id") %>%
    select(half_life, Schofield_et_al_K562) %>%
    na.omit() %>%
    cor(method = "spearman", use = "pairwise.complete.obs")
