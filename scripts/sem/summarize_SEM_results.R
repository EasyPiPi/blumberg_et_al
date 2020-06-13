library(tidyverse)

root_dir <- "~/Desktop/github_repo/blumberg_et_al"

sem_output_dir <- file.path(root_dir, "output/sem/output")
dir.create(sem_output_dir, showWarnings = FALSE, recursive = TRUE)
sem_fig_dir <- file.path(root_dir, "output/sem/figure")
dir.create(sem_fig_dir, showWarnings = FALSE, recursive = FALSE)

sem_df_names <- list.files(sem_output_dir, pattern = ".csv")

sem_dfs <- map(list.files(sem_output_dir, pattern = ".csv", full.names = TRUE),
    read_csv)
names(sem_dfs) <- str_split(sem_df_names, pattern = "_out", simplify = TRUE)[, 1]

create_plot_df <- function(df) {
    df <- df %>%
        mutate(
            covariate = str_replace(covariate, "_", " "),
            covariate = str_replace(covariate, "gc", "G+C"),
            covariate = str_replace(covariate, "exonJunDen", "spl.junc.den."),
            covariate = str_replace(covariate, "utr", "'UTR")
        ) %>%
        mutate(significance = case_when(
            pvalue < 0.0005 ~ "***",
            pvalue < 0.005 ~ "**",
            pvalue < 0.05 ~ "*",
            TRUE ~ ""
        )) %>%
        mutate(response = factor(response, levels = c("Transcription_rate", "Half_life")))
    return(df)
}

create_lineplot <- function(df, error = TRUE) {
    p <- df %>%
        ggplot(aes(x=covariate, y=est))

    if (error) {
        p <- p + geom_errorbar(data = df,
                      mapping = aes(ymin=ci.lower, ymax=ci.upper), width=.1)
    }
    p +
        geom_point(data = df) +
        geom_hline(yintercept=0, linetype="dashed", color = "grey", size=1) +
        geom_text(aes(x = covariate, y = ci.upper + 0.1, label = significance), size = 5) +
        labs(x = "", y = "Coefficient Value") +
        facet_wrap(~response, nrow = 2, labeller = labeller(response = response.lab)) +
        theme_classic(base_size = 15) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_rect(colour="white", fill="white"))


}

response.lab <- c("Transcription", "Half-life")
names(response.lab) <- c("Transcription_rate", "Half_life")

plots <- map(map(sem_dfs, create_plot_df), create_lineplot)
walk2(str_c(names(plots), ".png"),
      plots,
      ggsave, path = sem_fig_dir, width = 6, height = 6)

p <- create_lineplot(
    create_plot_df(sem_dfs$protein_coding_spliced_full_feature_no_replicates),
    error = FALSE)

ggsave(
    filename = "protein_coding_spliced_full_feature_no_replicates.png",
    plot = p, path = sem_fig_dir, width = 6, height = 6)
