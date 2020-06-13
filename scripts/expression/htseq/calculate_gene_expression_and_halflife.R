library(tidyverse)
library(broom)
library(MatchIt)

#### calculate gene expression using concat libraries and longest transcript annotation ####
root_dir <- "~/Desktop/github_repo/blumberg_et_al"

tx_features_dir <- file.path(root_dir, "data/features")
exp_dir <- file.path(root_dir, "data/htseq/concat/longest_tx")
length_dir <- file.path(root_dir, "data/length")
id_dir <- file.path(root_dir, "data/id_mapper")

hl_dir <- file.path(root_dir, "output/half_life")
dir.create(file.path(hl_dir, "table"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(hl_dir, "figure"), showWarnings = FALSE, recursive = TRUE)
exp_fig_dir <- file.path(root_dir, "output/expression/figure")
dir.create(exp_fig_dir, showWarnings = FALSE, recursive = TRUE)
exp_df_dir <- file.path(root_dir, "output/expression/table")
dir.create(exp_df_dir, showWarnings = FALSE, recursive = TRUE)

# read in id mapper and gene types
tx2gene <- read_csv(file.path(id_dir, "tx2gene.csv"))
genetype <- read_csv(file.path(id_dir, "genetype.csv"))

# read in read count dfs (using HTseq and longest transcript annotation)
exp_df_names <- list.files(exp_dir)

read_count_dfs <- str_split(exp_df_names, "\\.", simplify = T)
colnames(read_count_dfs) <- c("cell_line", "assay", "RNA_type", "resource", "gtf", "txt")

read_count_dfs <- as_tibble(read_count_dfs) %>%
    add_column(file_names = list.files(exp_dir, full.names = T)) %>%
    select(-txt) %>%
    mutate(read_count_df = map(file_names, read_delim, delim = "\t", col_names = c("ensembl_transcript_id", "read_count")))

# filtering rows which don't belong to genes
no_count_rownames <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")

read_count_dfs <- read_count_dfs %>%
    # move the no count rows to another col
    mutate(no_count_df = map(read_count_df, . %>% filter(ensembl_transcript_id %in% no_count_rownames))) %>%
    # remove the no count rows
    mutate(read_count_df = map(read_count_df, . %>% filter(!ensembl_transcript_id %in% no_count_rownames))) %>%
    mutate(assay = case_when(
        str_detect(read_count_dfs$gtf, "exbytx") ~ "RNAseq_ex",
        str_detect(read_count_dfs$gtf, "inbytx") ~ "RNAseq_in",
        TRUE ~ "PROseq"
    ))

# read in length dfs
len_df_names <- list.files(length_dir, full.names = T)

len_dfs <- tibble(file_names = len_df_names) %>%
    mutate(id = str_remove(basename(file_names), ".csv")) %>%
    mutate(len_df = map(file_names, read_csv)) %>%
    mutate(assay = case_when(
        id == "human_gene_len_bothendtrim500" ~ "PROseq",
        id == "human_in_len" ~ "RNAseq_in",
        id == "human_tx_len" ~ "RNAseq_ex"))

# calculate TPM
get_TPM <- function(read_count_df, len_df){
    left_join(read_count_df, len_df, by = "ensembl_transcript_id") %>%
        # 3th col is the gene length or tx length
        mutate(RPK = read_count / (.[[3]] / 1000)) %>%
        filter(RPK > 0, is.finite(RPK)) %>%
        mutate(TPM = RPK / (sum(RPK, na.rm = TRUE) / 1e6)) %>%
        select(ensembl_transcript_id, TPM)
}

TPM_dfs <- left_join(read_count_dfs, len_dfs, by = "assay") %>%
    mutate(TPM_df = map2(read_count_df, len_df, get_TPM)) %>%
    select(cell_line, assay, RNA_type, resource, gtf, TPM_df)

TPM_dfs <- rbind(TPM_dfs, TPM_dfs %>%
                     filter(resource == "2014NG") %>%
                     mutate(resource = "2014NG_polyA"))

TPM_dfs$group <- TPM_dfs$resource

TPM_dfs <- TPM_dfs %>%
    mutate(group = case_when(
        cell_line == "Hela" ~ "Hela_2017NG",
        ((cell_line == "K562") & (resource != "Amit") & (RNA_type == "total_RNA" | resource == "2014NG")) ~ "K562_2014NG_total",
        ((cell_line == "K562") & (resource != "Amit") & (RNA_type == "polyA_RNA" | resource == "2014NG_polyA")) ~ "K562_2014NG_polyA",
        ((cell_line == "K562") & (resource == "Amit")) ~ "K562_Amit"
    ))

TPM_melted_dfs <- TPM_dfs %>%
    select(group, assay, TPM_df) %>%
    spread(assay, TPM_df)

# # read published transcription rate
# trp_rate <- read_delim(
#     file.path(root_dir,
#               "data/transcription_rate/wachutka_et_al/elife-45056-supp3-v2"),
#     delim = "\t") %>%
#     rename(transcription_rate = `synthesis [1/cell/min]`) %>%
#     select(gene_id, transcription_rate) %>%
#     mutate(gene_id = map_chr(str_split(gene_id, "\\."), 1))
#
# df <- TPM_melted_dfs$PROseq[[4]]
# df <- df %>%
#     left_join(tx2gene, by = "ensembl_transcript_id") %>%
#     left_join(trp_rate, by = c("ensembl_gene_id" = "gene_id")) %>%
#     filter(TPM > 10) %>%
#     na.omit() %>%
#     select_if(is.numeric) %>%
#     mutate_all(log2)
# df %>%
#     ggplot(aes(x = TPM, y=transcription_rate)) +
#     geom_point()
# df %>% cor(method = "spearman", use = "pairwise.complete.obs")

#### Visualize PRO-seq vs. RNA-seq for different type of genes ####
create_plot_df <- function(proseq_df, rnaseq_df) {
    plot_df <- proseq_df %>%
        left_join(rnaseq_df, by = "ensembl_transcript_id",
                  suffix = c("_PROseq", "_RNAseq")) %>%
        left_join(tx2gene, by = "ensembl_transcript_id") %>%
        left_join(genetype, by = "ensembl_gene_id") %>%
        # filter not expressed genes
        filter(TPM_PROseq > 0, TPM_RNAseq > 0) %>%
        mutate_at(c("TPM_PROseq", "TPM_RNAseq"), log2) %>%
        filter(!is.na(biotype))
    return(plot_df)
}

create_TPM_scatter <- function(df, title, slope, corr,
                               xlab = "(PRO-seq TPM)",
                               ylab = "(RNA-seq TPM)",
                               alpha = 0.1, xlim = c(-8, 12), ylim = c(-8, 12)) {
    p <- df %>%
        ggplot(aes(x = TPM_PROseq, y = TPM_RNAseq)) +
        geom_point(alpha = alpha, color = "steelblue3") +
        geom_smooth(method ="lm", se = TRUE, color = "red") +
        coord_fixed(ratio = 1, xlim = xlim, ylim = ylim) +
        annotate(geom="text", x = -5, y = 10,
                 label = bquote(paste(italic("r"), " = ", .(round(corr, 2)))),
                 size = 5) +
        annotate(geom="text", x = -5, y = 8.5,
                 label = bquote(paste("slope = ",.(round(slope, 2)))),
                 size = 5) +
        annotate(geom="text", x = -5, y = 7,
                 label = bquote(paste(italic("n"), " = ", .(nrow(df)))),
                 size = 5) +
        labs(x = bquote(paste(log[2], .(xlab), sep = "")),
             y = bquote(paste(log[2], .(ylab), sep = "")),
             title = str_replace(title, "_", " ")) +
        theme_classic(base_size = 15) +
        theme(plot.title = element_text(hjust = 0.5))
    p
}

get_TPM_scatter <- function(plot_df, ...) {
    group_var <- enquos(...)

    plot_df %>%
        group_by(!!! group_var) %>%
        nest() %>%
        mutate(
            # build linear model for each gene biotype
            lm_models = map(data, ~ lm(TPM_RNAseq ~ TPM_PROseq, data = .)),
            # tidy lm results
            coefs = map(lm_models, ~ tidy(.)),
            # retrieve slopes
            slopes = map_dbl(coefs, function(df)
                df %>% filter(term == "TPM_PROseq") %>% pull(estimate))
        ) %>%
        mutate(corrs = map_dbl(data, ~ cor(.[c("TPM_PROseq", "TPM_RNAseq")])["TPM_PROseq", "TPM_RNAseq"])) %>%
        mutate(plots = pmap(list(data, biotype, slopes, corrs), create_TPM_scatter)) %>%
        select(!!! group_var, plots) %>%
        ungroup()
}

# plots for different gene categories
plot_dfs <- TPM_melted_dfs %>%
    mutate(plot_PR = map2(PROseq, RNAseq_ex, create_plot_df)) %>%
    mutate(subplot = map(plot_PR, get_TPM_scatter, biotype)) %>%
    mutate(plot_IE = map2(RNAseq_in, RNAseq_ex, create_plot_df))

plot_outputs <- plot_dfs %>%
    select(group, subplot) %>%
    unnest(cols = c(subplot)) %>%
    mutate(file_name = str_c(group, "_", biotype, ".png"))

walk2(plot_outputs$file_name, plot_outputs$plots,
      ggsave, path = exp_fig_dir, width = 5, height = 5)

# plot for spliced and non-spliced transcript
# process tx feature df
tx_features <-read_csv(file.path(tx_features_dir, "human_gene_features.csv"))
tx_spliced <- tx_features %>%
    mutate(spliced = ifelse(nexon > 1, "intron-containing", "intron-less")) %>%
    select(tx_name, spliced)

plot_dfs <- plot_dfs %>%
    mutate(plot_PR_spliced = map(plot_PR, ~ .x %>%
    left_join(tx_spliced, by = c(ensembl_transcript_id = "tx_name")))) %>%
    mutate(subplot_spliced = map(plot_PR_spliced, get_TPM_scatter, biotype, spliced))

plot_spliced_outputs <- plot_dfs %>%
    select(group, subplot_spliced) %>%
    unnest(cols = c(subplot_spliced)) %>%
    mutate(file_name = str_c(group, "_", biotype, "_", spliced, ".png"))

walk2(plot_spliced_outputs$file_name, plot_spliced_outputs$plots,
      ggsave, path = exp_fig_dir, width = 5, height = 5)

# plots for all TUs (PRO-seq and RNA-seq)
plot_all_TUs_PR <- plot_dfs %>%
    mutate(
        # build linear model for each gene biotype
        lm_models = map(plot_PR, ~ lm(TPM_RNAseq ~ TPM_PROseq, data = .)),
        # tidy lm results
        coefs = map(lm_models, ~ tidy(.)),
        # retrieve slopes
        slopes = map_dbl(coefs, function(df)
            df %>% filter(term == "TPM_PROseq") %>% pull(estimate))
    ) %>%
    mutate(corrs = map_dbl(plot_PR, ~ cor(.[c("TPM_PROseq", "TPM_RNAseq")])["TPM_PROseq", "TPM_RNAseq"])) %>%
    mutate(plots = pmap(list(plot_PR, "all TUs", slopes, corrs), create_TPM_scatter, xlim = c(-8, 12))) %>%
    select(group, plots)  %>%
    mutate(file_name = str_c(group, "_", "all_TUs_PR", ".png"))

walk2(plot_all_TUs_PR$file_name, plot_all_TUs_PR$plots,
      ggsave, path = exp_fig_dir, width = 5, height = 5)

# match mRNA and lincRNA PRO-seq expression in Amit's data
df <- plot_dfs %>% filter(group == "K562_Amit") %>% select(plot_PR) %>% pull()
df <- df[[1]] %>%
    filter(biotype %in% c("protein_coding", "lincRNA")) %>%
    mutate(biotype = ifelse(biotype == "protein_coding", 0, 1))

m.out <- matchit(biotype ~ TPM_PROseq, data = df,  method = "nearest")
m.data <- match.data(m.out)

m.data %>%
    mutate(biotype = ifelse(biotype == 0, "protein coding", "lincRNA")) %>%
    ggplot(aes(x = TPM_PROseq, y = TPM_RNAseq, col = biotype)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm") +
    labs(x = bquote(paste(log[2], "(PRO-seq TPM)", sep = "")),
         y = bquote(paste(log[2], "(RNA-seq TPM)", sep = ""))) +
    theme_classic(base_size = 15)

ggsave("K562_Amit_total_matched_PROseq_TUs_PR.png", path = exp_fig_dir, width = 7, height = 5)

# plots for all TUs (RNA-seq intronic and exonic reads)
plot_all_TUs_IE <- plot_dfs %>%
    mutate(
        # build linear model for each gene biotype
        lm_models = map(plot_IE, ~ lm(TPM_RNAseq ~ TPM_PROseq, data = .)),
        # tidy lm results
        coefs = map(lm_models, ~ tidy(.)),
        # retrieve slopes
        slopes = map_dbl(coefs, function(df)
            df %>% filter(term == "TPM_PROseq") %>% pull(estimate))
    ) %>%
    mutate(corrs = map_dbl(plot_IE, ~ cor(.[c("TPM_PROseq", "TPM_RNAseq")])["TPM_PROseq", "TPM_RNAseq"])) %>%
    mutate(plots = pmap(list(plot_IE, "all TUs", slopes, corrs, xlab = "(RNA-seq Intron TPM)",
                             ylab = "(RNA-seq Exon TPM)"), create_TPM_scatter, xlim = c(-12, 12), ylim = c(-12, 12))) %>%
    select(group, plots)  %>%
    mutate(file_name = str_c(group, "_", "all_TUs_IE", ".png"))

walk2(plot_all_TUs_IE$file_name, plot_all_TUs_IE$plots,
      ggsave, path = exp_fig_dir, width = 5, height = 5)

#### calculate half-life ####
get_half_life <- function(proseq_df, rnaseq_df){
    proseq_df %>%
        left_join(rnaseq_df, suffix = c("_PROseq", "_RNAseq"), by = "ensembl_transcript_id") %>%
        # filter by expression
        filter(TPM_PROseq >= 10, TPM_RNAseq > 1) %>%
        mutate(half_life = TPM_RNAseq / TPM_PROseq) %>%
        filter(!is.na(half_life))
}

hl_dfs <- TPM_melted_dfs %>%
    mutate(hl_PR = map2(PROseq, RNAseq_ex, get_half_life)) %>%
    mutate(hl_IE = map2(RNAseq_in, RNAseq_ex, get_half_life))

# Visualize gene expression level
# hl_dfs %>%
#     select(group, PROseq) %>%
#     unnest(PROseq) %>%
#     ggplot(aes(x = group, y = log2(TPM + 1))) +
#     geom_violin() +
#     geom_hline(yintercept = log2(11), color = "red")
#
# hl_dfs %>%
#     select(group, RNAseq_ex) %>%
#     unnest(RNAseq_ex) %>%
#     ggplot(aes(x = group, y = log2(TPM + 1))) +
#     geom_violin() +
#     geom_hline(yintercept = log2(2), color = "red")

# # write out tables
# walk2(hl_dfs$hl_PR, file.path(hl_dir, "table", str_c("half_life_PR_", hl_dfs$group, "_RNA.csv")), write_csv)
# walk2(hl_dfs$hl_IE, file.path(hl_dir, "table", str_c("half_life_IE_", hl_dfs$group, "_RNA.csv")), write_csv)

# exclude histone genes in polyA RNA-seq estimates which have not polyA tail
histone_gene <-
    read_delim(file.path(root_dir,
                         "data/gene_lists/histone_genes/HGNC_histone_20190404.txt"), delim = "\t")
histone_geneName <- histone_gene[["Ensembl gene ID"]]

create_hl_output <- function(hl_type = "hl_PR") {
    hl_type <- enquo(hl_type)
    hl_output_df <- hl_dfs %>%
        select(group, !!hl_type) %>%
        unnest(cols = c(!!hl_type)) %>%
        left_join(tx2gene, by = "ensembl_transcript_id") %>%
        select(group, ensembl_gene_id, half_life) %>%
        spread(group, half_life)
    hl_output_df <- hl_output_df %>%
        # filter histone genes in polyA RNA-seq libraries
        filter(ensembl_gene_id %in% histone_geneName) %>%
        mutate_at(vars(contains("polyA")), ~ NA) %>%
        rbind(hl_output_df %>% filter(!ensembl_gene_id %in% histone_geneName))
    return(hl_output_df)
}

hl_PR <- create_hl_output(hl_type = "hl_PR")
hl_IE <- create_hl_output(hl_type = "hl_IE")

hl_output <- hl_PR %>%
    full_join(hl_IE, by = "ensembl_gene_id", suffix = c("_PR", "_IE"))

write_csv(hl_output, file.path(hl_dir, "table", "half_life_concat_libraries.csv"))

# Visualization RNA half-life
# spliced vs non-spliced tx half-life
# only for Amit's data
df <- hl_dfs[hl_dfs$group == "K562_Amit", "hl_PR"] %>% unnest(cols = c(hl_PR))

df <- df %>%
    left_join(tx2gene, by = "ensembl_transcript_id") %>%
    left_join(genetype, by = "ensembl_gene_id")

df %>%
    left_join(tx_spliced, by = c(ensembl_transcript_id = "tx_name")) %>%
    filter(biotype %in% c("protein_coding", "lincRNA")) %>%
    mutate(log2half_life = log2(half_life)) %>%
    mutate(biotype = str_replace(biotype, "_", " ")) %>%
    ggplot(aes(x = biotype, y = log2half_life, col = spliced)) +
    geom_boxplot() +
    labs(x = "",
         y = bquote(paste(log[2], "(T"["1/2"] ^PR, ")", sep = "")),
         col = "Category") +
    theme_classic(base_size = 15)

ggsave("K562_Amit_halflife_PR_spliced_vs_nonspliced.png", path = file.path(hl_dir, "figure"), width = 7, height = 5)

# export expression for data generated by amit blumberg
df %>%
    na.omit() %>%
    write_csv(file.path(hl_dir, "table", "half_life_K562_Amit_total_RNA.csv"))

amit_proseq <- TPM_melted_dfs %>% filter(group == "K562_Amit") %>%
    select(PROseq) %>% unnest(cols = c(PROseq))
amit_rnaseq <- TPM_melted_dfs %>% filter(group == "K562_Amit") %>%
    select(RNAseq_ex) %>% unnest(cols = c(RNAseq_ex))

amit_exp <- amit_proseq %>%
    left_join(amit_rnaseq, by = "ensembl_transcript_id", suffix = c("_PROseq", "_RNAseq")) %>%
    left_join(tx2gene, by = "ensembl_transcript_id") %>%
    left_join(genetype, by = "ensembl_gene_id")

amit_exp %>%
    write_csv(file.path(exp_df_dir, "expression_K562_Amit_total_RNA.csv"))

# correct PRO-seq with elongation rate
elg_rate <- read_csv(file.path(root_dir, "data/elongation_rate/veloso_et_al/K562_elongation_rate.csv"))
colnames(elg_rate) <- c("ensembl_gene_id", "elongation_rate", "expression")
elg_rate <- elg_rate[elg_rate$elongation_rate > 0, ] %>% select(ensembl_gene_id, elongation_rate)

elg_rate_alex <- read_csv(file.path(root_dir,
                                    "data/elongation_rate/alex/elongation_rates.csv"))
elg_rate_alex <- elg_rate_alex %>%
    left_join(tx2gene, by = "ensembl_transcript_id") %>%
    select(ensembl_gene_id, everything(), -ensembl_transcript_id)

amit_exp <- amit_exp %>% select(ensembl_gene_id, biotype, TPM_PROseq, TPM_RNAseq) %>%
    left_join(elg_rate, by = "ensembl_gene_id") %>%
    left_join(elg_rate_alex, by = "ensembl_gene_id") %>%
    mutate(TPM_PROseq_crt1 = TPM_PROseq * 1000 / elongation_rate,
           TPM_PROseq_crt2 = TPM_PROseq * 1000 / Rate_60to120min)

amit_exp_crt1 <- amit_exp %>%
    select(TPM_PROseq, TPM_PROseq_crt1, TPM_RNAseq) %>%
    na.omit() %>%
    log2()

amit_exp_crt2 <- amit_exp %>%
    select(TPM_PROseq, TPM_PROseq_crt2, TPM_RNAseq) %>%
    na.omit() %>%
    log2()

amit_exp_crt1 %>% cor()
amit_exp_crt2 %>% cor()

make_scatter <- function(df, col1, col2, lab_x, lab_y, corr, n,
                         xlim = c(-5, 10), ylim = c(-5, 10)) {
    ggplot(df, aes_string(col1, col2)) +
        geom_point(alpha = 0.3, color = "steelblue3") +
        geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(x = bquote(paste(log[2], .(lab_x))),
             y = bquote(paste(log[2], .(lab_y)))) +
        coord_fixed(ratio = 0.8, xlim = xlim, ylim = ylim) +
        annotate(geom="text", x = -3, y = 9,
                 label = bquote(paste(italic("r"), " = ", .(round(corr, 2)))),
                 size = 5) +
        annotate(geom="text", x = -3, y = 8,
                 label = bquote(paste(italic("n"), " = ", .(n))),
                 size = 5) +
        theme_classic(base_size = 15)
}

make_scatter(amit_exp_crt1, "TPM_PROseq", "TPM_RNAseq", "(PRO-seq TPM)", "(RNA-seq TPM)", 0.55, 2204)
ggsave(file.path(exp_fig_dir, "K562_Amit_PROseq_vs_RNAseq.png"), width = 5, height = 5)

make_scatter(amit_exp_crt1, "TPM_PROseq_crt1", "TPM_RNAseq", "(PRO-seq TPM corrected by ER)", "(RNA-seq TPM)", 0.48, 2204)
ggsave(file.path(exp_fig_dir, "K562_Amit_PROseq_with_ER_correction_vs_RNAseq.png"), width = 5, height = 5)

make_scatter(amit_exp_crt2, "TPM_PROseq", "TPM_RNAseq", "(PRO-seq TPM)", "(RNA-seq TPM)", 0.87, 1543,
             xlim = c(-7, 10), ylim = c(-7, 10))
ggsave(file.path(exp_fig_dir, "K562_Amit_PROseq_2_vs_RNAseq.png"), width = 5, height = 5)

make_scatter(amit_exp_crt2, "TPM_PROseq_crt2", "TPM_RNAseq", "(PRO-seq TPM corrected by ER)", "(RNA-seq TPM)", 0.66,
             1543, xlim = c(-7, 10), ylim = c(-7, 10))
ggsave(file.path(exp_fig_dir, "K562_Amit_PROseq_with_ER_correction_2_vs_RNAseq.png"), width = 5, height = 5)


