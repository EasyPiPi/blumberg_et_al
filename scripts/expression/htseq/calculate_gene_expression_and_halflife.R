library(tidyverse)
library(broom)
library(MatchIt)
library(ggpubr)

#### calculate gene expression using concat libraries and dominant transcript annotation ####
root_dir <- "~/Desktop/github_repo/blumberg_et_al"

tx_features_dir <- file.path(root_dir, "data/features")
exp_dir <- file.path(root_dir, "data/htseq/concat/pc_linc_as")
# exp_dir <- file.path(root_dir, "data/htseq/concat/longest_tx")
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
tx2gene <- read_csv(file.path(id_dir, "tx2gene.csv"), col_types = cols())
genetype <- read_csv(file.path(id_dir, "genetype.csv"), col_types = cols())

# read in read count dfs (using HTseq and annotation of the dominant transcript)
exp_df_names <- list.files(exp_dir)

read_count_dfs <- str_split(exp_df_names, "\\.", simplify = T)
colnames(read_count_dfs) <- c("cell_line", "assay", "RNA_type", "resource", "gtf", "txt")

read_count_dfs <- as_tibble(read_count_dfs) %>%
    add_column(file_names = list.files(exp_dir, full.names = T)) %>%
    select(-txt) %>%
    mutate(read_count_df = map(file_names, read_delim,
                               delim = "\t", col_names = c("ensembl_transcript_id", "read_count"),
                               col_types = cols()))

# filtering rows which don't belong to genes
no_count_rownames <-
    c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")

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
    mutate(len_df = map(file_names, read_csv, col_types = cols())) %>%
    mutate(assay = case_when(
        id == "human_gene_len_bothendtrim500" ~ "PROseq",
        id == "human_in_len" ~ "RNAseq_in",
        id == "human_tx_len" ~ "RNAseq_ex"))

# calculate TPM
get_TPM <- function(read_count_df, len_df) {
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

TPM_dfs$group <- TPM_dfs$resource

TPM_dfs <- TPM_dfs %>%
    mutate(group = case_when(
        cell_line == "Hela" ~ "Hela_2017NG",
        ((cell_line == "K562") & (resource != "Amit") & (RNA_type == "polyA_RNA" | resource == "2014NG")) ~ "K562_2014NG",
        ((cell_line == "K562") & (resource == "Amit")) ~ "K562_Amit"
    ))

TPM_melted_dfs <- TPM_dfs %>%
    select(group, assay, TPM_df) %>%
    spread(assay, TPM_df)

#### Visualize PRO-seq vs. RNA-seq for different type of genes ####
create_plot_df <- function(proseq_df, rnaseq_df) {
    plot_df <- proseq_df %>%
        left_join(rnaseq_df, by = "ensembl_transcript_id",
                  suffix = c("_PROseq", "_RNAseq")) %>%
        left_join(tx2gene, by = "ensembl_transcript_id") %>%
        left_join(genetype, by = "ensembl_gene_id") %>%
        # filter not expressed genes
        filter(TPM_PROseq > 0, TPM_RNAseq > 0) %>%
        # filter genes without biotype
        filter(!is.na(biotype)) %>%
        # re-scaling for gene have non-zero expressions for both assays
        mutate(TPM_PROseq = TPM_PROseq / sum(TPM_PROseq, na.rm = TRUE) * 1e6,
               TPM_RNAseq = TPM_RNAseq / sum(TPM_RNAseq, na.rm = TRUE) * 1e6) %>%
        mutate_at(c("TPM_PROseq", "TPM_RNAseq"), log2)

    return(plot_df)
}

create_TPM_scatter <- function(df, title, slope, corr, corr_p,
                               xlab = "(PRO-seq TPM)",
                               ylab = "(RNA-seq TPM)",
                               alpha = 0.1, xlim = c(-8, 12), ylim = c(-8, 12)) {

    if (corr_p < 2.2e-16) {
        p_label <- bquote(paste(italic("p-value"), " < 2.2e-16"))
    } else {
        p_label <- bquote(paste(italic("p-value"), " = ", .(round(corr_p, 2))))
    }

    p <- df %>%
        # mutate_at(c("TPM_PROseq", "TPM_RNAseq"), log2) %>%
        ggplot(aes(x = TPM_PROseq, y = TPM_RNAseq)) +
        geom_point(alpha = alpha, color = "steelblue3") +
        # geom_smooth(method ="lm", se = TRUE, color = "red") +
        coord_fixed(ratio = 1, xlim = xlim, ylim = ylim) +
        annotate(geom="text", x = -5, y = 8.5,
                 # label = bquote(paste(italic("r"), " = ", .(round(corr, 2)))),
                 label = bquote(paste(italic("ρ"), " = ", .(round(corr, 2)))),
                 size = 5) +
        annotate(geom="text", x = -4, y = 7,
                 label = p_label,
                 size = 5) +
        # annotate(geom="text", x = -5, y = 8.5,
        #          label = bquote(paste("slope = ",.(round(slope, 2)))),
        #         size = 5) +
        annotate(geom="text", x = -5, y = 10,
                 label = bquote(paste("n", " = ", .(nrow(df)))),
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
        # mutate(corrs = map_dbl(data,
        #                       ~ cor(.[c("TPM_PROseq", "TPM_RNAseq")])["TPM_PROseq", "TPM_RNAseq"])) %>%
        mutate(corrs = map_dbl(data,
                               ~ cor(.[c("TPM_PROseq", "TPM_RNAseq")], method = "spearman")["TPM_PROseq", "TPM_RNAseq"]),
               corr_p = map_dbl(data,
                                ~ cor.test(.[["TPM_PROseq"]], .[["TPM_RNAseq"]], method = "spearman")$p.value)) %>%
        mutate(plots = pmap(list(data, biotype, slopes, corrs, corr_p), create_TPM_scatter)) %>%
        select(!!! group_var, plots) %>%
        ungroup()
}

# plots for different gene categories
suppressMessages(
    plot_dfs <- TPM_melted_dfs %>%
        mutate(plot_PR = map2(PROseq, RNAseq_ex, create_plot_df)) %>%
        mutate(subplot = map(plot_PR, get_TPM_scatter, biotype)) %>%
        mutate(plot_IE = map2(RNAseq_in, RNAseq_ex, create_plot_df))
)

plot_outputs <- plot_dfs %>%
    select(group, subplot) %>%
    unnest(cols = c(subplot)) %>%
    mutate(file_name = str_c(group, "_", biotype, ".png"))

suppressMessages(
    walk2(plot_outputs$file_name, plot_outputs$plots,
          ggsave, path = exp_fig_dir, width = 5, height = 5)
)
# plot for spliced and non-spliced transcript
# process tx feature df
tx_features <-read_csv(file.path(tx_features_dir, "human_gene_features.csv"), col_types = cols())
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

suppressMessages(
    walk2(plot_spliced_outputs$file_name, plot_spliced_outputs$plots,
          ggsave, path = exp_fig_dir, width = 5, height = 5)
)
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
    mutate(corrs = map_dbl(plot_PR, ~ cor(.[c("TPM_PROseq", "TPM_RNAseq")])["TPM_PROseq", "TPM_RNAseq"]),
           corr_p = map_dbl(plot_PR,
                            ~ cor.test(.[["TPM_PROseq"]], .[["TPM_RNAseq"]], method = "spearman")$p.value)) %>%
    mutate(plots = pmap(list(plot_PR, "all TUs", slopes, corrs, corr_p), create_TPM_scatter, xlim = c(-8, 12))) %>%
    select(group, plots)  %>%
    mutate(file_name = str_c(group, "_", "all_TUs_PR", ".png"))

suppressMessages(
    walk2(plot_all_TUs_PR$file_name, plot_all_TUs_PR$plots,
          ggsave, path = exp_fig_dir, width = 5, height = 5)
)

# match mRNA and lincRNA PRO-seq expression
# df <- plot_dfs %>% filter(group == "K562_Amit") %>% select(plot_PR) %>% pull()
df <- plot_dfs %>% filter(group == "K562_2014NG") %>% select(plot_PR) %>% pull()

df <- df[[1]] %>%
    filter(biotype %in% c("protein_coding", "lincRNA")) %>%
    mutate(biotype = ifelse(biotype == "protein_coding", 0, 1))

m.out <- matchit(biotype ~ TPM_PROseq, data = df,  method = "nearest")
m.data <- match.data(m.out)

p <- m.data %>%
    mutate(biotype = ifelse(biotype == 0, "protein coding", "lincRNA")) %>%
    ggplot(aes(x = TPM_PROseq, y = TPM_RNAseq, col = biotype)) +
    geom_point(alpha = 0.3) +
    # geom_smooth(method = "lm") +
    stat_cor(method = "spearman", cor.coef.name = "rho") +
    labs(x = bquote(paste(log[2], "(PRO-seq TPM)", sep = "")),
         y = bquote(paste(log[2], "(RNA-seq TPM)", sep = ""))) +
    theme_classic(base_size = 15)

# ggsave("K562_Amit_total_matched_PROseq_TUs_PR.png", path = exp_fig_dir, width = 7, height = 5)
suppressMessages(
    ggsave("K562_2014NG_total_matched_PROseq_TUs_PR.png", p, path = exp_fig_dir, width = 7, height = 5)
)

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
    mutate(corrs = map_dbl(plot_IE, ~ cor(.[c("TPM_PROseq", "TPM_RNAseq")])["TPM_PROseq", "TPM_RNAseq"]),
           corr_p = map_dbl(plot_PR,
                            ~ cor.test(.[["TPM_PROseq"]], .[["TPM_RNAseq"]], method = "spearman")$p.value)) %>%
    mutate(plots = pmap(list(plot_IE, "all TUs", slopes, corrs, corr_p, xlab = "(RNA-seq Intron TPM)",
                             ylab = "(RNA-seq Exon TPM)"), create_TPM_scatter, xlim = c(-12, 12), ylim = c(-12, 12))) %>%
    select(group, plots)  %>%
    mutate(file_name = str_c(group, "_", "all_TUs_IE", ".png"))

suppressMessages(
    walk2(plot_all_TUs_IE$file_name, plot_all_TUs_IE$plots,
          ggsave, path = exp_fig_dir, width = 5, height = 5)
)

#### calculate half-life ####
get_half_life <- function(proseq_df, rnaseq_df){
    proseq_df %>%
        left_join(rnaseq_df, suffix = c("_PROseq", "_RNAseq"), by = "ensembl_transcript_id") %>%
        # filter by expression
        filter(TPM_PROseq >= 10, TPM_RNAseq > 1) %>%
        # filter(TPM_PROseq >= 0, TPM_RNAseq > 0) %>%
        mutate(half_life = TPM_RNAseq / TPM_PROseq) %>%
        filter(!is.na(half_life))
}

hl_dfs <- TPM_melted_dfs %>%
    mutate(hl_PR = map2(PROseq, RNAseq_ex, get_half_life)) %>%
    mutate(hl_IE = map2(RNAseq_in, RNAseq_ex, get_half_life))

# exclude histone genes in polyA RNA-seq estimates (because the transcripts of these genes don't have polyA tails)
histone_gene <-
    read_delim(file.path(root_dir,
                         "data/gene_lists/histone_genes/HGNC_histone_20190404.txt"), delim = "\t",
               col_types = cols())
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
# only for 2014NG + ENCODE data
# df <- hl_dfs[hl_dfs$group == "K562_Amit", "hl_PR"] %>% unnest(cols = c(hl_PR))
df <- hl_dfs[hl_dfs$group == "K562_2014NG", "hl_PR"] %>% unnest(cols = c(hl_PR))

df <- df %>%
    left_join(tx2gene, by = "ensembl_transcript_id") %>%
    left_join(genetype, by = "ensembl_gene_id")

p <- df %>%
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

# ggsave("K562_Amit_halflife_PR_spliced_vs_nonspliced.png", path = file.path(hl_dir, "figure"), width = 7, height = 5)
suppressMessages(
    ggsave("K562_2014NG_halflife_PR_spliced_vs_nonspliced.png", p, path = file.path(hl_dir, "figure"), width = 7, height = 5)
)

# export expression for data from 2014NG + ENCODE
df <- df %>% na.omit()

df %>%
    # write_csv(file.path(hl_dir, "table", "half_life_K562_Amit_total_RNA.csv"))
    write_csv(file.path(hl_dir, "table", "half_life_K562_2014NG_polyA_RNA.csv"))

#### plot PRO-seq vs. RNA-seq ####
df <- df %>% mutate_if(is.numeric, log2)

p <- ggscatter(df, x = "TPM_PROseq", y = "half_life",
               facet.by = "biotype",
               color = "black", alpha = 0.2, # shape = 21, size = 3, # Points color, shape and size
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE, # Add confidence interval
               cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
               cor.coeff.args = list(method = "spearman", label.x = 3,
                                     label.sep = "\n", cor.coef.name = "rho"),
               xlab = expression("PRO-seq "*log[2]*"(TPM)"),
               ylab = expression(paste(log[2], "(T"["1/2"] ^ "PR", ")", sep = ""))
)

suppressMessages(
    ggsave(file.path(hl_dir, "figure", "K562_2014NG_polyA_RNA_PROseq_vs_PR_halflife.png"),
           plot = p, width = 10, height = 8)
)

# table(df$biotype)
# amit_proseq <- TPM_melted_dfs %>% filter(group == "K562_Amit") %>%
#     select(PROseq) %>% unnest(cols = c(PROseq))
# amit_rnaseq <- TPM_melted_dfs %>% filter(group == "K562_Amit") %>%
#     select(RNAseq_ex) %>% unnest(cols = c(RNAseq_ex))
#
# amit_exp <- amit_proseq %>%
#     left_join(amit_rnaseq, by = "ensembl_transcript_id", suffix = c("_PROseq", "_RNAseq")) %>%
#     left_join(tx2gene, by = "ensembl_transcript_id") %>%
#     left_join(genetype, by = "ensembl_gene_id")
#
# amit_exp %>%
#     write_csv(file.path(exp_df_dir, "expression_K562_Amit_total_RNA.csv"))

NG_proseq <- TPM_melted_dfs %>% filter(group == "K562_2014NG") %>%
    select(PROseq) %>% unnest(cols = c(PROseq))
ENCODE_rnaseq <- TPM_melted_dfs %>% filter(group == "K562_2014NG") %>%
    select(RNAseq_ex) %>% unnest(cols = c(RNAseq_ex))

NGENCODE_exp <- NG_proseq %>%
    left_join(ENCODE_rnaseq, by = "ensembl_transcript_id", suffix = c("_PROseq", "_RNAseq")) %>%
    left_join(tx2gene, by = "ensembl_transcript_id") %>%
    left_join(genetype, by = "ensembl_gene_id")

NGENCODE_exp %>%
    write_csv(file.path(exp_df_dir, "expression_K562_2014NG_polyA_RNA.csv"))

#### correct PRO-seq with elongation rate ####
# published data
elg_rate <- read_csv(file.path(root_dir, "data/elongation_rate/veloso_et_al/K562_elongation_rate.csv"),
                     col_types = cols())
colnames(elg_rate) <- c("ensembl_gene_id", "elongation_rate", "expression")
elg_rate <- elg_rate[elg_rate$elongation_rate > 0, ] %>% select(ensembl_gene_id, elongation_rate)

# unpublished data from Alex
elg_rate_alex <- read_csv(file.path(root_dir,
                                    "data/elongation_rate/alex/elongation_rates.csv"),
                          col_types = cols())
elg_rate_alex <- elg_rate_alex %>%
    left_join(tx2gene, by = "ensembl_transcript_id") %>%
    select(ensembl_gene_id, everything(), -ensembl_transcript_id)

# amit_exp <- amit_exp %>% select(ensembl_gene_id, biotype, TPM_PROseq, TPM_RNAseq) %>%
#     left_join(elg_rate, by = "ensembl_gene_id") %>%
#     left_join(elg_rate_alex, by = "ensembl_gene_id") %>%
#     mutate(TPM_PROseq_crt1 = TPM_PROseq * 1000 / elongation_rate,
#            TPM_PROseq_crt2 = TPM_PROseq * 1000 / Rate_60to120min)
#
# amit_exp_crt1 <- amit_exp %>%
#     select(TPM_PROseq, TPM_PROseq_crt1, TPM_RNAseq) %>%
#     na.omit() %>%
#     log2()
#
# amit_exp_crt2 <- amit_exp %>%
#     select(TPM_PROseq, TPM_PROseq_crt2, TPM_RNAseq) %>%
#     na.omit() %>%
#     log2()
#
# amit_exp_crt1 %>% cor()
# amit_exp_crt2 %>% cor()

NGENCODE_exp <- NGENCODE_exp %>% select(ensembl_gene_id, biotype, TPM_PROseq, TPM_RNAseq) %>%
    left_join(elg_rate, by = "ensembl_gene_id") %>%
    left_join(elg_rate_alex, by = "ensembl_gene_id") %>%
    mutate(TPM_PROseq_crt1 = TPM_PROseq * 1000 / elongation_rate,
           TPM_PROseq_crt2 = TPM_PROseq * 1000 / Rate_60to120min)

NGENCODE_exp_crt1 <- NGENCODE_exp %>%
    select(TPM_PROseq, TPM_PROseq_crt1, TPM_RNAseq) %>%
    na.omit() %>%
    log2()

NGENCODE_exp_crt2 <- NGENCODE_exp %>%
    select(TPM_PROseq, TPM_PROseq_crt2, TPM_RNAseq) %>%
    na.omit() %>%
    log2()

make_scatter <- function(df, col1, col2, lab_x, lab_y, corr, n,
                         xlim = c(0, 10), ylim = c(0, 10)) {
    ggplot(df, aes_string(col1, col2)) +
        geom_point(alpha = 0.3, color = "steelblue3") +
        # geom_smooth(method = "lm", se = TRUE, color = "red") +
        labs(x = bquote(paste(log[2], .(lab_x))),
             y = bquote(paste(log[2], .(lab_y)))) +
        coord_fixed(ratio = 0.8, xlim = xlim, ylim = ylim) +
        annotate(geom="text", x = 1, y = 9,
                 label = bquote(paste(italic("ρ"), " = ", .(round(corr, 2)))),
                 size = 5) +
        annotate(geom="text", x = 1, y = 8,
                 label = bquote(paste("n", " = ", .(n))),
                 size = 5) +
        theme_classic(base_size = 15)
}

# make_scatter(amit_exp_crt1, "TPM_PROseq", "TPM_RNAseq", "(PRO-seq TPM)", "(RNA-seq TPM)", 0.55, 2204)
# ggsave(file.path(exp_fig_dir, "K562_Amit_PROseq_vs_RNAseq.png"), width = 5, height = 5)
#
# make_scatter(amit_exp_crt1, "TPM_PROseq_crt1", "TPM_RNAseq", "(PRO-seq TPM corrected by ER)", "(RNA-seq TPM)", 0.48, 2204)
# ggsave(file.path(exp_fig_dir, "K562_Amit_PROseq_with_ER_correction_vs_RNAseq.png"), width = 5, height = 5)
#
# make_scatter(amit_exp_crt2, "TPM_PROseq", "TPM_RNAseq", "(PRO-seq TPM)", "(RNA-seq TPM)", 0.87, 1543,
#              xlim = c(-7, 10), ylim = c(-7, 10))
# ggsave(file.path(exp_fig_dir, "K562_Amit_PROseq_2_vs_RNAseq.png"), width = 5, height = 5)
#
# make_scatter(amit_exp_crt2, "TPM_PROseq_crt2", "TPM_RNAseq", "(PRO-seq TPM corrected by ER)", "(RNA-seq TPM)", 0.66,
#              1543, xlim = c(-7, 10), ylim = c(-7, 10))
# ggsave(file.path(exp_fig_dir, "K562_Amit_PROseq_with_ER_correction_2_vs_RNAseq.png"), width = 5, height = 5)

# corrected PRO-seq signal using published elongation rates
# NGENCODE_exp_crt1 %>% cor()

suppressMessages(
    p <- make_scatter(NGENCODE_exp_crt1, "TPM_PROseq", "TPM_RNAseq", "(PRO-seq TPM)", "(RNA-seq TPM)", 0.57, 1712)
)

suppressMessages(
    ggsave(file.path(exp_fig_dir, "K562_2014NG_PROseq_vs_RNAseq.png"), p, width = 5, height = 5)
)

suppressMessages(
    p <- make_scatter(NGENCODE_exp_crt1, "TPM_PROseq_crt1", "TPM_RNAseq",
                      "(PRO-seq TPM corrected by ER)", "(RNA-seq TPM)", 0.48, 1712)
)

suppressMessages(
    ggsave(file.path(exp_fig_dir, "K562_2014NG_PROseq_with_ER_correction_vs_RNAseq.png"), p, width = 5, height = 5)
)

# corrected using Alex data
# NGENCODE_exp_crt2 %>% cor()
# make_scatter(NGENCODE_exp_crt2, "TPM_PROseq", "TPM_RNAseq", "(PRO-seq TPM)", "(RNA-seq TPM)",
#             0.89, 794, xlim = c(-7, 10), ylim = c(-7, 10))
# ggsave(file.path(exp_fig_dir, "K562_2014NG_PROseq_2_vs_RNAseq.png"), width = 5, height = 5)

# make_scatter(NGENCODE_exp_crt2, "TPM_PROseq_crt2", "TPM_RNAseq",
#             "(PRO-seq TPM corrected by ER)", "(RNA-seq TPM)",
#             0.70, 794, xlim = c(-7, 10), ylim = c(-7, 10))
# ggsave(file.path(exp_fig_dir, "K562_2014NG_PROseq_with_ER_correction_2_vs_RNAseq.png"), width = 5, height = 5)
