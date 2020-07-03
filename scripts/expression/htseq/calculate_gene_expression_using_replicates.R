library(tidyverse)
library(corrplot)
library(skimr)

root_dir <- "~/Desktop/github_repo/blumberg_et_al"

# proseq_exp_dir <- file.path(root_dir, "data/htseq/K562_PROseq_Amit")
# rnaseq_exp_dir <- file.path(root_dir, "data/htseq/K562_RNAseq_Amit")

proseq_exp_dir <- file.path(root_dir, "data/htseq/NatGen2014")
rnaseq_exp_dir <- file.path(root_dir, "data/htseq/K562_RNAseq_ENCODE")

id_dir <- file.path(root_dir, "data/id_mapper")
length_dir <- file.path(root_dir, "data/length")
exp_fig_dir <- file.path(root_dir, "output/expression/figure")
dir.create(exp_fig_dir, showWarnings = FALSE, recursive = TRUE)
exp_df_dir <- file.path(root_dir, "output/expression/table")
dir.create(exp_df_dir, showWarnings = FALSE, recursive = TRUE)

# read in id mapper and gene types
tx2gene <- read_csv(file.path(id_dir, "tx2gene.csv"))
genetype <- read_csv(file.path(id_dir, "genetype.csv"))

# read in read count dfs
proseq_exp_df_names <- list.files(proseq_exp_dir)
rnaseq_exp_df_names <- list.files(rnaseq_exp_dir)

# proseq_exp_df_names <- str_subset(proseq_exp_df_names, "human_longest")
# rnaseq_exp_df_names <- str_subset(rnaseq_exp_df_names, "human_longest")

proseq_exp_df_names <- str_subset(proseq_exp_df_names, "pc_linc_as")
rnaseq_exp_df_names <- str_subset(rnaseq_exp_df_names, "pc_linc_as")

exp_df_names <- c(proseq_exp_df_names, rnaseq_exp_df_names)

read_count_dfs <- str_split(exp_df_names, "\\.", simplify = T)
colnames(read_count_dfs) <- c("RNA_type", "gtf", "region", "txt")

read_count_dfs <- as_tibble(read_count_dfs) %>%
    add_column(file_names = c(file.path(proseq_exp_dir, proseq_exp_df_names),
                              file.path(rnaseq_exp_dir, rnaseq_exp_df_names))) %>%
    select(-txt) %>%
    mutate(read_count_df = map(file_names, read_delim, delim = "\t", col_names = c("ensembl_transcript_id", "read_count"))) %>%
    mutate(assay = ifelse(str_detect(region, "transcript"), "PROseq", "RNAseq_ex"))

# filtering rows which don't belong to genes
no_count_rownames <- c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")

read_count_dfs <- read_count_dfs %>%
    # move the no count rows to another col
    mutate(no_count_df = map(read_count_df, . %>% filter(ensembl_transcript_id %in% no_count_rownames))) %>%
    # remove the no count rows
    mutate(read_count_df = map(read_count_df, . %>% filter(!ensembl_transcript_id %in% no_count_rownames)))

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
    inner_join(read_count_df, len_df, by = "ensembl_transcript_id") %>%
        # 3th col is the gene length or tx length
        mutate(RPK = read_count / (.[[3]] / 1000)) %>%
        mutate(TPM = RPK / (sum(RPK) / 1e6)) %>%
        select(ensembl_transcript_id, TPM)
}

TPM_dfs <- left_join(read_count_dfs, len_dfs, by = "assay") %>%
    mutate(TPM_df = map2(read_count_df, len_df, get_TPM)) %>%
    select(assay, RNA_type, region, gtf, TPM_df)

# join replicate TPM dfs into a single df
join_TPM_df <- function(df, exp_filter = 1) {
    df <- df %>% mutate(RNA_type = str_c(RNA_type, "_", assay)) %>%
        select(RNA_type, TPM_df) %>%
        unnest(cols = c(TPM_df)) %>%
        spread(key = RNA_type, value = TPM) %>%
        filter_if(is.numeric, all_vars(. > exp_filter)) %>%
        mutate_if(is.numeric, log2)
    return(df)
}

# filter > 1 TPM
log2TPM <- TPM_dfs %>% group_split(assay, gtf) %>% map(join_TPM_df)
names(log2TPM) <- TPM_dfs %>%
    select(assay, gtf) %>%
    mutate(file_name = str_c(assay, "_", gtf)) %>%
    select(file_name) %>% unique() %>% pull()

# output expression df for SEM
# log2TPM[[1]] %>% inner_join(log2TPM[[2]], by = "ensembl_transcript_id") %>%
#     write_csv(file.path(exp_df_dir, "expression_K562_Amit_total_RNA_log2TPM_1_replicates.csv"))
log2TPM[[1]] %>% inner_join(log2TPM[[2]], by = "ensembl_transcript_id") %>%
    write_csv(file.path(exp_df_dir, "expression_K562_2014NG_polyA_RNA_log2TPM_1_replicates.csv"))

# filter > 0 TPM
log2TPM <- TPM_dfs %>%
    group_split(assay, gtf) %>% map(join_TPM_df, exp_filter = 0)
names(log2TPM) <- TPM_dfs %>%
    select(assay, gtf) %>%
    mutate(file_name = str_c(assay, "_", gtf)) %>%
    select(file_name) %>% unique() %>% pull()

# output expression df for SEM
# log2TPM[[1]] %>% inner_join(log2TPM[[2]], by = "ensembl_transcript_id") %>%
#     write_csv(file.path(exp_df_dir, "expression_K562_Amit_total_RNA_log2TPM_0_replicates.csv"))

log2TPM[[1]] %>% inner_join(log2TPM[[2]], by = "ensembl_transcript_id") %>%
    write_csv(file.path(exp_df_dir, "expression_K562_2014NG_polyA_RNA_log2TPM_0_replicates.csv"))

# compute expression correlation for replicates
compute_corr <- function(df) {
    df <- df[2:ncol(df)] %>% cor()
    colnames(df) <- str_c(rep("Rep", ncol(df)), seq(1:ncol(df)), sep = " ")
    rownames(df) <- str_c(rep("Rep", ncol(df)), seq(1:ncol(df)), sep = " ")
    return(df)
}

corr_mx <- log2TPM %>% map(compute_corr)

# define handy functions
# corrplot
corrplot_partial <-
    partial(
        corrplot,
        method = "square",
        addrect = 2,
        tl.srt = 45,
        cl.lim = c(0, 1),
        addCoef.col = "white",
        tl.col = "black",
        diag = FALSE,
        tl.cex = 1.2, cl.cex = 1.5,
        number.digits = 4
    )
# save output figure
save_fig <- function(file_path, fig_fun) {
    png(file_path, width = 1000, height = 1000)
    fig_fun
    dev.off()
}

# output figures
for (i in 1:length(corr_mx)) {
    save_fig(file.path(
        exp_fig_dir,
        str_c(names(corr_mx[i]), "_rep_expression_correlation.png")
    ),
    corrplot_partial(corr_mx[[i]]))
}

