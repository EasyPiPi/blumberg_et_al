library(tidyverse)

root_dir <- "~/Desktop/github_repo/blumberg_et_al"

exp_dir <- file.path(root_dir, "output/expression/table")
id_dir <- file.path(root_dir, "data/id_mapper")
tx_features_dir <- file.path(root_dir, "data/features")
histone_dir <- file.path(root_dir, "data/histone")
sem_data_dir <- file.path(root_dir, "output/sem/data")
dir.create(sem_data_dir, showWarnings = FALSE, recursive = TRUE)

exp_df <- read_csv(file.path(exp_dir, "expression_K562_Amit_total_RNA_log2TPM_1_replicates.csv"))
feature_df <- read_csv(file.path(tx_features_dir, "human_gene_features.csv"))

gene_biotype_df <- read_csv(file.path(id_dir, "genetype.csv"))
tx2gene <- read_csv(file.path(id_dir, "tx2gene.csv"))

elg_df <- read_csv(file.path(root_dir, "data/elongation_rate/veloso_et_al/K562_elongation_rate.csv"))
elg_df <- elg_df %>% rename("ensembl_gene_id" = `Ensembl Gene ID`,
           "elongation_rate" = `K562 Elongation Rate (bp/min)`) %>%
    left_join(tx2gene, by = "ensembl_gene_id") %>%
    select("ensembl_transcript_id", "elongation_rate")

exp_feature_df <- exp_df %>%
    left_join(elg_df, by = "ensembl_transcript_id") %>%
    mutate_at(vars(contains("PROseq")), ~ . * elongation_rate / 1000) %>%
    na.omit() %>%
    select(-elongation_rate) %>%
    left_join(feature_df, by = c(ensembl_transcript_id = "tx_name")) %>%
    left_join(gene_biotype_df, by = c(gene_id = "ensembl_gene_id")) %>%
    mutate(spliced = ifelse(nexon > 1, "spliced", "non-spliced"))

# grep column names
pc_cols <- grep("5utr|cds|3utr|intron|exonJunDen", colnames(exp_feature_df), value = TRUE)
spliced_cols <- grep("_exon|intron|exonJunDen", colnames(exp_feature_df), value = TRUE)
nonspliced_cols <- grep("gc_gene|len_gene", colnames(exp_feature_df), value = TRUE)
exp_cols <- grep("ensembl_transcript_id|PROseq|RNAseq", colnames(exp_feature_df), value = TRUE)

exp_feature_gp <- exp_feature_df %>%
    group_nest(biotype, spliced)

# output dfs
exp_feature_df %>% filter(biotype == "protein_coding", spliced == "spliced") %>%
    select(c(exp_cols, pc_cols)) %>% na.omit() %>% mutate_at(pc_cols, scale) %>%
    write_csv(file.path(sem_data_dir, "protein_coding_spliced_full_feature_with_elongation_correction.csv"))

exp_feature_df %>% filter(biotype == "protein_coding", spliced == "spliced") %>%
    select(c(exp_cols, spliced_cols)) %>% na.omit() %>% mutate_at(spliced_cols, scale) %>%
    write_csv(file.path(sem_data_dir, "protein_coding_spliced_with_elongation_correction.csv"))

exp_feature_df %>% filter(biotype == "lincRNA", spliced == "spliced") %>%
    select(c(exp_cols, spliced_cols)) %>% na.omit() %>% mutate_at(spliced_cols, scale) %>%
    write_csv(file.path(sem_data_dir, "lincRNA_spliced_with_elongation_correction.csv"))

# exp_feature_df %>% filter(biotype == "protein_coding", spliced == "non-spliced") %>%
#     select(c(exp_cols, nonspliced_cols)) %>% na.omit() %>% mutate_at(nonspliced_cols, scale) %>%
#     write_csv(file.path(sem_data_dir, "protein_coding_nonspliced_with_elongation_correction.csv"))

# exp_feature_df %>% filter(biotype == "lincRNA", spliced == "non-spliced") %>%
#     select(c(exp_cols, nonspliced_cols)) %>% na.omit() %>% mutate_at(nonspliced_cols, scale) %>%
#     write_csv(file.path(sem_data_dir, "lincRNA_nonspliced_with_elongation_correction.csv"))

# no correction
exp_feature_no_corr_df <- exp_df %>%
    left_join(elg_df, by = "ensembl_transcript_id") %>%
    na.omit() %>%
    select(-elongation_rate) %>%
    left_join(feature_df, by = c(ensembl_transcript_id = "tx_name")) %>%
    left_join(gene_biotype_df, by = c(gene_id = "ensembl_gene_id")) %>%
    mutate(spliced = ifelse(nexon > 1, "spliced", "non-spliced"))

exp_feature_no_corr_df %>% filter(biotype == "protein_coding", spliced == "spliced") %>%
    select(c(exp_cols, pc_cols)) %>% na.omit() %>% mutate_at(pc_cols, scale) %>%
    write_csv(file.path(sem_data_dir, "protein_coding_spliced_full_feature_without_elongation_correction.csv"))

# SEM for histones
histone_dfs <-
    tibble(file_name = str_remove(list.files(histone_dir), ".txt"),
           histone_df = map(list.files(histone_dir, full.names = TRUE),
                            read_delim, delim = '\t')) %>%
    mutate(output_df =
               map(histone_df, ~ .x %>%
                       left_join(exp_feature_df %>% select("ensembl_transcript_id", contains("PROseq"), contains("RNAseq")),
                                 by = c("V1" = "ensembl_transcript_id")) %>%
                       na.omit()))

walk2(histone_dfs$output_df, file.path(sem_data_dir, paste0(histone_dfs$file_name, "_with_elongation_correction.csv")), write_csv)

histone_wo_correction_dfs <-
    tibble(file_name = str_remove(list.files(histone_dir), ".txt"),
           histone_df = map(list.files(histone_dir, full.names = TRUE),
                            read_delim, delim = '\t')) %>%
    mutate(output_df =
               map(histone_df, ~ .x %>%
                       left_join(exp_feature_no_corr_df %>% select("ensembl_transcript_id", contains("PROseq"), contains("RNAseq")),
                                 by = c("V1" = "ensembl_transcript_id")) %>%
                       na.omit()))

walk2(histone_wo_correction_dfs$output_df, file.path(sem_data_dir, paste0(histone_wo_correction_dfs$file_name, "_without_elongation_correction.csv")), write_csv)
