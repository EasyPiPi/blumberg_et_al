library(tidyverse)

root_dir <- path.expand("~/Desktop/github_repo/blumberg_et_al")

exp_dir <- file.path(root_dir, "output/expression/table")
id_dir <- file.path(root_dir, "data/id_mapper")
tx_features_dir <- file.path(root_dir, "data/features")
histone_dir <- file.path(root_dir, "data/histone")
sem_data_dir <- file.path(root_dir, "output/sem/data")
dir.create(sem_data_dir, showWarnings = FALSE, recursive = TRUE)

# exp_df <- read_csv(file.path(exp_dir, "expression_K562_Amit_total_RNA_log2TPM_1_replicates.csv"))
exp_df <- read_csv(file.path(exp_dir, "expression_K562_2014NG_polyA_RNA_log2TPM_1_replicates.csv"))
feature_df <- read_csv(file.path(tx_features_dir, "human_gene_features.csv"))

gene_biotype_df <- read_csv(file.path(id_dir, "genetype.csv"))
tx2gene <- read_csv(file.path(id_dir, "tx2gene.csv"))

exp_feature_df <- exp_df %>%
    #left_join(tx2gene, by = "ensembl_transcript_id") %>%
    left_join(feature_df, by = c(ensembl_transcript_id = "tx_name")) %>%
    # left_join(feature_gp_df, by = c(ensembl_gene_id = "gene_id")) %>%
    left_join(gene_biotype_df, by = c(gene_id = "ensembl_gene_id")) %>%
    # left_join(gene_biotype_df, by = "ensembl_gene_id") %>%
    # filter(biotype %in% c("protein_coding", "lincRNA")) %>%
    # filter(ensembl_transcript_id %in% prop_50_names) %>%
    mutate(spliced = ifelse(nexon > 1, "spliced", "non-spliced"))

# greo column names
pc_cols <- grep("5utr|cds|3utr|intron|exonJunDen", colnames(exp_feature_df), value = TRUE)
spliced_cols <- grep("_exon|intron|exonJunDen", colnames(exp_feature_df), value = TRUE)
nonspliced_cols <- grep("gc_gene|len_gene", colnames(exp_feature_df), value = TRUE)
exp_cols <- grep("ensembl_transcript_id|PROseq|RNAseq", colnames(exp_feature_df), value = TRUE)

exp_feature_gp <- exp_feature_df %>%
    group_nest(biotype, spliced)

# output dfs
exp_feature_df %>% filter(biotype == "protein_coding", spliced == "spliced") %>%
    select(c(exp_cols, pc_cols)) %>% na.omit() %>% mutate_at(pc_cols, scale) %>%
    write_csv(file.path(sem_data_dir, "protein_coding_spliced_full_feature.csv"))

exp_feature_df %>% filter(biotype == "protein_coding", spliced == "spliced") %>%
    select(c(exp_cols, spliced_cols)) %>% na.omit() %>% mutate_at(spliced_cols, scale) %>%
    write_csv(file.path(sem_data_dir, "protein_coding_spliced.csv"))

exp_feature_df %>% filter(biotype == "lincRNA", spliced == "spliced") %>%
    select(c(exp_cols, spliced_cols)) %>% na.omit() %>% mutate_at(spliced_cols, scale) %>%
    write_csv(file.path(sem_data_dir, "lincRNA_spliced.csv"))

exp_feature_df %>% filter(biotype == "protein_coding", spliced == "non-spliced") %>%
    select(c(exp_cols, nonspliced_cols)) %>% na.omit() %>% mutate_at(nonspliced_cols, scale) %>%
    write_csv(file.path(sem_data_dir, "protein_coding_nonspliced.csv"))

exp_feature_df %>% filter(biotype == "lincRNA", spliced == "non-spliced") %>%
    select(c(exp_cols, nonspliced_cols)) %>% na.omit() %>% mutate_at(nonspliced_cols, scale) %>%
    write_csv(file.path(sem_data_dir, "lincRNA_nonspliced.csv"))

# try not using replicates
# exp_feature_df %>% filter(biotype == "protein_coding", spliced == "spliced") %>%
#     select(c(exp_cols, pc_cols)) %>% na.omit() %>% mutate_at(pc_cols, scale) %>%
#     select(-c("K562_AB_2_PROseq", "total_RNA_2_RNAseq_ex",
#               "total_RNA_3_RNAseq_ex", "total_RNA_4_RNAseq_ex")) %>%
#     write_csv(file.path(sem_data_dir, "protein_coding_spliced_full_feature_no_replicates.csv"))

# keep a same set of genes which can be analyzed MLR with published studies
mlr_gene_set <-
    read_csv(file.path(root_dir, "output/mlr/MLR_hl_and_features.csv")) %>%
    select(X1) %>% pull()

exp_feature_df %>% filter(gene_id %in% mlr_gene_set) %>%
    filter(biotype == "protein_coding", spliced == "spliced") %>%
    select(c(exp_cols, pc_cols)) %>% na.omit() %>% mutate_at(pc_cols, scale) %>%
    write_csv(file.path(sem_data_dir, "protein_coding_spliced_full_feature_for_same_set_gene_in_MLR.csv"))

# SEM for histones
histone_dfs <-
    tibble(file_name = str_remove(list.files(histone_dir), ".txt"),
           histone_df = map(list.files(histone_dir, full.names = TRUE),
                            read_delim, delim = '\t')) %>%
    mutate(output_df =
               map(histone_df, ~ .x %>%
                       left_join(exp_df,by = "ensembl_transcript_id") %>%
                       na.omit()))

walk2(histone_dfs$output_df, file.path(sem_data_dir, paste0(histone_dfs$file_name, ".csv")), write_csv)
