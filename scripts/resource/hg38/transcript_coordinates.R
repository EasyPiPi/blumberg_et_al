library(GenomicFeatures)
library(rtracklayer)

root_dir <- path.expand("~/Desktop/github_repo/blumberg_et_al")

hg38_txdb <- makeTxDbFromBiomart()

hg38_tx <- transcripts(hg38_txdb)
hg38_tx$name <- hg38_tx$tx_name

dir.create(file.path(root_dir, "output/resource/hg38"), recursive = TRUE, showWarnings = FALSE)
export.bed(hg38_tx, file.path(root_dir, "output/resource/hg38/transcript_coordinates.bed"))
