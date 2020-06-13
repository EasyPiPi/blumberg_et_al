library(GenomicFeatures)

root_dir <- path.expand("~/Desktop/github_repo/blumberg_et_al")

hg19_txdb <- makeTxDbFromBiomart(host = "http://grch37.ensembl.org")
hg19_3UTR <- unlist(threeUTRsByTranscript(hg19_txdb, use.names = TRUE))
hg19_3UTR$ensembl_transcript_id <- names(hg19_3UTR)
names(hg19_3UTR) <- NULL

write.csv(as.data.frame(hg19_3UTR),
          file.path(root_dir, "output/resource/hg19/three_prime_utr_grng.csv"),
          row.names = FALSE)
