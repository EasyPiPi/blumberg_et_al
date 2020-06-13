library(tidyverse)

root_dir <- path.expand("~/Desktop/github_repo/blumberg_et_al")

elg_rate <- read_csv(file.path(root_dir, "data/elongation_rate/alex/elongation_rates.csv"))

pub_egl_rate <- read_csv(file.path(root_dir, "data/elongation_rate/veloso_et_al/K562_elongation_rate.csv"))

cv <- function(vec) {
    sd(vec, na.rm = TRUE) / mean(vec, na.rm = TRUE)
}

cv(elg_rate$Rate_60to120min)
cv(pub_egl_rate$`K562 Elongation Rate (bp/min)`)

summary(elg_rate$Rate_60to120min)
summary(pub_egl_rate$`K562 Elongation Rate (bp/min)`)

quantile(elg_rate$Rate_60to120min, probs = c(0.05, 0.1, 0.5, 0.9, 0.95))
quantile(pub_egl_rate$`K562 Elongation Rate (bp/min)`, probs = c(0.05, 0.1, 0.5, 0.9, 0.95))

