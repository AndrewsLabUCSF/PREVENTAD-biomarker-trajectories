
# LOAD PACKAGES -----------------------------------------------------------

install.packages("pacman")

pacman::p_load(tidyr, dplyr, ggplot2, gtsummary, gt, here, tidyverse, data.table)

DATA_PATHS <- list(
  output = list(
    figures = here("results", "figures"),
    tables = here("results", "tables"), 
    reports = here("results", "reports")
  ),
  data = list(
    cache = here("data", "cache"),
    processed = here("data", "processed")
  )
)
