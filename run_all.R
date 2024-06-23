## Run all

# Load libraries
library(knitr)
library(rmarkdown)

# Data analysis
source(here::here("code", "simulation.R"))
source(here::here("code", "summarise-results.R"))
source(here::here("code", "empirical.R"))


# Generate files
render("rmds/manuscript.Rmd", output_format = "all", output_file = "main-text.pdf", output_dir = "reports", envir = new.env())
render("rmds/figures.RmD", output_format = "all", output_file = "figs.pdf", output_dir = "reports", envir = new.env())
render("rmds/figures-supp.RmD", output_format = "all", output_file = "fig-supp.pdf", output_dir = "reports", envir = new.env())
render("rmds/tables.RmD", output_format = "all", output_file = "tabs.pdf", output_dir = "reports", envir = new.env())

## Zotero citation key format
# [auth:lower]_[Title:lower:skipwords:select,1,1]_[year]

# Rmarkdown template
# http://labrtorian.com/2019/08/26/rmarkdown-template-that-manages-academic-affiliations/
