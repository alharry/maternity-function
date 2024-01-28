## Run all

# Load libraries
library(knitr)
library(rmarkdown)

# Data analysis
source(here::here("code", "simulation.R")) # 48 hour run time on my computer
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

# TODO 
# - Bubley spiny dogfish - inconsistent definition (any with pups or large developing follicles)
# - Cotton Squalus mitsikuri uses pregnancy ogive as maternity ogive
# - Fujimani blue shark
# - Rochowski paper
# - Corro-Espinosa rhizoprionodon
# - mustelus henlei (file:///Users/avh/Downloads/2805-Article%20Text-420420528-1-10-20180628.pdf)
# - Squalus blainvillei (https://www-publish-csiro-au.libproxy.murdoch.edu.au/mf/pdf/MF19372)

# Notes
# - https://www.jstor.org/stable/2531207
# - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5978598/
# - https://www.researchgate.net/profile/Gunter_Maris/publication/237079179_On_Interpreting_the_Model_Parameters_for_the_Three_Parameter_Logistic_Model/links/5661554508aebae678aa7970.pdf
