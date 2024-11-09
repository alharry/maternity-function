# Setup

# Knitr options
knitr::opts_chunk$set(warning=FALSE, message=FALSE, echo = FALSE, dev = "pdf")

# Housekeeping
set.seed(3)

# Load libraries
libs<-c("MASS", "tidyverse", "knitcitations", "bibtex", "lubridate",
        "forcats", "here", "TMB", "broom", "rsample", "cowplot",
        "grid", "gridExtra", "xtable", "rsample")

lapply(libs, library, character.only = TRUE)

# Bibliography
cleanbib()
cite_options(citation_format = "pandoc", check.entries = FALSE)
