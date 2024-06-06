# require(devtools)
# use_readme_rmd()
# use_news_md()
# use_vignette("SamplingStrata")
# use_github_links()
# # use_travis()
# use_cran_badge()

# install.packages("tvthemes") # v1.1.0
library(tvthemes)

# devtools::install_version("pkgdown", version = "1.6.1")
# devtools::install_github("hadley/pkgdown")
# setwd("C:/Users/UTENTE/Google Drive/Genetic Algorithm/QGA")
setwd("D:\\Google Drive\\Genetic Algorithm\\QGA")
# 
# library(devtools)
# use_readme_rmd()
# use_vignette("QGA") 
library(pkgdown)
usethis::use_pkgdown()

# build_favicons(pkg = ".")
init_site(pkg = ".")
pkgdown::build_site()



