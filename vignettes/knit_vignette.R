library(knitr)
# remove file path vignettes
setwd("/Users/thchoi/Library/Mobile Documents/com~apple~CloudDocs/Drive/Packages/rankIC/")
replace = readLines("vignettes/rankIC.Rmd")
replace = gsub("<img src=\"vignettes/", "<img src=\"", replace)
fileConn = file("vignettes/rankIC.Rmd")
writeLines(replace, fileConn)
close(fileConn)

