#remotes::install_github("jhelvy/xaringanBuilder")
#remotes::install_github('rstudio/chromote')
#install.packages('pdftools')

library(xaringanBuilder)

basedir = "/Users/clyde/Dropbox/sta601/Course-Website"

convert = function(file) {
  input=paste0(file, ".Rmd")
  handout = paste0(file, "-handout.pdf")
  slides = paste0(file, ".pdf")
  build_pdf(input, output = handout,complex_slides=TRUE)
  build_pdf(input, output=slides, complex_slides=TRUE, partial_slides = TRUE)
}

setwd("/Users/clyde/Dropbox/sta601/Course-Website/static/slides")

# Lectures

#convert("00-course-overview")
#convert("01-basics-of-bayes")
#convert("02-loss-functions")
#convert("03-normal-predictive-distributions")
#convert("04-predictive-checks")
#convert("05-hypothesis-testing")
#convert("06-hypothesis-testing-cont")
#convert("07-hierarchical-models")
#convert("08-metropolis-diagnostics")
#convert("09-adaptive-metropolis")
#convert("10-gibbs")
#convert("11-gibbs-data-augmentation")
#convert("12-multiple-testing")
#convert("13-Bayes-multiple-testing")
#convert("14-model-selection")
#convert("15-BMA")
#convert("18-random-effects")
#convert("19-mixed-effects")
#convert("20-missing-data")
convert("21-HMC")
setwd(basedir)
