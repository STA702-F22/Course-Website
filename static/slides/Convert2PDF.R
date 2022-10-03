#remotes::install_github("jhelvy/renderthis")
#remotes::install_github('rstudio/chromote')
#install.packages('pdftools')

library(renderthis)

basedir = "/Users/clyde/Dropbox/sta702/Course-Website"

convert = function(file) {
  input=paste0(file, ".Rmd")
  handout = paste0(file, "-handout.pdf")
  slides = paste0(file, ".pdf")
  to_pdf(input, to = handout,complex_slides=TRUE)
  to_pdf(input, to=slides, complex_slides=TRUE, partial_slides = TRUE)
}

setwd("/Users/clyde/Dropbox/sta702/Course-Website/static/slides")

# Lectures

#convert("00-course-overview")
#convert("01-basics-of-bayes")
#convert("02-loss-functions")
#convert("03-normal-predictive-distributions")
#convert("04-predictive-checks")
#convert("05-hierarchical-models")
#convert("06-metropolis")
#convert("07-adaptive-metropolis")
#convert("08-gibbs")
#convert("09-gibbs-data-augmentation")
convert("10-hypothesis-testing")
#convert("11-hypothesis-testing-cont")
#convert("12-multiple-testing")
#convert("13-Bayes-multiple-testing")
#convert("14-model-selection")
#convert("15-BMA")
#convert("18-random-effects")
#convert("19-mixed-effects")
#convert("20-missing-data")
#convert("21-HMC")
setwd(basedir)

