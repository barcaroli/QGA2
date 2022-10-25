prepare_data <-function() {
  # Prepare data for fitness evaluation
  data(iris)
  iris$id <- c(1:nrow(iris))
  iris$dom <- 1
  library(SamplingStrata)
  frame <- buildFrameDF(
    df = iris,
    id = "id",
    domainvalue = "dom",
    X = c("id"),
    Y = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
  )
  cv <- as.data.frame(list(
    DOM = "DOM1",
    CV1 = 0.03,
    CV2 = 0.03,
    CV3 = 0.03,
    CV4 = 0.03,
    domainvalue = 1
  ))
  
  nvar = length(grep("Y",colnames(frame)))
  source("aggrStrata.R")
  source("aggrStrata2.R")
  atomic_strata = aggrStrata2(dataset=frame,
                              model=NULL,
                              vett=c(1:150),
                              dominio=1)
}