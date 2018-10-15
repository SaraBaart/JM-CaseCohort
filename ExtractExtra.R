ExtractExtra <- function(survObject, lmeObject, timeVar){
  Time <- survObject$y[, 1]
  TermsX <- lmeObject$terms
  formYz <- formula(lmeObject$modelStruct$reStruct[[1]]) 
  data <- lmeObject$data[all.vars(TermsX)]
  
  mfX <- model.frame(TermsX, data = data)
  mfZ <- model.frame(terms(formYz), data = data)
  TermsZ <- attr(mfZ, "terms")
  TermsT <- survObject$terms
  
  Extra1 <- list(Time = Time, timeVar = timeVar, TermsX = TermsX, TermsZ = TermsZ, TermsT = TermsT)
  return(Extra1)
}
