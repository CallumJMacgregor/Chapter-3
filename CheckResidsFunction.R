chkres <- function(model) {
  require(RVAideMemoire)
  sresid <- resid(model, type = "deviance")
  hist(sresid)
  fitted.glmm <- fitted(model, level=1)        # Extract the fitted (predicted) values
  plot(sresid ~ fitted.glmm)                   # Check for homoscedasticity
  plot(model)                                  # gives a heteroscedasticity plot #fans out, not good
  require(arm)
  binnedplot(fitted(model),resid(model))
}




chkres.PQL <- function(model) {
  require(RVAideMemoire)
  plot(model27a)
  sresid <- resid(model27a, type = "pearson")
  hist(sresid)
  fitted.glmm <- fitted(model27a, level=1)        # Extract the fitted (predicted) values
  plot(sresid ~ fitted.glmm)                   # Check for homoscedasticity
  
}