chkres <- function(model) {
  require(RVAideMemoire)
  p1 <- plotresid(model)
  sresid <- resid(model, type = "deviance")
  hist(sresid)
  fitted.glmm <- fitted(model, level=1)        # Extract the fitted (predicted) values
  plot(sresid ~ fitted.glmm)                   # Check for homoscedasticity
  plot(model)                                  # gives a heteroscedasticity plot #fans out, not good
  require(arm)
  binnedplot(fitted(model),resid(model))
}