getInitialEstimates <- function(timeVar, d01, X, robQuantile = NULL){
  tempFrame <- data.frame(timeVar, d01)
  marginalFit <- survfit(Surv(timeVar, d01) ~ 1, data = tempFrame)
  
  quantileFun <- stepfun(1 - marginalFit$surv[!duplicated(marginalFit$surv)],
                         c(0, marginalFit$time[!duplicated(marginalFit$surv)]))
  
  if(is.null(robQuantile)){
    robQuantile <- (1 - min(marginalFit$surv))/2
  }
  intercept <- log(- log(1 - robQuantile) / quantileFun(robQuantile))
  coef <- coxrobust::coxr(Surv(timeVar, d01) ~ ., data.frame(timeVar, d01, X), trunc = 0.80)$coef
  c(intercept, coef)
}