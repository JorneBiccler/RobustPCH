setwd("~/Scripts/")
source("./cvPCH.R")

library(survival)
library(coxrobust)
library(MASS)
set.seed(87913)

data <- Aids2[sample(1:nrow(Aids2), 1000), ]
data$time <- data$death - data$diag + 1
data$d01 <- data$status == "D"

coxr(Surv(time, 1 - d01) ~ age + sex , data)

breaks <- c(getEvalTimes(data$time[data$d01 == 1], 
                         nEvalTimes = 9, extremeTails = c(0, 1),
                         includeZero = T, removeLargest = T), Inf)

modMat <- model.matrix(Surv(time,  d01) ~ age + sex, data = data)[, - 1]

evalTimes <- getEvalTimes(data$time[data$ d01 == 1], nEvalTimes = 31,
                          includeZero = F,
                          extremeTails = c(0, 1), removeLargest = F)

MLEsolution <- cvMLE(timeVar = data$time,
                     d01 = data$d01,
                     breaks = breaks,
                     X = modMat,
                     nLambda = 100,
                     logLambdaRatio = 5,
                     nFolds = 10,
                     cvCriterion = "logLik",
                     nCores = 10, repeatedCV = 10)

censoringFit <- survfit(Surv(time, 1 -  d01) ~ 1, data)
cProbEvalTimes <- predictSurvProb(censoringFit, newdata = data, times =  evalTimes)

## pretty inneficient way of doing this but ok ...
cProbTimeVar <- diag(predictSurvProb(censoringFit, newdata = data, 
                                     times =  data$time))[rank(data$time)]
cProbTimeVar[cProbTimeVar == 0] <- min(cProbTimeVar[cProbTimeVar != 0])

invisible(capture.output(BrierSolution <- cvBrier(timeVar = data$time,
                                                  d01 = data$d01, 
                                                  breaks = breaks, 
                                                  X = modMat,
                                                  cProbTimeVar = cProbTimeVar,
                                                  cProbEvalTimes = cProbEvalTimes,
                                                  evalTimes = evalTimes,
                                                  nLambda = 100,
                                                  logLambdaRatio = 5, 
                                                  nFolds = 10,
                                                  nCores = 10,
                                                  cvCriterion = "brier",
                                                  repeatedCV = 10,
                                                  trim = 1)))

saveRDS(BrierSolution, "./Applications/BrierSolutiontimeAids.rds")
saveRDS(MLEsolution, "./Applications/MLESolutiontimeAids.rds")

coefNames <- c("Baseline","Age","Gender: male")

nCov <- ncol(modMat)
breaks[length(breaks)] <- 1000
transMat <- paramTransMat(nIntervals = length(breaks) - 1, nParameters = nCov + 1)

nEstInterval <- length(breaks) - 1
nCoef <- length(coefNames)

MLE <- getParam(MLEsolution, nIntervals = length(breaks) - 1, nParameters = length(coefNames))
Brier <- getParam(BrierSolution, nIntervals = length(breaks) - 1, nParameters = length(coefNames))

## 
breaks[length(breaks)] <- max(data$time) + 1
pchResidsBrier <- rep(NA, nrow(data))
pchResidsMLE <- rep(NA, nrow(data))
for(i in 1:nrow(data)){
  predictedIBrier <- predictPCH(parameters = Brier, 
                                X = modMat, 
                                breaks = breaks, 
                                evalTimes = c(0, data$time[i]))
  predictedIMLE <- predictPCH(parameters = MLE, 
                              X = modMat, 
                              breaks = breaks, 
                              evalTimes = c(0, data$time[i]))
  pchResidsBrier[i] <- qnorm(predictedIBrier[i, 2] /
                               (2 - data$d01[i]))
  pchResidsMLE[i] <- qnorm(predictedIMLE[i, 2] /
                             (2 - data$d01[i]))
}

hist(pchResidsMLE)
hist(pchResidsBrier)

data[which(abs(pchResidsBrier) > 3), ]
data[which(abs(pchResidsMLE) > 3), ]

cleanData <- data[!(abs(pchResidsBrier) > 3), ]
modMatClean <- model.matrix(Surv(time,  d01) ~ age + sex, data = cleanData)[, - 1]

MLEsolutionClean <- cvMLE(timeVar = cleanData$time,
                     d01 = cleanData$d01,
                     breaks = breaks,
                     X = modMatClean,
                     nLambda = 100,
                     logLambdaRatio = 5,
                     nFolds = 5,
                     cvCriterion = "logLik",
                     nCores = 10, repeatedCV = 10)

MLEClean <- getParam(MLEsolutionClean, nIntervals = length(breaks) - 1, nParameters = length(coefNames))

breaks[length(breaks)] <- 1000
PlotFrame <- data.frame(start = rep(breaks[- length(breaks)], each = 3),
                        stop = rep(breaks[-1], each = 3),
                        estimate = c(Brier, MLE, MLEClean),
                        Method = factor(rep(c("BSLE", "MLE", "MLE, outliers\nremoved"), each = length(Brier)),
                                        levels = c("MLE, outliers\nremoved", "MLE", "BSLE")),
                        Variable = factor(rep(coefNames, length(Brier) / length(coefNames)),
                                          levels = c(coefNames)))

coefPlot <- ggplot(PlotFrame, aes(x = start/ 365.24, y = estimate, col = Method)) + 
  geom_segment(aes(xend = stop/ 365.24, yend = estimate), size = 1) +
  facet_wrap( ~ Variable, scales = "free_y", ncol =  2) +
  xlab("Years since diagnosis") +
  ylab("Estimated coefficient") +
  theme_bw() +
  scale_color_brewer(palette = "Set2", guide = guide_legend(reverse=TRUE))

ggsave("./Figures/coefPlotAIDS.eps", coefPlot, width = 6, height = 4)
