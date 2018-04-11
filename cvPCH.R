library(survival)
source("./penalizedPCH.R")
source("./initialEstimates.R")

lambdaMaxBrier <- function(timeVar, d01, X, breaks, init, evalTimes, 
                           cProbEvalTimes, cProbTimeVar, penalType = "totalVar"){
  ## calculate the smalles lambda values for which the model
  ## has all the differences equal to zero, this is a good
  ## value from which to start a solution path.
  
  nParam <- (ncol(X) + 1)
  
  # create some objects needed in the gradient
  lastIndex <- firstIndexIncInterval(timeVar, breaks)
  nIntervals <- length(breaks) - 1 
  exY <- extendedY(evalTimes = evalTimes, 
                   timeVar = timeVar)
  IPWWeights <- exY^(1 - d01)
  IPWWeights <- IPWWeights / pmax(as.numeric(cProbEvalTimes), cProbTimeVar)
  timeAtRiskMat <- timeAtRisk(evalTimes, breaks, nrowRep = length(timeVar), 1)
  
  transformMat <- paramTransMat(nIntervals,
                                nParam)
  X2 <- cbind(1, X)
  largeX <- X2[rep(1:length(timeVar), length(evalTimes)), ]
  
  ## get the gradient of the brier score at the initial solution
  brierGrad <- brierGradient(IPWWeights = IPWWeights, 
                             timeAtRiskMat = timeAtRiskMat, 
                             exY = exY, 
                             evalTimes = evalTimes, 
                             transformMat = transformMat, 
                             largeX = largeX,
                             X = X, 
                             X2 = X2,
                             parameters = init, 
                             breaks =  breaks)
  
  ## get the lambda value
  if(penalType == "totalVar"){
    lambdaMax <- 0
    for(i in 1:nParam){
      paramIndices <- i + (1:(nIntervals - 1)) * nParam
      gradientNorm <- sqrt(sum(brierGrad[paramIndices]^2))
      if(gradientNorm > lambdaMax){
        lambdaMax <- gradientNorm
      }
    }
  } else if (penalType == "fusedLike"){
    lambdaMax <- max(abs(brierGrad[- (1:nParam)]))
  }
  lambdaMax
}

library(pec)
library(robustHD)
library(doParallel)
library(doRNG)

cvBrier <- function(timeVar, d01, X, breaks, evalTimes,
                    nLambda = 20, logLambdaRatio = 5, nFolds = 10,
                    cProbEvalTimes = NULL, cProbTimeVar = NULL,
                    cvCriterion = "brier", penalType = "totalVar", 
                    nCores, repeatedCV = 1, trim = 0.8){
  ## funtion that selects the best lambda value by using a cross validated
  ## Brier score as evaluation criterion.
  
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  ## transform the robust estimates into the parameterization used in the
  ## rest of the functions
  X <- robStandardize(X, fallback = T)
  
  # tempDat <- data.frame(timeVar, X)
  # aftMod <- TML.censored(log(timeVar) ~ ., delta = d01, data = tempDat)
  # initSol <- as.numeric(- aftMod$th1 / aftMod$v1)
  
  initSol <- getInitialEstimates(timeVar, d01, X)
  
  init <- penalizedBrier(timeVar, d01, X, evalTimes, cProbEvalTimes, cProbTimeVar, 
                         breaks = c(0, Inf), init = initSol, lambda = 0, penalType = penalType)
  
  init <- c(init, rep(0, (length(breaks) - 2) * (ncol(X) + 1)))
  ## get the maximum lambda value of interest and create a sequence of lambda values
  logMaxLambda <- - log(lambdaMaxBrier(timeVar = timeVar, d01 = d01, 
                                       breaks = breaks, X = X,
                                       evalTimes = evalTimes,
                                       cProbEvalTimes = cProbEvalTimes,
                                       cProbTimeVar = cProbTimeVar,
                                       penalType = penalType,
                                       init = init))
  
  lambdaSeq <- exp(- seq(logMaxLambda * (1 - 1 / nLambda), 
                         logMaxLambda * logLambdaRatio, 
                         length.out = nLambda))
  
  # sample CV indices
  randSample <- sample(nrow(X))
  CVIndicesList <- split(randSample, ceiling(seq_along(randSample)/(length(randSample)/nFolds)))
  
  transformMat <- paramTransMat(nIntervals = length(breaks) - 1,
                                nParameters = ncol(X) + 1)
  
  # create a vector in which to save the brier scores
  cvScores <- rep(0, nLambda)
  
  if(evalTimes[1] == 0){
    evalTimes <- evalTimes[- 1]
    cProbEvalTimes <- cProbEvalTimes[, - 1]
  }
  
  # perform the CV steps
  for(l in 1:repeatedCV){
    cvResults <- foreach(i = 1:nFolds) %dorng% {
      
      source("./penalizedPCH.R")
      
      trTimeVar <- timeVar[- CVIndicesList[[i]]]
      testTimeVar <- timeVar[CVIndicesList[[i]]]
      
      trd01 <- d01[- CVIndicesList[[i]]]
      testd01 <- d01[CVIndicesList[[i]]]
      
      trX <- X[- CVIndicesList[[i]], ]
      testX <- X[CVIndicesList[[i]], ]
      
      foldPath <- penalizedPathBrier(timeVar = trTimeVar,
                                     d01 = trd01,
                                     init = init,
                                     X = trX, 
                                     breaks = breaks,
                                     evalTimes = evalTimes,
                                     cProbEvalTimes = cProbEvalTimes[- CVIndicesList[[i]], ],
                                     cProbTimeVar = cProbTimeVar[- CVIndicesList[[i]]],
                                     lambdaSeq = lambdaSeq,
                                     penalType = penalType)
      
      predicted <- list()
      cvScores <- rep(0, nLambda)
      
      for(j in 1:nLambda){
        if(cvCriterion == "brier"){
          predicted[[j]] <- predictPCH(transformMat %*% foldPath[j, ],
                                       breaks, testX, evalTimes)
        } else if(cvCriterion == "logLik") {
          cvScores[j] <- - sum(getLogLik(timeVar = testTimeVar, 
                                         d01 = testd01,  
                                         breaks = breaks, 
                                         X = testX,
                                         parameters = transformMat %*% foldPath[j, ]))
        }
      }
      if(cvCriterion == "brier"){
        predicted
      } else if (cvCriterion == "logLik"){
        cvScores
      }
    }
    
    if(cvCriterion == "brier"){
      predicted <- lapply(1:nLambda, function(x) matrix(0, nrow = nrow(X), ncol = length(evalTimes)))
      for(j in 1:nLambda){
        for(i in 1:nFolds){
          predicted[[j]][CVIndicesList[[i]], ] <- cvResults[[i]][[j]]
        }
        cvScores[j] <- cvScores[j] + getTrimmedBrier(cProbEvalTimes, cProbTimeVar, timeVar, d01,
                                                    evalTimes, predicted[[j]], trim)
      }
      
      tempFrame <- data.frame(timeVar, d01, X)
      # cvScores <- cvScores + sapply(predicted, function(x) ibs(pec(list(cbind(1, x)),
      #                                                   Surv(timeVar, d01) ~ 1,
      #                                                   data = tempFrame,
      #                                                   times = c(0,evalTimes),
      #                                                   reference = F,
      #                                                   exact = F,
      #                                                   verbose = F,
      #                                                   maxtime = max(evalTimes))))
    } else if (cvCriterion == "logLik") {
      cvScores <- cvScores + rowSums(do.call(cbind, cvResults))
    }
  }
  
  stopCluster(cl)
  
  ## calculate the solution of the best lambda value
  ## the intial parameter values are the ones obtained in the
  ## last CV step and for the corresponding lambda value
  finalPath <- penalizedPathBrier(timeVar = timeVar,
                                  d01 = d01,
                                  init = init,
                                  X = X, 
                                  breaks = breaks,
                                  evalTimes = evalTimes,
                                  cProbEvalTimes = cProbEvalTimes,
                                  cProbTimeVar = cProbTimeVar,
                                  lambdaSeq = lambdaSeq,
                                  penalType = penalType)
  
  list(optSol = finalPath[which.min(cvScores), ], 
       usedLambdas = lambdaSeq, 
       cvScores = cvScores,
       finalPath = finalPath,
       scale = attr(X, "scale"),
       center = attr(X, "center"),
       initBeta = initSol[1])
}


lambdaMaxMLE <- function(timeVar, d01, X, breaks, penalType = "totalVar"){
  ## calculate the smalles lambda values for which the model
  ## has all the differences equal to zero, this is a good
  ## value from which to start a solution path.
  
  nParam <- (ncol(X) + 1)
  
  # create some objects needed for the function that calculates the gradient
  R <- timeAtRisk(timeVar, breaks, 1, 1)
  
  lastIndex <- firstIndexIncInterval(timeVar, breaks)
  O <- matrix(0, ncol = ncol(R), nrow = nrow(R))
  O[cbind(1:length(timeVar), lastIndex)] <- d01
  
  X2 <- cbind(1, X)
  
  nIntervals <- length(breaks) - 1 
  transformMat <- paramTransMat(nIntervals,
                                nParam)
  
  ## fit a standard constant hazards model (exponential distribution).
  tempFrame <- data.frame(timeVar, d01, X)
  initSol <- - coef(survreg(Surv(timeVar, d01) ~. , data = tempFrame, dist = "exponential"))
  
  ## calculate the gradient of the MLE
  mleGrad <- mleGradient(O, R, X2, breaks, transformMat,
                         c(initSol, rep(0, (nParam) * (length(breaks) - 2))))
  ## get the lambda value
  if(penalType == "totalVar"){
    lambdaMax <- 0
    for(i in 1:nParam){
      paramIndices <- i + (1:(nIntervals - 1)) * nParam
      gradientNorm <- sqrt(sum(mleGrad[paramIndices]^2))
      if(gradientNorm > lambdaMax){
        lambdaMax <- gradientNorm
      }
    }
  } else if (penalType == "fusedLike"){
    lambdaMax <- max(abs(mleGrad[- (1:nParam)]))
  }
  lambdaMax
}

cvMLE <- function(timeVar, d01, X, breaks,
                  nLambda = 20, logLambdaRatio = 5, nFolds = 10,
                  evalTimes = NULL,
                  cvCriterion = "brier", penalType = "totalVar",
                  nCores, repeatedCV = 1){
  
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  X <- standardize(X)
  
  ## get the maximum lambda value of interest and create a sequence of lambda values
  logMaxLambda <- - log(lambdaMaxMLE(timeVar = timeVar, d01 = d01, 
                                     breaks = breaks, X = X, 
                                     penalType = penalType))
  lambdaSeq <- exp(- seq(logMaxLambda * (1 - 1 / nLambda), 
                         logMaxLambda * logLambdaRatio, 
                         length.out = nLambda))
  
  # sample CV indices
  randSample <- sample(nrow(X))
  CVIndicesList <- split(randSample, ceiling(seq_along(randSample)/(length(randSample)/nFolds))) 
  
  transformMat <- paramTransMat(nIntervals = length(breaks) - 1,
                                nParameters = ncol(X) + 1)
  # create a vector in which to save the cv performance scores
  cvScores <- rep(0, nLambda)

  # perform the CV steps
  for(l in 1:repeatedCV){
    cvResults <- foreach(i = 1:nFolds) %dorng% {
      source("./penalizedPCH.R")
      
      trTimeVar <- timeVar[- CVIndicesList[[i]]]
      testTimeVar <- timeVar[CVIndicesList[[i]]]
      
      trd01 <- d01[- CVIndicesList[[i]]]
      testd01 <- d01[CVIndicesList[[i]]]
      
      trX <- X[- CVIndicesList[[i]], ]
      testX <- X[CVIndicesList[[i]], ]
      
      foldPath <- penalizedPathMLE(timeVar = trTimeVar, d01 = trd01, 
                                   X = trX, breaks = breaks,
                                   lambdaSeq = lambdaSeq, 
                                   penalType = penalType)
      predicted <- list()
      cvScores <- rep(0, nLambda)
      
      for(j in 1:nLambda){
        if(cvCriterion == "brier"){
          predicted[[j]] <- predictPCH(transformMat %*% foldPath[j, ],
                                       breaks, testX, evalTimes)
        } else if(cvCriterion == "logLik") {
          cvScores[j] <- - sum(getLogLik(timeVar = testTimeVar, 
                                         d01 = testd01,  
                                         breaks = breaks, 
                                         X = testX,
                                         parameters = transformMat %*% foldPath[j, ]))
        }
      }
      
      if(cvCriterion == "brier"){
        predicted
      } else if (cvCriterion == "logLik"){
        cvScores
      }
      
    }
    
    if(cvCriterion == "brier"){
      predicted <- lapply(1:nLambda, function(x) matrix(0, nrow = nrow(X), ncol = length(evalTimes)))
      for(j in 1:nLambda){
        for(i in 1:nFolds){
          predicted[[j]][CVIndicesList[[i]], ] <- cvResults[[i]][[j]]
        }
      }
      tempFrame <- data.frame(timeVar, d01, X)
      
      cvScores <- sapply(predicted, function(x) ibs(pec(list(cbind(1, x)), 
                                                        Surv(timeVar, d01) ~ 1, 
                                                        data = tempFrame,
                                                        times = c(0,evalTimes),
                                                        reference = F,
                                                        exact = F,
                                                        verbose = F,
                                                        maxtime = max(evalTimes))))
    } else if (cvCriterion == "logLik") {
      cvScores <- cvScores + rowSums(do.call(cbind, cvResults))
    }
  }
  
  stopCluster(cl)
  
  ## calculate the solution of the best lambda value
  ## the intial parameter values are the ones obtained in the
  ## last CV step and for the corresponding lambda value
  finalPath <- penalizedPathMLE(timeVar = timeVar, d01 = d01, 
                                X = X, breaks = breaks,
                                lambdaSeq = lambdaSeq, 
                                penalType = penalType)

  list(optSol = finalPath[which.min(cvScores), ],
       usedLambdas = lambdaSeq, 
       cvScores = cvScores,
       finalPath = finalPath,
       scale = attr(X, "scale"),
       center = attr(X, "center"))
}

## function to remove the scaling
getParam <- function(solution, nIntervals, nParameters){
  transMat <- paramTransMat(nIntervals, nParameters)
  
  estPar <- as.numeric(transMat %*% solution[[1]])
  tempMat <- matrix(estPar, nrow = nParameters)
  tempMat[- 1, ] <- tempMat[- 1, ] / solution$scale
  tempMat[1, ] <- tempMat[1, ] - colSums(tempMat[- 1, ] * solution$center)
  
  estPar <- as.numeric(tempMat)
}

## function to calculate evaluation / cutpoint time-points
getEvalTimes <- function(eventTimes, nEvalTimes = 17, 
                         includeZero = T, extremeTails = c(0, 0.7), removeLargest){
  ## function that returns nEvaltimes equally spaced quantiles of eventTimes
  ## if includeZero is true an extra 0 is appended to the beginning of the vector
  estQuantiles <- quantile(eventTimes, seq(extremeTails[1], extremeTails[2], length.out =  nEvalTimes + 2))
  if(includeZero) estQuantiles <- c(0, estQuantiles[- c(1, nEvalTimes + 2)])
  else estQuantiles <- estQuantiles[- c(1, nEvalTimes + 2)]
  
  duplicatedEntries <- duplicated(estQuantiles)
  estQuantiles <- estQuantiles[!duplicatedEntries]
  
  if(removeLargest){
    estQuantiles <- estQuantiles[estQuantiles != max(eventTimes)]
  }
  estQuantiles
}