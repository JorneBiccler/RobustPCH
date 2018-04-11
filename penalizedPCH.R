## load the necessary packages
library(survival)
#########################
### GENERAL FUNCTIONS ###
#########################

firstIndexIncInterval <- function(evalTimes, breaks){
  ## function that returns the last index of the breaks vector
  ## for which the evalTimes value is larger than the corresponding breaks value
  tempMatrix <- matrix(breaks,
                       nrow = length(evalTimes),
                       ncol = length(breaks),
                       byrow = T)
  rowSums(tempMatrix <= evalTimes)
}

organizeParameters <- function(nCoef, nIntervals, parameters){
  ## function that returns a list containing a vector of intercepts
  ## and a matrix of coefficients in which each row corresponds with
  ## an interval and contains the coefficients of the covariates.
  coef <- matrix(parameters[- seq(from = 1, to = length(parameters), by = nCoef + 1)],
                 ncol = nCoef, byrow = T)
  intercepts <- parameters[seq(from = 1, to = length(parameters), by = nCoef + 1)]
  list(coef = coef, intercepts = intercepts)
}


predictPCH <- function(parameters, breaks, X, evalTimes){
  ## Function to calculate the predicted survival of new observations
  ## at the times provided.
  ## ARGS:
  ##   parameters: a vector containing the intercepts and coefficients,
  ##               the structure should be as follows: first the first intercept
  ##               then the coefficients corresponding to the first interval
  ##               after which the same is done for the second interval etc.
  ##   breaks: the breakpoints 
  ##   X: covariate matrix
  ##   evalTimes: times for which predictions should be returned
  
  nIntervals <- length(breaks) - 1
  nCoef <- ncol(X)
  organizedParam <- organizeParameters(nCoef, nIntervals, parameters)
  
  intercept <- organizedParam$intercepts
  coef <- organizedParam$coef
  
  ## case of an exponential model (i.e. breaks = c(0, Inf) or length(breaks) == 2)
  if(length(breaks) == 2){
    estCumHaz <- exp(intercept[1] + as.numeric(X %*% coef[1, ])) %*% t(evalTimes)
    
    ## jump out of the function and return the survival estimates
    predicted <- exp(- estCumHaz)
    return(predicted)
  }
  
  ## fill a matrix with the hazard contribution of the intervalls that are fully
  ## covered at the evaluation times


  if(length(breaks) == 3){
    
    fullHaz <- rep(diff(breaks[- (nIntervals + 1)]) * 
                     exp(intercept[- nIntervals]), each = nrow(X)) * 
      exp(X %*% coef[- (nIntervals), ])
    
    ## if there are only 2 intervals, the result of the apply is a vector,
    ## to make everything consistent this is turned into the correct matrix
    fullHaz <- as.matrix(fullHaz)
  } else {
    fullHaz <- rep(diff(breaks[- (nIntervals + 1)]) * 
                     exp(intercept[- nIntervals]), each = nrow(X)) * 
      exp(X %*% t(coef[- (nIntervals), ]))
    ## sum the different chunks into a cumulative hazards
    fullHaz <- t(apply(fullHaz, 1, cumsum))
  }
  
  ## create some help-vectors to show to which interval the different evaluation times belong
  intPlus1 <- firstIndexIncInterval(evalTimes, breaks)
  zeroIntervalsCovered <- sum(intPlus1 == 1)
  fullyCoveredIntervals <- intPlus1[intPlus1 != 1]
  
  ## calculate the cumulative hazard
  estCumHaz <- cbind(matrix(0, 
                            nrow = nrow(X),
                            ncol = zeroIntervalsCovered), 
                     fullHaz[, fullyCoveredIntervals - 1]) +
    (rep((evalTimes - breaks[intPlus1]) * exp(intercept[intPlus1]), each = nrow(X)) * 
       exp(X %*% t(coef[intPlus1, ])))
  
  ## return the survival estimates
  predicted <- exp(- estCumHaz)
  predicted
}

timeAtRisk <- function(evalTimes, breaks, ncolRep = 1, nrowRep = 1){
  ## function that creates a matrix with as components the time in the interval.
  ## the rows correspond with the evalTimes the columns with the intervals
  ## each row is repeated nrowRep times (usually the number of observations)
  ## and each column is repeated ncolRepTimes (usually the number of covariates)
  
  nInterval <- length(breaks) - 1
  diffMatrix <- matrix(diff(breaks),
                       byrow = T,
                       nrow = length(evalTimes),
                       ncol = nInterval)
  diffMatrix[matrix(breaks[-1], nrow = length(evalTimes), ncol = nInterval, byrow = T) > evalTimes] <- 0
  
  colIndices <- firstIndexIncInterval(evalTimes, breaks)
  diffMatrix[cbind(1:length(evalTimes), colIndices)] <- evalTimes - breaks[colIndices]
  
  if(length(breaks) == 2){
    diffMatrix <- as.matrix(diffMatrix)
    diffMatrix <- as.matrix(diffMatrix[rep(1:nrow(diffMatrix), each = nrowRep), ])
    return(diffMatrix)
  }
  
  diffMatrix <- diffMatrix[rep(1:nrow(diffMatrix), each = nrowRep), ]
  diffMatrix[, rep(1:ncol(diffMatrix), each = ncolRep)]
}

extendedY <- function(evalTimes, timeVar){
  ## function that transforms the time variable into 
  ## a new outcome variable with length(evaltimes) * length(timeVar) entries
  ## such that the entries are as specified in the least squares equations.
  
  rep(timeVar, length(evalTimes)) > rep(evalTimes, each = length(timeVar))
}

paramTransMat <- function(nIntervals, nParameters){
  ## returns a matrix to go from the differences parameterization to
  ## the standard parameterization
  matrixOfOnes <- kronecker(matrix(1, ncol = nIntervals, nrow = nIntervals),
                            diag(nParameters))
  matrixOfOnes[upper.tri(matrixOfOnes)] <- 0
  matrixOfOnes
}

proximalStep <- function(curSol, nParam, gradient, stepsize, lambda, penalType = "totalVar"){
  ## Perform a proximal gradient descent step corresponding to a group lasso of differences penalty
  ## the curSol is the current solution (in differences representation),
  ## the gradient is the gradient (in differences representation),
  ## the stepsize is the stepzise that should be taken during the proximal gradient descent step
  ## lambda is the penalization parameter
  
  nInterval <- length(curSol) / nParam
  newSol <- curSol - stepsize * gradient
  
  if(lambda != 0 & penalType == "totalVar"){
    for(i in 1:nParam){
      paramIndices <- i + (1:(nInterval - 1)) * nParam
      newSolNorm <- sqrt(sum(newSol[paramIndices]^2))
      if(newSolNorm <= lambda * stepsize){
        newSol[paramIndices] <- 0
      } else {
        newSol[paramIndices] <- newSol[paramIndices] - 
          lambda * stepsize * (newSol[paramIndices]) / newSolNorm
      }
    }
  } else if(lambda != 0 & penalType == "fusedLike"){
    newSol[- (1:nParam)] <- ifelse(newSol[- (1:nParam)] <= lambda * stepsize,
                                   0,
                                   newSol[- (1:nParam)] - 
                                     lambda * stepsize * (newSol[- (1:nParam)]) / 
                                     abs(newSol[- (1:nParam)]))
  }
  newSol
}

getLogLik <- function(timeVar, d01, breaks, X, parameters){
  ## function that returns the log-likelihood contribution of
  ## each observation.
  
  ## order the observations according to the time vaiable
  ## this is done to make some computations easier
  ## possible speed increase: rewrite the predictPCH fuction
  ## to deal with non-sorted evalTimes
  obsOrder <- order(timeVar)
  ordTimeVar <- timeVar[obsOrder]
  ordD01 <- d01[obsOrder]
  ordX <- X[obsOrder, ]
  
  ## get the predicted survival at the time times
  predicted <- diag(predictPCH(parameters, breaks, X = ordX, evalTimes = ordTimeVar))
  ## calculate the hazards rates at the time times
  indices <- firstIndexIncInterval(ordTimeVar, breaks)
  organizedParam <- organizeParameters(ncol(X), nIntervals, parameters)
  
  haz <- exp(organizedParam$intercept[indices] + 
               rowSums(organizedParam$coef[indices, ] * ordX))
  
  ## calculate the log likelihood contributions and return them in the original ordering
  result <- rep(NA, length(timeVar))
  result[obsOrder] <- log(haz^ordD01 * predicted)
  result
}

getTrimmedBrier <- function(cProbEvalTimes, cProbTimeVar, timeVar, d01, 
                            evalTimes, predictedSurvival, trim){
  ## returns a trimmed mean of the Brier residuals
  ## the trimming is performed one-sided
  n <- length(timeVar)
  
  exY <- extendedY(evalTimes = evalTimes,
                   timeVar = timeVar)
  
  IPWWeights <- exY^(1 - d01)
  
  IPWWeights <- IPWWeights / pmax(as.numeric(cProbEvalTimes), cProbTimeVar)
  
  ## this uses the composite midpoint rule.
  ## note: doesn't matter if 0 is included in evaltimes, the difference is going to be zero
  ## and hence will have weigh t= 0 -> no contribution to the score function
  ## normalize the weights
  
  brierResid <- IPWWeights * (exY - predictedSurvival) ^ 2
  brierResid <- matrix(brierResid, nrow = n)

  brierResid <- rowSums(brierResid)
  
  sortedBrierResid <- sort(brierResid)
  trimedBrierResid <- sortedBrierResid[1:floor(trim * length(sortedBrierResid))]
  
  mean(trimedBrierResid)
}

#############################################
### PCH MODEL WITH A PENALIZED BRIER LOSS ###
#############################################

brierGradient <- function(IPWWeights, timeAtRiskMat, exY, evalTimes, transformMat, 
                          largeX, X, X2, parameters, breaks){
  ## returns the gradient of the Brier score
  ## note: doesn't matter if 0 is included in evaltimes, the difference is going to be zero
  ## and hence will have weigh t= 0 -> no contribution to the score function
  evalTimesWeights <- diff(c(0, evalTimes))
  evalTimesWeights <- evalTimesWeights / mean(evalTimesWeights)
  
  ## transform the values to the original parameterizaiton
  oPar <- transformMat %*% parameters
  ## get the predicted survival at the evaluation times
  predictedSurvival <- as.numeric(predictPCH(oPar, breaks, X, evalTimes))
  
  ## get the gradient of the original parameterization
  gradient <- NULL
  
  ctVec <- rep(evalTimesWeights, each = nrow(X)) * IPWWeights * (exY - predictedSurvival) * predictedSurvival
  
  for(i in 1:(length(breaks) - 1)){
    gradient <- c(gradient, t(ctVec * timeAtRiskMat[, i] *
                              as.numeric(exp(X2 %*% oPar[1:(ncol(X2)) + (i - 1) * ncol(X2)]))) %*% largeX)
  }
  
  ## transform the gradient to the gradient of the differences parameterization
  ## and return the vector
  t(transformMat) %*% gradient / (length(IPWWeights))
}

penalizedBrier <- function(timeVar, d01, X, evalTimes, 
                           cProbEvalTimes, cProbTimeVar,
                           breaks, init, 
                           lambda,
                           penalType = "totalVar",
                           epsilon = 10^(-6),
                           stepSize = 0.01,
                           bbSteps = T,
                           maxStep = 10^(5),
                           minStep = 10^(-5)){
  ## Returns the penalized Brier solution.
  
  # init is the initial solution (in the differences parameterization)
  # perform the first step and set-up some vectors
  oldSol <- init
  initStep <- stepSize
  
  ## calculate some data that's used in the gradient
  exY <- extendedY(evalTimes = evalTimes, 
                   timeVar = timeVar)
  IPWWeights <- exY^(1 - d01)
  
  IPWWeights <- IPWWeights / pmax(as.numeric(cProbEvalTimes), cProbTimeVar)
  timeAtRiskMat <- timeAtRisk(evalTimes, breaks, nrowRep = length(timeVar), 1)

  transformMat <- paramTransMat(nIntervals = length(breaks) - 1,
                                nParameters = ncol(X) + 1)
  X2 <- cbind(1, X)
  largeX <- X2[rep(1:length(timeVar), length(evalTimes)), ]
  
  ## get the gradient of the brier score at the initial solution
  gradientOld <- brierGradient(IPWWeights = IPWWeights, 
                               timeAtRiskMat = timeAtRiskMat, 
                               exY = exY, 
                               evalTimes = evalTimes, 
                               transformMat = transformMat, 
                               largeX = largeX,
                               X = X, 
                               X2 = X2,
                               parameters = oldSol, 
                               breaks =  breaks)
  
  ## perform a the proximal step
  newSol <- proximalStep(oldSol,
                         ncol(X) + 1,
                         gradientOld,
                         stepSize,
                         lambda,
                         penalType)
  
  gradientCur <- gradientOld
  solDif <- oldSol - newSol
  i <- 0
  ## keep performing the proximal gradient steps until a convergence
  ## criterion is met
  while(sum(abs(solDif)) / sum(abs(oldSol)) > epsilon | i < 1){
    #  Barzilai-Borwein steps
    gradientCur <- brierGradient(IPWWeights = IPWWeights, 
                                 timeAtRiskMat = timeAtRiskMat, 
                                 exY = exY, 
                                 evalTimes = evalTimes, 
                                 transformMat = transformMat, 
                                 largeX = largeX,
                                 X = X,
                                 X2 = X2, 
                                 parameters = newSol, 
                                 breaks =  breaks)
    
    if(sum(solDif^2, na.rm = T) == 0){
      break()
    }
    
    if(bbSteps){
      gradientDif <- gradientCur - gradientOld
      gradientOld <- gradientCur
      stepSize <- sum(gradientDif * solDif) / sum(gradientDif^2)
      if(is.na(stepSize) | stepSize < 0){
        stepSize <- initStep
      } else if(stepSize < minStep) {
        stepSize <- minStep
      } else if(stepSize > maxStep){
        stepSize <- maxStep
      }
    }
    oldSol <- newSol
    newSol <- proximalStep(newSol,
                           ncol(X) + 1,
                           gradientCur,
                           stepSize,
                           lambda, 
                           penalType)
    solDif <- newSol - oldSol
    i <- i + 1
  }
  newSol
}

penalizedPathBrier <- function(timeVar, d01, X, evalTimes, 
                               cProbEvalTimes, cProbTimeVar,
                               init,
                               breaks, lambdaSeq, penalType = "totalVar", ...){
  ## calculate a path of solutions, this helps with convergence issues
  ## by using warm-starts and a "smart" initial solution
  ## the result is a matrix in which each row contains esitmated coefficients
  ## for the corresponding lambda value

  # create a matrix in which the coefficients are saved
  solutions <- matrix(0, nrow = length(lambdaSeq), ncol = (ncol(X) + 1) * (length(breaks) - 1))

  ## get the first solution using the robust estimates as intial values
  solutions[1, ] <- penalizedBrier(timeVar = timeVar, d01 = d01, X = X, cProbEvalTimes = cProbEvalTimes,
                                   cProbTimeVar = cProbTimeVar, evalTimes = evalTimes, breaks = breaks,
                                   init = init, 
                                   lambda = lambdaSeq[1],
                                   penalType = penalType)

  ## calculate the rest of the path by using warm starts
  ## i.e. use the solution of the previous lambda value
  ## as starting value 
  for(i in 2:length(lambdaSeq)){
    solutions[i, ] <- penalizedBrier(timeVar = timeVar, d01 = d01, X = X, cProbEvalTimes = cProbEvalTimes,
                                     cProbTimeVar = cProbTimeVar, evalTimes = evalTimes, breaks = breaks,
                                     init = solutions[i - 1, ], lambda = lambdaSeq[i],
                                     penalType = penalType)
  }
  solutions
}

####################################
### PENALIZED PCH MODEL WITH MLE ###
####################################

mleGradient <- function(O, R, X2, breaks, transformMat, parameters){
  ## O: matrix with 1 if the event happens within the interval (column) for observation i (rows)
  ## R: time at risk matrix within each interval (rows: observations, columsn: intervals)
  ## X2: covariate matrix with the first column consisting of ones (the intercept)
  ## rest: see other comments somewhere else
  
  ## function that calculates the gradient of the log likelihood
  
  ## calculate the parameters in the standard parameterization
  originalParameters <- transformMat %*% parameters
  
  # get the gradient by using the original parameterization
  gradient <- NULL
  for(j in 1:(length(breaks) - 1)){
    gradient <- c(gradient,
                  t(O[, j]) %*% X2 - 
                    t(R[, j]) %*% (as.numeric(exp(X2 %*% originalParameters[1:ncol(X2) + (j - 1) * ncol(X2)])) 
                                   * X2))
  }
  
  # transform the gradient into the gradient of the differences parameterization
  - 1 * t(transformMat) %*% gradient / (nrow(X2))
}

penalizedMLE <- function(timeVar, d01, X, 
                         breaks, init, 
                         lambda, penalType = "totalVar",
                         epsilon = 10^(-6),
                         stepSize = 0.01,
                         maxStep = 10^(5),
                         minStep = 10^(-5),
                         bbSteps = T){
  # init is the initial solution (in the differences parameterization)
  # perform the first step and set-up some vectors
  oldSol <- init
  initStep <- stepSize
  
  # create some objects needed in the gradient
  R <- timeAtRisk(timeVar, breaks, 1, 1)
  
  lastIndex <- firstIndexIncInterval(timeVar, breaks)
  
  O <- matrix(0, ncol = ncol(R), nrow = nrow(R))
  O[cbind(1:length(timeVar), lastIndex)] <- d01
  
  X2 <- cbind(1, X)
  
  nIntervals <- length(breaks) - 1
  nParameters <- (ncol(X) + 1)
  transformMat <- paramTransMat(nIntervals,
                                nParameters)
  
  ## get the gradient at the provided intial solution
  gradientOld <- mleGradient(O, R, X2, breaks, transformMat, oldSol)
  
  if(any(is.na(gradientOld))){
    return(oldSol)
  }
    
  ## perform the proximal gradient step
  newSol <- proximalStep(oldSol,
                         ncol(X) + 1,
                         gradientOld,
                         stepSize,
                         lambda,
                         penalType)
  ## let's force the entrance of the loop:
  solDif <- oldSol
  
  ## perform the proximal gradient steps until the convergence
  ## criterion is met
  while(sum(abs(solDif)) / sum(abs(oldSol)) > epsilon){

    gradientCur <- mleGradient(O, R, X2, breaks, transformMat, newSol)
    if(any(is.na(gradientCur))){
      return(oldSol)
    }
    solDif <- newSol - oldSol
    gradientDif <- gradientCur - gradientOld
    
    if(sum(gradientDif^2, na.rm = T) == 0){
      break()
    }
    
    #  Barzilai-Borwein steps    
    if(bbSteps){
      gradientOld <- gradientCur
      stepSize <- sum(gradientDif * solDif) / sum(gradientDif^2)
      if(is.na(stepSize) | stepSize < 0){
        stepSize <- initStep
      } else if(stepSize < minStep) {
        stepSize <- minStep
      } else if(stepSize > maxStep){
        stepSize <- maxStep
      }
    }
    oldSol <- newSol
    newSol <- proximalStep(newSol,
                           ncol(X) + 1,
                           gradientCur,
                           stepSize,
                           lambda,
                           penalType)
  }
  newSol
}

penalizedPathMLE <- function(timeVar, d01, X, breaks, 
                             lambdaSeq, penalType = "totalVar"){
  ## create a path of penalized MLE's.
  ## This tends to be more efficient than calculating them one by one
  ## since warm starts etc. can be used. The result is a matrix with rows
  ## corresponding to the solutions of different lambda values
  
  ## create a matrix in which the path will be saved
  solutions <- matrix(0, nrow = length(lambdaSeq), ncol = (ncol(X) + 1) * (length(breaks) - 1))
  
  ## use the standards survreg package to obtain estimates of a constant hazards (exponential) model
  ## this function is more stable and faster than using the proximal gradient method.
  tempFrame <- data.frame(timeVar, d01, X)

  ## try to fit an expontial model, the values are used as initial solution
  ## this corresponds with a model that has a very large penalty term
  maxIterInit <- 1000
  initRegMod <- survreg(Surv(timeVar, d01) ~. , data = tempFrame, 
                        dist = "exponential", maxiter = maxIterInit)
  initSol <- - coef(initRegMod)
  
  
  ## if even an exponential model cannot be fitted, there is no real use in trying to fit
  ## models that include breaks, hence the solution of the exponential model si returned.
  if(initRegMod$iter == maxIterInit){
    return(matrix(initSol, 
                  nrow = length(lambdaSeq), 
                  ncol = (ncol(X) + 1) * (length(breaks) - 1),
                  byrow = T))
  }

  ## use the MLE of the expenential model as initial values for the first solution
  solutions[1, ] <- penalizedMLE(timeVar = timeVar, d01 = d01, X = X, breaks = breaks,
                                 init = c(initSol, rep(0, (ncol(X) + 1)* (length(breaks) - 2))), 
                                 lambda = lambdaSeq[1], penalType = penalType)
  
  ## use warm-starts to calculate the solutions for the other lambda values
  for(i in 2:length(lambdaSeq)){
    solutions[i, ] <- penalizedMLE(timeVar = timeVar, d01 = d01, X = X, breaks = breaks,
                                   init = solutions[i - 1, ], 
                                   lambda = lambdaSeq[i], penalType = penalType)
  }
  solutions
}
