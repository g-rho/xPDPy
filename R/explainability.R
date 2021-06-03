#' @importFrom colorRamps blue2red
#' @importFrom GISTools add.alpha
#' @importFrom graphics axis
#' @importFrom graphics layout
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics text
#' @importFrom graphics legend
#' @importFrom stats aggregate
#' @importFrom stats predict
#' @importFrom dplyr sample_n
#' @importFrom dplyr group_by
#' @importFrom dplyr group_indices
#' @importFrom dplyr left_join
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom magrittr %>%




#' @title Explainability and matchplot for PDP
#'
#' @description Computes explainability and matchplots for a partial dependence function.
#'
#' @param model       A model with corresponding predict function that returns numeric values.
#' @param x           Data frame.
#' @param vnames      Character vector of the variable set for which the patial dependence function is to be computed.
#' @param viz         Logical specifying whether a matchplot should be created.
#' @param parallel    Logical specifying whether computation should be parallel.
#' @param sample.frac fraction-size for sampling of x.
#' @param pfunction   User generated predict function with arguments \code{pfunction(model, newdata)}.
#' @param ...         Further arguments to be passed to the \code{predict()} method of the \code{model}.
#'
#' @return Numeric value with the explainability of the partial dependence function for the variables specified in \code{vnames}.
#'
#' @export
#'
#' @examples
#' library(pdp)
#' library(randomForest)
#' data(boston)
#' set.seed(42)
#' boston.rf <- randomForest(cmedv ~ ., data = boston)
#' xpy(boston.rf, boston, c("lstat", "rm"))
#'
#' # example 2: user-defined predict function for classification
#' library(klaR)
#' data("GermanCredit")
#'
#' library(randomForest)
#' rfGC <- randomForest(credit_risk ~., data = GermanCredit)
#'
#' # define predict function that returns continuous value (here: prob(yes))
#' pfunc = function(model, x) predict(model, x, type = "prob")[,2]
#' xpy(rfGC, GermanCredit, vnames = "duration", parallel = FALSE, pfunction = pfunc)

#'
#' @author \email{gero.szepannek@@web.de}
#'
#' @references \itemize{
#'     \item Szepannek, G. (2019): How Much Can We See? A Note on Quantifying Explainability of Machine Learning Models,
#'           \href{https://arxiv.org/abs/1910.13376}{\emph{arXiv:1910.13376 [stat.ML]}}.
#'   }
#'
#' @rdname xpy
xpy <- function(model, x, vnames, viz = TRUE, parallel = TRUE, sample.frac = 1, pfunction = NULL, ...){

  # number of unique observations to give to processors each step
  ssize <- 20

  # Sampling
  x <- sample_n(x, round(nrow(x)*sample.frac))

  # required packages to be loaded in parallel R-threads
  # this should be extracted from variable model
  pckgs = c("randomForest")

  # checking for predict_function parameter
  if (is.null(pfunction)) {
    # predict_function not specified
    pfunction <- predict
  }

  # adding indices for unique vnames-values
  x <- x %>% group_by(x[,names(x) %in% vnames])
  x$ind <- x %>% group_indices()

  # ignore duplicates in xs
  xs <- x[!duplicated(x[,"ind"]), names(x) %in% append("ind", vnames)]

  # keep column ind for join
  xrest <- x[, !names(x) %in% append("ind", vnames), drop = FALSE]

  # checking if parallel computation should be used
  if (parallel) {
    # calculate number of steps of size ssize
    steps <- ceiling(nrow(xs)/ssize)

    # create cluster of all available cores on host
    cl <- makeCluster(detectCores())    # uses parallel::detectCores() + ::makeCluster

    # register cluster for foreach
    registerDoParallel(cl)

    # parallel loop 1 to steps, results are concatenated ('c'),
    # load pckgs in every thread
    pdps <- foreach(i = 1:steps, .combine = 'c',.packages = pckgs) %dopar% {  # uses foreach::foreach

      # beginning of subset
      ifrom <- (i-1)*ssize +1

      # end of subset
      ito <- min(i*ssize, nrow(xs))

      # cartesian product
      xx <- merge(xs[ifrom:ito,], xrest, by = NULL)

      # predict
      xx$yhat <- pfunction(model, xx, ...)

      # compute average prediction
      return(pred = tapply(xx$yhat, xx$ind, mean))
    }

    stopCluster(cl)

    # convert results to dataframe
    pdps <- stack(pdps)

    x$pred <-  pfunction(model, x, ...)
    avpred <- mean(x$pred)

    # merge results
    preds <- merge(pdps, x[, c("pred","ind"), drop = FALSE])
  }
  else {
    # parallel = FALSE
    pdps <- NULL
    from <- seq(1, nrow(xs), by = ssize)
    for (i in from){
      xx <- merge(xs[i:min((i+ssize-1),nrow(xs)),], xrest , by.x = NULL, by.y = NULL)
      ID <- xx[,"ind"]
      xx$yhat <- pfunction(model, xx, ...)
      pdps <- rbind(pdps, aggregate(xx$yhat, list(ID), mean))
    }
    x$pred <- pfunction(model, x, ...)
    avpred <- mean(x$pred)
    preds <- merge(pdps, x[, c("pred","ind"), drop = FALSE], by.x = "Group.1", by.y = "ind")
    colnames(preds)[2] <- "values"
  }

  # plotting
  if(viz){
    rnge <- range(c(range(preds$pred), range(preds$values)))
    plot(preds$values, preds$pred, xlim = rnge, ylim = rnge, xlab = "PDP", ylab = "Prediction", pch = 4, main = "PDP vs. Predictions")
    lines(rnge, rnge, lty = "dotted")
  }

  ASE <- mean((preds$pred - preds$values)^2)
  ASE0 <- mean((preds$pred - avpred)^2)
  xty <- 1 - ASE / ASE0
  xty
}


#' @title Forward variable selection for PDP explanation
#'
#' @description Computes forward variable selection for partial dependence function based on explainability.
#'
#' @param model       A model with corresponding predict function that returns numeric values.
#' @param x           Data frame.
#' @param target      Character specifying the name of the (numeric) target variable (must be contained in data).
#' @param parallel    Logical specifying whether computation should be parallel.
#' @param sample.frac fraction-size for sampling of x.
#' @param pfunction   User generated predict function with arguments \code{pfunction(model, newdata)}.
#' @param ...         Further arguments to be passed to the \code{predict()} method of the \code{model}.
#'
#' @return Object of class \code{vsexp} containing the following elements:
#' @return \item{selection.order}{Vector of variable names in order of their entrance to the PD function during the variable selection process.}
#' @return \item{explainability}{Explainabilities after the different iterations of the variable selection.}
#' @return \item{details}{A data frame containing the explainabilities for all variables (rows) if they were entered in the model at each step (columns).}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(pdp)
#' library(randomForest)
#' data(boston)
#' set.seed(42)
#' boston.rf <- randomForest(cmedv ~ ., data = boston)
#' vs <- fw.xpy(boston.rf, boston, "cmedv")
#' vs
#'
#' plot(vs)
#'}
#'
#' @author \email{gero.szepannek@@web.de}
#'
#' @references \itemize{
#'     \item Szepannek, G. (2019): How Much Can We See? A Note on Quantifying Explainability of Machine Learning Models,
#'           \href{https://arxiv.org/abs/1910.13376}{\emph{arXiv:1910.13376 [stat.ML]}}.
#'   }
#'
#' @rdname fw.xpy
fw.xpy <- function(model, x, target, parallel = TRUE, sample.frac = 1, pfunction = NULL, ...){
  
  # number of unique observations to give to processors each step
  ssize <- 20
  
  # Sampling
  x <- sample_n(x, round(nrow(x)*sample.frac))
  
  # required packages to be loaded in parallel R-threads
  # this should be extracted from variable model
  pckgs = c("randomForest")
  
  # checking for predict_function parameter
  if (is.null(pfunction)) {
    # predict_function not specified
    pfunction <- predict
  }
  
  # Initialization and selection of first variable
  n <- 1
  cat("Step", n, "\n")
  
  sel   <- NULL
  trace <- NULL
  nms <- nms.full <- names(x)[-which(names(x) == target)]
  xpys <- rep(NA, length(nms))
  names(xpys) <- nms
  
  x$pred <-  pfunction(model, x, ...)
  avpred <- mean(x$pred)
  
  if(parallel) {
    cl <- makeCluster(detectCores())
    registerDoParallel(cl)
  }
  
  # helper function
  fw.xpy.help <- function(model, xh, vnames, ssize, pckgs, cl, parallel, pfunction, ...){
    xh <- xh %>% group_by(xh[,names(xh) %in% vnames])
    xh$ind <- xh %>% group_indices()
    xs <- xh[!duplicated(xh[,"ind"]), names(xh) %in% append("ind", vnames), drop = FALSE]
    xrest <- xh[, !names(xh) %in% append("ind", vnames), drop = FALSE]
    
    if(parallel) {
      steps <- ceiling(nrow(xs)/ssize)
      pdps <- foreach(i = 1:steps, .combine = 'c',.packages = pckgs) %dopar% {
        ifrom <- (i-1)*ssize +1
        ito <- min(i*ssize, nrow(xs))
        xx <- merge(xs[ifrom:ito,], xrest, by = NULL)
        xx$yhat <- pfunction(model, subset(xx, select = -c(ind)), ...)
        return(pred = tapply(xx$yhat, xx$ind, mean))
      }
      pdps <- stack(pdps)
      
      # this could be done via JOIN
      preds <- merge(pdps, xh[, c("pred","ind"), drop = FALSE])
    }
    else {
      # parallel = FALSE
      pdps <- NULL
      from <- seq(1, nrow(xs), by = ssize)
      for (i in from){
        xx <- merge(xs[i:min((i+ssize-1),nrow(xs)),], xrest , by.x = NULL, by.y = NULL)
        ID <- xx[,"ind"]
        xx$yhat <- pfunction(model, xx, ...)
        pdps <- rbind(pdps, aggregate(xx$yhat, list(ID), mean))
      }
      preds <- merge(pdps, xh[, c("pred","ind"), drop = FALSE], by.x = "Group.1", by.y = "ind")
      colnames(preds)[2] <- "values"
    }
    ASE <- mean((preds$pred - preds$values)^2)
    ASE0 <- mean((preds$pred - avpred)^2)
    xty <- 1 - ASE / ASE0
    xty
  }
  
  
  for(v in nms) xpys[which(names(xpys) == v)] <- fw.xpy.help(model, x, v, viz = F, ssize, pckgs, cl, parallel, pfunction,  ...)
  
  sel   <- c(sel, which.max(xpys))
  trace <- c(trace, max(xpys, na.rm = T))
  print(xpys)
  cat("\n", nms.full[sel], max(xpys, na.rm = T), "\n\n")
  
  # forward selection variables such that explainability is maximized
  while(length(nms) > 1){
    n <- n + 1
    cat("Step", n, "\n")
    
    nms <- nms.full[-sel]
    xpys <- cbind(xpys, NA)
    for(v in nms) xpys[which(rownames(xpys) == v), ncol(xpys)] <- fw.xpy.help(model, x, c(names(sel), v), viz = F, ssize, pckgs, cl, parallel, pfunction, ...)
    
    sel <- c(sel, which.max(xpys[,ncol(xpys)]))
    colnames(xpys) <- c(paste("Step", 1:n))
    trace <- c(trace, max(xpys, na.rm = T))
    
    print(xpys)
    cat("\n", nms.full[sel], max(xpys[,ncol(xpys)], na.rm=T), "\n\n")
  }
  
  if(parallel) {
    stopCluster(cl)
  }
  
  res <- list(selection.order = sel, explainability = trace, details = xpys)
  class(res) <- "vsexp"
  return(res)
}

#' @export
plot.vsexp <- function(x, ...){
  plot(0:length(x$explainability), c(0,x$explainability), type = "l", xaxt = "n", xlab = "", ylim = c(0, 1), ylab = "explainability")
  axis(1, at = 1:length(x$selection.order), labels = names(x$selection.order), las = 2)
}


#' @export
print.vsexp <- function(x, ...) print(cbind(x$selection.order, x$explainability))

###################previous version of fw.xpy()
###################
#' fw.xpy <- function(model, x, target, parallel = TRUE, sample.frac = 1, pfunction = NULL, ...){
#'   
#'   
#'   
#'   # Initialization and selection of first variable
#'   n <- 1
#'   cat("Step", n, "\n")
#'   
#'   sel   <- NULL
#'   trace <- NULL
#'   nms <- nms.full <- names(x)[-which(names(x) == target)]
#'   xpys <- rep(NA, length(nms))
#'   names(xpys) <- nms
#'   
#'   for(v in nms) xpys[which(names(xpys) == v)] <- xpy(model, x, v, viz = F, ...)
#'   
#'   sel   <- c(sel, which.max(xpys))
#'   trace <- c(trace, max(xpys, na.rm = T))
#'   print(xpys)
#'   cat("\n", nms.full[sel], max(xpys, na.rm = T), "\n\n")
#'   
#'   # forward selection variables such that explainability is maximized
#'   while(length(nms) > 1){
#'     n <- n + 1
#'     cat("Step", n, "\n")
#'     
#'     nms <- nms.full[-sel]
#'     xpys <- cbind(xpys, NA)
#'     for(v in nms) xpys[which(rownames(xpys) == v), ncol(xpys)] <- xpy(model, x, c(names(sel), v), viz = F, ...)
#'     
#'     sel <- c(sel, which.max(xpys[,ncol(xpys)]))
#'     colnames(xpys) <- c(paste("Step", 1:n))
#'     trace <- c(trace, max(xpys, na.rm = T))
#'     
#'     print(xpys)
#'     cat("\n", nms.full[sel], max(xpys[,ncol(xpys)], na.rm=T), "\n\n")
#'   }
#'   
#'   res <- list(selection.order = sel, explainability = trace, details = xpys)
#'   class(res) <- "vsexp"
#'   return(res)
#' }
#' 
#' 
#' #' @export
#' plot.vsexp <- function(x, ...){
#'   plot(0:length(x$explainability), c(0,x$explainability), type = "l", xaxt = "n", xlab = "", ylim = c(0, 1), ylab = "explainability")
#'   axis(1, at = 1:length(x$selection.order), labels = names(x$selection.order), las = 2)
#' }
#' 
#' 
#' #' @export
#' print.vsexp <- function(x, ...) print(cbind(x$selection.order, x$explainability))



#' @title Explanation gap visualization
#'
#' @description Visualization of 2D PDP vs. unexplained residual predictions.
#'
#' @param model   A model with corresponding predict function that returns numeric values.
#' @param x       Data frame.
#' @param vnames  Character vector of the variable set for which the patial dependence function is to be computed.
#' @param type    Character, either \code{"pdp"}, \code{"gap"} or \code{"both"} specifying the meaning of the colours
#' of scatter: either the value of the PD function or the differences between PD and the model's predcitions.
#' In case of \code{"both"} two plots are created.
#' @param depth   Integer specifiying the number colours in the heat map.
#' @param alpha   Numeric value for alpha blending of the points in the scatter plot.
#' @param right   Position where to place the legend relative to the range of the x axis.
#' @param top     Position where to place the legend relative to the range of the y axis.
#' @param digits  Nuber of digits for rounding in the legend.
#' @param parallel    Logical specifying whether computation should be parallel.
#' @param sample.frac fraction-size for sampling of x.
#' @param pfunction   User generated predict function with arguments \code{pfunction(model, newdata)}.
#' @param ...     Further arguments to be passed to the \code{predict()} method of the \code{model}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(pdp)
#' library(randomForest)
#' data(boston)
#' set.seed(42)
#' boston.rf <- randomForest(cmedv ~ ., data = boston)
#' pdp2d(boston.rf, boston, c("lstat", "rm"), type = "both")
#' }
#'
#' @author \email{gero.szepannek@@web.de}
#'
#' @references \itemize{
#'     \item Szepannek, G. (2019): How Much Can We See? A Note on Quantifying Explainability of Machine Learning Models,
#'           \href{https://arxiv.org/abs/1910.13376}{\emph{arXiv:1910.13376 [stat.ML]}}.
#'   }
#'
#' @rdname pdp2d
pdp2d <- function(model, x, vnames, type = "pdp", depth = 21, alpha = 2/3, right = 0.8, top = 0.95, digits = 1, parallel = TRUE, sample.frac = 1, pfunction = NULL, ...){
  if(sum(names(x) %in% vnames) != 2) stop("You should use 2 Variables in order to compute scatterplots!")
  
  # number of unique observations to give to processors each step
  ssize <- 20
  
  # Sampling
  x <- sample_n(x, round(nrow(x)*sample.frac))
  
  # required packages to be loaded in parallel R-threads
  # this should be extracted from variable model
  pckgs = c("randomForest")
  
  # checking for predict_function parameter
  if (is.null(pfunction)) {
    # predict_function not specified
    pfunction <- predict
  }
  
  # adding indices for unique vnames-values
  x <- x %>% group_by(x[,names(x) %in% vnames])
  x$ind <- x %>% group_indices()
  
  # ignore duplicates in xs
  xs <- x[!duplicated(x[,"ind"]), names(x) %in% append("ind", vnames)]
  
  # keep column ind for join
  xrest <- x[, !names(x) %in% append("ind", vnames), drop = FALSE]
  
  # checking if parallel computation should be used
  if (parallel) {
    # calculate number of steps of size ssize
    steps <- ceiling(nrow(xs)/ssize)
    
    # create cluster of all available cores on host
    cl <- makeCluster(detectCores())    # uses parallel::detectCores() + ::makeCluster
    
    # register cluster for foreach
    registerDoParallel(cl)
    
    # parallel loop 1 to steps, results are concatenated ('c'),
    # load pckgs in every thread
    pdps <- foreach(i = 1:steps, .combine = 'c',.packages = pckgs) %dopar% {  # uses foreach::foreach
      
      # beginning of subset
      ifrom <- (i-1)*ssize +1
      
      # end of subset
      ito <- min(i*ssize, nrow(xs))
      
      # cartesian product
      xx <- merge(xs[ifrom:ito,], xrest, by = NULL)
      
      # predict
      xx$yhat <- pfunction(model, xx, ...)
      
      # compute average prediction
      return(pred = tapply(xx$yhat, xx$ind, mean))
    }
    
    stopCluster(cl)
    
    # convert results to dataframe
    pdps <- stack(pdps)
    
    pdps$ind <- as.numeric(as.character(pdps$ind))
    
  }
  else {
    # parallel = FALSE
    pdps <- NULL
    from <- seq(1, nrow(xs), by = ssize)
    for (i in from){
      xx <- merge(xs[i:min((i+ssize-1),nrow(xs)),], xrest , by.x = NULL, by.y = NULL)
      ID <- xx[,"ind"]
      xx$yhat <- pfunction(model, xx, ...)
      pdps <- rbind(pdps, aggregate(xx$yhat, list(ID), mean))
    }
    colnames(pdps)[1] <- "ind"
    colnames(pdps)[2] <- "values"
    
  }
  
  # merge results
  pdps <- left_join(x[,"ind", drop = FALSE], pdps, by = c("ind" = "ind"))
  
  
  pred <- pfunction(model, x, ...)
  
  if(type == "both") par(mfrow = c(1,2))
  if(type != "gap"){
    minmax <- range(pdps$values)
    cols <- seq(minmax[1], minmax[2], length.out = depth)
    
    cols <- sapply(pdps$values, function(z) which.min(abs(z-cols)))
    cols <- colorRamps::blue2red(depth)[cols]
    cols <- GISTools::add.alpha(cols, alpha)
    plot(x[[vnames[1]]], x[[vnames[2]]], pch = 16, col = cols, xlab = vnames[1], ylab = vnames[2], main = "PD(x)")
    
    # legend
    vals   <- seq(minmax[1], minmax[2], length.out = depth)
    cols   <- colorRamps::blue2red(depth)
    steps  <- seq(1,depth, length.out = 5)
    xpos   <- range(xs[,2]) %*% c(1-right, right)
    ypos   <- range(xs[,1]) %*% c(1-top, top)
    legend(xpos, ypos, round(vals[steps],1),  fill = cols[steps], bty = "n", cex = 0.7)
    
  }
  if(type != "pdp"){
    diffs <- pdps$values - pred
    cols <- seq(-max(abs(range(diffs))), max(abs(range(diffs))), length.out = depth)
    cols <- sapply(diffs, function(z) which.min(abs(z-cols)))
    cols <- colorRamps::blue2red(depth)[cols]
    cols <- GISTools::add.alpha(cols, alpha)
    plot(x[[vnames[1]]], x[[vnames[2]]], pch = 16, col = cols, xlab = vnames[1], ylab = vnames[2], main = "PD(x) - f(x)")
    
    # legend
    vals   <- seq(-max(abs(range(diffs))), max(abs(range(diffs))), length.out = depth)
    cols   <- colorRamps::blue2red(depth)
    steps  <- seq(1,depth, length.out = 5)
    xpos   <- range(xs[,2]) %*% c(1-right, right)
    ypos   <- range(xs[,1]) %*% c(1-top, top)
    legend(xpos, ypos, round(vals[steps],1),  fill = cols[steps], bty = "n", cex = 0.7)
    
  }
}

#' @title Scatterplot matrix of 2D partial dependence plots
#'
#' @description Creates a scatterplot matrix of 2D partial dependence plots.
#'
#' @param model   A model with corresponding predict function that returns numeric values.
#' @param x       Data frame.
#' @param vnames  Character vector of the variable set for which the patial dependence function is to be computed.
#' @param depth   Integer specifiying the number colours in the heat map.
#' @param alpha   Numeric value for alpha blending of the points in the scatter plot.
#' @param sample.frac fraction-size for sampling of x.
#' @param pfunction   User generated predict function with arguments \code{pfunction(model, newdata)}.
#' @param ...     Further arguments to be passed to the \code{predict()} method of the \code{model}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(pdp)
#' library(randomForest)
#' data(boston)
#' set.seed(42)
#' boston.rf <- randomForest(cmedv ~ ., data = boston)
#' PDPmatrix(boston.rf, boston, vnames = c("lstat", "rm", "lon", "nox"))
#' }
#'
#' @author \email{gero.szepannek@@web.de}
#'
#' @references \itemize{
#'     \item Szepannek, G. (2019): How Much Can We See? A Note on Quantifying Explainability of Machine Learning Models,
#'           \href{https://arxiv.org/abs/1910.13376}{\emph{arXiv:1910.13376 [stat.ML]}}.
#'   }
#'
#' @rdname PDPmatrix
PDPmatrix <- function(model, x, vnames, depth = 21, alpha = 2/3, sample.frac = 1, pfunction = NULL, ...){
  n <- length(vnames)
  mat <- matrix(0,n,n)
  k <- 1
  for(i in 1:n){
    for(j in i:n){
      mat[i,j] <- k
      k <- k+1
    }
  }
  layout(mat, rep(2,n), rep(2,n))
  
  # number of unique observations to give to processors each step
  ssize <- 20
  
  # Sampling
  x <- sample_n(x, round(nrow(x)*sample.frac))
  
  # checking for predict_function parameter
  if (is.null(pfunction)) {
    # predict_function not specified
    pfunction <- predict
  }
  
  pred <- pfunction(model, x, ...)
  minmax <- range(pred)
  colrefs <- seq(minmax[1], minmax[2], length.out = depth)
  
  for(i in 1:(n-1)){
    
    plot(5, 5, type="n", axes=FALSE, ann=FALSE, xlim=c(0, 10), ylim = c(0,10))
    text(5,5, vnames[i], cex = 2)
    
    for(j in (i+1):n){
      
      x <- x %>% group_by(x[,names(x) %in% vnames[c(i,j)], drop = FALSE])
      x$ind <- x %>% group_indices()
      xs <- x[!duplicated(x[,"ind"]), names(x) %in% append("ind", vnames[c(i,j)])]
      xrest <- x[, !names(x) %in% append("ind", vnames[c(i,j)]), drop = FALSE]
      
      pdps <- NULL
      from <- seq(1, nrow(xs), by = ssize)
      for (k in from){
        xx <- merge(xs[k:min((k+ssize-1),nrow(xs)),], xrest , by.x = NULL, by.y = NULL)
        ID <- xx[,"ind"]
        xx$yhat <- pfunction(model, xx, ...)
        pdps <- rbind(pdps, aggregate(xx$yhat, list(ID), mean))
      }
      colnames(pdps)[1] <- "ind"
      colnames(pdps)[2] <- "values"
      
      pdps <- left_join(x[,"ind", drop = FALSE], pdps, by = c("ind" = "ind"))
      
      
      cols <- sapply(pdps$values, function(z) which.min(abs(z-colrefs)))
      cols <- colorRamps::blue2red(depth)[cols]
      cols <- GISTools::add.alpha(cols, alpha)
      par(mar=c(1,1,1,1))
      plot(x[[vnames[j]]], x[[vnames[i]]], pch = 16, col = cols, xlab = vnames[1], ylab = vnames[2], main = "")
    }
  }
  
  plot(5, 5, type="n", axes=FALSE, ann=FALSE, xlim=c(0, 10), ylim = c(0,10))
  text(5,5, vnames[n], cex = 2)
}