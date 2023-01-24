#Functions from the cisTopic package (https://github.com/aertslab/cisTopic) modified
#to run without connection to X11 display

library(plyr)

selectModelModified <- function (object, pathFolder, select = NULL, type = "derivative", keepBinaryMatrix = TRUE,
          keepModels = TRUE, ...)
{
  if (as.vector(class(object)) == "cisTopic") {
    models <- object@models
    if (length(models) < 1) {
      stop("Please, run runCGSModels() or runWarpLDAModels() first.")
    }
  }
  else if (is.list(object)) {
    models <- object
  }
  else {
    stop("Incorrect input. Is it a cisTopic object or a list of models?")
  }
  if (length(models) < 2 && type == "derivative") {
    stop("You need at least 3 models to use the derivative method.")
  }
  if (is.null(object@calc.params[["runWarpModels"]]) && type ==
      "derivative") {
    print("Are these CGS models? Please, use type=\"maximum\"")
  }
  loglikelihood <- sapply(seq_along(models), FUN = function(i) models[[i]]$log.likelihood[2,
                                                                                          ncol(models[[i]]$log.likelihood)])
  topics <- sapply(seq_along(models), FUN = function(i) nrow(models[[i]]$topics))
  object.log.lik <- data.frame(topics = topics, LL = loglikelihood)
  object.log.lik <- object.log.lik[order(object.log.lik$topics),
  ]
  ll <- object.log.lik$LL
  if (length(models) > 2) {
    object.log.lik$first_derivative <- c(-Inf, (diff(ll)/diff(topics)))
    object.log.lik$second_derivative <- c(-Inf, -Inf, diff(object.log.lik$first_derivative)[-1]/diff(topics[-1]))
    if (!is.null(models[[1]]$perplexity)) {
      object.log.lik$perplexity <- sapply(seq_along(models),
                                          FUN = function(i) models[[i]]$perplexity)
    }
    else {
      if (type == "perplexity") {
        stop("The perplexity option is only available for WarpLDA models.")
      }
    }
  }
  pdf(file = file.path(pathFolder, 'modelSelection.pdf'))
  par(bty = "n")
  plot(object.log.lik$topics, object.log.lik$LL, xlab = "Number of topics",
       ylab = "log P(D|M,T)", type = "o", pch = 16, col = "black",
       main = "Model selection")
  if (is.null(select)) {
    if (type == "maximum") {
      points(object.log.lik$topics[which(object.log.lik$LL ==
                                           max(object.log.lik$LL))], max(object.log.lik$LL),
             pch = 4, col = "red", lwd = 7)
      title(sub = paste("Best model:", object.log.lik$topics[which(object.log.lik$LL ==
                                                                     max(object.log.lik$LL))], "topics"))
       
      selected.model <- models[[which(object.log.lik$LL ==
                                        max(object.log.lik$LL))]]
    }
    else if (type == "derivative" && length(models) > 2) {
      points(object.log.lik$topics[which(object.log.lik$second_derivative ==
                                           max(object.log.lik$second_derivative))], object.log.lik$LL[which(object.log.lik$second_derivative ==
                                                                                                              max(object.log.lik$second_derivative))], pch = 4,
             col = "red", lwd = 7)
      title(sub = paste("Best model:", object.log.lik$topics[which(object.log.lik$second_derivative ==
                                                                     max(object.log.lik$second_derivative))], "topics"))
       
      selected.model <- models[[which(object.log.lik$second_derivative ==
                                        max(object.log.lik$second_derivative))]]
    }
    else if (type == "perplexity" & !is.null(models[[1]]$perplexity)) {
      points(object.log.lik$topics[which(object.log.lik$perplexity ==
                                           min(object.log.lik$perplexity))], object.log.lik$LL[which(object.log.lik$perplexity ==
                                                                                                       min(object.log.lik$perplexity))], pch = 4, col = "red",
             lwd = 7)
      title(sub = paste("Best model:", object.log.lik$topics[which(object.log.lik$perplexity ==
                                                                     min(object.log.lik$perplexity))], "topics"))
       
      selected.model <- models[[which(object.log.lik$perplexity ==
                                        min(object.log.lik$perplexity))]]
    }
  }
  else {
    points(object.log.lik$topics[which(object.log.lik$topics ==
                                         select)], object.log.lik$LL[which(object.log.lik$topics ==
                                                                             select)], pch = 4, col = "red", lwd = 7)
    title(sub = paste("Selected model:", object.log.lik$topics[which(object.log.lik$topics ==
                                                                       select)], "topics"))
     
    selected.model <- models[[which(object.log.lik$topics ==
                                      select)]]
  }
  if (length(models) > 2) {
    plot(object.log.lik$topics[-c(1, 2)], object.log.lik$second_derivative[-c(1,
                                                                              2)], xlab = "Number of topics", ylab = "Second derivative on the log-likelihood curve",
         type = "o", pch = 16, col = "black", main = "Model selection")
    if (is.null(select)) {
      if (type == "maximum") {
        points(object.log.lik$topics[which(object.log.lik$LL ==
                                             max(object.log.lik$LL))], object.log.lik$second_derivative[which(object.log.lik$LL ==
                                                                                                                max(object.log.lik$LL))], pch = 4, col = "red",
               lwd = 7)
        title(sub = paste("Best model:", object.log.lik$topics[which(object.log.lik$LL ==
                                                                       max(object.log.lik$LL))], "topics"))
         
      }
      else if (type == "derivative") {
        points(object.log.lik$topics[which(object.log.lik$second_derivative ==
                                             max(object.log.lik$second_derivative))], max(object.log.lik$second_derivative),
               pch = 4, col = "red", lwd = 7)
        title(sub = paste("Best model:", object.log.lik$topics[which(object.log.lik$second_derivative ==
                                                                       max(object.log.lik$second_derivative))], "topics"))
         
      }
      else if (type == "perplexity" & !is.null(models[[1]]$perplexity)) {
        points(object.log.lik$topics[which(object.log.lik$perplexity ==
                                             min(object.log.lik$perplexity))], object.log.lik$second_derivative[which(object.log.lik$perplexity ==
                                                                                                                        min(object.log.lik$perplexity))], pch = 4,
               col = "red", lwd = 7)
        title(sub = paste("Best model:", object.log.lik$topics[which(object.log.lik$perplexity ==
                                                                       min(object.log.lik$perplexity))], "topics"))
         
      }
    }
    else {
      points(object.log.lik$topics[which(object.log.lik$topics ==
                                           select)], object.log.lik$second_derivative[which(object.log.lik$topics ==
                                                                                              select)], pch = 4, col = "red", lwd = 7)
      title(sub = paste("Selected model:", object.log.lik$topics[which(object.log.lik$topics ==
                                                                         select)], "topics"))
       
      selected.model <- models[[which(object.log.lik$topics ==
                                        select)]]
    }
  }
  if (!is.null(models[[1]]$perplexity)) {
    plot(object.log.lik$topics, object.log.lik$perplexity,
         xlab = "Number of topics", ylab = "Perplexity", type = "o",
         pch = 16, col = "black", main = "Model selection")
    if (is.null(select)) {
      if (type == "maximum") {
        points(object.log.lik$topics[which(object.log.lik$LL ==
                                             max(object.log.lik$LL))], object.log.lik$perplexity[which(object.log.lik$LL ==
                                                                                                         max(object.log.lik$LL))], pch = 4, col = "red",
               lwd = 7)
        title(sub = paste("Best model:", object.log.lik$topics[which(object.log.lik$LL ==
                                                                       max(object.log.lik$LL))], "topics"))
         
      }
      else if (type == "derivative" && length(models) >
               2) {
        points(object.log.lik$topics[which(object.log.lik$second_derivative ==
                                             max(object.log.lik$second_derivative))], object.log.lik$perplexity[which(object.log.lik$second_derivative ==
                                                                                                                        max(object.log.lik$second_derivative))], pch = 4,
               col = "red", lwd = 7)
        title(sub = paste("Best model:", object.log.lik$topics[which(object.log.lik$second_derivative ==
                                                                       max(object.log.lik$second_derivative))], "topics"))
          
      }
      else if (type == "perplexity") {
        points(object.log.lik$topics[which(object.log.lik$perplexity ==
                                             min(object.log.lik$perplexity))], min(object.log.lik$perplexity),
               pch = 4, col = "red", lwd = 7)
        title(sub = paste("Best model:", object.log.lik$topics[which(object.log.lik$perplexity ==
                                                                       min(object.log.lik$perplexity))], "topics"))
         
        selected.model <- models[[which(object.log.lik$perplexity ==
                                          min(object.log.lik$perplexity))]]
      }
    }
    else {
      points(object.log.lik$topics[which(object.log.lik$topics ==
                                           select)], object.log.lik$perplexity[which(object.log.lik$topics ==
                                                                                       select)], pch = 4, col = "red", lwd = 7)
      title(sub = paste("Selected model:", object.log.lik$topics[which(object.log.lik$topics ==
                                                                         select)], "topics"))
       
      selected.model <- models[[which(object.log.lik$topics ==
                                        select)]]
    }
  }
  dev.off()
  if (as.vector(class(object)) == "cisTopic") {
    object@log.lik <- object.log.lik
    object@selected.model <- selected.model
    if (keepBinaryMatrix != TRUE) {
      object@binary.count.matrix <- NULL
    }
    if (keepModels != TRUE) {
      object@models <- list()
    }
    return(object)
  }
  else if (is.list(object)) {
    return(selected.model)
  }
}

binarizecisTopicsModified <- function (object, pathFolder, method = "GammaFit", thrP = 0.99, plot = FALSE,
          cutoffs = NULL)
{
  scores <- .getScores(object)
  object.binarized.cisTopics <- list()
  pdf(file = file.path(pathFolder, 'binarizedCisTopics.pdf'))
  par.opts <- par()
  if (method == "GammaFit") {
    if (!"fitdistrplus" %in% installed.packages()) {
      stop("Please, install fitdistrplus: \n install.packages(\"fitdistrplus\")")
    }
    else {
      require(fitdistrplus)
    }
    if (is.null(thrP)) {
      stop("If GammaFit is selected as method, a probability threshold for cutoff in the distribution must be provided.")
    }
    for (i in 1:ncol(scores)) {
      distr <- suppressWarnings(fitdist(scores[, i], "gamma",
                                        method = "mme"))
      cutoff <- as.numeric(unlist(quantile(distr, probs = thrP))[1])
      if (plot) {
       
        hist(scores[, i], breaks = 100, prob = TRUE,
             main = colnames(scores)[i], col = adjustcolor("dodgerblue",
                                                                       alpha.f = 0.8), xlab = "Score")
        curve(dgamma(x, rate = as.numeric(unlist(distr[1]$estimate[2])),
                     shape = as.numeric(unlist(distr[1]$estimate[1]))),
              add = TRUE, col = "magenta", lwd = 2)
        abline(v = cutoff, lty = 3, lwd = 2, col = "grey")
        title(sub = paste("Regions with score > ", signif(cutoff,
                                                          2), sep = ""))
        
      }
      object.binarized.cisTopics[[colnames(scores)[i]]] <- scores[which(scores[,
                                                                               i] > cutoff), i, drop = FALSE]
    }
  }
  if (method == "Predefined") {
    if (is.null(cutoffs)) {
      stop("If Predefined is selected as method, the threshold cutoffs must be provided.")
    }
    if (length(cutoffs) == 1) {
      cutoffs <- rep(cutoffs, ncol(scores))
    }
    else if (length(cutoffs) != ncol(scores)) {
      stop("Thresholds for all topics are not provided. Check the length of the cutoff vector.")
    }
    for (i in 1:length(cutoffs)) {
      ranking <- scores[order(scores[, i], decreasing = TRUE),
                        i]
      cutoffs[i] <- ranking[cutoffs[i] + 1]
    }
    if (plot == TRUE) {
      
      
      
      for (i in 1:ncol(scores)) {
        par(mfrow = c(1, 1))
        par(las = 1)
        hist(scores[, i], breaks = 100, prob = TRUE,
             main = colnames(scores)[i], col = adjustcolor("dodgerblue",
                                                                       alpha.f = 0.8), xlab = "Score")
        abline(v = cutoffs[i], lty = 3, lwd = 2, col = "grey")
        title(sub = paste("Regions with score > ", signif(cutoffs[i],
                                                          2), sep = ""))
      }
      
    }
    object.binarized.cisTopics <- llply(1:ncol(scores), function(i) scores[which(scores[,
                                                                                        i] > cutoffs[i]), i, drop = FALSE])
  }
  object.binarized.cisTopics <- llply(1:ncol(scores), function(i) object.binarized.cisTopics[[i]][order(object.binarized.cisTopics[[i]][,
                                                                                                                                        1], decreasing = TRUE), , drop = FALSE])
  names(object.binarized.cisTopics) <- colnames(scores)
  object@binarized.cisTopics <- object.binarized.cisTopics
  
  
  par(mfrow = c(1, 1))
  par(las = 2)
  barplot(sapply(object.binarized.cisTopics, nrow), col = adjustcolor("dodgerblue",
                                                                                  alpha.f = 0.8), main = "Number of regions selected per topic")
  
  suppressWarnings(par(par.opts))
  dev.off()
  
  return(object)
}

#Helper Function
.getScores <- function(
    object
){
  scores <- object@region.data[, grep('Scores_Topic', colnames(object@region.data))]
  
  # Check info
  if (ncol(scores) < 1){
    stop('Please, run getRegionsScores() first.')
  }
  
  colnames(scores) <- paste('Topic', 1:ncol(scores), sep='')
  return(scores)
}