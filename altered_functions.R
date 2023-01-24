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
  pdf(file = file.path(pathFolder, 'model_selection.pdf'))
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
