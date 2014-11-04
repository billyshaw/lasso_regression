# This code uses lasso2's l1ce function to do L1-shrinkage on the ratios data
#
# In order to choose the best bound, it implements a cross-validation with folds = 5, 
# find the bound with the least mean residual, and use that bound for the best L1-shrinkage model
#
# Function best.model returns the l1ce model, and prints out the list of best predictors for the model
#
# This code uses GBAA0052 as example 


library("lasso2")
library("glmnet")

# load ratios table
load("baa.ratios.rda")
ratios.rows <- row.names(ratios)

# load and find predictors list
predictors.table <- read.table("baa.TFs.txt")
predictors.list <- as.character(predictors.table[,1]) # list of predictors

# function to split subsets for cross validation
cv.folds <- function(n, folds = 10) {
  split(sample(1:n), rep(1:folds, length = n))
}

# function to predict for cross validation
predict.from.l1ce <- function( l1ce.model, x) {
  x <- as.matrix(x)
  if (class(l1ce.model) != "l1ce") {
    stop("input class is not l1ce")
    return(FALSE)
  }
  coeff <- l1ce.model$coefficients
  n <- dim(x)[1]
  tmp.mat <- t(cbind(rep(1,n), x)) * coeff
  y.hat <- apply(tmp.mat, 2, sum)
  invisible(y.hat)
}

best.model <- function (ratios, predictors.list, gene) {
  
  # data validation: check if all predictors in list exist
  predictors.list <- predictors.list[predictors.list%in%rownames(ratios)]
  
  # data validation: check if gene is own predictor
  if (typeof(gene) == "character") {
    predictors.list <- predictors.list[!predictors.list==gene]
  } 
   
  # create a subset of ratios with only predictors
  predictors.ratios <- t(ratios[predictors.list,])
  
  # find the gene info passed into the function
  row <- ratios[gene,] 
  
  bounds <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cv.bounds <- NULL
  
  for (b in 1:length(bounds)) {
    
    # cross validation with folds = 5
    cv.subsets <- cv.folds(nrow(predictors.ratios), folds = 5)
    cv.rss <- numeric(5)
    
    l1ce.model <- l1ce(row ~ predictors.ratios, sweep.out=~1, standardize = TRUE, absolute.t = FALSE)
    
    # loop through each subset
    for (i in 1:5) {
      y.hat <- predict.from.l1ce(l1ce.model, predictors.ratios[cv.subsets[[i]],] )
      # print (row[cv.subsets[[i]]] - y.hat)
      cv.rss[i] <- mean( (row[cv.subsets[[i]]] - y.hat ) ** 2)
    }
    
    cv = mean(cv.rss)
    cv.bounds <- rbind(cv.bounds, cv)
    
  }
  
  best.idx <- match(min(cv.bounds), cv.bounds)
  best.bound <- bounds[best.idx]
  
  best.model <- l1ce.model <- l1ce(row ~ predictors.ratios, sweep.out=~1, standardize = TRUE, absolute.t = FALSE, bound = best.bound)
  
  cat("The best model for GBAA0052 is L1-shrinkage (l1ce) model with the bound")
  print(best.bound)
  cat("\nThe best predictors for the model are as follows\n ")
  best.models.list <- which(abs(best.model$coefficients) > 0)
  print(best.models.list)
  
  return(best.model)
  
}

best.model(ratios, predictors.list, "GBAA0052")
