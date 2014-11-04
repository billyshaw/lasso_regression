This code uses lasso2's l1ce function to do <b>L1-shrinkage</b> on the ratios data

In order to choose the best bound, it implements a cross-validation with folds = 5, 
find the bound with the least mean residual, and use that bound for the best L1-shrinkage model

Function best.model returns the l1ce model, and prints out the list of best predictors for the model

This code uses GBAA0052 as example 

Made for NYU Computing with Large Data Sets (CSCI-480-001) with <a href ="http://bonneaulab.bio.nyu.edu/">Richard Bonneau</a>