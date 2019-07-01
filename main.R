####################################################################
# MAIN SCRIPT TO RUN FOR EMPIRICAL BAYES DECISION TREE EXPERIMENTS #
####################################################################

# Variables and screen cleaning:
graphics.off();cat("\014");rm(list=ls());options(warn=-1)

# Directory path tuning and libraries loading:
path <- "C:/..."
setwd(path)
library(rpart); library(rpart.plot); library(caret); library(xlsx); 
library(ggpubr); library(MASS); source("myLib.R")

# INPUTS:
datasetNames <- c("banknote", "mammo", "pima", "ticTacToe") # from c("banana", "bankLoan", "banknote", "mammo", "pima", "ticTacToe")
cps          <- c(0, 0.001, 0.005, 0.01, 0.05) # seq(from = 0, to = 1, by = 0.5)
nRuns        <- 100
nFold        <- 4
granUncEval  <- 2 # granularity of the search in the unertain predictive probability interval [p1Inf, p1Sup], i.e. number of considered points
inputs       <- list(datasetNames=datasetNames, cps=cps, nRuns=nRuns, nFold=nFold)

# RUNS:
logLossResults    <- list(); brierResults      <- list(); aucResults     <- list()
logLossInfResults <- list(); logLossSupResults <- list()
brierInfResults   <- list(); brierSupResults   <- list()
aucInfResults     <- list(); aucSupResults     <- list()
for (datasetName in datasetNames){
  logLossResults[[datasetName]]    <- data.frame(); brierResults[[datasetName]]      <- data.frame(); aucResults[[datasetName]] <- data.frame()
  logLossInfResults[[datasetName]] <- data.frame(); logLossSupResults[[datasetName]] <- data.frame()
  brierInfResults[[datasetName]]   <- data.frame(); brierSupResults[[datasetName]]   <- data.frame()
  aucInfResults[[datasetName]]     <- data.frame(); aucSupResults[[datasetName]]     <- data.frame()
}

iRun=1
for (iRun in 1 : nRuns){
  cat("\n\nrun", iRun, ": ")
  logLossResults1run    <- data.frame(); brierResults1run      <- data.frame(); aucResults1run <- data.frame()
  logLossInfResults1run <- data.frame(); logLossSupResults1run <- data.frame() 
  brierInfResults1run   <- data.frame(); brierSupResults1run   <- data.frame() 
  aucInfResults1run     <- data.frame(); aucSupResults1run     <- data.frame()

  datasetName=datasetNames[1]
  for (datasetName in datasetNames){
    cat("\n", datasetName, ": ")
    dataset              <- read.csv(paste0("./datasets/", datasetName, ".csv"), stringsAsFactors=F, sep=";")
    dataset$class        <- standardizeLabels(dataset$class)
    dataset              <- dataset[sample(nrow(dataset)), ]
    classLabels          <- unique(dataset$class)
    folds                <- createFolds(dataset$class, nFold, list = T, returnTrain = F)
    logLossResults1CV    <- data.frame(); brierResults1CV    <- data.frame(); aucResults1CV    <- data.frame()
    logLossInfResults1CV <- data.frame(); brierInfResults1CV <- data.frame(); aucInfResults1CV <- data.frame()
    logLossSupResults1CV <- data.frame(); brierSupResults1CV <- data.frame(); aucSupResults1CV <- data.frame()
    fold="Fold1"
    for (fold in names(folds)){
      cat(fold, "")
      
      # Division dataset -> trainSet + validationSet:
      validationSet <- dataset[folds[[fold]], ]
      trainSet      <- dataset[setdiff(1:nrow(dataset), folds[[fold]]), ]
      
      # Trees learning:
      trees <- learnTrees(cps, trainSet)
      
      # Probabilist predictions computing:
      probaPreds <- predictProb(trees, validationSet)
      
      # Evaluation of the probabilist predictions:
      trueLabels  <- validationSet$class
      res_logLoss <- c(); res_brier <- c(); res_auc <- c()
      res_logLossInf <- c(); res_logLossSup <- c()
      res_brierInf   <- c(); res_brierSup   <- c()
      res_aucInf     <- c(); res_aucSup     <- c()
      pred=probaPreds[[1]]
      for (pred in probaPreds){
        p1          <- pred[, 2] # probabilities of class label "1"
        logLoss     <- MLmetrics::LogLoss(p1, trueLabels)
        brier       <- DescTools::BrierScore(trueLabels, p1)
        auc         <- MLmetrics::AUC(p1, trueLabels)
        res_logLoss <- c(res_logLoss, logLoss); res_brier <- c(res_brier, brier); res_auc <- c(res_auc, auc)
      }
      # Min and max values computation for each metric:
      cp=cps[1]
      for (cp in c(as.character(cps), "pruned")){
        pred       <- probaPreds[[cp]]
        p1         <- pred[, 2]
        predEB     <- probaPreds[[paste0(cp, "EB")]]
        p1EB       <- predEB[, 2]
        p1Inf      <- pmax(p1 - abs(p1 - p1EB)/2, 0)
        p1Sup      <- pmin(p1 + abs(p1 - p1EB)/2, 1)
        logLossInf <- c(); logLossSup <- c(); brierInf <- c(); brierSup <- c(); aucInf <- c(); aucSup <- c()
        i=1
        for (i in 1 : length(p1Inf)){
          p1Inf_i    <- p1Inf[i]
          p1Sup_i    <- p1Sup[i]
          grid       <- seq(from = p1Inf_i, to = p1Sup_i, by = (p1Sup_i - p1Inf_i)/granUncEval)
          logLossVec <- c(); brierVec <- c(); aucVec <- c()
          for (p in grid){
            logLossVec <- c(logLossVec, MLmetrics::LogLoss(p, trueLabels[i]))
            brierVec   <- c(brierVec,   DescTools::BrierScore(trueLabels[i], p))
          }
          logLossInf <- c(logLossInf, min(logLossVec)); logLossSup <- c(logLossSup, max(logLossVec))
          brierInf   <- c(brierInf, min(brierVec));     brierSup   <- c(brierSup,   max(brierVec))
        }
        aucInf <- MLmetrics::AUC(p1Inf, trueLabels); aucSup <- MLmetrics::AUC(p1Sup, trueLabels) # AUC has to be computed on vectors
        
        # In case of decreasing evaluation functions:
        minLogLoss <- min(mean(logLossInf), mean(logLossSup)); maxLogLoss <- max(mean(logLossInf), mean(logLossSup))
        minBrier   <- min(mean(brierInf),  mean(brierSup));    maxBrier   <- max(mean(brierInf), mean(brierSup))
        minAUC     <- min(aucInf, aucSup);                     maxAUC     <- max(aucInf, aucSup)
        
        res_logLossInf <- c(res_logLossInf, minLogLoss); res_logLossSup <- c(res_logLossSup, maxLogLoss)
        res_brierInf   <- c(res_brierInf, minBrier);     res_brierSup   <- c(res_brierSup, maxBrier)
        res_aucInf     <- c(res_aucInf, minAUC);         res_aucSup     <- c(res_aucSup, maxAUC)
      }
      logLossResults1CV <- rbind(logLossResults1CV, res_logLoss)
      brierResults1CV   <- rbind(brierResults1CV,   res_brier)
      aucResults1CV     <- rbind(aucResults1CV,     res_auc)
      
      logLossInfResults1CV <- rbind(logLossInfResults1CV, res_logLossInf); logLossSupResults1CV <- rbind(logLossSupResults1CV, res_logLossSup)
      brierInfResults1CV   <- rbind(brierInfResults1CV,   res_brierInf);   brierSupResults1CV   <- rbind(brierSupResults1CV,   res_brierSup)
      aucInfResults1CV     <- rbind(aucInfResults1CV,     res_aucInf);     aucSupResults1CV     <- rbind(aucSupResults1CV,     res_aucSup)

    }
    names(logLossResults1CV)        <- names(trees)
    names(brierResults1CV)          <- names(trees)
    names(brierResults1CV)          <- names(trees)
    row.names(logLossResults1CV)    <- names(folds); row.names(brierResults1CV) <- names(folds); row.names(aucResults1CV) <- names(folds)
    names(logLossInfResults1CV)     <- c(as.character(cps), "pruned") ; names(logLossSupResults1CV) <- c(as.character(cps), "pruned")
    names(brierInfResults1CV)       <- c(as.character(cps), "pruned") ; names(brierSupResults1CV)   <- c(as.character(cps), "pruned")
    names(aucInfResults1CV)         <- c(as.character(cps), "pruned") ; names(aucSupResults1CV)     <- c(as.character(cps), "pruned")
    row.names(logLossInfResults1CV) <- names(folds); row.names(logLossSupResults1CV) <- names(folds); 
    row.names(brierInfResults1CV)   <- names(folds); row.names(brierSupResults1CV)   <- names(folds)
    row.names(aucInfResults1CV)     <- names(folds); row.names(aucSupResults1CV)     <- names(folds)
    logLossResults1run              <- rbind(logLossResults1run, colMeans(logLossResults1CV))
    brierResults1run                <- rbind(brierResults1run,   colMeans(brierResults1CV))
    aucResults1run                  <- rbind(aucResults1run,     colMeans(aucResults1CV))
    
    logLossInfResults1run <- rbind(logLossInfResults1run, colMeans(logLossInfResults1CV))
    logLossSupResults1run <- rbind(logLossSupResults1run, colMeans(logLossSupResults1CV))
    brierInfResults1run   <- rbind(brierInfResults1run,   colMeans(brierInfResults1CV))
    brierSupResults1run   <- rbind(brierSupResults1run,   colMeans(brierSupResults1CV))
    aucInfResults1run     <- rbind(aucInfResults1run,     colMeans(aucInfResults1CV))
    aucSupResults1run     <- rbind(aucSupResults1run,     colMeans(aucSupResults1CV))
    
    logLossResults[[datasetName]]   <- rbind(logLossResults[[datasetName]], colMeans(logLossResults1CV))
    brierResults[[datasetName]]     <- rbind(brierResults[[datasetName]],   colMeans(brierResults1CV))
    aucResults[[datasetName]]       <- rbind(aucResults[[datasetName]],     colMeans(aucResults1run))
    
    logLossInfResults[[datasetName]] <- rbind(logLossInfResults[[datasetName]], colMeans(logLossInfResults1CV))
    logLossSupResults[[datasetName]] <- rbind(logLossSupResults[[datasetName]], colMeans(logLossSupResults1CV))
    brierInfResults[[datasetName]]   <- rbind(brierInfResults[[datasetName]],   colMeans(brierInfResults1CV))
    brierSupResults[[datasetName]]   <- rbind(brierSupResults[[datasetName]],   colMeans(brierSupResults1CV))
    aucInfResults[[datasetName]]     <- rbind(aucInfResults[[datasetName]],     colMeans(aucInfResults1CV))
    aucSupResults[[datasetName]]     <- rbind(aucSupResults[[datasetName]],     colMeans(aucSupResults1CV))
  }
  names(logLossResults1run)        <- names(trees); names(logLossResults1run) <- names(trees); names(aucResults1run) <- names(trees)
  names(logLossInfResults1run)     <- c(as.character(cps), "pruned"); names(logLossSupResults1run) <- c(as.character(cps), "pruned")
  names(brierInfResults1run)       <- c(as.character(cps), "pruned"); names(brierSupResults1run)   <- c(as.character(cps), "pruned")
  names(aucInfResults1run)         <- c(as.character(cps), "pruned"); names(aucSupResults1run)     <- c(as.character(cps), "pruned")
  row.names(logLossResults1run)    <- datasetNames; row.names(brierResults1run) <- datasetNames; row.names(aucResults1run) <- datasetNames
  row.names(logLossInfResults1run) <- datasetNames; row.names(logLossSupResults1run) <- datasetNames
  row.names(brierInfResults1run)   <- datasetNames; row.names(brierSupResults1run)   <- datasetNames
  row.names(aucInfResults1run)     <- datasetNames; row.names(aucSupResults1run)     <- datasetNames
}
for (datasetName in datasetNames){
  names(logLossResults[[datasetName]])        <- names(trees)
  names(logLossInfResults[[datasetName]])     <- c(as.character(cps), "pruned")
  names(logLossSupResults[[datasetName]])     <- c(as.character(cps), "pruned")
  names(brierResults[[datasetName]])          <- names(trees)
  names(brierInfResults[[datasetName]])       <- c(as.character(cps), "pruned")
  names(brierSupResults[[datasetName]])       <- c(as.character(cps), "pruned")
  names(aucResults[[datasetName]])            <- names(trees)
  names(aucInfResults[[datasetName]])         <- c(as.character(cps), "pruned")
  names(aucSupResults[[datasetName]])         <- c(as.character(cps), "pruned")
  row.names(logLossResults[[datasetName]])    <- paste("run", 1:nRuns)
  row.names(logLossInfResults[[datasetName]]) <- paste("run", 1:nRuns)
  row.names(logLossSupResults[[datasetName]]) <- paste("run", 1:nRuns)
  row.names(brierResults[[datasetName]])      <- paste("run", 1:nRuns)
  row.names(brierInfResults[[datasetName]])   <- paste("run", 1:nRuns)
  row.names(brierSupResults[[datasetName]])   <- paste("run", 1:nRuns)
  row.names(aucResults[[datasetName]])        <- paste("run", 1:nRuns)
  row.names(aucInfResults[[datasetName]])     <- paste("run", 1:nRuns)
  row.names(aucSupResults[[datasetName]])     <- paste("run", 1:nRuns)
}
finalResults <- list(logLoss=logLossResults, logLossInf=logLossInfResults, logLossSup=logLossSupResults, 
                     brier=brierResults,     brierInf=brierInfResults,     brierSup=brierSupResults, 
                     AUC=aucResults,         AUCInf=aucInfResults,         AUCSup=aucSupResults)

boxplots <- plotResults(inputs, finalResults)
ggarrange(boxplots$logLoss$banana, boxplots$logLoss$bankLoan, boxplots$logLoss$banknote, 
          boxplots$logLoss$mammo, boxplots$logLoss$pima, boxplots$logLoss$ticTacToe, ncol = 3, nrow = 2)

# p+ geom_errorbar(aes(ymin = c(0.35, 0.35, 0.25, 0.25, 0.33, 0.33, 0.31, 0.31, 0.37, 0.37, 0.37, 0.37, 0.32, 0.32, 0.31, 0.31),
#                      ymax = c(0.45, 0.45, 0.35, 0.35, 0.43, 0.43, 0.41, 0.41, 0.47, 0.47, 0.47, 0.47, 0.42, 0.42, 0.41, 0.41),
#                  width = 0.2))



