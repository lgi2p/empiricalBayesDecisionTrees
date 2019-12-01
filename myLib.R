#####################
# Functions library #
#####################

######################################################################################################
pruneRpart <- function(tree){
  opt        <- which.min(tree$cptable[,"xerror"])
  cpOpt      <- tree$cptable[opt, "CP"]
  prunedTree <- prune(tree, cp = cpOpt)
  return(prunedTree)
}
######################################################################################################
standardizeLabels <- function(classLabels) {
  labels                                <- unique(classLabels)
  classLabels[classLabels == labels[1]] <- "aa" # uniformisation
  classLabels[classLabels == labels[2]] <- "bb" # of
  classLabels[classLabels == "aa"]      <- 1    # labels
  classLabels[classLabels == "bb"]      <- 0    #
  classLabels                           <- as.numeric(classLabels)
}
######################################################################################################
learnTrees <- function(cps, trainSet){
  
  learntTrees    <- list()
  for (cp in cps){
    tree   <- rpart(class ~ ., data = trainSet, method = "class", cp = cp)
    treeEB <- EBcorrect(tree)
    learntTrees[[as.character(cp)]]               <- tree
    learntTrees[[paste0(as.character(cp), "EB")]] <- treeEB
  }
  
  largeTree          <- rpart(class ~ ., data = trainSet, cp = 0, method = "class")
  prunedTree         <- pruneRpart(largeTree)
  learntTrees$pruned <- prunedTree
  
  largeTree            <- rpart(class ~ ., data = trainSet, cp = 0, method = "class")
  largeEBTree          <- EBcorrect(largeTree)
  prunedTree           <- pruneRpart(largeTree)
  learntTrees$prunedEB <- prunedTree
  
  return(learntTrees)
}
######################################################################################################
predictProb <- function(trees, validationSet){
  predictions <- list()
  for (modelName in names(trees)){
    predictions[[modelName]] <- predict(trees[[modelName]], validationSet)
  }
  return(predictions)
}
######################################################################################################
plotResults <- function(inputs, results){
  
  evalIndexes  <- names(results)[!grepl("nf", names(results)) & !grepl("up", names(results))]
  # uncertEvalIndexes  <- names(results)[grepl("nf", names(results))|grepl("up", names(results))]
  datasetNames <- names(results[[1]])
  models       <- names(results[[1]][[1]])
  boxplots     <- list()
  for (evalIndex in evalIndexes){
    for (datasetName in datasetNames){
      preciseResults <- results[[evalIndex]][[datasetName]]
      uncertEval <- rbind(colMeans(results[[paste0(evalIndex, "Inf")]][[datasetName]]),
                          colMeans(results[[paste0(evalIndex, "Sup")]][[datasetName]]))
      row.names(uncertEval) <- c("inf", "sup")
      uncertEval <- as.data.frame(uncertEval)
      p <- myBoxplot(preciseResults, evalIndex, datasetName) + ggtitle(datasetName) + labs(x="", y = "", cex=0.8)+
      labs(x="", y = "", cex=0.8) 
      if (inputs$uncertainEval){
        for (model in names(uncertEval)){
          p <- p + geom_segment(aes(x = model, y = uncertEval['inf', model],
                                    xend = model, yend = uncertEval['sup', model]), color = 'blue')
        }
        p <- p +
          geom_segment(aes(x = "0", y = uncertEval['inf', "0"],
                           xend = "0", yend = uncertEval['sup', "0"]), color = 'blue') +
          geom_segment(aes(x = "0.005", y = uncertEval['inf', "0.005"],
                           xend = "0.005", yend = uncertEval['sup', "0.005"]), color = 'blue') +
          geom_segment(aes(x = "0.01", y = uncertEval['inf', "0.01"],
                           xend = "0.01", yend = uncertEval['sup', "0.01"]), color = 'blue') +
          geom_segment(aes(x = "pruned", y = uncertEval['inf', "pruned"],
                           xend = "pruned", yend = uncertEval['sup', "pruned"]), color = 'blue')
      }
      boxplots[[evalIndex]][[datasetName]] <- p
    }
  }
  
  for (iEvalIndex in 1 : length(evalIndexes)){
    evalIndex <- evalIndexes[iEvalIndex]
    str <- paste0("boxplots", evalIndex, " <- ggarrange(")
    for (iDataset in 1 : length(datasetNames)){
      str <- paste0(str, "boxplots[[", iEvalIndex, "]][[", iDataset, "]],")
    }
    if (length(datasetNames) <= 2){
      str <- paste0(str, "ncol=2, nrow = 1)")
    } else {
      if (length(datasetNames) <= 4){
        str <- paste0(str, "ncol=2, nrow = 2)")
      } else {
        if (length(datasetNames) <= 6){
          str <- paste0(str, "ncol=2, nrow = 3)")
        } else if (length(datasetNames) <= 9){
          str <- paste0(str, "ncol=3, nrow = 3)")
        }
      }
    }
    str <- paste0(str, "; boxplots", evalIndex, "<- annotate_figure(boxplots", evalIndex, ", top = text_grob('", evalIndex,
                  "', color = 'blueviolet', face = 'bold', size = 20))")
    
    eval(parse(text = str))
    # eval(parse(text = str2))
    eval(parse(text = paste0("ggsave(paste0('./results/boxplots", evalIndex, 
                             "', gsub(':', '_', paste0(Sys.time())), '.png'), boxplots",
                             evalIndex, ", width = 23, height = 10, units = 'cm')")))
  }
  
  # Mean results (1 result per eval index and per dataset) computation:
  meanResults <- list()
  for (evalIndex in names(results)){
    meanResults1index <- c()
    for (datasetName in datasetNames){
      meanResults1index1dataset <- colMeans(results[[evalIndex]][[datasetName]])
      meanResults1index         <- rbind(meanResults1index, meanResults1index1dataset)
    }
    row.names(meanResults1index) <- datasetNames
    meanResults[[evalIndex]]     <- meanResults1index
  }
  row.names(meanResults1index) <- datasetNames
  meanResults[[evalIndex]]     <- meanResults1index
  
  # Total results (on 1 dataframe) formatting:
  completeResults <- c()
  for (evalIndex in evalIndexes){
    completeResults1index <- c()
    for (datasetName in datasetNames){
      df <- results[[evalIndex]][[datasetName]]
      completeResults1index1dataset <- c()
      for (iCol in 1 : ncol(df)){
        df2 <- data.frame(model = rep(names(df)[iCol], nrow(df)), dataset = datasetName,
                          index = evalIndex, value = df[, iCol])
        completeResults1index1dataset <- rbind(completeResults1index1dataset, df2)
      }
      completeResults1index <- rbind(completeResults1index, completeResults1index1dataset)
    }
    completeResults <- rbind(completeResults, completeResults1index)
  }
  
  # Xls file (containing 1 sheet per model):
  xlsFileName <- paste0("results/results ", gsub(":", "_", paste0(Sys.time())), ".xlsx")
  for (evalIndex in names(results)){
    write.xlsx(meanResults[[evalIndex]], file = xlsFileName,
               sheetName = paste0("MEAN ", evalIndex), append = T)
  }
  write.xlsx(completeResults, file = xlsFileName, sheetName = "Total results", append = T, row.names=F)
  
  return(boxplots)
}
######################################################################################################
resSep <- function(res){
  
  namesDT        <- setdiff(names(res), names(res)[grepl("EB", names(res))])
  resDT          <- res[, namesDT]
  resDT          <- data.frame(cbind('dt', t(colMeans(resDT))))
  names(resDT)   <- c("type", namesDT)
  resEBdt        <- res[, grepl("EB", names(res))]
  resEBdt        <- data.frame(cbind('EBdt', t(colMeans(resEBdt))))
  names(resEBdt) <- c("type", namesDT)
  res2           <- rbind(resDT, resEBdt)
  return(res2)
}
######################################################################################################
myBoxplot <- function(df, evalIndex, datasetName){
  
  df2 <- data.frame(model=rep(names(df)[1], nrow(df)), eval = df[[names(df)[1]]])
  for (j in 2 : ncol(df)){
    df2 <- rbind(df2, data.frame(model=rep(names(df)[j], nrow(df)), eval = df[[names(df)[j]]]))
  }
  
  df3      <- df2
  df3$cp   <- unlist(strsplit(as.character(df3$model), "EB"))
  df3$type <- NA
  for (i in 1 : nrow(df3)){
    if (grepl("EB", as.character(df3$model[i]))){
      df3$type[i] <- "EB"
    } else {
      df3$type[i] <- "original"
    }
  }
  # df3$type <- c("original", "EB")[grep("EB", as.character(df3$model))]
  df3 <- subset(df3, select = - model)
  
  boxplot <- ggplot(df3, aes(x=cp, y=eval, color=type)) + 
    # geom_violin() + 
    geom_boxplot(width=0.1) +
    # ggtitle(datasetName) +
    # labs(x="", y = "aaa", cex=0.8) +
    bgcolor("#BFD5E3") +
    theme(plot.title = element_text(color="#0072B2", size=15),
          plot.subtitle = element_text(color="#56B4E9", size=7),
          axis.text.x = element_text(color="#993333", size=11),
          axis.text.y = element_text(color="#993333", size=11))
  return(boxplot)
}
######################################################################################################
EBcorrect <- function(tree, target = "global"){
  if (nrow(tree$frame) == 1){ # initial node (i.e. naive tree) case
    EBtree <- tree
  } else {
    df <- data.frame(tree$frame$yval2[tree$frame$var == "<leaf>", 2 : 3])
    names(df) <- c("n0", "n1")
    df$n      <- rowSums(df)
    df$p0     <- df$n0/df$n
    df$p1     <- df$n1/df$n
    df$p1[df$p1 == 0] <- 10^(-3) # bug when p1==0 fpr fitting the beta distribution
    df$p1[df$p1 == 1] <- 1-10^(-3) # bug when p1==0 fpr fitting the beta distribution
    means     <- df$p1
    means2    <- c()                              # generation of the
    for (i in 1 : nrow(df)){                      # corresponding
      means2 <- c(means2, rep(df$p1[i], df$n[i])) # proportions
    }                                             # sample.
    if (target == "global"){
      model <- fitdistr(means2, dbeta, start = list(shape1 = 1, shape2 = 10))
    } else if (target == "onLeaves"){
      model <- fitdistr(means, dbeta, start = list(shape1 = 1, shape2 = 10))
    }
    alpha0  <- model$estimate[1]                                          
    beta0   <- model$estimate[2]                                             
    EBmeans <- (df$n1 + alpha0)/(df$n + alpha0 + beta0)
    df$EBp0 <- 1 - EBmeans 
    df$EBp1 <- EBmeans
    EBtree  <- tree
    EBtree$frame$yval2[tree$frame$var == "<leaf>", 4 : 5] <- as.matrix(df[, c("EBp0", "EBp1")])
  }
  
  return(EBtree)
}
######################################################################################################

