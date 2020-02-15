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
EBcorrect <- function(tree, target = "global"){
  if (nrow(tree$frame) == 1){ # initial node (i.e. naive tree) case
    EBtree <- tree
  } else {
    df <- data.frame(tree$frame$yval2[tree$frame$var == "<leaf>", 2 : 3])
    names(df) <- c("n0", "n1")
    df$n      <- rowSums(df)
    df$p0     <- df$n0/df$n
    df$p1     <- df$n1/df$n
    df$p1[df$p1 == 0] <- 10^(-3) # bug when p1==0 for fitting the beta distribution
    df$p1[df$p1 == 1] <- 1-10^(-3) # bug when p1==0 for fitting the beta distribution
    means     <- df$p1
    means2    <- c()                              # GENERATION OF
    for (i in 1 : nrow(df)){                      # THE ARTIFICIAL
      means2 <- c(means2, rep(df$p1[i], df$n[i])) # PROPORTIONS
    }                                             # SAMPLE
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
    EBtree$frame$yval2[tree$frame$var == "<leaf>", 4 : 5] <- as.matrix(df[, c("EBp0", "EBp1")]) # replacement with corrected
  }
  
  return(EBtree)
}
######################################################################################################
saveXLS <- function(finalResults){

  evalIndexes  <- names(finalResults)[!grepl("nf", names(finalResults)) & !grepl("up", names(finalResults))]
  
  # Mean results (1 result per eval index and per dataset) computation:
  meanResults <- list()
  for (evalIndex in names(finalResults)){
    meanResults1index <- c()
    for (datasetName in datasetNames){
      meanResults1index1dataset <- colMeans(finalResults[[evalIndex]][[datasetName]])
      meanResults1index         <- rbind(meanResults1index, meanResults1index1dataset)
    }
    row.names(meanResults1index) <- datasetNames
    meanResults[[evalIndex]]     <- meanResults1index
  }
  row.names(meanResults1index) <- datasetNames
  meanResults[[evalIndex]]     <- as.data.frame(meanResults1index)
  
  # Total results (on 1 dataframe) formatting:
  completeResults <- data.frame()
  for (evalIndex in evalIndexes){
    completeResults1index <- c()
    for (datasetName in datasetNames){
      df <- finalResults[[evalIndex]][[datasetName]]
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
  for (evalIndex in names(finalResults)){
    eval(parse(text = paste0(evalIndex, " <- as.data.frame(meanResults[[evalIndex]])")))
  }
  WriteXLS(c(names(finalResults), "completeResults"), xlsFileName)

}
######################################################################################################
myBoxplot <- function(df, evalIndex, datasetName){
  
  df2 <- data.frame(model = rep(names(df)[1], nrow(df)), eval = df[[names(df)[1]]])
  for (j in 2 : ncol(df)){
    df2 <- rbind(df2, data.frame(model = rep(names(df)[j], nrow(df)), eval = df[[names(df)[j]]]))
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
  df3 <- subset(df3, select = - model)
  cps <- setdiff(unique(df3$cp), "pruned")
  df3$cp <- factor(df3$cp, levels = c(cps[order(cps, decreasing = T)], "pruned"))
  
  boxplot <- ggplot(df3, aes(x=cp, y=eval, fill=type)) + 
    geom_boxplot(width=0.3) +
  theme(axis.text.x = element_text(color="#993333", size=11),
        axis.text.y = element_text(color="#993333", size=11))
  return(boxplot)
}
######################################################################################################
generate_boxplots <- function(inputs, finalResults){
  evalIndexes  <- names(finalResults)[!grepl("nf", names(finalResults)) & !grepl("up", names(finalResults))]
  datasetNames <- names(finalResults[[1]])
  models       <- names(finalResults[[1]][[1]])
  boxplots     <- list()
  evalIndex=evalIndexes[1]
  datasetName=datasetNames[1]
  
  for (evalIndex in evalIndexes){
    for (datasetName in datasetNames){
      local({
        preciseResults <- finalResults[[evalIndex]][[datasetName]]
        uncertEval <- rbind(colMeans(finalResults[[paste0(evalIndex, "Inf")]][[datasetName]]),
                            colMeans(finalResults[[paste0(evalIndex, "Sup")]][[datasetName]]))
        row.names(uncertEval) <- c("inf", "sup")
        uncertEval <- as.data.frame(uncertEval)
        p <- myBoxplot(preciseResults, evalIndex, datasetName) + ggtitle(datasetName) + labs(y = evalIndex, cex = 0.8)
        if (inputs$uncertainEval){
          str <- "p <- p"
          for (cp in c(as.character(inputs$cps), "pruned")){
            new_segt <- paste0("+ geom_segment(aes(x = '", cp, "', y = uncertEval['inf', '", 
                               cp, "'], xend = '", cp, "', yend = uncertEval['sup', '", cp, 
                               "'], linetype = 'uncertain\nevaluation'), color = '#9900FF', size=1)")
            str <- paste0(str, new_segt)
          }
          eval(parse(text = str))
        }
        p <- p + theme(legend.title=element_blank()) + scale_fill_manual(values=c("#FF3300", "#99CCFF"))
        
        boxplot_title <- paste0("boxplots_",inputs$nRuns, "x", inputs$nFold, '_foldCV_', datasetName, "_", evalIndex, "_", Sys.time(), '.png')
        boxplot_title <- gsub(' ', '_', boxplot_title)
        boxplot_title <- gsub(':', '_', boxplot_title)
        boxplot_title <- gsub('-', '_', boxplot_title)
        ggsave(filename = paste0('./results/', boxplot_title), plot = p)
  
        boxplots[[evalIndex]][[datasetName]] <<- p
      })
    }
  }
    
  return(boxplots)
}
######################################################################################################

