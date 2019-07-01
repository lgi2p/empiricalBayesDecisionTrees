######################
# EB pruning example #
######################

# Variables and screen cleaning:
graphics.off();cat("\014");rm(list=ls());options(warn=-1)

# Directory path and libraries tuning:
path <- "C:/..."
setwd(path)
source("myLib.R"); library(MASS); library(rpart); library(rpart.plot)

# 1. Empirical Bayes (EB) simple example:
df        <- data.frame(n1=c(2, 20, 200, 40), n=c(5, 40, 500, 50)) # generation of a decision 
means     <- df$n1/df$n                                            # tree's class effective
means2    <- c()                              # generation of the
for (i in 1 : length(means)){                 # corresponding
  means2 <- c(means2, rep(means[i], df$n[i])) # proportions
}                                             # sample.
model         <- fitdistr(means2, dbeta,             # set 'means' to tend to  #
                          start = list(shape1 = 1,   # leave mean or 'means2'  # 
                                       shape2 = 10)) # to tend to global mean. # Empirical Bayes (EB)
alpha0        <- model$estimate[1]                                             # corrections (of frequentist 
beta0         <- model$estimate[2]                                             # estimated probabilities)
EBmeans       <- (df$n1 + alpha0)/(df$n + alpha0 + beta0)                      #
df$p1         <- means                                                         #
df$EBp1       <- EBmeans                                  
row.names(df) <- paste("leaf", 1 : 4)
m <- data.frame(global = sum(df$n1)/sum(df$n),
                leave  = sum(df$n1/df$n)/nrow(df),
                estim  = alpha0/(alpha0+beta0))
row.names(m) = "mean"
print.data.frame(df, digits = 3)
print.data.frame(m, digits = 3)

# 2. EB tree correction:
datasetName   <- "banana"
dataset       <- read.csv(paste0("./datasets/", datasetName, ".csv"), stringsAsFactors=F, sep=";")
dataset$class <- standardizeLabels(dataset$class)
tree          <- rpart(class ~ ., dataset, method = "class", cp=0.15)
EBtree <- tree
# EBtree <- EBcorrect(tree)
preds         <- predict(tree, dataset)
print(head(dataset))
prp(tree, extra = 1)
print(head(preds))
print(tree$frame$yval2)
df        <- data.frame(tree$frame$yval2[tree$frame$var == "<leaf>", 2 : 3])
names(df) <- c("n0", "n1")
df$n      <- rowSums(df)
df$p0     <- df$n0/df$n
df$p1     <- df$n1/df$n
print(df)

means     <- df$p1
means2    <- c()                              # generation of the
for (i in 1 : nrow(df)){                      # corresponding
  means2 <- c(means2, rep(df$p1[i], df$n[i])) # proportions
}                                             # sample.
model         <- fitdistr(means2, dbeta,             # set 'means' to tend to  #
                          start = list(shape1 = 1,   # leave mean or 'means2'  #
                                       shape2 = 10)) # to tend to global mean. # Empirical Bayes (EB)
alpha0        <- model$estimate[1]                                             # corrections (of frequentist
beta0         <- model$estimate[2]                                             # estimated probabilities)
EBmeans       <- (df$n1 + alpha0)/(df$n + alpha0 + beta0)                      #                                                       #
df$EBp0       <- 1 - EBmeans
df$EBp1       <- EBmeans
print(df)

EBtree$frame$yval2[tree$frame$var == "<leaf>", 4 : 5] <- as.matrix(df[, c("EBp0", "EBp1")]) # CORRECTION
EBpreds <- predict(EBtree, dataset)
prp(EBtree, extra = 1)
print(head(preds))
print(head(EBpreds))
p1        <- preds[, 2]   # probabilities of
EBp1      <- EBpreds[, 2] # class label "1"
logLoss   <- MLmetrics::LogLoss(p1, dataset$class)
EBlogLoss <- MLmetrics::LogLoss(EBp1, dataset$class)
cat("logLoss =", logLoss, ", EBlogLoss =", EBlogLoss)

