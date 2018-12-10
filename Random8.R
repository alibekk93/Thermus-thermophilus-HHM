### Random secondary structure HMM (3 states)

# Load packages
library(HMM)

# Load data and add a fold number for cross-validation
data.full <- data.frame(read.csv2('E:/Project/Bioinformatics project/ss-filtered.csv'))
n5 <- round(length(data.full$ID)/5)
data.full <- data.full[sample(1:nrow(data.full)), ]
data.full$Fold <- 1
data.full$Fold[1:n5] <- 2
data.full$Fold[n5+1:n5] <- 3
data.full$Fold[n5*2+1:n5] <- 4
data.full$Fold[n5*3+1:n5] <- 5
s.structures <- c('-', 'S', 'E', 'T', 'H', 'B', 'G', 'I')
a.acids <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y')

tosingles <- function(data)
{
  return(unlist(strsplit(data, '')))
}

# Function to:
### split data into Train / Test sets
runHMM <- function(data, testfold)
{
  train.set <- subset(data, subset = (Fold != testfold))
  test.set <- subset(data, subset = (Fold == testfold))
  
  ### derive HMM and get probability matrixes
  hmm <<- initHMM(s.structures, a.acids, transProbs = matrix(c(rep(0.125, 64)), 8))
  #hmm <- baumWelch(hmm, AAs.train, 10)
  
  ### get a predicted amino acid chain for the test set
  Viterbi.results.table <- data.frame('Real SS' = test.set$SS_seq, 'Predicted SS' = c(rep(0, length(test.set$SS_seq))), 'Length' = c(rep(0, length(test.set$SS_seq))), 'Trues' = c(rep(0, length(test.set$SS_seq))), 'Falses' = c(rep(0, length(test.set$SS_seq))), 'Accuracy' = c(rep(0, length(test.set$SS_seq))))
  
  for(l in 1:length(test.set$AA_seq))
  {
    prediction <- viterbi(hmm, tosingles(toString(test.set$AA_seq[l])))
    result <- ''
    for(c in prediction)
    {
      result <- paste(result, c, sep = '')
    }
    Viterbi.results.table$Predicted.SS[l] <- result
    Viterbi.results.table$Length[l] <- length(prediction)
    t <- 0
    f <- 0
    for(i in 1:length(prediction))
    {
      if(tosingles(toString(test.set$SS_seq[l]))[i] == prediction[i])
      {
        t <- t + 1
      }
      else if(tosingles(toString(test.set$SS_seq[l]))[i] != prediction[i])
      {
        f <- f + 1
      }
    }
    Viterbi.results.table$Trues[l] <- t
    Viterbi.results.table$Falses[l] <- f
    Viterbi.results.table$Accuracy[l] <- t/i
  }
  return(Viterbi.results.table)
}

results <- c(0, 0, 0, 0, 0)

for(i in 1:5)
{
  filename <- paste('E:/Project/Bioinformatics project/', i, '.csv', sep="")
  result <- runHMM(data.full, i)
  results[i] <- mean(result$Accuracy)
  write.csv(result, filename)
  hist(result$Accuracy, breaks = 50)
}

print(matrix(c(results, mean(results)), dimnames = list(c('1', '2', '3', '4', '5', 'Mean'), 'Accuracy')))

