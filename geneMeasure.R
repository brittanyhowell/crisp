setwd("~/Desktop/Petra Palenikova/")

CtrlGrowthEffect <- read.csv(header = TRUE, "dataSUR/ControlGrowth.csv")
CtrlNoCut <- read.csv(header = TRUE, "dataSUR/ControlNoCutting.csv")
CtrlNoGrowthEffect <- read.csv(header = TRUE, "dataSUR/ControlNoGrowth.csv") # Combine 5 gRNAs into one. (Randomly)
JakSTAT <- read.csv(header = TRUE, "dataSUR/JakStatPathway.csv")

sub.CNGE <- head(CtrlNoGrowthEffect)

rand.sub.CNGE <- sub.CNGE[sample(1:nrow(sub.CNGE)), ]

# test, lets make a new data frame
test <- data.frame(replicate(10,sample(0:2,10,rep=TRUE)))

# And now we can, try to make a new table from the sums of this one. 
# First thing, sum two things together.
# Row sum?


## make a list of numbers in groups of 5. 

# determine how many rows you need
numRows <- dim(CtrlGrowthEffect)[1]
numIndex <- ceiling(numRows/5)

# Make the index: make a list of numbers in groups of five, order the list, and then truncate it to the exact number of records required.  
index <- rep(1:numIndex, 5)
index <- index[order(index)]
index <- index[1:numRows]

# Concatenate the five names of the gRNAs
gRNA.names <- data.frame(CtrlGrowthEffect$Kosuke_Id)
gRNA.names$index <- index

aggregatedNames <- aggregate(x=gRNA.names$CtrlGrowthEffect.Kosuke_Id, by = list(index), FUN=print)
## Okay so this is unfinished, but I am going to return to it if for whatever reason, I need to identify which gRNAs have been aggregated

# Turn numbers into summable numbers  
indx <- sapply(CtrlGrowthEffect, is.factor)
CtrlGrowthEffect[indx] <- lapply(CtrlGrowthEffect[indx], function(x) as.numeric(as.character(x)))

## Combine numbers on the basis of their index 
aggregatedCtrls <- aggregate(x=CtrlGrowthEffect, by = list(index), FUN=sum)


#apply(matrix, rows (1) or columns (2), function)
testSum <- apply(test,2,sum)

# Compute difference between Pool & Sample 


#+- combine replicates after stage 2


new_test <- data.frame()
t(sapply(1:3,function(x){
  tosum <- test_index[which(test_index$index == x),]
  summed <- colSums(tosum)
}))

