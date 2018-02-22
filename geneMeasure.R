setwd("~/Desktop/Petra Palenikova/")

CtrlGrowthEffect <- read.csv(header = TRUE, "dataSUR/ControlGrowth.csv")
CtrlNoCut <- read.csv(header = TRUE, "dataSUR/ControlNoCutting.csv")
CtrlNoGrowthEffect <- read.csv(header = TRUE, "dataSUR/ControlNoGrowth.csv") # Combine 5 gRNAs into one. (Randomly)
JakSTAT <- read.csv(header = TRUE, "dataSUR/JakStatPathway.csv")


rand.CtrlNoGrowthEffect <- CtrlNoGrowthEffect[sample(1:nrow(CtrlNoGrowthEffect)), ]

# And now we can, try to make a new table from the sums of this one. 
# First thing, sum two things together.
# Row sum?


## make a list of numbers in groups of 5. 

# determine how many rows you need
numRows <- dim(rand.CtrlNoGrowthEffect)[1]
numIndex <- ceiling(numRows/5)

# Make the index: make a list of numbers in groups of five, order the list, and then truncate it to the exact number of records required.  
index <- rep(1:numIndex, 5)
index <- index[order(index)]
index <- index[1:numRows]

# # Concatenate the five names of the gRNAs
# gRNA.names <- data.frame(CtrlGrowthEffect$Kosuke_Id)
# gRNA.names$index <- index
# 
# aggregatedNames <- aggregate(x=gRNA.names$CtrlGrowthEffect.Kosuke_Id, by = list(index), FUN=print)
# ## Okay so this is unfinished, but I am going to return to it if for whatever reason, I need to identify which gRNAs have been aggregated

# Turn numbers into summable numbers  
indx <- sapply(rand.CtrlNoGrowthEffect, is.factor)
rand.CtrlNoGrowthEffect[indx] <- lapply(rand.CtrlNoGrowthEffect[indx], function(x) as.numeric(as.character(x)))

## Combine numbers on the basis of their index 
aggregatedCtrls <- aggregate(x=rand.CtrlNoGrowthEffect, by = list(index), FUN=sum)

# Compute Y=log2(x+16)-log2(CTRL+16) for all guides, in all samples. 
# The CTRL is the Jak STAT pool number. 16 is a Leo chosen number, used to make the fold changes seem a little less dramatic. This calculation is the difference between the control and the effect of the gRNA - so it is the measure of how much that gene is needed for the cell to prolifterate.

# Compute difference between Pool & Sample 


#+- combine replicates after stage 2
