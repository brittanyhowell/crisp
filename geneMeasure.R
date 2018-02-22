setwd("~/Desktop/Petra Palenikova/")

CtrlGrowthEffect <- read.csv(header = TRUE, "dataSUR/ControlGrowth.csv") 
# Control, cutting gRNAs targeting genes that DO affect cell fitness (aka essential genes). These gRNAs will kill cells (i.e. drop out) no matter what is the effect on the STAT5 reporte

CtrlNoCut <- read.csv(header = TRUE, "dataSUR/ControlNoCutting.csv") 
# Control, non-cutting gRNAs (completely ineffective) 

CtrlNoGrowthEffect <- read.csv(header = TRUE, "dataSUR/ControlNoGrowth.csv") 
# Control, cutting gRNAs targeting genes that do not affect cell fitness

JakSTAT <- read.csv(header = TRUE, "dataSUR/JakStatPathway.csv")
# gRNAs targeting genes belonging to or associated to Jak/STAT pathway. Part of these gRNAs may be essential and kill cells, no matter what is the effect on the STAT5 reporter. Other gRNAs may have little or no “direct” effect on cell fitness but if neomycin is applied the can increase/decrease it depending on their positive/negative effect on the reporter by increasing/decreasing the expression of the neoR gene.



# Mix up the order of the cutting, but not affecting the growth controls
rand.CtrlNoGrowthEffect <- CtrlNoGrowthEffect[sample(1:nrow(CtrlNoGrowthEffect)), ]

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


        # Dummy df 
        test.df <- data.frame(replicate(5,sample(0:1,5,rep=TRUE)))
        colnames(test.df) <- c("pool", " 1", " 2", " 3", " 4")
        
        # Remove the pool column
        pool <- test.df$pool
        pool <- pool*10

        test.df = subset(test.df, select = -c(pool) )
        
        
JakSTAT.pool <- aggregatedCtrls$JakSTATpool
aggregatedCtrls.countOnly <- subset(aggregatedCtrls, select = -c(Kosuke_Id, Group.1, JakSTATpool))

FoldChange <- as.data.frame(apply(aggregatedCtrls.countOnly, 2, function(x) log2(x+16)-log2(JakSTAT.pool+16) ))






# Compute difference between Pool & Sample 


#+- combine replicates after stage 2
