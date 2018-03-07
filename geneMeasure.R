setwd("~/Desktop/Petra Palenikova/")

library(ggplot2)
library(reshape)
library(tidyr)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)

CtrlGrowthEffect.file <- read.csv(header = TRUE, "dataSUR/ControlGrowth.csv") 
# Control, cutting gRNAs targeting genes that DO affect cell fitness (aka essential genes). These gRNAs will kill cells (i.e. drop out) no matter what is the effect on the STAT5 reporter

CtrlNoCut.file <- read.csv(header = TRUE, "dataSUR/ControlNoCutting.csv") 
# Control, non-cutting gRNAs (completely ineffective) 

CtrlNoGrowthEffect.file <- read.csv(header = TRUE, "dataSUR/ControlNoGrowth.csv") 
# Control, cutting gRNAs targeting genes that do not affect cell fitness

JakSTAT.file <- read.csv(header = TRUE, "dataSUR/JakStatPathway.csv")
# gRNAs targeting genes belonging to or associated to Jak/STAT pathway. Part of these gRNAs may be essential and kill cells, no matter what is the effect on the STAT5 reporter. Other gRNAs may have little or no “direct” effect on cell fitness but if neomycin is applied the can increase/decrease it depending on their positive/negative effect on the reporter by increasing/decreasing the expression of the neoR gene.



# Mix up the order of the cutting, but not affecting the growth controls
rand.CtrlNoGrowthEffect <- CtrlNoGrowthEffect.file[sample(1:nrow(CtrlNoGrowthEffect.file)), ]

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
CtrlNoGrowthEffect.agg <- aggregate(x=rand.CtrlNoGrowthEffect, by = list(index), FUN=sum)
CtrlNoGrowthEffect.agg <- subset(CtrlNoGrowthEffect.agg, select = -c(Group.1))

# Compute Y=log2(x+16)-log2(CTRL+16) for all guides, in all samples. 
# The CTRL is the Jak STAT pool number. 16 is a Leo chosen number, used to make the fold changes seem a little less dramatic. This calculation is the difference between the control and the effect of the gRNA - so it is the measure of how much that gene is needed for the cell to prolifterate.


        # Dummy df 
        # test.df <- data.frame(replicate(5,sample(0:1,5,rep=TRUE)))
        # colnames(test.df) <- c("pool", " 1", " 2", " 3", " 4")
        
        # Remove the pool column
        # pool <- test.df$pool
        # pool <- pool*10

        # test.df = subset(test.df, select = -c(pool) )
        
        
# aggregatedCtrls.pool <- aggregatedCtrls$JakSTATpool
# aggregatedCtrls.countOnly <- subset(aggregatedCtrls, select = -c(Kosuke_Id, Group.1, JakSTATpool))

# FoldChange.ctrlNoGrowthEffect <- as.data.frame(apply(aggregatedCtrls.countOnly, 2, function(x) log2(x+16)-log2(JakSTAT.pool+16) ))



##
# Make the gRNA name the label of the row

CtrlGrowthEffect <- CtrlGrowthEffect.file[,-1]
rownames(CtrlGrowthEffect) <- CtrlGrowthEffect.file[,1]

CtrlNoCut <- CtrlNoCut.file[,-1]
rownames(CtrlNoCut) <- CtrlNoCut.file[,1]

## Note cannot make NAs rownames, and gRNAs were lost in the aggregation of 5gRNAs.
CtrlNoGrowthEffect <- CtrlNoGrowthEffect.agg[,-1]
# rownames(CtrlNoGrowthEffect) <- CtrlNoGrowthEffect.agg[,1]

JakSTAT <- JakSTAT.file[,-1]
rownames(JakSTAT) <- JakSTAT.file[,1]


## Get the JakSTAT Pool as a different vector and slice it out


CtrlGrowthEffect.pool <- CtrlGrowthEffect$JakSTATpool
CtrlGrowthEffect <- subset(CtrlGrowthEffect, select = -c(JakSTATpool))

CtrlNoCut.pool <- CtrlNoCut$JakSTATpool
CtrlNoCut <- subset(CtrlNoCut, select = -c(JakSTATpool))

CtrlNoGrowthEffect.pool <- CtrlNoGrowthEffect$JakSTATpool
CtrlNoGrowthEffect <- subset(CtrlNoGrowthEffect, select = -c(JakSTATpool))

JakSTAT.pool <- JakSTAT$JakSTATpool
JakSTAT <- subset(JakSTAT, select = -c(JakSTATpool))

## Calculate foldchange for all samples

# make sure this actually works
# test.df <- data.frame(replicate(5,sample(0:1,5,rep=TRUE)))
# test.col <- data.frame(replicate(1,sample(0:3,5,rep=TRUE)))
# colnames(test.col) <- "Cat"
# testApply <- as.data.frame(apply(test.df, 2, function(x) x+test.col ))


rawFoldChange.CtrlGrowthEffect <- as.data.frame(apply(CtrlGrowthEffect, 2, function(x) log2(x+16)-log2(CtrlGrowthEffect.pool+16) ))
rawFoldChange.CtrlNoCut <- as.data.frame(apply(CtrlNoCut, 2, function(x) log2(x+16)-log2(CtrlNoCut.pool+16) ))
rawFoldChange.CtrlNoGrowthEffect <- as.data.frame(apply(CtrlNoGrowthEffect, 2, function(x) log2(x+16)-log2(CtrlNoGrowthEffect.pool+16) ))
rawFoldChange.JakSTAT <- as.data.frame(apply(JakSTAT, 2, function(x) log2(x+16)-log2(JakSTAT.pool+16) ))

## Finding the median across guides.
medCalc.CtrlGrowthEffect <- as.data.frame(apply(rawFoldChange.CtrlGrowthEffect, 1, median ))
colnames(medCalc.CtrlGrowthEffect) <- "median"
med.CtrlGrowthEffect <- medCalc.CtrlGrowthEffect$median

medCalc.CtrlNoCut <- as.data.frame(apply(rawFoldChange.CtrlNoCut, 1, median  ))
colnames(medCalc.CtrlNoCut) <- "median"
med.CtrlNoCut <- medCalc.CtrlNoCut$median

medCalc.CtrlNoGrowthEffect <- as.data.frame(apply(rawFoldChange.CtrlNoGrowthEffect, 1, median  ))
colnames(medCalc.CtrlNoGrowthEffect) <- "median"
med.CtrlNoGrowthEffect <- medCalc.CtrlNoGrowthEffect$median

medCalc.JakSTAT <- as.data.frame(apply(rawFoldChange.JakSTAT,  1, median  ))
colnames(medCalc.JakSTAT) <- "median"
med.JakSTAT <- medCalc.JakSTAT$median



## Subtract the median from the foldchange
FC.CtrlGrowthEffect <- as.data.frame(apply(rawFoldChange.CtrlGrowthEffect, 2, function(x) x - med.CtrlGrowthEffect))
FC.CtrlNoCut <- as.data.frame(apply(rawFoldChange.CtrlNoCut, 2, function(x) x - med.CtrlNoCut))
FC.CtrlNoGrowthEffect <- as.data.frame(apply(rawFoldChange.CtrlNoGrowthEffect, 2, function(x) x - med.CtrlNoGrowthEffect))
FC.JakSTAT <- as.data.frame(apply(rawFoldChange.JakSTAT, 2, function(x) x - med.JakSTAT))


## Testing the median subtraction works
# test.df <- data.frame(replicate(5,sample(0:10,5,rep=TRUE)))
# test.vec <- data.frame(replicate(1,sample(1:5,5,rep=TRUE)))
# colnames(test.df) <- c("first", "second", "third", "fourth", "fifth")
# colnames(test.vec) <- "vec"
# combine <- as.data.frame(apply(test.df, 2, function(x) x-test.vec))

### Okay choice here. Combine all of the genes or just the exons?

## Subset data to make 
# SUR1.FC.JakSTAR <- subset(FC.JakSTAT, grepl(, SUR1, )
sub4dman <- as.numeric(head(FC.JakSTAT, n = 5))

# subset(combine.FC.JakSTAT, select = -c(gene1, gene2, other))


# HEY LOOK A PLOT
ggplot(data=gene.FC.CtrlGrowthEffect, aes(x=SUR1.250x.A.DPI7, y=SUR1.250x.A.DPI22)) +
  theme_bw() +
    geom_point() +
    geom_abline(intercept = 0, slope = 1)
# There are indeed differences in foldchange across that there sample. Back to the plan. 
       


## Calculting median(noGrowthGuides)

# Combine all data into a single column
oneCol.CtrlNoGrowthEffect <- melt(CtrlNoGrowthEffect,var='remain')
CtrlNoGrowthEffect.all.med <- median(oneCol.CtrlNoGrowthEffect$value)



## Test combine columns
# test.df <- data.frame(replicate(5,sample(0:10,5,rep=TRUE)))
# colnames(test.df) <- c("first", "second", "third", "fourth", "fifth")
# melt(test.df,var='remain')

# Compute difference between Pool & Sample 


#+- combine replicates after stage 2



## Extract the gene name, make gene and exon new columns
    
    # # Test on a smol df
    # smol <- head(FC.CtrlGrowthEffect, n=4)
    # smol <- subset(smol, select = c(gRNA, SUR7.50x.C.DPI7,SUR7.50x.C.DPI15))
    
    # sep.smol <- separate(smol, gRNA, c("gene1", "gene2", "other", "exon"))
    # combine.smol <- unite(sep.smol, "gene", c("gene1", "gene2"), sep = "_", remove = FALSE)
    # new.smol <- subset(combine.smol, select = -c(gene1, gene2, other))
    
    # Make gRNA name a column again
    FC.CtrlGrowthEffect$gRNA <- rownames(FC.CtrlGrowthEffect)
    FC.CtrlNoCut$gRNA <- rownames(FC.CtrlNoCut)
    FC.CtrlNoGrowthEffect$gRNA <- rownames(FC.CtrlNoGrowthEffect)
    FC.JakSTAT$gRNA <- rownames(FC.JakSTAT)
    
    # Separate gRNA into gene and exon labels
    sep.FC.CtrlGrowthEffect <- separate(FC.CtrlGrowthEffect, gRNA, c("gene1", "gene2", "other", "exon"))
    sep.FC.CtrlNoCut <- separate(FC.CtrlNoCut, gRNA, c("gene1", "gene2", "other", "exon"))
    sep.FC.CtrlNoGrowthEffect <- separate(FC.CtrlNoGrowthEffect, gRNA, c("gene1", "gene2", "other", "exon"))
    sep.FC.JakSTAT <- separate(FC.JakSTAT, gRNA, c("gene1", "gene2", "other", "exon"))
    
    # Combine gene name columns
    combine.FC.CtrlGrowthEffect <- unite(sep.FC.CtrlGrowthEffect, "gene", c("gene1", "gene2"), sep = "_", remove = FALSE)
    combine.FC.CtrlNoCut <- unite(sep.FC.CtrlNoCut, "gene", c("gene1", "gene2"), sep = "_", remove = FALSE)
    combine.FC.CtrlNoGrowthEffect <- unite(sep.FC.CtrlNoGrowthEffect, "gene", c("gene1", "gene2"), sep = "_", remove = FALSE)
    combine.FC.JakSTAT <- unite(sep.FC.JakSTAT, "gene", c("gene1", "gene2"), sep = "_", remove = FALSE)
    
    # Remove the unnecessary columns
    gene.FC.CtrlGrowthEffect <- subset(combine.FC.CtrlGrowthEffect, select = -c(gene1, gene2, other))
    gene.FC.CtrlNoCut <- subset(combine.FC.CtrlNoCut, select = -c(gene1, gene2, other))
    gene.FC.CtrlNoGrowthEffect <- subset(combine.FC.CtrlNoGrowthEffect, select = -c(gene1, gene2, other))
    gene.FC.JakSTAT <- subset(combine.FC.JakSTAT, select = -c(gene1, gene2, other))
