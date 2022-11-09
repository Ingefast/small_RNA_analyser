## ---------------------------
##
## sRNA.ordinatio_analyser.r
##
## Checking replicability and similarity between replicates and conditions using sRNA
## normalized expression on genes, TEs
##
## Author: Juan Santos, SLU, 2020 December
##
## ---------------------------

## Cleans the workspace and disable scientific notation 
rm(list = ls())
options(scipen=999)

## set working directory
setwd("/example/output")

## required R libraries

library(vegan)
library(MASS)


## input table

data1 <- read.table("/example/output/total_table.size_24.genes.txt", header = TRUE, row.names = 1)

head(data1)
dim(data1)
summary(data1)

# selecting columns
data1 <- data1[ , c("mutantA_rep1_rpm", "mutantA_rep2_rpm", "mutantB_rep1_rpm", "mutantB_rep2_rpm", "wt_rep1_rpm", "wt_rep2_rpm")]
head(data1)
dim(data1)

#transposing the table
data1 <- t(data1)
head(data1)
dim(data1)

#Removes columns with only zeroes
data1 <- data1[, colSums(data1 != 0, na.rm = TRUE) > 0]

dim(data1)


##########################################################################
################  CHOOSING ORDINATION METHOD         #####################
##########################################################################

## making the ordination with PCA
#ord <- rda(data1)

## making the ordination with DCA
#ord <- decorana(data1)

## making the ordination with NMDS
ord <- metaMDS(data1,autotransform = FALSE, zerodist = "ignore")



##########################################################################
#########################    THE PLOTTING      ###########################
##########################################################################

plot(ord, display = "sites", type="n", main = "NMDS ordination based on\n 24nt sRNA RPM values in genes")

#Defines a series of symbols
vpch<- c(15, 17)

#Defines a series of colors
vcols <- c("gold", "gold", "red", "red", "blue", "blue" )

#Defines the names of the samples for plot identity
ftypes<-c("mutantA_rep1_rpm", "mutantA_rep2_rpm", "mutantB_rep1_rpm", "mutantB_rep2_rpm", "wt_rep1_rpm", "wt_rep2_rpm")

#Shows the points labeled after Plot identity
points(ord, col = rep(vcols,each=nrow(data1)/length(unique(ftypes))), cex = 1.5, pch = rep(vpch,each=nrow(data1)/length(unique(ftypes))))

#creates a vector for the category names to be shown in the legend
unique(ftypes)->leg.txt

#creates a legend inside the diagram
legend("topleft", cex=0.6, ncol=2, inset=.02, title="Samples", legend=leg.txt, col=vcols, pch=vpch, horiz=FALSE)




