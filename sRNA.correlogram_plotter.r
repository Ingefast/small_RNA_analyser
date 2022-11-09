## ---------------------------
##
## sRNA.correlogram_plotter.r
##
## Checking replicability and correlation between replicates across and within conditions
## in sRNA genomic tracks (genes, TEs) or bedGraph files.
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
library(corrgram)
library(lattice)

## input table
total <- read.table("/example/total_table.size_24.genes.txt", header = T)

#title and name of the final figure file
sample_name <- c("Correlogram of 24 nt sRNA expression in genes (RPM)")

total <- total[c(-1)] 
head(total)
total[total == 0] <- NA
head(total)

#this removes outliers
#total <- ifelse(total <= -5, NA, total)
#total <- ifelse(total <= -5, NA, total)

# selecting columns
total <- total[ , c(7:12)]

#saving to png file
png(paste(sample_name,".png", sep=""), width=13, height=13, units="in", res=150)

#adds the color scale
Lab.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")

#Tells the panel.cor function to consider only pairwise complete obs otherwise you #get no corr coeff

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use="pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}

#plotting of only the RPM values, not the raw counts
pairs(main = sample_name, total[ , c(1:6)], xlim=c(0, 0.05), ylim=c(0, 0.05), lower.panel = function(...) smoothScatter(..., colramp = Lab.palette, nrpoints = 0, add = TRUE), upper.panel = panel.cor, diag.panel = panel.minmax)

dev.off()
