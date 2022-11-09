#Juan Santos, August 2021
## ---------------------------
##
## sRNA.size_distribution_plotter.r
##
## Plot sRNA size distribution of reads mapping into genes and TEs
##
## Author: Juan Santos, SLU, 2020 December
##
## ---------------------------

## Cleans the workspace and disable scientific notation 
rm(list = ls())
options(scipen=999)

## set working directory
setwd("/example/output")


## read table with number of mapped reads in each sample as a baseline to make a RPM normalization

baseline <- read.table("read_n_baseline.txt", header=TRUE) 

head(baseline)

## reads output files one by one

#first-rep
wt_rep1 <- read.table("/example/output/wt_rep1.size_1830.te.reads.txt", check.names = FALSE, header=FALSE) 

wt_rep1$len <- wt_rep1$V3 - wt_rep1$V2 
wt_rep1_count <- data.frame(table(wt_rep1$len))
wt_rep1_count$rpm <-  round ((wt_rep1_count$Freq * 1000000 ) / baseline[ baseline$sample_name == "wt_rep1", 2], digits = 0 )


#second-rep
wt_rep2 <- read.table("/Users/jsantos/Desktop/ingefast_github/sRNA_github/example/output/wt_rep2.size_1830.te.reads.txt", check.names = FALSE, header=FALSE) 

wt_rep2$len <- wt_rep2$V3 - wt_rep2$V2 
wt_rep2_count <- data.frame(table(wt_rep2$len))
wt_rep2_count$rpm <-  round ((wt_rep2_count$Freq * 1000000 ) / baseline[ baseline$sample_name == "wt_rep2", 2], digits = 0 )


#first-rep
mutantA_rep1 <- read.table("/Users/jsantos/Desktop/ingefast_github/sRNA_github/example/output/mutantA_rep1.size_1830.te.reads.txt", check.names = FALSE, header=FALSE) 

mutantA_rep1$len <- mutantA_rep1$V3 - mutantA_rep1$V2 
mutantA_rep1_count <- data.frame(table(mutantA_rep1$len))
mutantA_rep1_count$rpm <-  round ((mutantA_rep1_count$Freq * 1000000 ) / baseline[ baseline$sample_name == "mutantA_rep1", 2], digits = 0 )

#second-rep
mutantA_rep2 <- read.table("/Users/jsantos/Desktop/ingefast_github/sRNA_github/example/output/mutantA_rep2.size_1830.te.reads.txt", check.names = FALSE, header=FALSE) 

mutantA_rep2$len <- mutantA_rep2$V3 - mutantA_rep2$V2 
mutantA_rep2_count <- data.frame(table(mutantA_rep2$len))
mutantA_rep2_count$rpm <-  round ((mutantA_rep2_count$Freq * 1000000 ) / baseline[ baseline$sample_name == "mutantA_rep2", 2], digits = 0 )

#first-rep
mutantB_rep1 <- read.table("/Users/jsantos/Desktop/ingefast_github/sRNA_github/example/output/mutantB_rep1.size_1830.te.reads.txt", check.names = FALSE, header=FALSE) 

mutantB_rep1$len <- mutantB_rep1$V3 - mutantB_rep1$V2 
mutantB_rep1_count <- data.frame(table(mutantB_rep1$len))
mutantB_rep1_count$rpm <-  round ((mutantB_rep1_count$Freq * 1000000 ) / baseline[ baseline$sample_name == "mutantB_rep1", 2], digits = 0 )

#second-rep
mutantB_rep2 <- read.table("/Users/jsantos/Desktop/ingefast_github/sRNA_github/example/output/mutantB_rep2.size_1830.te.reads.txt", check.names = FALSE, header=FALSE) 

mutantB_rep2$len <- mutantB_rep2$V3 - mutantB_rep2$V2 
mutantB_rep2_count <- data.frame(table(mutantB_rep2$len))
mutantB_rep2_count$rpm <-  round ((mutantB_rep2_count$Freq * 1000000 ) / baseline[ baseline$sample_name == "mutantB_rep2", 2], digits = 0 )


########################################################################################
########################################################################################
###################      the making of a total table    ################################
########################################################################################
########################################################################################


total_frame <- data.frame( c(seq( from = 18 , to = 30 , by = 1) ))

colnames(total_frame) <- c("id")
head(total_frame)
dim(total_frame)

MMM_wt_rep1 <- merge(total_frame, wt_rep1_count, by.x = "id", by.y = "Var1", all.x = TRUE)
dim(MMM_wt_rep1)
head(MMM_wt_rep1)

MMM_wt_rep2 <- merge(total_frame, wt_rep2_count, by.x = "id", by.y = "Var1", all.x = TRUE)
dim(MMM_wt_rep2)
head(MMM_wt_rep2)


MMM_mutantA_rep1 <- merge(total_frame, mutantA_rep1_count, by.x = "id", by.y = "Var1", all.x = TRUE)
dim(MMM_mutantA_rep1)
head(MMM_mutantA_rep1)

MMM_mutantA_rep2 <- merge(total_frame, mutantA_rep2_count, by.x = "id", by.y = "Var1", all.x = TRUE)
dim(MMM_mutantA_rep2)
head(MMM_mutantA_rep2)


MMM_mutantB_rep1 <- merge(total_frame, mutantB_rep1_count, by.x = "id", by.y = "Var1", all.x = TRUE)
dim(MMM_mutantB_rep1)
head(MMM_mutantB_rep1)

MMM_mutantB_rep2 <- merge(total_frame, mutantB_rep2_count, by.x = "id", by.y = "Var1", all.x = TRUE)
dim(MMM_mutantB_rep2)
head(MMM_mutantB_rep2)


total_norm <- data.frame( MMM_wt_rep1$id, MMM_wt_rep1$rpm, MMM_wt_rep2$rpm, MMM_mutantA_rep1$rpm, MMM_mutantA_rep2$rpm, MMM_mutantB_rep1$rpm, MMM_mutantB_rep2$rpm)

colnames(total_norm) <- c( "size", "wt_rep1", "wt_rep2", "mutantA_rep1", "mutantA_rep2", "mutantB_rep1", "mutantB_rep2")



##################################################
##################################################
################## PLOTTING ##################
##################################################
##################################################


## add extra space to right margin of plot within frame
par(mar=c(5, 5, 8, 5) + 0.5)

plot(total_norm$size, total_norm$wt_rep1, main = c(paste("sRNA size distribution over TEs"), paste("n of categories =", nrow(total_norm))), type = "o", col= 'red', bg = "white", lwd = 3, pch = 21, ylim = c(0, 50000), xlim = c(18, 30), axes = FALSE, xlab = "", ylab = "")
mtext("RPM", side=2, col="black", line=4, cex=1) 
axis(2, at=c(0, 50000), col = "black", col.axis = "black", las = 1, cex = 1)

## Draw the X axis
axis(1, at=c(total_norm$size), col = "black", col.axis = "black", labels = T)
mtext("sRNA size (nt)", side = 1, line = 2.5, cex = 1)  

## Allow a second plot on the same graph
par(new=TRUE)
plot(total_norm$size, total_norm$wt_rep2, ylim = c(0, 50000), xlim = c(18, 30), xlab = "", ylab = "", type = "o", col= 'orange', bg = "white", lwd = 3,  pch = 21, axes = FALSE)

## Allow a third plot on the same graph
par(new=TRUE)
plot(total_norm$size, total_norm$mutantA_rep1, ylim = c(0, 50000), xlim = c(18, 30), xlab = "", ylab = "", type = "o", col= 'darkolivegreen4', bg = "white", lwd = 3,  pch = 21, axes = FALSE)

## Allow a fourth plot on the same graph
par(new=TRUE)
plot(total_norm$size, total_norm$mutantA_rep2, ylim = c(0, 50000), xlim = c(18, 30), xlab = "", ylab = "", type = "o", col= 'darkolivegreen1', bg = "white", lwd = 3,  pch = 21, axes = FALSE)

## Allow a fifth plot on the same graph
par(new=TRUE)
plot(total_norm$size, total_norm$mutantB_rep1, ylim = c(0, 50000), xlim = c(18, 30), xlab = "", ylab = "", type = "o", col= 'blue', bg = "white", lwd = 3,  pch = 21, axes = FALSE)

## Allow a sixth plot on the same graph
par(new=TRUE)
plot(total_norm$size, total_norm$mutantB_rep2, ylim = c(0, 50000), xlim = c(18, 30), xlab = "", ylab = "", type = "o", col= 'cyan', bg = "white", lwd = 3,  pch = 21, axes = FALSE)

legend( "topright", c( "wt_rep1", "wt_rep2", "mutantA_rep1", "mutantA_rep2", "mutantB_rep1", "mutantB_rep2"), title.col = "black", col = c('red', 'orange', 'darkolivegreen4', 'darkolivegreen1', 'blue', 'cyan'), inset=0.01, cex = 1, text.col = c('red', 'orange', 'darkolivegreen4', 'darkolivegreen1', 'blue', 'cyan'), lty = c(1), lwd=c(5,5,5,5), pch = c(-1), bty="n")





