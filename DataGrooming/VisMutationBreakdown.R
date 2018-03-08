# Title     : Plot information associated with Mutation Breakdown generated from MetaDataTableMaker
# Objective : See Title
# Created by: schencro
# Created on: 3/6/18

library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)

# args <- commandArgs(trailingOnly = TRUE)
# wd <- args[1]
# setwd(wd)
# data <- args[2]
# df <- read.csv(data, header=TRUE, sep="\t", stringsAsFactors=F)

# For Development
setwd("/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/DataTables/")
df <- read.csv('MutationBreakdown.txt', sep='\t', header=T, stringsAsFactors = F)
df$Cancer <- as.factor(df$Cancer)
df$TotalSNVs <- df$TotalSNVs/1000000

toPlot <- data.frame(matrix(NA, ncol=length(colnames(df)), nrow=0))
toPlotDecile <- toPlot <- data.frame(matrix(NA, ncol=2, nrow=0))
colnames(toPlotDecile) <- c("Cancer","TotalSNVs")
for(i in seq(1,length(unique(df$Cancer)))){
  data = subset(df, df$Cancer==unique(df$Cancer)[i])
  data = data %>% filter(TotalSNVs < quantile(data$TotalSNVs, 0.9))
  data = data %>% filter(TotalSNVs > quantile(data$TotalSNVs, 0.1))
  intervalVariants = quantile(data$TotalCodingRegion, probs=1:10/10)
  dataOut <- data.frame(Cancer=rep(unique(df$Cancer)[i],length(intervalVariants)), TotalSNVs=intervalVariants)
  toPlot <- rbind(toPlot, data)
  toPlotDecile <- rbind(toPlotDecile, dataOut)
}

maxForEach <- aggregate(TotalSNVs ~ Cancer, data = toPlot, mean)
colnames(maxForEach) <- colnames(toPlotDecile)
out <- factor(maxForEach, levels=unique(maxForEach$Cancer[order(maxForEach$TotalSNVs)]), ordered=TRUE)
customOrder <- levels(out)
p1 <- ggplot(toPlot, aes(x=Cancer, y=TotalSNVs, color="lightred", group=interaction(toPlot$TotalSNVs,toPlot$Cancer))) +  geom_point(size=0.4,position=position_dodge(width=0.6), show.legend=F) +
  scale_x_discrete(limits=customOrder) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid.minor.y = element_blank()) +
  # stat_summary(fun.y=mean, geom="point", shape="\U2014", size=3, colour="red") +
  ylab("SNVs per Mb") + xlab("")


toPlot <- data.frame(matrix(NA, ncol=length(colnames(df)), nrow=0))
toPlotDecile <- toPlot <- data.frame(matrix(NA, ncol=2, nrow=0))
colnames(toPlotDecile) <- c("Cancer","TotalNonCodingRegion")
for(i in seq(1,length(unique(df$Cancer)))){
  data = subset(df, df$Cancer==unique(df$Cancer)[i])
  data = data %>% filter(TotalNonCodingRegion < quantile(data$TotalNonCodingRegion, 0.9))
  data = data %>% filter(TotalNonCodingRegion > quantile(data$TotalNonCodingRegion, 0.1))
  intervalVariants = quantile(data$TotalCodingRegion, probs=1:10/10)
  dataOut <- data.frame(Cancer=rep(unique(df$Cancer)[i],length(intervalVariants)), TotalSNVs=intervalVariants)
  toPlot <- rbind(toPlot, data)
  toPlotDecile <- rbind(toPlotDecile, dataOut)
}

dfTypes <- data.frame(Cancer=toPlot$Cancer, TotalCodingRegion = toPlot$TotalCodingRegion, TotalNonCodingRegion=toPlot$TotalNonCodingRegion)
dfTypes <- melt(dfTypes, id.vars=c("Cancer"))
colnames(dfTypes) <- c("Cancer","Region","SNVs")
dfTypes$SNVs <- dfTypes$SNVs/1000000

p2 <- ggplot(dfTypes, aes(x=Cancer, y=SNVs, color=Region, group=interaction(SNVs,Region))) +  geom_point(size=0.4, position = position_dodge(width=0.75)) +
  scale_x_discrete(limits=customOrder) + theme_minimal() + scale_color_brewer(labels = c("Coding Region", "Non-Coding Region"), palette="Set2") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position ="top", legend.title=element_blank(), panel.grid.minor.y = element_blank()) +
  ylab("SNVs per Mb") + xlab("")

grid.arrange(p1,p2)
