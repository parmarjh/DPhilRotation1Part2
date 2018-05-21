library(reshape2)
library(ggplot2)
# df <- read.csv('/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/MutagenesisExtractionAndProcess/Predictions/Diffs.txt', header=F,sep='\t', stringsAsFactors = F)
# 
# setwd('/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/MutagenesisExtractionAndProcess/Predictions/')
# 
# df <- df[order(-df$V4),]
# 
# lowest <- df[seq(length(df$V1)-1000,length(df$V1),1),]
# colnames(lowest) <- c("Site","WT","SNV","Diff")
# lowest$Site <- as.factor((lowest$Site))
# lowest <- lowest[order(lowest$WT),]
# dfLowest = melt(lowest, id.vars=c('Site','Diff'))
# 
# saveRDS(dfLowest, file='Diff.lowest.rds')
# 
# highest <- df[seq(1,1000),]
# colnames(highest) <- c("Site","WT","SNV","Diff")
# highest$Site <- as.factor((highest$Site))
# highest <- highest[order(highest$WT),]
# dfhighest = melt(highest, id.vars=c('Site','Diff'))
# 
# saveRDS(dfhighest, file='Diff.highest.rds')
# 
# dfOut <- rbind(dfLowest, dfhighest)
# 
# saveRDS(dfOut, file='Diff.highestANDlowest.rds')



dfHighLow <- readRDS('Diff.highestANDlowest.rds')

h <- ggplot(dfHighLow, aes(x=reorder(Site, -value),y=value, colour=variable, group=Site)) + 
  geom_line(color="black", size=0.1) + geom_point(size=0.5) +
  scale_color_manual(values=c('black','red')) + 
  theme_minimal() + xlab("Genomic Position") + ylab("Predicted Accessibility") + coord_flip() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y=element_blank()) + 
  scale_y_continuous(limits = c(0,1.0))

write.csv(dfHighLow, file='./DiffsHighLow', row.names=FALSE)


vaf_v_diff <- read.csv('/Users/schencro/Desktop/Oxford/Rotation_1/PCAWGArm/MutagenesisExtractionAndProcess/Predictions/COAD.matches.txt',header=F,stringsAsFactors = F, sep='\t')
colnames(vaf_v_diff) <- c('Site','WT','SNV','Diff','Sample','VAF','Class')
vaf_v_diff$VAF

ggplot(vaf_v_diff, aes(x=VAF,y=Diff)) + geom_point()

toKMeans = data.frame(Diff=vaf_v_diff$Diff, VAF=vaf_v_diff$VAF)
fit <- kmeans(toKMeans, centers = 4)

toKMeans$Cluster = as.factor(fit$cluster)
pos <- subset(toKMeans, toKMeans$Diff>0)
pos$Val <- rep("Pos.Diff",length(pos$Diff))
neg <- subset(toKMeans, toKMeans$Diff<0)
neg$Val <- rep("Neg.Diff",length(neg$Diff))

toKMeans <- rbind(pos, neg)

ggplot(toKMeans, aes(x=VAF,y=Diff, colour=Cluster)) + geom_point() + theme_bw() + xlab("Variant Allele Frequency") +
  ylab('Predicted Difference in Accessibility') + scale_color_brewer(palette = "Set2")

ggplot(toKMeans, aes(VAF, fill=Val, colour=Val, alpha=0.05)) + geom_density() + scale_color_brewer(palette = "Set1") +  scale_fill_brewer(palette = "Set1") +
  xlab("Variant Allele Frequency") +
  ylab('Density') + theme_bw()
  


