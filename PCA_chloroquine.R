# 1. Tell R which packages you need to use

library(edgeR)
library(limma)
library(WGCNA)
library(ggplot2)
library(readr)
library(rafalib)
library(devtools)
library(plyr)

# 2. Upload data *replace ###### with path to file "all_patient_count_subset.txt"

vivax_counts <- read_delim("~/XXXXXX/all_patient_count_subset.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)

# 3. Convert counts to TPM (transcripts per million)

vivax_counts_temp <- vivax_counts[,1:4]
Gene_size <- c(vivax_counts[,4])
vivax_counts_temp[,c(5:44)] <- 1000*vivax_counts[, c(5:44)]/Gene_size
#Transpose twice to divide each element by column sum
vivax_counts_temp[,5:44] <- t(1000000*t(vivax_counts_temp[,5:44])/rowSums(t(vivax_counts_temp[,5:44]))) 

# 4. Remove genes with low expression

keep <- rowSums(cpm(vivax_counts_temp[,5:44])>4) >= 10
vivax_counts_temp <- vivax_counts_temp[keep,]

colnames(vivax_counts_temp) <- colnames(vivax_counts)
Gene_size <- c(vivax_counts[,4])

vivax_counts <- vivax_counts_temp

# 4. Calculate principle coordinates

pcadata <- data.frame(vivax_counts[,5:44])

pcadata2 <- t(pcadata)
pcadata2 <- pcadata2[,colSums(pcadata2)>0]
pcaOut <- prcomp(pcadata2, scale = T, center = F)
summary(pcaOut)

explained <- as.list((pcaOut$sdev)^2 / sum(pcaOut$sdev^2))
names(explained) <- colnames(pcaOut$x)
print(plot(as.numeric(explained[1:10]), main = "baseName"))
PC <- as.data.frame(pcaOut$x)

# 5. Write file with all PC values
write.table(PC, na="NA", file="20160614_CQ_PC")

# 6. Create figures

PC$label <- rownames(pcadata2)
PC$sample <- c("A","A","B","B","C","C","D","D","E","E","F","F","G","G","H","H","I","I","J","J","K","K","L","L","M","M","N","N","O","O","P","P","Q","Q","S","S","T","T","U","U")
PC$otherLabel <- c("Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ","Before CQ","After CQ")

# 7. PCA of all samples without labels
# (Figure 3B)

ggplot(PC, aes(x = PC1, y = PC2, label=sample, group = sample)) + 
  geom_point(size = 3, aes(colour=as.factor(otherLabel))) + 
  xlab(paste("PC1", " (", round(as.numeric(explained["PC1"]) * 100, 2), "%)", sep = "")) +
  ylab(paste("PC2", " (", round(as.numeric(explained["PC2"]) * 100, 2), "%)", sep = "")) +
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_color_manual(values = c( "#FF3300" , "#0000CC")) +
  scale_fill_manual(values = c( "#FF3300" , "#0000CC")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "PCA of samples before and after CQ")) +
  theme(legend.position=c(0.2,0.8))

# 8. PCA of all samples with labels
# (Unpublished)

ggplot(PC, aes(x = PC1, y = PC2, label=sample, group = sample)) + 
  geom_point(size = 2, aes(colour=as.factor(otherLabel))) + 
  geom_segment(data=allPC, arrow = arrow(angle = 15, length = unit(0.15, "inches"), type = "closed"), size = 1, aes(x=start_PC1, y=start_PC2, xend=end_PC1, yend=end_PC2)) +
  geom_text( nudge_y = 2) +
  xlab(paste("PC1", " (", round(as.numeric(explained["PC1"]) * 100, 2), "%)", sep = "")) +
  ylab(paste("PC2", " (", round(as.numeric(explained["PC2"]) * 100, 2), "%)", sep = "")) +
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_color_manual(values = c("#0000CC", "#FF3300" )) +
  scale_fill_manual(values = c("#0000CC", "#FF3300")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "PCA of samples before and after CQ")) +
  theme(legend.position="none")

# 9. Create dataframe to point arrows in correct direction (Before CQ -> After Chloroquine)

PC$otherLabel <- c("0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8","0","8")
testsp <- split(PC, PC$otherLabel)
testsp
df0 <- testsp[[1]]
colnames(df0) <- paste("start", colnames(df0), sep = "_")
colnames(df0)[42] <- "sample"
df8 <- testsp[[2]]
colnames(df8) <- paste("end", colnames(df8), sep = "_")
colnames(df8)[42] <- "sample"
allPC <- merge(df0, df8, by = "sample")

# 10. PCA of all samples with arrows (Before CQ -> After Chloroquine)
# (Supplemental Figure 6)

ggplot(PC, aes(x = PC1, y = PC2, label=sample, group = sample)) + 
  geom_point(size = 2, aes(colour=as.factor(otherLabel))) + 
  geom_segment(data=allPC, arrow = arrow(angle = 15, length = unit(0.15, "inches"), type = "closed"), size = 1, aes(x=start_PC1, y=start_PC2, xend=end_PC1, yend=end_PC2)) +
  geom_text(aes(label=sample, colour=otherLabel), nudge_y = 2) +
  xlab(paste("PC1", " (", round(as.numeric(explained["PC1"]) * 100, 2), "%)", sep = "")) +
  ylab(paste("PC2", " (", round(as.numeric(explained["PC2"]) * 100, 2), "%)", sep = "")) +
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_color_manual(values = c("black", "white" )) +
  scale_fill_manual(values = c("#0000CC", "#FF3300")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "PCA of samples before and after CQ")) +
  theme(legend.position="none")

# 11. PCA of all samples with faded arrows (Before CQ -> After Chloroquine)
# (Unpublished)

ggplot(PC, aes(x = PC1, y = PC2, label=sample, group = sample)) + 
  geom_point(size = 2, aes(colour=as.factor(otherLabel))) + 
  geom_segment(data=allPC, arrow = arrow(angle = 15, length = unit(0.15, "inches"), type = "closed"), size = 1, alpha = 0.5, aes(x=start_PC1, y=start_PC2, xend=end_PC1, yend=end_PC2)) +
  geom_text(nudge_y = 2) +
  xlab(paste("PC1", " (", round(as.numeric(explained["PC1"]) * 100, 2), "%)", sep = "")) +
  ylab(paste("PC2", " (", round(as.numeric(explained["PC2"]) * 100, 2), "%)", sep = "")) +
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  scale_color_manual(values = c("black", "black" )) +
  scale_fill_manual(values = c("#0000CC", "#FF3300")) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "PCA of samples before and after CQ")) +
  theme(legend.position="none")

summary(pcat)

