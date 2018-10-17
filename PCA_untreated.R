# 1. Tell R which packages you need to use

library(edgeR)
library(limma)
library(WGCNA)
library(ggplot2)
library(readr)
library(rafalib)
library(devtools)
library(plyr)

# 2. Upload data *replace ###### with path to file "all_patient_count.txt"

vivax_counts <- read_delim("~/XXXXXX/all_patient_count.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)

# 3. Remove unrelated data
RelapseDrops <- c("AR_1","BR_1","CR_1","DR_1","ER_1","FR_1","KR_0","KR_1", "X59")
vivax_counts <- vivax_counts[ , !(names(vivax_counts) %in% RelapseDrops)]

CqDrops <- c("A_1","B_1","C_1","D_1","E_1","F_1","G_1","H_1","I_1","J_1","K_1","KR_1","L_1","M_1","N_1","O_1","P_1","Q_1","S_1","T_1","U_1", "X59")
vivax_counts <- vivax_counts[ , !(names(vivax_counts) %in% CqDrops)]

# 3. Convert counts to TPM (transcripts per million)

vivax_counts_temp <- vivax_counts[,1:4]
Gene_size <- c(vivax_counts[,4])
vivax_counts_temp[,c(5:30)] <- 1000*vivax_counts[, c(5:30)]/Gene_size
vivax_counts_temp[,5:30] <- t(1000000*t(vivax_counts_temp[,5:30])/rowSums(t(vivax_counts_temp[,5:30])))  #Transpose twice to divide each element by column sum

# 4. Remove genes with low expression

keep <- rowSums(cpm(vivax_counts_temp[,5:30])>20) >= 6
vivax_counts_temp <- vivax_counts_temp[keep,]

colnames(vivax_counts_temp) <- colnames(vivax_counts)
Gene_size <- c(vivax_counts[,4])

vivax_counts <- vivax_counts_temp

# 5. Calculate principle coordinates

pcadata <- data.frame(vivax_counts[,5:30])
pcadata2 <- t(pcadata)
pcadata2 <- pcadata2[,colSums(pcadata2)>0]
pcaOut <- prcomp(pcadata2, scale = T, center = F)
summary(pcaOut)

explained <- as.list((pcaOut$sdev)^2 / sum(pcaOut$sdev^2))
names(explained) <- colnames(pcaOut$x)
print(plot(as.numeric(explained[1:10]), main = "baseName"))
PC <- as.data.frame(pcaOut$x)

# 6. Write file with all PC values

write.table(PC, na="NA", file="20160614_Samples_PC")

# 7. Create figures

PC$label <- rownames(pcadata2)
PC$otherLabel <- c("0","R","0","R","0","R","0","R","0","R","0","R","0","0","0","0","0","R","0","0","0","0","0","0","0","0")
PC$sample <- c("A","AR","B","BR","C","CR","D","DR","E","ER","F","FR","G","H","I","J","K","L","M","N","O","P","Q","S","T","U")

# 8. PCA of all samples untreated
# (Supplemental Figure 2)

ggplot(PC, aes(x = PC1, y = PC2, label=sample)) + 
  geom_point(size = 2, aes()) + 
  #geom_point(size = 2, aes(colour=as.factor(otherLabel))) + 
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
  labs(list(title = "PCA of samples")) +
  theme(legend.position=c(.9, .9))

summary(pcat)
