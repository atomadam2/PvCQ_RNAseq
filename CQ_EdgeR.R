# 1. Tell R which packages you need to use

library(edgeR)
library(limma)
library(ggplot2)
library(psych)
library(statmod)

install.packages("devtools")
library(devtools) # get from CRAN with install.packages("devtools")
install_github("ririzarr/rafalib")
library(rafalib)

# 2. Upload data *replace ###### with path to file "all_patient_count_subset.txt"

vivax_counts <- read_delim("~/XXXXXX/all_patient_count_subset.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)

# 3. Remove unrelated data

Drops <- c("X47", "KR_0", "KR_1")
vivax_counts <- vivax_counts[ , !(names(vivax_counts) %in% Drops)]

# 4. Remove genes with low expression

vivax_counts <- vivax_counts[-which(apply(vivax_counts[,5:44], 1, function(x)sum(x=="0"))>=30),]

# 5. Preliminary MDS plot of data
y <- DGEList(counts=vivax_counts[,5:44], genes=vivax_counts[,1:4])
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]

nrow(y)
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
y$samples
plotMDS(y, col=c(rep("black",1), rep("red",1)))

# 6. Create Design Matrix

Drug <- factor(c(rep(c(0,8),20)))
Sample <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20))
data.frame(Sample=colnames(y),Sample,Drug)
design <- model.matrix(~Sample+Drug)
rownames(design) <- colnames(y)
design

# 7. Estimate Dispersion

y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion
plotBCV(y)

# 8. Calculate Differentially Expressed Genes

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)
colnames(design)
o <- order(lrt$table$PValue)
cpm(y)[o[1:10],]
summary(de <- decideTestsDGE(lrt))
detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1,1), col="blue")

# 9. Write file with all differentially expressed genes
# (Supplementary Table 4)
write.table(topTags(lrt, n=25000), file="20160614_CQ_DE")

# 10. Volcano plot
# (Figure 3C)

df_allhits <- data.frame(topTags(lrt, n=15000))
df_allhits$threshold <- as.factor(abs( df_allhits$FDR < 0.1))
#df_allhits$threshold <- as.factor(abs(df_allhits$logFC) > 2 & df_allhits$FDR < 0.05)
volcano1 <- ggplot(df_allhits, aes(logFC, -log10(PValue), colour=threshold))
volcano1 + geom_point(size=1)+
  #xlim(-10, 10) +                          # x axis
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  #  geom_vline(aes(x = 0))+
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
  labs(list(title = "Volcano Plot of Gene Expression Changes", x = "log2 fold change", y = "-log10 p-value"))

