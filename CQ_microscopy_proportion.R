# 1. Tell R which packages you need to use

library(ggplot2)
library(reshape2)
library(readr)

# 2. Upload data *replace ###### with path to file "parasitemia_reads.txt"

Para_Prop <- read_delim("~/XXXXXX/parasitemia_reads.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# 3. Load melted data into variable p for ggplot

finaldat <- melt(Para_Prop)
p <- ggplot(finaldat, aes(variable, value))

colnames(Para_Prop) <- c("Sample","Microscopy","Proportion")

# 4. Boxplot of effect of chloroquine on parasites as measured by microscopy and proportion of reads
# (Figure 3A)

p+geom_boxplot(outlier.shape = NA) +
  theme_bw() + 
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  ylim(0, 100) +
  scale_fill_grey(start = 0.8, end = 1)+
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "Chloroquine effect on Parasitemia and RNA", y="Percent")) +
  geom_jitter( alpha=0.9, 
               position=position_jitter(w=0.1,h=0.1)) 

# 5. Scatterplot of effect of chloroquine on parasites as measured by microscopy compared to proportion of reads
# (Supplementary Figure 4)

colnames(Para_Prop) <- c("Sample","Microscopy","Proportion")

Para_Prop$sample <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","S","T","U")
ggplot(Para_Prop, aes(x = Microscopy, y = Proportion, label=sample)) + 
  geom_point(size = 2) + 
  #geom_text( nudge_y = 8, ) +
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) + 
  scale_x_continuous(name="Microcopy", limits=c(0, 100)) +
  scale_y_continuous(name="Proportion of Reads", limits=c(0, 100))+
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "Chloroquine effect on Parasitemia and RNA")) +
  theme(legend.position=c(.8, .2)) 


