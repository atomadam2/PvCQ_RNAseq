# 1. Tell R which packages you need to use

library(ggplot2)

# 2. Upload data *replace ###### with path to file "Pv48_Pv47.txt"

Gametocyte_Scatter <- read_delim("~/XXXXXX/Pv48_Pv47.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# 3. Scatterplot of all Pv48 vs Pv47 expression in all samples
# (Figure 2A)

ggplot(Gametocyte_Scatter, aes(x = Pv48, y = Pv47, label=Name)) + 
  labs(x = "Pvs48/45 (TPM)", y = "Pv47 (TPM)") +
  geom_point(size = 2) + 
  #geom_text( nudge_y = 8) +
  theme_bw() +                             # theme with white background
  theme(                                 # eliminates background, gridlines, and chart border
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(axis.line = element_line(color = 'black', #draws x and y axis line
                                 size = 1, 
                                 linetype = "solid"))  +
  theme(axis.text.x = element_text(face="bold", color="#000000", size=14),
        axis.text.y = element_text(face="bold", color="#000000", size=14),
        legend.text = element_text( size = 14)) +
  theme(axis.title=element_text(size=20), plot.title = element_text(size=24), 
        legend.title=element_text(size=16)) + 
  labs(list(title = "Pv47 vs Pv48/45")) +
  theme(legend.position=c(.8, .2)) 
