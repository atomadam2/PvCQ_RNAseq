# 1. Install necessary packages
# Leave these commented out (#) unless you need to install or reinstall

install.packages("installr")
install.packages("ggplot2")
install.packages("rlang")
install.packages("ggdendro")
install.packages("heatmaply")
install.packages("plotly")
install.packages("viridis")
install.packages("viridisLite")
install.packages("install_github")
install_github('hadley/ggplot2')
install.packages("pvclust")

# 2. Tell R which packages you need to use

library(ggplot2)
library(ggdendro)
library(plyr)
library(reshape2)
library(scales)
library(viridis)
library(viridisLite)
library(readr)
library(heatmaply)
library(pvclust)


# 3. Upload data *replace ###### with path to file "Gametocyte_Corr.txt"

targets <- read_delim("~/XXXXXX/Gametocyte_Corr.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# 4. Simplify data

targets.a <- targets[,c(1:22)]
targets.l <- (targets.a[,2:22])
row.names(targets.l) <- targets.a$X1

# 5. pvclust for significance of 2 major groups

result <- pvclust(data=targets.l, method.hclust = "single", method.dist = "euclidean")
plot(result)

# 6. Create dendrogram using hclust

hc <- hclust(dist(targets.l), "ave")
ggdendrogram(hc, rotate = FALSE, size = 2)

d <- dist(x=targets.a[,2:3], method = "euclidean")
targets.t  <- hclust(d, method = "single")
dendcompletem <- as.dendrogram(targets.t)

targets.m <- melt(targets.a)

x <- as.matrix(scale(targets.l))
dd.col <- as.dendrogram(hclust(dist(targets.l)))
dd.row <- as.dendrogram(hclust(dist(t(targets.l))))
dx <- dendro_data(dd.row)
dy <- dendro_data(dd.col)

ggdend <- function(df) {
  ggplot() +
    geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
    labs(x = "", y = "") + theme_minimal() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

# 7. x and y axis dendrogram

px <- ggdend(dx$segments)
py <- ggdend(dy$segments) + coord_flip()
px
py

# 8. Heatmap
#(Figure 2B)

col.ord <- order.dendrogram(dd.col)
row.ord <- order.dendrogram(dd.row)
xx <- scale(targets.l)[col.ord, row.ord]
xxb <- data.matrix(targets.l[col.ord, row.ord])
xx_names <- attr(xx, "dimnames")
df <- as.data.frame(xxb)
colnames(df) <- xx_names[[2]]
df$gene <- xx_names[[1]]
df$gene <- with(df, factor(gene, levels=gene, ordered=TRUE))
mdf <- reshape2::melt(df, id.vars="gene")
p <- ggplot(mdf, aes(x = variable, y = gene)) + geom_tile(aes(fill = value)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_fill_gradientn(colours = c("white", "lightpink",  "red"), values = c(-1,0,1))
p
