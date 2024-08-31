#Below are codes for making the overlap column to show genes overlapping IL-15 signature.
library(ggplot2)
df <- read.table("OverlapCol.txt", header = T)

ggplot(df, aes(x = X, y = Genes, fill = YorN)) +
    geom_tile() +
    scale_fill_gradient2(low = "white", high = "black")

#Manually export (save as image) into a png file with Width = 900 and Height = 9000.