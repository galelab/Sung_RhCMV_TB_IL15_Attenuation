#Make a proportional Venn diagram.
library(eulerr)

fit <- euler(c(IL15 = 55, TB6Ag = 0, pp71 = 0, "IL15&TB6Ag" = 2, "IL15&pp71" = 60, "TB6Ag&pp71" = 0, "IL15&TB6Ag&pp71" = 69))
plot(fit, fills = c("lightskyblue1", "palegreen3", "khaki"), edges = c("cornflowerblue", "palegreen4", "khaki4"), legend = FALSE, labels = identical(legend, FALSE))

#Export (save as image) into a png file with Width = 1200 and Height = 1200.