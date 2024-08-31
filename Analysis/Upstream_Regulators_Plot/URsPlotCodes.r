#Below are codes to make the plot for the IPA upstream regulators analysis results.
library(ggplot2)
library(RColorBrewer)
data <- read.csv("TopSharedURs.csv")

#Plot:
ggplot(data, aes(x = Abbre, y = reorder(UpstreamRegulators, abs(Zscore)), color = Zscore, size = Pvalue)) +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.15)), axis.title = element_text(size=rel(1.15))) +
  scale_x_discrete(limits=c("T_3", "T_7", "T_14", "T_98", "T_101", "T_105", "T_112", "T_560", "T_743", "T_746", "T_750", "T_757", "T_1177", "T_1191", "T_1219", "P_3", "P_7", "P_14", "P_98", "P_101", "P_105", "P_112", "P_560", "P_743", "P_746", "P_750", "P_757", "P_1177", "P_1191", "P_1219")) +
  xlab("Vaccine_TimePoints") +
  ylab("Upstream Regulators") +
  ggtitle("Top Upstream Regulators") +
  theme(plot.title = element_text(hjust=0.5, face = "bold")) +
  scale_color_gradient2(name = "Activation \n z-score", midpoint = 0, high = "red3", mid = "white", low = "blue3") +
  scale_size(name = "Significance \n -log10(p-value)") +
  theme(legend.title = element_text(size=rel(1.15), hjust=0.5, face="bold"))

#Save as png file (Width = 800 and Height = 400) and pdf file (8.00 and 4.00 inches).