#R codes below were modified using Leanne Whitmore's codes.
library(ggplot2)
library(data.table)

#Theme:
theme_minimal_LW <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.3, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

#Plot:
df <- read.csv("CBC_Data.csv", header=T)
targetfile <- read.csv("Metadata.csv", header=T)
Animal <- unique(targetfile$animal)
df <- df[df$RM %in% Animal, ]
dfsub <- df[,c("Day", "RM", "BasoPCNT", "EosinoPCNT", "LymphoPCNT", "MonoPCNT", "NeutroPCNT")]
dfmelt <- data.table::melt(dfsub, id.vars=c("Day", "RM"))

vaccine <- c()
for (i in dfmelt$RM) {
    x <- targetfile[targetfile$animal==i, "Vaccine"][1]
    vaccine <- c(vaccine, x)
}
dfmelt$Vaccine <- vaccine
dfmelt$Day <- as.character(dfmelt$Day)
dfmelt$Day <- factor(dfmelt$Day, levels=c("0", "14", "98", "101", "105", 
    "112", "560", "742", "746", "757", "1177", "1191", "1219"))
ggplot(dfmelt, aes(x=Day, y=value, group=RM, color=Vaccine)) +
    theme_minimal_LW() + facet_wrap(~variable, scales = "free_y",ncol=3) +
    scale_color_manual(values=c("68_1_RhCMV_TB_6Ag"="#000000", "68_1_RhCMV_pp71_TB_6Ag"="#FF0000")) +
    stat_summary(aes(y = value, group = Vaccine), fun = "mean", geom = "line") +
    stat_summary(aes(y = value, group = Vaccine), fun.data = "mean_se", 
                    geom = "errorbar", width=0.3) +
    theme(strip.text.x = element_text(size = 8)) +
    theme(legend.position = "top", legend.direction = "horizontal",
        strip.background = element_rect(colour = "black", fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)) + labs(x="Day Post-Vaccination", y="CBC in Percentage (%)")

#Save as a 8 inch x 5 inch pdf plot (also a 800x500 png plot).
#ggsave("cbc.png", width=7.5, height=7, dpi=300).