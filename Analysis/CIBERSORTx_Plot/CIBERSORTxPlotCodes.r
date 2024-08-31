#R codes below were modified using Leanne Whitmore's codes.
library(ggplot2)
library(data.table)

#Functions.
load_cibersort_data <- function(filename) {
    df <- read.csv(filename,
        header = TRUE, as.is = TRUE,
        check.names = FALSE, stringsAsFactors = FALSE
    )

    df <- df[,1:(length(colnames(df))-3)]
}

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

calculate_stats <- function(df, cell, time, group, fileoutput) {
    gr1 <- df[(df$variable == cell & df$Time_Point == "D0" & df$Vaccine == group), "value"]
    gr2 <- df[(df$variable == cell & df$Time_Point == time & df$Vaccine == group), "value"]

    if ((length(gr1) == 0) || (length(gr2) == 0)) {
        message("STATUS: gr1 has ", length(gr1), " samples and gr2 has ", length(gr2), " samples")
    } else {
        wil <- wilcox.test(gr2, gr1)
        ttest <- t.test(gr2, gr1)
        write(
            paste0(
                cell, ",", time, ",", mean(gr2), ",", mean(gr1), ",",
                wil$statistic, ",", wil$p.value, ",", group
            ),
            file.path(paste0(fileoutput, "wilcox.csv")),
            append = TRUE
        )
        write(
            paste0(
                cell, ",", time, ",", mean(gr2), ",", mean(gr1), ",",
                ttest$statistic, ",", ttest$p.value, ",", group
            ),
            file.path(paste0(fileoutput, "ttest.csv")),
            append = TRUE
        )
    }
}

#Load data.
target <- read.csv("Metadata.csv", sep = ",", row.names=1,
        as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
cibfile <- "CIBERSORTx_Results.csv"
cib <- load_cibersort_data(cibfile)
bcells <- rowSums(cib[, c("B cells naive", "B cells memory")])
CD4cells <- rowSums(cib[, c("T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated")])
Nkcells <- rowSums(cib[, c("NK cells resting", "NK cells activated")])
dendriticcells <- rowSums(cib[, c("Dendritic cells resting", "Dendritic cells activated")])
dfsub <- cbind(CD4cells, cib[, "T cells CD8"], Nkcells, bcells, cib[, "Monocytes"], dendriticcells, cib[, "Neutrophils"])
rownames(dfsub) <- cib$Mixture
colnames(dfsub) <- c("CD4 cells", "CD8 T cells", "NK cells", 
    "B cells", "Monocytes", "Dendritic cells", "Neutrophils")
totaltarget <- merge(dfsub, target[, c("Time_Point", "Vaccine")],
     by.x="row.names", by.y="row.names")
totaltargetmelt <- reshape2::melt(totaltarget)

#Figure.
totaltargetmelt$Time_Point <- factor(totaltargetmelt$Time_Point, levels = c("D0", "D3", "D7", "D14", 
    "D98", "D101", "D105", "D112", "D560", "D743", "D746", "D750", "D757", "D1177", "D1191", "D1219"))
totaltargetmelt$value <- totaltargetmelt$value*100
ggplot(totaltargetmelt, aes(x=Time_Point, y=value, group = Vaccine, color=Vaccine)) +
    facet_wrap(~variable, scales = "free_y", ncol=4) + theme_minimal_LW() + theme(strip.text.x = element_text(size = 5)) +
    theme(axis.text.y=element_text(size=6), axis.text.x=element_text(size=6)) + labs(y="mean % abundance") +
    stat_summary(aes(y = value, group = Vaccine), fun = "mean", geom = "line") +
    stat_summary(aes(y = value, group = Vaccine), fun.data = "mean_se", 
                    geom = "errorbar", width=0.3) +
    scale_color_manual(values=c("68_1_RhCMV_TB_6Ag"="black", "68_1_RhCMV_pp71_TB_6Ag"="red")) +
    theme(axis.title.x=element_text(size=8), axis.title.y=element_text(size=8)) +
    theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_text(size=7),
        legend.text = element_text(size=7),strip.background = element_rect(colour = "black", fill = "white"),
        strip.text.x = element_text(size = 8), axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7))
#Save as a 8 inch x 5 inch pdf plot (also a 800x500 png plot).
#ggsave("CibersortFig.png", width = 7.5, height = 4.5)
#ggsave("CibersortFig.pdf", width = 7.5, height = 4.5)
#ggsave("CibersortFig.svg", width = 7.5, height = 4.5)
fileoutput <- "stats4"
write(
    "cell,time,mean gr2,BLN,wilcox.statistic,pvalue,group",
    paste0(fileoutput, "wilcox.csv"))
write(
    "cell,time,mean gr2,BLN,t.statistic,pvalue,group",
    paste0(fileoutput, "ttest.csv")
)

for (group in unique(totaltargetmelt$Vaccine)) {
    for (cell in unique(totaltargetmelt$variable)) {
        for (time in unique(totaltargetmelt$Time_Point)) {
            if (time!="D0") {
                calculate_stats(totaltargetmelt, cell = cell, 
                    time = time, group = group, fileoutput = fileoutput)
            }
        }
    }
}

dfttest <- read.csv("stats4ttest.csv", header = T, row.names = NULL)
dfttest <- dfttest[!is.na(dfttest$pvalue), ]
adjp <- c()
for (g in unique(dfttest$group)) {
    dfsub <- dfttest[dfttest$group==g, ]
    adj <- p.adjust(dfsub$pvalue, method = "BH")
    adjp <- c(adjp, adj)
}
dfttest$adj.fdr.pval <- adjp
write.csv(dfttest, "stats4ttestadj.csv", row.names = F, quote = F)

dfwil <- read.csv("stats4wilcox.csv", header = T, row.names = NULL)
dfwil <- dfwil[!is.na(dfwil$pvalue), ]
adjp <- c()
for (g in unique(dfwil$group)) {
    dfsub <- dfwil[dfwil$group == g, ]
    adj <- p.adjust(dfsub$pvalue, method = "BH")
    adjp <- c(adjp, adj)
}
dfwil$adj.fdr.pval <- adjp
write.csv(dfwil, "stats4wilcoxadj.csv", row.names = F, quote = F)