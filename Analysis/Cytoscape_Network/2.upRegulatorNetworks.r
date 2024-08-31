d7tf <- c("IMMUNOGLOBULIN (complex)", "INTERFERON ALPHA (family)", "IL4", "IL15", "STING1", "STAT3", "TBX21", "TNF", "IL18", "IL2", "IL21", "STAT1", "IFNG")
df <- read.table("D7_Ag_upstreamRegulators_LW.txt", sep="\t", header=T)
genes <- df$overlapGenes
totaltargets <- c()
write("networkClusterD7.sif", append=F)
for (i in d7tf) {
    x = strsplit(df[df$Upstream.Regulator==i, "Target.Molecules.in.Dataset"], "\\,")
    for (g in x[[1]]) {
        totaltargets <- c(totaltargets, g)
        write(paste(g, "pd", i,sep = "\t"), "networkClusterD7.sif", append = "T")
    }
}

d750tf <- c("IMMUNOGLOBULIN (complex)", "INTERFERON ALPHA (family)", "IL4", "IL15", "STAT3", "TBX21", "TNF", "IL18", "IL2", "IL21", "STAT1", "IFNG")
df <- read.table("D750_Ag_upstreamRegulators_LW.txt", sep = "\t", header = T)
genes <- df$overlapGenes
write("networkClusterD750.sif", append = F)
for (i in d7tf) {
    x <- strsplit(df[df$Upstream.Regulator == i, "Target.Molecules.in.Dataset"], "\\,")
    for (g in x[[1]]) {
        totaltargets <- c(totaltargets, g)
        write(paste(g, "pd", i, sep = "\t"), "networkClusterD750.sif", append = "T")
    }
}

totaltargets <- unique(totaltargets)
lfcs <- read.csv("LFC_Padj_HGNC_Sig.csv", header=TRUE)
upstreamregs <- unique(c(d7tf, d750tf))
overlap <- intersect(lfcs$HGNC.symbol, c(upstreamregs, totaltargets))
lfcssub <- lfcs[lfcs$HGNC.symbol %in% overlap, ]
write.csv(lfcssub, "LFCs_cytoscape_subset.csv")

df <- read.table("networkClusterD7.sif", sep="\t")
d7x <- table(df$V1)
total <- length(d7x)
d7x1 <- d7x[d7x > 2]
print(d7x1)
totalabove1 <- length(d7x1)
message("Of ", total, " target genes ", totalabove1, " have more than one upstream regulator: percentage is ", totalabove1/total )

df <- read.table("networkClusterD750.sif", sep = "\t")
d750x <- table(df$V1)
total <- length(x)
d750x1 <- d750x[d750x > 2]
print(d750x1)
totalabove1 <- length(d750x1)
message("Of ", total, " target genes ", totalabove1, " have more than one upstream regulator: percentage is ", totalabove1 / total)