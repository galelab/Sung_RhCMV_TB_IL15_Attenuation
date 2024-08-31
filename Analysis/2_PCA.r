#s2 feature reduction.
#This function does feature reduction and visualization on normalized voom counts.

library(edgeR)
library(factoextra)
library(umap)
library(Rtsne)
library(randomcoloR)
library(Polychrome)
library(ggplot2)

s2_feature_reduction <- function(countfile, targetfile,
                                 target_columns = c(4, 6),
                                 pca = TRUE, UMAP = FALSE, tsne = FALSE,
                                 base_file_name = "vnorm.png",
                                 figres = 300) {
    pdf(NULL)
    files <- loadfiles(count_file = countfile, target_file = targetfile)
    results_path <- generate_folder("s2_feature_reduction_results")

    #print("STATUS: Running MDS feature reduction")
    #visualize_feature_reduction_data(files$counts,
    #                                 files$targets[,target_columns[1]],
    #                                 results_path, base_file_name,
    #                                 figres)

    if (all.equal(rownames(files$targets), colnames(files$counts)) != TRUE) {
        print("WARNING: order of samples in target and count file is not the same, this needs fixing before code can proceed")
    } else {
        if (pca == TRUE) {
            print("STATUS: Running PCA feature reduction")
            pca_fun(files$counts, files$target,
                    results_path, base_file_name,
                    target_columns, 300, 3)
        }
        if (UMAP == TRUE) {
            print("STATUS: Running UMAP feature reduction")
            umap_fun(files$counts, files$target,
                     results_path, base_file_name,
                     target_columns, figres)
        }
        if (tsne == TRUE) {
            print("STATUS: Running tSNE feature reduction")
            tsne_fun(files$counts, files$target,
                     results_path, base_file_name,
                     target_columns, figres)
        }
    }
}

visualize_feature_reduction_data <- function(data,
                                            labels,
                                            results_path,
                                            base_file_name,
                                            figres = 100) {
    #MDS (multidimensional scaling) uses log fold changes between genes as distances.
    MDS <- plotMDS(data, gene.selection = "pairwise", cex = 0.8)
    minx <- min(MDS$x)
    maxx <- max(MDS$x)
    miny <- min(MDS$y)
    maxy <- max(MDS$y)
    png(file.path(results_path,
        paste0("2.mds_", base_file_name)),
        res = figres)
    plot(MDS$x, MDS$y, cex = 1, xlim = c(minx - 1, maxx + 1),
         ylim = c(miny - 1, maxy + 1),
         xlab = paste0(MDS$axislabel, " 1"),
         ylab = paste0(MDS$axislabel, " 2"), frame = FALSE)
    text(MDS$x, MDS$y, labels, cex = 0.6, pos = 4)
    dev.off()
}

pca_fun <- function(exprs, labels, results_path,
                    base_file_name, target_columns,
                    figres = 100, size = 1, pca=FALSE, legend = "right") {
    # Run PCA/SVD reduction
    if (isFALSE(pca)) {
        pca <- prcomp(t(exprs))
    }
    E <- get_eig(pca)
    cx <- sweep(t(exprs), 2, colMeans(t(exprs)), "-")
    sv <- svd(cx)


    visualize_pca(
        file.path(results_path, paste0("2.svd_", base_file_name)),
        sv$u, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, E, size, legend
    )
    visualize_pca(
        file.path(results_path, paste0("2.pca_", base_file_name)),
        pca$x, labels[, target_columns[1]],
        labels[, target_columns[2]],
        figres, E, size, legend
    )
    visualize_scree_plot(
        file.path(
            results_path,
            paste0("scree_", base_file_name)
        ), pca, figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    is_pc1_0 <- loadingscores$PC1 > 0
    is_pc2_0 <- loadingscores$PC2 > 0

    loadingscores <- loadingscores[is_pc1_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC1)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc1", base_file_name, ".txt")),
        loadingscores["PC1"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc2_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC2)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc2", base_file_name, ".txt")),
        loadingscores["PC2"], figres
    )
    return(pca)
}

visualize_pca <- function(plot_file, PCA, class1, class2, figres, E, size, legend) {
    # Visualize PCA  results
    library(Polychrome)
    minx <- min(PCA[, 1])
    maxx <- max(PCA[, 1])
    miny <- min(PCA[, 2])
    maxy <- max(PCA[, 2])

    #P36 <- colorRampPalette(c("greenyellow", "darkgreen"))(length(levels(factor(class2))))
    #P36 <- distinctColorPalette(length(levels(factor(class2))))
    #P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
    P36 <- createPalette(length(levels(factor(class2))), c("#d73131", "#1437ff"))
    Vaccine_levels = c("68_1_RhCMV_TB_6Ag", "68_1_RhCMV_pp71_TB_6Ag")
    qplot(PCA[, 1], PCA[, 2], color = factor(class2, levels = Vaccine_levels), shape = factor(class1), size = I(size)) +
        theme_Publication() +
        theme(legend.title = element_blank()) +
        xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
        ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
        theme(legend.position = legend) +
        scale_color_manual(values = as.character(P36)) +
        scale_fill_manual(values = as.character(P36)) +
        scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
       
    ggsave(plot_file, width = 8, height = 6, units = "in", dpi = figres)
}

visualize_scree_plot <- function(plot_file, PCA, figres) {
    #Visualize principle component variation results.
    scree.plot <- fviz_eig(PCA, addlabels = TRUE, hjust = -0.3)
    png(plot_file, res = figres)
    print(scree.plot)
    dev.off()
}

save_loading_scores <- function(write_file, df, figres) {
    #Save list of genes that have a positive effect on variation of principle.
    #Component 1 and 2 sorted from most influential.
    write.table(df, file = write_file)
}

theme_Publication <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
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
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face = "italic"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

umap_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns,
                     figres = 100, size = 1, legend = "right") {
    #Run default parameters of umap.
    U <- umap(t(exprs))
    visualize_umap(file.path(results_path, paste0("2.umap_", base_file_name)),
                   U$layout, labels[,target_columns[1]],
                   labels[,target_columns[2]], figres, size, legend)
}

visualize_umap <- function(plot_file, U, class1, class2, figres, size, legend) {
    #Visualize umap reduction.
    minx <- min(U[,1])
    maxx <- max(U[,1])
    miny <- min(U[,2])
    maxy <- max(U[,2])
    P28 <- createPalette(length(levels(factor(class2))), distinctColorPalette(28))
    qplot(U[,1], U[,2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
        theme_Publication() + theme(legend.title = element_blank()) +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        scale_color_manual(values = as.character(P28)) +
        scale_fill_manual(values = as.character(P28)) +
        xlim(minx, maxx) + ylim(miny, maxy) +
        theme(legend.position = legend) +
        scale_shape_manual(values = seq(1, length(levels(factor(class1)))))

    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = figres)
}

tsne_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns, figres = 100) {
    #Run default parameters of tsne.
    T <- Rtsne(t(exprs), perplexity = 1)
    visualize_tSNE(file.path(results_path, paste0("2.tsne_", base_file_name)),
                   T$Y, labels[,target_columns[1]],
                   labels[,target_columns[2]], figres)
}

visualize_tSNE <- function(plot_file, U, class1, class2, figres) {
    #Visualize tsne reduction.
    minx <- min(U[,1])
    maxx <- max(U[,1])
    miny <- min(U[,2])
    maxy <- max(U[,2])
    png(plot_file, res = figres)
    par(mar = c(5, 4, 2, 4), xpd = TRUE)
    plot(U[,1], U[,2], frame = FALSE,
         ylim = c(miny - 1, maxy + 1), xlim = c(minx - 1, maxx + 1),
         pch = as.numeric(as.factor(class1)),
         col = as.numeric(as.factor(class2)),
         xlab = "Dim 1", ylab = "Dim 2")
    legend("topright", inset = c(-0.25, -0.1), bty = "n",
           pch = as.numeric(levels(as.factor(as.numeric(as.factor(class1))))),
           legend = levels(as.factor(class1)))
    legend("bottomright", inset = c(-0.25, 0), bty = "n", pch = "-",
           col = levels(as.factor(as.numeric(as.factor(class2)))),
           legend = c(levels(as.factor(class2))))
    dev.off()
}

#Loadfiles: This loads count and target file info (internal function only).
loadfiles <- function(count_file, target_file) {
    #Load in count data and target information from csv.
    if (typeof(count_file) == "character") {
        counts <- read.table(count_file, header = TRUE,
                             sep = "\t", row.names = 1,
                             as.is = TRUE, check.names = FALSE)
        counts <- counts[,order(names(counts))]
    } else if ((isFALSE(count_file)) || (is.null(count_file))) {
        counts <- FALSE
    }

    if (typeof(target_file) == "character") {
        targets <- read.table(target_file, header = TRUE,
                              sep = ",", row.names = 1,
                              as.is = TRUE, check.names = FALSE)
        targets <- targets[order(rownames(targets)),]
    } else if ((isFALSE(target_file)) || (is.null(count_file))) {
        targets <- FALSE
    }

    results <- list("counts" = counts, "targets" = targets)
    return(results)
}

#Generate folder: This generates new folder (internal function only).
generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}