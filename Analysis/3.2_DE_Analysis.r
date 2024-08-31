#' s3 differential gene analysis
#'
#' This function runs differential gene analysis (linear modeling on normalized counts)
#' @param targetfile target file.
#' @param folder2voomobjects folder where necessary voom R objects are stored (should be results folder from step 1) defualt is s1_norm_raw_counts_results
#' @param results_folder User can specify name of output folder they want results to go in (default is s3_DE_results/)
#' @param gene_conversion_file file with alternative gene names (we usually convert Rhesus Ensembl genes to HGNC) 
#' @param blocking_column column to account sampling from the same animal multiple times (needs to be the same as whatever was specified in step s1_normalize_raw_counts)
#' @param matrixfile text file outlining different contrasts to do in DE analysis 
#' @param pvalue parameter for determining significantly DE genes (default 0.05)
#' @param logfoldchange parameter for determining significantly DE genes (default 1.5)
#' @param reset_design_matrix parameter for user to reset design matrix if they want to do this differently than they did in step 1 (defualt is FALSE)
#' @param target_columns vector of columns from target file to use as treatment group, this needs to be set if reset_design_matrix is TRUE
#' @param batch_column column from target file to be used as batch effect if reseting design matrix
#' @keywords differential expression linear modeling
#' @export
#' @import edgeR
#' @import limma
#' @import gplots
#' @import ggplot2
#' @import data.table
#' @import stringr
#' @examples
#' s3_DE_analysis(countfile="./s1_norm_raw_counts_results/1.norm_matrix.txt", targetfile="./p1_modified_count_matrix_results/target_file.csv", gene_conversion_file="rhesus2human.csv", blocking_column=2, matrixfile="./MATRIX.txt")

library(edgeR)
library(stringr)
library(data.table)
library(ggplot2)

s3_DE_analysis <- function(countfile, targetfile,
                           folder2voomobjects="s1_norm_raw_counts_results",
                           results_folder=FALSE,  gene_conversion_file=FALSE,
                           blocking_column=FALSE, matrixfile=FALSE, pvalue=0.05,
                           logfoldchange=1.5, reset_design_matrix=FALSE,
                           target_columns=FALSE, batch_column=FALSE) {
    if (typeof(results_folder) == "logical") {
        results_path <- generate_folder("s3_DE_results")
        unlink("./s3_DE_results/*")
    } else {
        results_path <- generate_folder(results_folder)
        unlink(paste0(results_folder, "/*"))
    }

    files <- loadfiles(count_file = FALSE, target_file = targetfile)

    if (file.exists(paste0(folder2voomobjects, "/1.voomobject.rds"))) {
        eset_voom <- readRDS(paste0(folder2voomobjects, "/1.voomobject.rds"))
    } else {
        print("WARNING: Could not find voom file...make sure you are in correct working directory, correctly pointing to where voom objects are stored or repeat step (s1_normalize_raw_counts)")
    }

    if (isTRUE(reset_design_matrix)) {
        eset_voom$design <- ""
        if (typeof(batch_column) == "logical") {
            counter <- 0
            treats <- list()
            if (length(target_columns) > 1) {
                for (target in target_columns) {
                    counter <- counter + 1
                    tr <- factor(files$targets[, target],
                        levels = unique(files$targets[, target]))
                    treats[[paste0("tr", counter)]] <- tr
                }
                if (length(target_columns) == 2) {
                    tr1 <- treats$tr1
                    tr2 <- treats$tr2
                    design <- model.matrix(~ 0 + tr1:tr2)
                } else if (length(target_columns) == 3) {
                    tr1 <- treats$tr1
                    tr2 <- treats$tr2
                    tr3 <- treats$tr3
                    design <- model.matrix(~ 0 + tr1:tr2:tr3)
                } else {
                    print("WARNING: not set up to work with more than 3 treatments")
                }
                rownames(design) <- colnames(eset_voom$E)
                colnames(design) <- make.names(colnames(design))
                design <- design[, colnames(design)[order(tolower(colnames(design[, ])))]]
                design <- design[, colSums(design) > 0]
            } else {
                tr <- factor(files$targets[, target_columns],
                    levels = unique(files$targets[, target_columns])
                )
                design <- model.matrix(~ 0 + tr)
                rownames(design) <- colnames(eset_voom$E)
                colnames(design) <- make.names(colnames(design))
            }
        } else {
            counter <- 0
            treats <- list()
            if (length(target_columns) > 1) {
                for (target in target_columns) {
                    counter <- counter + 1
                    tr <- factor(files$targets[, target],
                        levels = unique(files$targets[, target]))
                    treats[[paste0("tr", counter)]] <- tr
                }
                Xid <- factor(files$targets[, batch_column],
                    levels = unique(files$targets[, batch_column])
                )

                if (length(target_columns) == 2) {
                    tr1 <- treats$tr1
                    tr2 <- treats$tr2
                    design <- model.matrix(~ 0 + tr1:tr2 + batch)
                } else if (length(target_columns) == 3) {
                    tr1 <- treats$tr1
                    tr2 <- treats$tr2
                    tr3 <- treats$tr3
                    design <- model.matrix(~ 0 + tr1:tr2:tr3 + Xid)
                } else {
                    print ("WARNING: not set up to work with more than 3 treatments")
                }

                rownames(design) <- colnames(eset_voom$E)
                colnames(design) <- make.names(colnames(design))
                design <- design[, colnames(design)[order(tolower(colnames(design[, ])))]]
                design <- design[, colSums(design) > 0]
                excludeAll <- nonEstimable(design)
                if ("tr1" %in% excludeAll) {
                    return("interactions term non estimable")
                }
                design <- design[, !colnames(design) %in% excludeAll]
                if (!is.fullrank(design)) {
                    return("not full rank")
                }

            } else {

                tr <- factor(files$targets[, target_columns],
                    levels = unique(files$targets[, target_columns])
                )
                Xid <- factor(files$targets[, batch_column],
                    levels = unique(files$targets[, batch_column])
                )
                design <- model.matrix(~ 0 + tr + Xid)

                rownames(design) <- colnames(eset_voom$E)
                colnames(design) <- make.names(colnames(design))
                design <- design[, colnames(design)[order(tolower(colnames(design[, ])))]]
                design <- design[, colSums(design) > 0]
                excludeAll <- nonEstimable(design)
                if ("tr" %in% excludeAll) {
                    return("interactions term non estimable")
                }
                design <- design[, !colnames(design) %in% excludeAll]
                if (!is.fullrank(design)) {
                    return("not full rank")
                }

            }
        }
    } else {
        if (file.exists(paste0(folder2voomobjects, "/1.designobject.rds"))) {
            design <- readRDS(paste0(folder2voomobjects, "/1.designobject.rds"))
        } else {
            print ("WARNING: Could not find design file...make sure you are in correct working directory, correctly pointing to where voom objects are stored or repeat step (s1_normalize_raw_counts)")
        }
    }

    if (file.exists(paste0(folder2voomobjects, "/1.corfit.rds"))) {
        corfit <- readRDS(paste0(folder2voomobjects, "/1.corfit.rds"))
    } else {
        print ("WARNING: Could not find corfit file...not necessary for analysis but just warning for the user")
    }

    if (exists("design") == FALSE | exists("eset_voom") == FALSE) {
        print ("WARNING: could not find design and voom model... make sure user is in the correct working directory")
    } else {
        if (typeof(blocking_column) != "logical") {
            BLOCKID  <- factor(files$targets[, blocking_column],
                               levels = unique(files$targets[, blocking_column]))

            if (exists("corfit") == TRUE) {
                V.lmfit <- lmFit(eset_voom, design = design,
                             block = BLOCKID,
                             correlation = corfit$consensus)
            } else {
                V.lmfit <- lmFit(eset_voom, design = design, block = BLOCKID)
            }
        } else {
            if (exists("corfit") == TRUE) {
                V.lmfit <- lmFit(eset_voom, design = design, correlation = corfit$consensus)
            } else {
                V.lmfit <- lmFit(eset_voom, design = design)
            }
        }
        if (typeof(matrixfile) == "character") {
            print("STATUS: getting DE genes...")

            matrix_contrast <- scan(matrixfile,  character(), quote = "")
            contmatrix <- makeContrasts(contrasts = matrix_contrast,
                                         levels = design)

            fit <- contrasts.fit(V.lmfit, contrasts = as.matrix(contmatrix))
            ##############EBAYES PARAMETERS###########
            # robust: should the estimation of df prior and var prior be robustified for 
            ## outlier variance (defualt is FALSE, we set to TRUE)
            # trend: should an intensity be allowed for prior variance (defualt is that prior variance is consistent FALSE)
            # portion: portion of genes which are differentially expressed (defualt is 0.01)
            # stdev.coef.lim: lower and upper limits for standard deviation of LFC for differentially expressed genes 
            ## defualt is c(0.1, 4)
            # winsor.tail.p: gives the left and right tail proportions of x to winsor (only used when robust=TRUE)
            ## (defualt is c(0.5, 0.1))
            ##Winsorizing - the transformation of statistic by limiting extreme values in the statistical data to reduce 
            ## the effect of outliers. Set outliers to a specifc percentile - a 90% winsoriation would see all the data below the 5tj
            ## percentile set to the 5th percentile and data above the 95th percentile set to the 95th percentile
            fit <- eBayes(fit, robust = TRUE, trend = TRUE)
            results <- decideTests(fit, lfc = round(log2(logfoldchange), digits = 2),
                                   method = "separate",
                                   adjust.method = "BH", p.value = pvalue)
            ###############################################################

            #########################OUTPUT FILES##########################

            # OUTPUT ALL fit DATA
            write.fit(fit, file = file.path(results_path, "3.All_data.csv"),
                      digits = 3, method = "separate", adjust = "BH", sep = " ,")

            #OUTPUT ALL LFCS
            write.csv(fit$coefficients,
                        file = file.path(results_path, "3.All_LFC.csv"))

            #OUTPUT ALL LFCS AFTER TRANSLATING GENE NAMES
            convert2HGNC(
                gene_conversion_file,
                "3.All_LFC.csv",
                "3.All_LFC_HGNC.csv",
                results_path
            )

            # OUTPUT ALL PVALS
            write.csv(fit$p.value, file = file.path(results_path,
                                                      "3.All_Pvalues.csv"))

            # OUTPUT ALL PVALS AFTER TRANSLATING GENE NAMES
            convert2HGNC(
                gene_conversion_file,
                "3.All_Pvalues.csv",
                "3.All_Pvalues_HGNC.csv",
                results_path
            )
            # OUTPUT ALL ADJUSTED PVALS
            for (i in 1:ncol(fit$p.value)) {
                 fit$p.value[, i] <- p.adjust(fit$p.value[, i], method = "BH")
            }
            write.csv(fit$p.value, file = file.path(results_path, "3.All_Pvalues_adj.csv"), quote = F)


            # OUTPUT ALL PVALS AFTER TRANSLATING GENE NAMES
            convert2HGNC(
                gene_conversion_file,
                "3.All_Pvalues_adj.csv",
                "3.All_Pvalues_adj_HGNC.csv",
                results_path
            )
 
            # OUTPUT ALL t VALUES
            write.csv(fit$t, file = file.path(results_path, "3.All_tvalues.csv"))

            # OUTPUT ALL t VALUES AFTER GENE TRANSLATION
            convert2HGNC(
                gene_conversion_file,
                "3.All_tvalues.csv",
                "3.All_tvalues_HGNC.csv",
                results_path
            )
            ####################################################################
            if (typeof(gene_conversion_file) == "character") {

                rhesus2human <- read.csv(file = gene_conversion_file, header = TRUE,
                                        stringsAsFactors = FALSE)
            }
            ##################Pull out signifcantly expressed genes#############
            ##SIGNIFICANT LOGVALUES
            #Extract results of differential expression
            #LogFold change is the coefficients
            dataMatrix <- fit$coefficients
            sigMask <- dataMatrix * (results**2) # 1 if significant, 0 otherwise

            # filter for significant genes
            ExpressMatrixLFC <- subset(dataMatrix,
                                       rowSums(sigMask) != 0)
            sigMask <- subset(sigMask, rowSums(sigMask) != 0)

            write.csv(ExpressMatrixLFC, file = file.path(results_path,
                                        "3.Significant_separate_LFC.csv"))
            #temp = ExpressMatrixLFC * sigMask
            convert2HGNC(gene_conversion_file,
                        "3.Significant_separate_LFC.csv",
                        "3.Significant_separate_LFC_HGNC_AV.csv",
                         results_path)
            ##SIGNIFICANT T values
            #Extract results of differential expression
            #LogFold change is the coefficients
            dataMatrix <- fit$t
            sigMask <- dataMatrix * (results**2) # 1 if significant, 0 otherwise

            #filter for significant genes
            ExpressMatrixtvalue <- subset(dataMatrix,
                                          rowSums(sigMask) != 0)

            sigMask <- subset(sigMask, rowSums(sigMask) != 0)
            write.csv(ExpressMatrixtvalue,
                      file = file.path(results_path,
                             "3.Significant_separate_tvalues.csv"))

            convert2HGNC(gene_conversion_file,
                         "3.Significant_separate_tvalues.csv",
                         "3.Significant_separate_tvalues_HGNC_AV.csv",
                         results_path)
            ##SIGNIFICANT P values
            #Extract results of differential expression
            #LogFold change is the coefficients
            dataMatrix <- fit$p.value
            sigMask <- dataMatrix * (results**2) # 1 if significant, 0 otherwise
            #filter for significant genes
            ExpressMatrixPvalue <- subset(dataMatrix, rowSums(sigMask) != 0)
            sigMask <- subset(sigMask, rowSums(sigMask) != 0)
            write.csv(ExpressMatrixPvalue,
                      file = file.path(results_path,
                             "3.Significant_separate_Pvaluesadj.csv"))
            convert2HGNC(gene_conversion_file,
                         "3.Significant_separate_Pvaluesadj.csv",
                         "3.Significant_separate_Pvaluesadj_HGNC_AV.csv",
                         results_path)

            results_path2 <- generate_folder("s3_DE_results/enrichfiles")
            unlink("./s3_DE_results/enrichfiles/*")

            for (i in colnames(fit$coefficients)) {
                #filter for significant genes
                ExpressMatrixLFC1 <- subset(fit$coefficients[, i],
                                            results[, i] != 0)
                allgenes <- fit$coefficients[, i]
                if (typeof(gene_conversion_file) == "character") {
                    if (length(ExpressMatrixLFC1) > 0) {
                        sig_HGNC <- merge(rhesus2human, ExpressMatrixLFC1,
                                          by.x = "Gene.stable.ID",
                                          by.y = "row.names",
                                          all.X = T, all.Y = T)
                        sig_HGNC <- sig_HGNC[ , !(names(sig_HGNC) %in% c("Gene.stable.ID"))]
                        dimensions <- dim(sig_HGNC)
                        if  (dimensions[1] == 0) {
                            print (paste0("WARNING: no sucessful translations from Ensembl to HGNCs so using ensembl IDs for comparison ", i))
                            write.table(ExpressMatrixLFC1,
                                        file = file.path(results_path2,
                                                         paste0(i, "_sig.rnk")),
                                        row.names = FALSE, col.names = FALSE,
                                        sep = "\t", quote = FALSE)
                            write.table(ExpressMatrixLFC1,
                                        file = file.path(results_path2,
                                                paste0(i, "_sig4GSEA.rnk")),
                                        row.names = FALSE, col.names = FALSE,
                                        sep = "\t", quote = FALSE)
                        }
                        else {
                            sig_HGNC <- avereps(sig_HGNC,
                                                ID = sig_HGNC$HGNC.symbol)
                            write.table(sig_HGNC,
                                        file = file.path(results_path2,
                                                         paste0(i, "_sig.rnk")),
                                         row.names = FALSE, col.names = FALSE,
                                         sep = "\t", quote = FALSE)
                            sig_HGNCnodup <- sig_HGNC[!duplicated(sig_HGNC[,2]), ]
                            write.table(sig_HGNCnodup,
                                        file = file.path(results_path2,
                                                paste0(i, "_sig4GSEA.rnk")),
                                        row.names = FALSE, col.names = FALSE,
                                        sep = "\t", quote = FALSE)
                        }
                    } else {
                        print (paste0("WARNING: comparison ", i, " only has 0 significantly different gene therefore not generating enrichment file"))
                    }
                    all_HGNC <- merge(rhesus2human, allgenes,
                                      by.x = "Gene.stable.ID",
                                      by.y = "row.names",
                                      all.X = T, all.Y = T)
                    all_HGNC <- all_HGNC[ , !(names(all_HGNC) %in% c("Gene.stable.ID"))]
                    all_HGNC <- avereps(all_HGNC, ID = all_HGNC$HGNC.symbol)
                    ###This file may have duplicates
                    write.table(all_HGNC, file = file.path(results_path2,
                                                 paste0(i, "_all.rnk")),
                                row.names = FALSE, col.names = FALSE,
                                sep = "\t", quote = FALSE)
                    all_HGNCnodup <- all_HGNC[!duplicated(all_HGNC[,2]), ]
                    ###No Duplicates in this file
                    write.table(all_HGNCnodup, file = file.path(results_path2,
                                               paste0(i,"_all4GSEA.rnk")),
                               row.names = FALSE, col.names = FALSE,
                               sep = "\t", quote = FALSE)

                } else {
                    if (length(significantgenes) > 0) {
                        write.table(significantgenes,
                                    file = file.path(results_path2,
                                            paste0(i, "_sig.rnk")),
                                    col.names = FALSE, sep = "\t",
                                    quote = FALSE)
                    } else {
                        print(paste0("WARNING: comparison ", i, " only has 0 significantly different gene therefore not generating enrichment file"))                        
                    }
                    write.table(allgenes, file = file.path(results_path2,
                                                 paste0(i, "_all.rnk")),
                                col.names = FALSE, sep = "\t", quote = FALSE)
                }
            }
            new_colnames <- c()
            for (i in colnames(ExpressMatrixLFC)) {
                i <- str_remove_all(i, "tr")
                new_colnames <- c(new_colnames, i)
            }
            colnames(ExpressMatrixLFC) <- new_colnames

            if (dim(ExpressMatrixLFC)[1] > 0) {
                hm_results      <- vizualize_DE_genes_HM(ExpressMatrixLFC,
                                                        file.path(results_path,
                                                        "3.heatmap_djn.png"))
                global_modules  <- hm_results$modules
                write.csv(global_modules,
                        file = file.path(results_path, "3.modules.csv"))
                global_modulesM <- as.matrix(global_modules)
                GM_HGNC         <- merge(rhesus2human, global_modulesM,
                                        by.x = "Gene.stable.ID",
                                        by.y = "row.names",
                                        all.X = T, all.Y = T)
                write.csv(GM_HGNC, file = file.path(results_path,
                                                "3.modules_HGNC.csv"))

                clustermatrix  <- hm_results$clustermatrix
                #invert row order
                clustermatrix   <- clustermatrix[order(nrow(clustermatrix):1), ]
                write.csv(clustermatrix, file = file.path(results_path,
                                                        "3.Clustered_LFC.csv"),
                        quote = FALSE)

                colnames(results) <- new_colnames
                vizualize_DE_genes_bp(results, file.path(results_path,
                                            "3.barplot_NumDEgenes.png"))
                dotplot4bulkheatmap(hm_results$modules, clustermatrix, ExpressMatrixPvalue, breaks = c(15, 30), file.path(results_path, "3.dotplot.png"))

            } else {
                print ("WARNING: NO SIGNIFICANT GENES SO NOT GENERATE FIGURES OR EXTRA FILES")
            }

            SigModules <- read.csv("./s3_DE_results/3.modules.csv")
            Results <- results[SigModules$X,]
            SigHM(Results, file.path(results_path, "3.SigHM.png"))
            write.csv(Results, file = file.path(results_path,
                                                "3.Results.csv"))

        } else {
            print("WARNING: need to specify matrix file")
        }
    }
    return(hm_results)
}

SigHM <- function(Results, plot_file){
    png(plot_file, width = 9, height = 6, units = "in", res = 300)
    Rdf <- reshape2::melt(Results, c("x", "y"), value.name = "z")
    breaks = as.numeric(Rdf$x[c(15.5, 30.5)])
    ggplot(Rdf, aes(x=y,y=x,fill=z))+
        geom_tile()+
        scale_fill_gradient2(low = "mediumblue",
        mid = "white",
        high = "red")+
        geom_vline(xintercept = breaks+.5,
        linetype=1, colour="black")
    ggsave(plot_file, dpi = 300)
    dev.off()
}

convert2HGNC <- function(gene_conversion_file, input_file,
                        output_file, results_path) {

    if (typeof(gene_conversion_file) == "character") {

        rhesus2human <- read.csv(file = gene_conversion_file, header = TRUE,
                                 stringsAsFactors = FALSE)
        DE_HGNC_LFC <- read.csv(file.path(results_path, input_file), header = T,
                                row.names = 1, check.names = FALSE, sep = ",")
        if (dim(DE_HGNC_LFC)[1] > 0) {
            DE_HGNC_LFC <- merge(rhesus2human, DE_HGNC_LFC,
                                by.x = "Gene.stable.ID", by.y = "row.names")
            DE_HGNC_LFC <- avereps(DE_HGNC_LFC, ID = DE_HGNC_LFC$HGNC.symbol)
            write.csv(DE_HGNC_LFC, file = file.path(results_path, output_file))
        } else {
            print("WARNING NO SIGNIFICANT GENES SO CAN'T GENERATE HGNC FILE")
        }
    } else {
        print("WARNING: need to specify conversion file to convert Ensembls to HGNCs")
    }
}

vizualize_DE_genes_HM <- function(data, plot_file) {
    print("STATUS: Generating heatmap of DE genes...")
    png(plot_file, width = 8, height = 10, units = "in", res = 300)
    global_modules <- heatmap.F.4(data, cutoff = 1, distmethod = "pearson",
                                  clustermethod = "ward.D2", clusterdim = "row")
    dev.off()
    return(global_modules)
}

vizualize_DE_genes_bp <- function(results, plot_file) {
    print("STATUS: Generating bar plot of number of DE genes...")
    png(plot_file, width = 8, height = 10, units = "in", res = 300)
    results_t <- t(summary(results))
    results_t <- results_t[, -2]

    for (i in 1:(length(row.names(results_t)))) {
        results_t[i, 1] <- results_t[i, 1] * -1
    }

    DE <- as.data.frame(results_t)
    DE <- setnames(DE, old = c("Var1", "Var2", "Freq"),
                   new = c("Time_Point", "group", "DE_genes"))

    #Create plot
    ggplot(DE, aes(x = Time_Point, y = DE_genes, fill = group,
           label = DE$DE_genes)) +
           geom_bar(stat = "identity", position = "identity") +
    # geom_text(size = 5, position = position_stack(vjust = 0) )+
    scale_fill_manual(values = c("#ff4d4d", "#9d9dff"), limits = c("Up", "Down")) +
    ylab("Number of Differentially Expressed Genes") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggsave(plot_file, dpi = 300)
}

theme_Publicationdot <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size)
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
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
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

dotplot4bulkheatmap <- function(clusters, clustermatrix, pvals, breaks=NULL, figurename) {
    df <- data.frame(matrix(ncol = length(colnames(clustermatrix)), 
        nrow = length(unique(clusters))))
    colnames(df) <- colnames(clustermatrix)
    rownames(df) <-unique(clusters)
    dfpval <- data.frame(matrix(ncol = length(colnames(clustermatrix)), 
        nrow = length(unique(clusters))))
    colnames(dfpval) <- colnames(clustermatrix)
    rownames(dfpval) <- unique(clusters)
    colnames(pvals) <- colnames(clustermatrix)
    for (cluster in unique(clusters)) {
        genes = which(clusters==cluster)
        for (col in colnames(clustermatrix)) {
            tmp = median(clustermatrix[names(genes), col])
            df[cluster, col] <- tmp
            countpvals =sum(pvals[names(genes), col] < 0.05)
            dfpval[cluster, col] <- countpvals
        }
    }
    dfmelt <- reshape2::melt(as.matrix(df))
    colnames(dfmelt)[3] <- "LFC"
    dfpvalsmelt <- reshape2::melt(as.matrix(dfpval))
    colnames(dfpvalsmelt)[3] <- "pval.count"
    if (all.equal(dfmelt$Var2, dfpvalsmelt$Var2)==TRUE) {
        dftotal <- cbind(dfmelt, dfpvalsmelt[,3])
        colnames(dftotal)[4] <- "pval.count"
    } else {
        stop("WARNING: DATA NOT IN THE SAME ORDER")
    }
    myPalette <- colorRampPalette(c("blue", "skyblue", "white", "orange", "red"))(100)
    sc <- scale_fill_gradientn(
        colours = myPalette,
        values = scales::rescale(c(
        min(dftotal$LFC), 0,
        0, max(dftotal$LFC)
        )), 
    )
    dftotal$LFC <- as.numeric(dftotal$LFC)
    dftotal$pval.count <- as.numeric(dftotal$pval.count)

    pl <- ggplot(data = dftotal, aes(
        x = Var2,
        y = Var1,
        size = pval.count
    )) + geom_point(aes(fill=LFC),
       colour="black",pch=21) + theme_Publicationdot() + labs(fill="Median\nLFC") +
       scale_size_continuous(range = c(1, 8)) + geom_vline(xintercept = breaks+.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
    sc
    ggsave(figurename, height=8, width=15, dpi=300)
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