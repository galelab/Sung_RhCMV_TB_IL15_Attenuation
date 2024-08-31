#s1 normalize raw counts function following Gale lab protocol.
#This function allows the normalization of raw counts following a general protocol developed by previous members of the Gale lab.

library(edgeR)

s1_normalize_raw_counts <- function(countfile, targetfile,
                                    gene_conversion_file = FALSE,
                                    target_column = FALSE, batch_column = FALSE,
                                    blocking_column = FALSE,
                                    visualize_data = TRUE,
                                    filter_genes_below_counts = 0,
                                    filter_method = "sum",
                                    norm_method = "none",
                                    results_folder = "s1_norm_raw_counts_results",
                                    figres = 100) {

    #Read in files.
    print("STATUS: loading files")
    files <- loadfiles(count_file = countfile,
                      target_file = targetfile)
    DE_DF <- DGEList(counts = files$counts)

    #Filter out genes with low counts.
    print("STATUS: filtering out genes with low counts")
    if (filter_method == "sum") {
        A <- rowSums(DE_DF$counts)
        isexpr <- A >= filter_genes_below_counts
        DE_DF <- DE_DF[isexpr,]
    }
    else if (filter_method == "mean") {
        A <- rowMeans(DE_DF$counts)
        isexpr <- A >= filter_genes_below_counts
        DE_DF <- DE_DF[isexpr,]
    }

    #Normalize via TMM.
    print("STATUS: getting normalizing factors (method TMM)")
    DE_DF_fl <- calcNormFactors(DE_DF)

    #Get biological coefficients of variation.
    #print("STATUS: getting biological coefficients of variation (takes time...)")
    #count_matrix_flv <- biological_coefficients_variation(DE_DF_fl)

    #Set up model design.
    print("STATUS: setting up model design")

    if (length(files$targets[,target_column]) != length(colnames(DE_DF))) {
        print("WARNING: different number of treatments and column names in count file")
        print(paste0("Length of treatments: ", length(files$targets[,target_column])))
        print(paste0("Length of column names in count/normalized matrix: ", length(colnames(DE_DF))))
    }
    else if (all.equal(rownames(files$targets), colnames(files$counts)) != TRUE) {
        print("WARNING: order of samples in target file does not match order in count file (needs fixing before we can proceed)")
        #files$counts <- files$counts[,rownames(files$targets)]
        #if (all.equal(rownames(files$targets), colnames(files$counts)) == TRUE) {
            #print("WARNING: order of samples has been corrected!")
        #}
    }
    else {
        results_path <- generate_folder(results_folder)
        unlink(paste0(results_folder, "/*"))
        factors <- list()

        if (typeof(batch_column) == "logical") {
            treatment   <- factor(files$targets[,target_column],
                           levels = unique(files$targets[,target_column]))
            design      <- model.matrix(~0 + treatment)
            rownames(design) <- colnames(DE_DF$counts)
            colnames(design) <- make.names(colnames(design))
            design <- design[,colnames(design)[order(tolower(colnames(design[,])))]]
        } else {
            treatment   <- factor(files$targets[,target_column],
                           levels = unique(files$targets[,target_column]))
            batch       <- factor(files$targets[,batch_column],
                           levels = unique(files$targets[,batch_column]))
            design      <- model.matrix(~0 + treatment + batch)
            rownames(design) <- colnames(DE_DF$counts)
            colnames(design) <- make.names(colnames(design))
            design <- design[,colnames(design)[order(tolower(colnames(design[,])))]]
        }

        if (is.fullrank(design) == TRUE & is.null(nonEstimable(design))) {
            #colnames(design) <- levels(CLASS1)
            if (blocking_column != FALSE) {
                BLOCKID <- factor(files$targets[,blocking_column],
                                  levels = unique(files$targets[,blocking_column]))
                corfit <- duplicateCorrelation(DE_DF_fl$counts, design,
                                               block = BLOCKID)
            }

            #Run voom.
            print("STATUS: running voom")
            png(file.path(results_path, "1.voomplot.png"), res = figres)
            #Transform count data to log2-counts per million.
            V.CPM <- voom(DE_DF_fl, normalize.method = norm_method, design = design, plot = TRUE, span = 0.1)
            dev.off()

            #Save normalized counts and design variable used for linear modeling later.
            orig.cols    <- colnames(files$counts)
            #orig.cols    <- append(orig.cols, "Name", after = 0)
            orig.rows    <- rownames(V.CPM$E)
            write.table(data.frame(V.CPM$E), sep = "\t",
                       row.names = orig.rows, col.names = orig.cols,
                       file = file.path(results_path, "1.norm_matrix.txt"))
            norm_matrix  <- V.CPM$E

            if (typeof(gene_conversion_file) == "character") {
                rhesus2human <- read.csv(file = gene_conversion_file,
                                         header = TRUE,
                                         stringsAsFactors = FALSE)
                nm_hgnc      <- merge(rhesus2human, norm_matrix,
                                      by.x = "Gene.stable.ID", by.y = "row.names")
                nm_hgnc      <- avereps(nm_hgnc, ID = nm_hgnc$Gene.stable.ID)
                write.table(nm_hgnc, sep = "\t",
                            file = file.path(results_path, "1.norm_matrix_HGNC.txt"))
            }

            saveRDS(V.CPM, file.path(results_path, "1.voomobject.rds"))

            if (blocking_column != FALSE) {
                saveRDS(corfit, file.path(results_path, "1.corfit.rds"))
            }

            saveRDS(design, file = file.path(results_path, "1.designobject.rds"))

            if (visualize_data == TRUE) {
                print("STATUS: generating figures")
                visualize_counts(files$counts, V.CPM$E, files$targets,
                                 figres = figres,
                                 results_path = results_path)
            }

            #results_norm <- list("norm_exprs_voom" = V.CPM, "design" = design)

            #return(results_norm)

        } else {
            print("WARNING: error with design matrix... rethink how it is being set up")
            print(paste0("WARNING: full rank check is ", is.fullrank(design)))
            print(paste0("WARNING: nonEstimatable check is ", is.null(nonEstimable(design))))
        }
    }
}

#OTHER FUNCTIONS used by normalize_raw_counts.

#Visualize counts.
visualize_counts <- function(countsmatrix, norm_exprs, labels,
                             figres = 100, results_path) {
    #Generate figures for counts.
    print("STATUS: generating log2 boxplot of counts")
    generate_boxplots(log2(countsmatrix + 1), labels[,1],
                      file.path(results_path, "1.boxplot_raw_count_matrix.png"),
                      figres, maintitle = "Raw Count Matrix",
                      ylabtitle = "Log2 Expression")

    print("STATUS: generating boxplot of normalized voom counts")
    generate_boxplots(norm_exprs, labels[,1],
                      file.path(results_path, "1.boxplot_vnorm_matrix.png"),
                      figres, maintitle = "Normalized Count Matrix",
                      ylabtitle = "Voom Normalized Expression")

    print("STATUS: generating density plot of all log counts")
    generate_density_plot(log2(countsmatrix + 1), labels[,1],
                          file.path(results_path,
                                    "1.densities_raw_log_counts.png"),
                          figres)

    print("STATUS: generating density plot of raw counts")
    generate_density_plot(countsmatrix, labels[,1],
                          file.path(results_path, "1.densities_raw_counts.png"),
                          figres)

    print("STATUS: generating density plot of normalized voom counts")
    generate_density_plot(norm_exprs, labels[,1],
                          file.path(results_path,
                                    "1.densities_vnorm_matrix.png"),
                          figres)

    #print("STATUS: generating biological variation vs abundance")
    #png(file.path(results_path, "1.biologicalcoefficientvariation_raw.png"),
    #    res = figres)
    ##par(mar = c(1,1,1,1))
    #plotBCV(count_matrix_flv, cex = 0.4,
    #        main = "Biological Coefficient of Variation (BCV) vs Abundance")
    #dev.off()
}

#Biological coefficients of variation.
biological_coefficients_variation <- function(count_matrix_fl) {
    count_matrix_flv <- estimateCommonDisp(count_matrix_fl,
                                           verbose = TRUE)  #Print the BCV value.
    count_matrix_flv <- estimateTrendedDisp(count_matrix_flv)
    count_matrix_flv <- estimateTagwiseDisp(count_matrix_flv)
    return(count_matrix_flv)
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

#Generate boxplots.
generate_boxplots <- function(data, labels, filename, figres, maintitle, ylabtitle) {
    png(filename, res = figres)
    #par(mar = c(1,1,1,1))
    minvalue <- min(data)
    maxvalue <- max(data)
    boxplot(data, labels = labels, ylim = c(minvalue - 1, maxvalue + 1),
            ylab = ylabtitle, main = maintitle, cex.axis = .6, las = 2,
            frame = FALSE)
    dev.off()
}

#Generate density plot.
generate_density_plot <- function(data, labels, filename, figres) {
    png(filename, res = figres)
    par(xpd = TRUE)
    if (length(labels) > 10) {
        plotDensities(data, legend = FALSE)
    } else {
        plotDensities(data, legend = "topright",
                      inset = c(-0.2,0), levels(labels))
    }
    dev.off()
}