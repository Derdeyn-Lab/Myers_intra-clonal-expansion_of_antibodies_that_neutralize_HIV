############## --Libraries--##############

# setwd("/mnt/home/lwhitmo/projects/D001_clonotypes_RUp16/WhitmoreAnalysis/2.AnalysisWhitmoreProcessedData")
options(future.globals.maxSize= 11377049600)
library(Seurat)
library(cowplot)
library(ggplot2)
library(stringr)
library(dplyr)
library(scRepertoire)
library(Polychrome)
library(xlsx)
library(SingleCellExperiment)
library(Ibex)
library(viridis)
library(biomaRt)
library(data.table)
library(ggVennDiagram)
library(glmnet)

############## --Folders Path--##############

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

############## --Functions--################

FeaturePlotSingleUcellscore <- function(obj, feature, metadata_column=NULL, ...) {
    all_cells <- colnames(obj)
    if (!is.null(metadata_column)) {
        print("getting groups")
        groups <- levels(obj@meta.data[, metadata_column])
        print(groups)
    }
    # the minimal and maximal of the value to make the legend scale the same.
    minimal <- min(obj@meta.data[, feature])
    maximal <- max(obj@meta.data[, feature])
    print(minimal)
    print(maximal)
    ps <- list()
    if (!is.null(metadata_column)) {
        print("STATUS: Multiple groups")
        for (group in groups) {
            print(group)
            subset_indx <- obj@meta.data[, metadata_column] == group
            subset_cells <- all_cells[subset_indx]
            p <- FeaturePlot(obj, features = feature, raster = FALSE, cells = subset_cells, ...) +
                # scale_color_viridis_c(limits = c(minimal, maximal), direction = 1) +
            scale_color_viridis_c(values=scales::rescale(x = c(minimal, 0.1, maximal)), to = c(0, 1), 
                    from = c(minimal, maximal), direction = 1) +
                ggtitle(group) +
                theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 12)) +
                labs(x = "UMAP 1", y = "UMAP 2") +
                theme(plot.title = element_text(size = 11), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8)) +
                theme(axis.ticks = element_blank(), axis.text = element_blank())

            ps[[group]] <- p
        }
    } else {
        p <- FeaturePlot(obj, features = feature, cells = all_cells, ...) +
            # scale_color_viridis_c(limits = c(minimal, maximal), direction = 1) +
            scale_color_viridis_c(values = scales::rescale(
                x = c(0, 0.1, maximal), to = c(0, 1),
                from = c(0, maximal)
            ), direction = 1) +
            theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 12)) +
            labs(x = "UMAP 1", y = "UMAP 2") +
            theme(plot.title = element_text(size = 12), axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8)) +
            theme(axis.ticks = element_blank(), axis.text = element_blank())
        ps[["all"]] <- p
    }
    return(ps)
}

FeaturePlotSingle <- function(obj, features, metadata_column, ...) {
    # the minimal and maximal of the value to make the legend scale the same.
    library(patchwork)
    ps <- list()
    minimal <- min(obj[["RNA"]]@data[features, ])
    maximal <- max(obj[["RNA"]]@data[features, ])
    for (feature in features) {
        # minimal <- min(obj[["RNA"]]@data[feature, ])
        # maximal <- max(obj[["RNA"]]@data[feature, ])
        p <- FeaturePlot(obj, features = feature,  ...) +
            scale_color_gradient(limits = c(minimal, maximal)) +
            theme(plot.title = element_text(size = 8),
                  axis.title.x = element_text(size = 6), axis.title.y = element_text(size = 6),
                  axis.ticks = element_blank(), axis.text = element_blank()) +
                  ggtitle(feature) + 
                 labs(x = "UMAP 1", y = "UMAP 2")
        ps[[feature]] <- p
    }
    return(ps)
}

theme_Publication <- function(base_size = 14, base_family = "arial") {
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
            legend.key.size = unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

initial_SC_Seurat <- function(data_dir, sampleID_1, sampleID_2,tissue="pbmc",ndims=10) {
    results_path = file.path("RUp16_results", "initial_seurat_GEX") 
    generate_folder(results_path)
    results_path = file.path("RUp16_results", "initial_seurat_GEX", sampleID_1) 
    generate_folder(results_path)
    message("STATUS: processing sample ", sampleID_1)
    ## -- Load data -- ##
    soh5 <- Read10X_h5(file.path(data_dir, "count","sample_filtered_feature_bc_matrix.h5"))
    soh5GEX <- soh5$'Gene Expression'
    so <- CreateSeuratObject(counts = soh5GEX, min.cells = 3, min.features = 200, project = sampleID_2)
    v = intersect(rownames(so), mtgenes)
    so[["percent.mt"]] <- PercentageFeatureSet(so, features = v)
    so@meta.data$tissue <- rep(tissue, ncol(so))
    splitdata = str_split(so@meta.data$orig.ident, "_")
    animals <- c()
    time <- c()
    for (i in splitdata) {
        animals <- c(animals, i[1])
        time <- c(time, i[2])
    }
    so@meta.data$animal <- animals
    so@meta.data$time <- time
    ## -- Generate filter figures --##
    VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(file.path(results_path, paste0(sampleID_2, "_VnPlotMt.png")), dpi = 300)
    so <- subset(so, subset = nFeature_RNA > 200 & percent.mt < 20)
    VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(file.path(results_path, paste0(sampleID_2, "_VnPlot_Filtered.png")), dpi = 300)
    ## -- set meta data  -- ##
    message("STATUS: Number of cells ", length(colnames(so)))
    if (length(colnames(so)) > 100) {
        print("Normalize")
        so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
        so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
        ##-- Identify the 10 most highly variable genes --##
        top10 <- head(VariableFeatures(so), 10)
        ## -- plot variable features with and without labels -- ##
        plot1 <- VariableFeaturePlot(so)
        LabelPoints(plot = plot1, points = top10, repel = TRUE)
        ggsave(file.path(results_path, paste0(sampleID_2, "_HighlyVariableGenes.png")), bg="white", dpi = 300)
        ## --scale data --##
        all.genes <- rownames(so)
        so <- ScaleData(so, features = all.genes)
        ## -- Feature reduction --##
        so <- RunPCA(so, features = VariableFeatures(object = so))
        ElbowPlot(so, ndims = 50)
        ggsave(file.path(results_path, paste0(sampleID_2, "_elbowplot.png")), bg = "white", dpi = 300)
        VizDimLoadings(so, dims = 1:2, reduction = "pca")
        ggsave(file.path(results_path, paste0(sampleID_2, "_PCA_Loadings.png")), bg = "white", dpi = 300)
        DimPlot(so, reduction = "pca")
        ggsave(file.path(results_path, paste0(sampleID_2, "_PCA.png")), bg = "white", dpi = 300)
        so <- RunUMAP(so, dims = 1:ndims)
        DimPlot(so, reduction = "umap")
        ggsave(file.path(results_path, paste0(sampleID_2, "_UMAP.png")), bg = "white", dpi = 300)
        so <- FindNeighbors(so, reduction="pca", dims = 1:ndims)
        so <- FindClusters(so, resolution = 0.5) #Default values ued 
        DimPlot(so, reduction = "umap")
        ggsave(file.path(results_path, paste0(sampleID_2, "_clusters.png")), bg = "white", dpi = 300)
        for (cluster in unique(Idents(so))) {
            cluster.markers <- FindMarkers(so, ident.1 = cluster, min.pct = 0.25)
            cluster.markers[["gene_name"]] <- rownames(cluster.markers)
            cluster.markers <- cluster.markers %>%
                filter(p_val_adj <= 0.05) %>%
                dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
            cluster.markers <- cluster.markers[order(-cluster.markers$avg_log2FC), ]
            rownames(cluster.markers) <- cluster.markers$gene_name
            cluster.markers$gene_name <- NULL
            x <- head(cluster.markers, n = 6)
            VlnPlot(so, pt.size = 0.1, features = rownames(x))
            ggsave(file.path(results_path, paste0(cluster, "_cluster_markers.png")), height = 6, width = 6, dpi = 300)
            write.csv(cluster.markers, file.path(results_path, paste0(cluster, "_cluster_markers.csv")))
        }
        return(list("seurat"=so))
    } else {
        message("STATUS: Not enough cells for sample ", sampleID_2)
        return(NULL)
    }
}

get_DE_between_conditions <- function(ident_1, ident_2, compare,
                                      combined,
                                      result_folder, fontsize = 9, height = 8,
                                      foldchange = 0.26,
                                      pval = 0.05, percent_cells = NULL, pctthreshdefault=0.01) {
    message("STATUS: getting DEs...")
    ### Parameters for FindMarkers
    #### test.use: mast (default is wilcox)
    #### min.pct: 0.1(default: 0.1) filter out genes (features) that are detected at less than 10 percent frequency in cells in ident_1 or ident_2
    #### logfc: 0 (default is .25) logfc must be higher than 0 (I set this at 0 because I filter this out at a later stage - this was primarily done to
    ####                          understand how the tool (Seurat) works
    #### min.cells.feature: 3 (default is 3) minimum number of cells expressing the feature in at least one of the two groups (similar to min.pct)
    #### min.cells.group: 3 (defualt is 3) minimum number of cells in the group
    #### max.cells.per.ident: Inf (default Inf-means no down sampling) Down sample each identity class (cluster of cells) to a max number
    #### min.dif.pct: Inf (default Inf) Only test genes that show minimum difference in the fraction of detection between the two identities (cluster)
    #### recorrect_umi = FALSE to use the existing corrected counts
    if (file.exists(file.path(result_folder, paste0("DEgenes_full_", ident_1, "_", ident_2, "_", compare, ".txt")))) {
        message("STATUS: DE has already been done filtering from full file ")
        DEgenes <- read.table(file.path(result_folder, paste0("DEgenes_full_", ident_1, "_", ident_2, "_", compare, ".txt")), header = TRUE, row.names = 1)
    } else {
        DEgenes <- FindMarkers(combined,
            assay = "RNA",
            slot="data", 
            ident.1 = ident_1,
            ident.2 = ident_2,
            test.use = "wilcox",
            logfc.threshold = 0, 
            min.pct = pctthreshdefault
            # recorrect_umi = FALSE
        )
        write.table(data.frame(DEgenes),
            file.path(
                result_folder,
                paste0(
                    "DEgenes_full_",
                    ident_1,
                    "_", ident_2, "_", compare,
                    ".txt"
                )
            ),
            sep = "\t", quote = FALSE
        )
    }

    DEgenes[["gene_name"]] <- rownames(DEgenes)

    DE_sig_final <- DEgenes %>%
        filter(avg_log2FC >= foldchange | avg_log2FC <= -foldchange) %>%
        dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj) ## Note log2 of 1.2 = 0.26
    DE_sig_final <- DE_sig_final %>%
        filter(p_val_adj <= pval) %>%
        dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
    DE_sig_final <- DE_sig_final %>%
        filter(pct.1 >= pctthreshdefault |pct.2 >= pctthreshdefault) %>%
        dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)

    if (!is.null(percent_cells)) {
        DE_sig_final <- DE_sig_final %>%
            filter(pct.1 >= percent_cells | pct.2 >= percent_cells) %>%
            dplyr::select(gene_name, p_val, avg_log2FC, pct.1, pct.2, p_val_adj)
    }
    rownames(DE_sig_final) <- DE_sig_final$gene_name

    DE_sig_final$gene_name <- NULL
    DE_sig_final <- DE_sig_final[order(-DE_sig_final$avg_log2FC), ]
    if (is.null(percent_cells)) {
        write.table(data.frame(DE_sig_final),
            file.path(
                result_folder,
                paste0(
                    "DEgenes_sig_",
                    ident_1, "_", ident_2, "_", compare,
                    ".txt"
                )
            ),
            sep = "\t", quote = FALSE
        )
    } else {
        write.table(data.frame(DE_sig_final),
            file.path(
                result_folder,
                paste0(
                    "DEgenes_sig_",
                    ident_1, "_", ident_2, "_", compare,
                    "_", percent_cells, ".txt"
                )
            ),
            sep = "\t", quote = FALSE
        )
    }

    if (length(rownames(DE_sig_final)) > 1) {
        message("STATUS: Making BarPlot")
        pl <- SC_DE_barplot(DE_sig_final, ident_1, ident_2, fontsize = fontsize)
        if (is.null(percent_cells)) {
            ggsave(file.path(
                result_folder, paste0(
                    "DEgenes_sig_", ident_1, "_",
                    ident_2, "_", compare, "barPlot.png"
                )
            ),
            width = 5, height = height, units = "in", dpi = 300
            )
        } else {
            ggsave(file.path(
                result_folder, paste0(
                    "DEgenes_sig_", ident_1, "_",
                    ident_2, "_", compare, "_", percent_cells, "barPlot.png"
                )
            ),
            width = 5, height = height, units = "in", dpi = 300
            )
        }
    }
    return(DE_sig_final)
}

run_gsea_and_ora <- function(finalrankedgenes, gmt.file, universe, region, results_folder, GSEA=TRUE) {
    library(fgsea)
    pathways <- gmtPathways(gmt.file)
    if (isTRUE(GSEA)) {
        fgseaRes <- fgsea(
            pathways = pathways,
            stats = finalrankedgenes,
            minSize = 15,
            maxSize = 1000
        )
        message("done gsea")
        fgseaRes <- fgseaRes[order(pval), ]
        fgseaResSig <- fgseaRes[fgseaRes$padj <= 0.05, ]
        fgseaResdf <- as.data.frame(fgseaRes)
        fwrite(as.data.frame(fgseaResdf), file.path(results_folder, paste0(region, "_GSEA_allResults.csv")))
        fwrite(as.data.frame(fgseaResSig), file.path(results_folder, paste0(region, "_GSEA_allResultsSig.csv")))
        message("done writing gsea")
    }
    foraRes <- fora(pathways, finalrankedgenes, universe, minSize = 5, maxSize = Inf)
    foraRes <- foraRes[order(pval), ]
    foraResSig <- foraRes[foraRes$padj <= 0.05, ]
    fwrite(as.data.frame(foraRes), file.path(results_folder, paste0(region, "_ORA_allResults.csv")))
    fwrite(as.data.frame(foraResSig), file.path(results_folder,paste0(region, "_ORA_allResultsSig.csv")))
    message("done ora")

    # return(list("ORA" = foraRes, "sigORA" = foraResSig, "sigGSEA" =fgseaResSig, "GSEA" = fgseaRes))
}
 stop()

############## --Load Data for VDJB RUp16 without wk36d1 remove group RUp16_PBMC_wk42_96--##############
results_folder <- "RUp16_results"
generate_folder(results_folder)
RUPsamples = c("010_RUp16_weeks_18-20-48", "010_RUp16_weeks_26-36-48", "O10_RUp16_Wk80-81","010_RUp16_Wk87", "010-RUp16-Efficacy_Wk83LN_Rsf16_Wk50",
    "010-RUp16-Efficacy_Wk93","010-RUp16-Efficacy_Wk83LN_Rsf16_Wk50", "RUp16_PBMC_wk26_36_RSf16_PBMC_wk48")
#RUp16_PBMC_wk42_96
contig.list <- list()
for (f in RUPsamples) {
    dirs <- list.dirs(file.path("/mnt/home/lwhitmo/projects/D001_clonotypes_RUp16/preprocessing/NCBIAnno_v103/", f, "outs/per_sample_outs/"), full.names = FALSE, recursive = FALSE)
    for (d in dirs) { 
        contigs <-  read.csv(file.path("/mnt/home/lwhitmo/projects/D001_clonotypes_RUp16/preprocessing/NCBIAnno_v103/", f, "outs/per_sample_outs/", d, "vdj_b", "filtered_contig_annotations.csv"))
        ID = str_remove_all(d, "-Efficacy")
        ID = str_remove_all(ID, "010-")
        ID = str_remove_all(ID, "010_")
        ID = str_remove_all(ID, "_LN")
        ID = str_replace_all(ID, "Rup", "RUp")
        if ((isTRUE(str_detect(ID, "Rup"))) || (isTRUE(str_detect(ID, "RUp")))) { 
            print(ID)
            if ((isTRUE(str_detect(f, "RUp16_PBMC_wk26_36_RSf16_PBMC_wk48"))) & (isTRUE(str_detect(ID, "WK26"))))  {
                contig.list[[paste0(ID, "v2")]] <- contigs
            } else if (isTRUE(str_detect(ID, "RUp16_WK36_D1"))) {
                print(paste0("skipping ", ID))
            } else {
                contig.list[[ID]] <- contigs
            }
        }
    }
}

results_folder <- file.path("RUp16_results", "VDJAnalysisNoWeek36D1_no_42-96")
generate_folder(results_folder)

saveRDS(contig.list, file.path(results_folder, "contiglist.rds"))
results_folder <- file.path("RUp16_results", "VDJAnalysisNoWeek36D1_no_42-96")
contig.list <- readRDS(file.path(results_folder, "contiglist.rds"))

## Final clone calls Threshold 0.85
results_folder <- file.path("RUp16_results", "VDJAnalysisNoWeek36D1_no_42-96")
generate_folder(results_folder)

combined.BCR <- combineBCR( contig.list, 
    samples = names(contig.list), 
    removeNA=TRUE, threshold =  0.85)
combined.ex.BCR <- combineExpandedBCR( contig.list, 
    samples = names(contig.list), 
    removeNA=TRUE, threshold =  0.85,
    filterNonproductive = TRUE)

totalscRep <- c()
for (s in names(combined.BCR)) {
    totalscRep <- rbind(totalscRep, combined.BCR[[s]])
}
totalscRep <- as.data.frame(totalscRep)
newids <- unique(totalscRep$CTstrict)
newidsdf <- as.data.frame(newids) 
newidsdf$LW_clones <- c(1:length(newids))

for (s in names(combined.BCR)) {
    tmp <- c()
    for (i in combined.BCR[[s]]$CTstrict) {
        tmp <- c(tmp, newidsdf[newidsdf$newids==i, "LW_clones"])
    }
    combined.BCR[[s]]$LW_clones <- tmp
}

clonalHomeostasis(combined.BCR, cloneCall = "LW_clones",cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded =
    1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5))
ggsave(file.path(results_folder, "ClonalHomeostasis.png"),width=7, height=4, dpi=300)

clonalQuant(combined.BCR, 
            cloneCall="LW_clones", 
            chain = "both", 
            scale = TRUE) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5))
ggsave(file.path(results_folder, "ClonalQuant.png"))

pl <- clonalCompare(combined.BCR, 
            top.clones = 50,
            #   relabel.clones = TRUE,
            cloneCall="LW_clones",  
            graph = "alluvial") + theme(legend.position = "none",
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5))
ggsave(file.path(results_folder, "ClonalCompare.png"),  width=4, height=4, dpi=300)

lwcoi <- c("684", "695", "1303","386","385","599", "659", "692","951","611","762", "651")

for (i in lwcoi) {
    pl <- clonalCompare(combined.BCR, 
                top.clones = 30,
                  highlight.clones= c(i),
                #   relabel.clones = TRUE,
                cloneCall="LW_clones",  
                graph = "alluvial") + theme(legend.position = "none") +
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5))
    ggsave(file.path(results_folder, paste0(i, "ClonalCompare.png")))
}
lwcoi <- c("684", "695", "1303","386","385","599")
pl <- clonalCompare(combined.BCR, 
    top.clones = 30,
        highlight.clones= lwcoi,
    #    relabel.clones = TRUE,
    cloneCall="LW_clones",  graph = "alluvial") + 
    theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5))
ggsave(file.path(results_folder, "ClonalCompare_allcoi_LWcall.png"),  width=4, height=4, dpi=300)
ggsave(file.path(results_folder, "ClonalCompare_allcoi_LWcall.pdf"),  width=4, height=4, dpi=300)

pl <- clonalCompare(combined.BCR, 
    top.clones = 30,
        highlight.clones= lwcoi,
    #    relabel.clones = TRUE,
    cloneCall="LW_clones",  graph = "alluvial") + 
    theme(, 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5))
ggsave(file.path(results_folder, "ClonalCompare_allcoi_LWcallLEGEND.png"),  width=8, height=4, dpi=300)
ggsave(file.path(results_folder, "ClonalCompare_allcoi_LWcallLEGEND.pdf"),  width=8, height=4, dpi=300)

saveRDS(combined.BCR, file.path(results_folder, "combinedBCR.rds"))
results_folder <- file.path("RUp16_results", "VDJAnalysisNoWeek36D1_no_42-96")
combined.BCR <- readRDS(file.path(results_folder, "combinedBCR.rds"))
# creaate pivot 
countlist <- list()
for (n in names(combined.BCR)) {
    df <- combined.BCR[[n]]
    df1 <- df %>% dplyr::count(LW_clones)
    df1 <- df1[!is.na(df1$LW_clones), ]
    # rownames(df1) <- as.character(df1$BK_clones)
    # df1$BK_clones <- NULL
    colnames(df1) <- c("LW_clones",n)
    countlist[[n]] <- as.data.frame(df1)
}
merged_df <- Reduce(function(x, y) merge(x, y, by = "LW_clones",all = TRUE), countlist)

merged_df_count <- merged_df
rownames(merged_df_count) <- as.character(merged_df_count$LW_clones)
merged_df_count$LW_clones <- NULL
merged_df_count[merged_df_count > 0] <- 1
merged_df_count[is.na(merged_df_count)] <- 0
A <- rowSums(merged_df_count)
A <- as.data.frame(A)
colnames(A) <- c("Week_count")
A$LW_clones <- rownames(A)
A <- A[rev(order(A$Week_count)),]
A <- A[A$Week_count>2, ]
A$LW_clones <- factor(A$LW_clones, levels=c(A$LW_clones))
ggplot(A, aes(x=LW_clones, y=Week_count)) + geom_col() +
theme_Publication() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3.5))
ggsave(file.path(results_folder, "bar_weekcount.png"), width=5.5, height = 3, dpi=300)

##Total Merged doc for Cindy and August 
total <- c()
totalex <- c()
for (i in names(combined.BCR)) {
    total <- rbind(total, combined.BCR[[i]])
    totalex <- rbind(totalex, combined.ex.BCR[[i]])
}
totalcom <- merge(total, totalex[,c("barcode", "CTnt","CTaa")], by="barcode")
colnames(totalcom)[15:16] <- c("CTnt.full","CTaa.full")
write.csv(total, file.path(results_folder, "merged_final.csv"), quote=F, row.names=F)

############## --Combined VDJB time points for RUp16 and run final analysis without wk36d1 remove group RUp16_PBMC_wk42_96--##############
results_folder <- file.path("RUp16_results", "VDJAnalysisNoWeek36D1_no_42-96")
combined.BCR <- readRDS(file.path(results_folder, "combinedBCR.rds"))

results_folder <- file.path("RUp16_results", "VDJAnalysisNoWeek36D1_no_42-96_final")
generate_folder(results_folder)
combined.BCR.comb <- combined.BCR
combined.BCR.comb[['RUp16_Wk26']] <- rbind(combined.BCR.comb[['RUp16_Wk26']], combined.BCR.comb[['RUp16_WK26v2']])
combined.BCR.comb[['RUp16_WK26v2']] <- NULL


clonalHomeostasis(combined.BCR.comb, cloneCall = "LW_clones",cloneSize = c(Rare = 0.001, Small = 0.02, Medium = 0.1, Large = 0.3, Hyperexpanded =
    1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5), axis.title.x=element_blank(), legend.text = element_text(size=8)) +
    scale_x_discrete(labels = c("Wk18", "Wk20", "Wk26", "Wk36","Wk48","Wk49","Wk80","Wk81","LN-Wk83", "Wk87", "Wk93"))
ggsave(file.path(results_folder, "Figure1A.png"),width=4.5, height=4, dpi=300)
ggsave(file.path(results_folder, "Figure1A.pdf"),width=4.5, height=4, dpi=300)

lwcoi <- c("684", "695", "1303","386","385","599")
pl <- clonalCompare(combined.BCR.comb, 
    top.clones = 30,
        highlight.clones= lwcoi, palette="viridis",
    #    relabel.clones = TRUE,
    cloneCall="LW_clones",  graph = "alluvial") + 
    theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5)) +
    scale_x_discrete(labels = c("Wk18", "Wk20", "Wk26", "Wk36","Wk48","Wk49","Wk80","Wk81","LN-Wk83", "Wk87", "Wk93"))
ggsave(file.path(results_folder, "ClonalCompare_allcoi_LWcall.png"),  width=3, height=3, dpi=300)
ggsave(file.path(results_folder, "ClonalCompare_allcoi_LWcall.pdf"),  width=3, height=3, dpi=300)

pl <- clonalCompare(combined.BCR.comb, 
    top.clones = 30,
        highlight.clones= lwcoi, palette="viridis",
    #    relabel.clones = TRUE,
    cloneCall="LW_clones",  graph = "alluvial") + 
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5))
ggsave(file.path(results_folder, "ClonalCompare_allcoi_LWcallLEGEND.png"),  width=8, height=4, dpi=300)
ggsave(file.path(results_folder, "ClonalCompare_allcoi_LWcallLEGEND.pdf"),  width=8, height=4, dpi=300)

pl <- clonalCompare(combined.BCR.comb, 
    top.clones = 30,
        highlight.clones= lwcoi,
    #    relabel.clones = TRUE,
    cloneCall="LW_clones",  graph = "alluvial") + 
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5))
ggsave(file.path(results_folder, "ClonalCompare_allcoi_LWcallLEGEND.pdf"),  width=8, height=4, dpi=300)

saveRDS(combined.BCR.comb, file.path(results_folder, "combinedBCRcomb.rds"))
results_folder <- file.path("RUp16_results", "VDJAnalysisNoWeek36D1_no_42-96_final")
combined.BCR.comb <- readRDS(file.path(results_folder, "combinedBCRcomb.rds"))
dfcellnum <- c()
for (n in names(combined.BCR.comb)) {
    dfcellnum <- rbind(dfcellnum, c(n, nrow(combined.BCR.comb[[n]])))
}
dfcellnum <- as.data.frame(dfcellnum)
colnames(dfcellnum) <- c("sample", "cell.number")
dfcellnum$cell.number <- as.numeric(dfcellnum$cell.number)
dfcellnum$sample <- factor(dfcellnum$sample, levels=dfcellnum$sample )
ggplot(dfcellnum, aes(x=sample, y=cell.number)) +
  geom_col() + theme_Publication() + theme(panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5),
    axis.text.y = element_text(size=7.5), axis.title.x=element_blank(), axis.title.y=element_text(size=9)) + labs(y="VDJ sequences") +
     scale_x_discrete(labels = c("Wk18", "Wk20", "Wk26", "Wk36","Wk48","Wk49","Wk80","Wk81","LN-Wk83", "Wk87", "Wk93"))
ggsave(file.path(results_folder, "Figure1B.png"), width=3, height = 3,units="in",dpi=300)
ggsave(file.path(results_folder, "Figure1B.pdf"), width=3, height = 3,units="in",dpi=300)

# creaate pivot 
countlist <- list()
for (n in names(combined.BCR.comb)) {
    df <- combined.BCR.comb[[n]]
    df1 <- df %>% dplyr::count(LW_clones)
    df1 <- df1[!is.na(df1$LW_clones), ]
    # rownames(df1) <- as.character(df1$BK_clones)
    # df1$BK_clones <- NULL
    colnames(df1) <- c("LW_clones",n)
    countlist[[n]] <- as.data.frame(df1)
}
merged_df <- Reduce(function(x, y) merge(x, y, by = "LW_clones",all = TRUE), countlist)

merged_df_count <- merged_df
rownames(merged_df_count) <- as.character(merged_df_count$LW_clones)
merged_df_count$LW_clones <- NULL
merged_df_count[merged_df_count > 0] <- 1
merged_df_count[is.na(merged_df_count)] <- 0
A <- rowSums(merged_df_count)
A <- as.data.frame(A)
colnames(A) <- c("Week_count")
A$LW_clones <- rownames(A)
A <- A[rev(order(A$Week_count)),]
A <- A[A$Week_count>2, ]
A$LW_clones <- factor(A$LW_clones, levels=c(A$LW_clones))

periods <- list("SOSIP 1"=c("RUp16_Wk18", "RUp16_Wk20"),
"SOSIP 2"=c("RUp16_Wk26", "RUp16_Wk36"),
"SOSIP 3"=c("RUp16_Wk48", "RUp16_Wk49"),
"SOSIP 4"=c("RUp16_Wk80", "RUp16_Wk81","RUp16_Wk83"),
"Longterm"=c("RUp16_Wk87","RUp16_Wk93"))
countperiods <- list()

for (i in names(periods)) {
    print(i)
    print(periods[[i]])
    df1 <- merged_df_count[, periods[[i]]]
    if (length(periods[[i]])>1){
        print(i)
        A<-rowSums(df1)
        A[A>0]<- 1
    } else { 
        A <- df1
    }

    A <- as.data.frame(A)
    colnames(A) <- c(i)
    print(dim(A))
    countperiods[[i]] <- A
}
combined_df <- do.call(cbind, countperiods)
Atotal <- rowSums(combined_df)
Atotal <- as.data.frame(Atotal)
colnames(Atotal) <- c("SOSIP_count")
Atotal$LW_clones <- rownames(Atotal)
Atotal <- Atotal[rev(order(Atotal$SOSIP_count)),]
Atotal <- Atotal[Atotal$SOSIP_count>2, ]
Atotal$LW_clones <- factor(Atotal$LW_clones, levels=c(Atotal$LW_clones))
LWclones_SOSIP <- Atotal$LW_clones
ggplot(Atotal, aes(x=LW_clones, y=SOSIP_count)) + geom_col() + labs(y="SOSIP count") + 
theme_Publication(base_size = 10) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 3.5),
    axis.title.x=element_blank(),panel.grid.major = element_blank())
ggsave(file.path(results_folder, "SupplFig3.png"), width=5.5, height = 3, dpi=300)
ggsave(file.path(results_folder, "SupplFig3.pdf"), width=5.5, height = 3, dpi=300)

cells <- c()
for (n in names(combined.BCR.comb)) {
    cells <- c(cells,combined.BCR.comb[[n]][combined.BCR.comb[[n]]$LW_clones %in% as.character(Atotal$LW_clones), "barcode"])
}
write.csv(Atotal, file.path(results_folder, "bar_SOSIP.csv"))
write.csv(cells, file.path(results_folder, "bar_SOSIP_cells.csv"))

dflength <-c()
for (n in names(combined.BCR.comb)) {
    df <- combined.BCR.comb[[n]]
    cdr3heavy <- df$cdr3_aa1
    cdr3light <- df$cdr3_aa2
    full.length <- c()
    heavy.length <- c() 
    light.length <- c() 
    for (i in cdr3heavy) {
        heavy.length <- c(heavy.length, nchar(i))
    }
    for (i in cdr3light) {
        light.length <- c(light.length, nchar(i))
    }

    dftmp <- cbind(df$LW_clones, heavy.length)
    dftmp <- cbind(dftmp, light.length)

    dflength <- rbind(dflength, dftmp)
}
dflength <- as.data.frame(dflength)
full.length <- rowSums(dflength[, c("heavy.length", "light.length")])
dflength$full.length <- full.length 
colnames(dflength) <- c("LW_clones", "heavy.length", "light.length","full.length")
dfmelt <-  reshape2::melt(dflength,id.vars="LW_clones")
p <- ggplot(dfmelt, aes(x=value)) + 
  geom_density() + theme_Publication() + facet_wrap(~variable, ncol=3, scales = "free")
ggsave(file.path(results_folder, "cdr3length.png"), width=7, height=3,dpi=300)
dfmeltsub <- dfmelt[dfmelt$variable!="full.length",]
dfmeltsub$variable <- as.character(dfmeltsub$variable)
dfmeltsub$variable[dfmeltsub$variable=="heavy.length"] <- "heavy"
dfmeltsub$variable[dfmeltsub$variable=="light.length"] <- "light"

dummyheavy <- dfmeltsub[dfmeltsub$variable=="heavy", ] %>%
  group_by(variable) %>%
  summarize(median = median(value))
p <- ggplot(dfmeltsub[dfmeltsub$variable=="heavy", ], aes(x=value)) +
  geom_density() + theme_Publication() + theme(panel.grid.major = element_blank(), axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5,size=7), 
  axis.text.y=element_text(size=8), plot.title = element_text(size=10), axis.title = element_text(size=8)) +
  geom_vline(data = dummyheavy, aes(xintercept = median), linetype="dashed", color="blue") +# facet_wrap(~variable, ncol=1, scale="free") +
  labs(title="heavy", x="CDR3 aa length") + scale_x_continuous(breaks = seq(min(dfmeltsub[dfmeltsub$variable=="heavy", "value"]), max(dfmeltsub[dfmeltsub$variable=="heavy", "value"]), by = 1))
ggsave(file.path(results_folder, "SupplFig2C.png"), width=3, height=3,dpi=300)
ggsave(file.path(results_folder, "SupplFig2C.pdf"), width=3, height=3,dpi=300)
dummylight <- dfmeltsub[dfmeltsub$variable=="light", ] %>%
  group_by(variable) %>%
  summarize(median = median(value))
p <- ggplot(dfmeltsub[dfmeltsub$variable=="light", ], aes(x=value)) +
  geom_density() + theme_Publication() + theme(panel.grid.major = element_blank(), axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5,size=8), 
  axis.text.y=element_text(size=8), plot.title = element_text(size=10), axis.title = element_text(size=8)) +
  geom_vline(data = dummylight, aes(xintercept = median), linetype="dashed", color="blue") +# facet_wrap(~variable, ncol=1, scale="free") +
  labs(title="light", x="CDR3 aa length") + scale_x_continuous(breaks = seq(min(dfmeltsub[dfmeltsub$variable=="light", "value"]), max(dfmeltsub[dfmeltsub$variable=="light", "value"]), by = 1))
ggsave(file.path(results_folder, "SupplFig2D.png"), width=3, height=3,dpi=300)
ggsave(file.path(results_folder, "SupplFig2D.pdf"), width=3, height=3,dpi=300)

############## --Get heavy light chain comninations --##############
results_folder <- file.path("RUp16_results","heavylightchaincombinations")
generate_folder(results_folder)
df <- read.csv("./RUp16_results/VDJAnalysisNoWeek36D1_no_42-96/merged_final.csv")
sarray <- strsplit(df$CTstrict, "\\.")
LWclones <- df$LW_clones
dfhl <- c()
counter=1
for (s in sarray) {
    x <- strsplit(s[[3]], "\\_")
    dfhl <- rbind(dfhl, c(x[[1]][1], s[[5]], LWclones[counter]))
    counter=counter+1
}
dfhl <- as.data.frame(dfhl)
colnames(dfhl) <- c("heavy", "light", "LW_clones")

dfhlcount <- dfhl %>% group_by(heavy, LW_clones) %>%
    dplyr::count()
dfhlcount <- as.data.frame(dfhlcount)
dfhlcount <- dfhl %>% group_by(heavy) %>%
    dplyr::count()
dfhlcount <- as.data.frame(dfhlcount)

dfhlcount <- dfhlcount[rev(order(dfhlcount$n)),]
dfhlcount$per <- (dfhlcount$n/sum(dfhlcount$n))*100
dfhlcountsub <- head(dfhlcount, 20)
dfhlcountsub$comb <- factor(dfhlcountsub$heavy, levels=dfhlcountsub$heavy)

dfheavy <- dfhlcountsub
dfhlcount <- dfhl %>% group_by(light, LW_clones) %>%
    dplyr::count()
dfhlcount <- as.data.frame(dfhlcount)
dfhlcount <- dfhl %>% group_by(light) %>%
    dplyr::count()
dfhlcount <- as.data.frame(dfhlcount)

dfhlcount <- dfhlcount[rev(order(dfhlcount$n)),]
dfhlcount$per <- (dfhlcount$n/sum(dfhlcount$n))*100
dfhlcountsub <- head(dfhlcount, 20)
dfhlcountsub$comb <- factor(dfhlcountsub$light, levels=dfhlcountsub$light)
dflight <- dfhlcountsub

dfheavy$chain  <- rep("heavy", nrow(dfheavy))
dflight$chain <- rep("light", nrow(dflight))
colnames(dfheavy)[1] <- "gene"
colnames(dflight)[1] <- "gene"
df <- rbind(dfheavy, dflight)
ggplot(df, aes(x=comb, y=per)) + geom_bar(stat = "identity",  fill="black") + facet_wrap(~chain, ncol=2, scale="free_x") + ylim(0,15) +
    theme_Publication()+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.title.x=element_blank()) + labs(y="percentage") + labs(y="percentage") + theme(panel.grid.major = element_blank(), 
    strip.background = element_rect(colour = "white", fill = "white"))
ggsave(file.path(results_folder, "SupplFig2A-B.png"),width=7, height=4, dpi=300)
ggsave(file.path(results_folder, "SupplFig2A-B.pdf"),width=7, height=4, dpi=300)
