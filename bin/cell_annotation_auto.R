#!/usr/bin/env Rscript

# ============================================================
# GSE138852 snRNA-seq
# Cell type annotation
#
# Input:
#   merged_qc_filtered.rds
#
# Output:
#   annotated.rds
#   celltype_umap.pdf
#   celltype_markers.csv
#   celltype_annotation.csv
#   module_scores.csv
#
# Main cell types:
#   Microglia
#   Astrocyte
#   Neuron
#   Oligodendrocyte
#   OPC
#   Endothelial
# ============================================================


suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
})


# ============================================================
# 1. Parse command line arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {

    idx <- which(args == flag)

    if (length(idx) == 0) {
        return(default)
    }

    if (idx[1] >= length(args)) {
        stop(
            paste0(
                "Missing value for argument: ",
                flag
            )
        )
    }

    args[idx[1] + 1]
}


input_file  <- get_arg("--input")
output_file <- get_arg("--output")
umap_file   <- get_arg("--umap")
marker_file <- get_arg("--markers")
annot_file  <- get_arg("--annotation")


if (is.null(input_file)) {
    stop("Missing --input")
}

if (is.null(output_file)) {
    output_file <- "annotated.rds"
}

if (is.null(umap_file)) {
    umap_file <- "celltype_umap.pdf"
}

if (is.null(marker_file)) {
    marker_file <- "celltype_markers.csv"
}

if (is.null(annot_file)) {
    annot_file <- "celltype_annotation.csv"
}


# ============================================================
# 2. Load Seurat object
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Cell type annotation\n")
cat("============================================================\n")

cat("Input:", input_file, "\n")

obj <- readRDS(input_file)

cat("Cells:", ncol(obj), "\n")
cat("Features:", nrow(obj), "\n")


# ============================================================
# 3. Check required metadata
# ============================================================

required_meta <- c(
    "library",
    "condition"
)

missing_meta <- setdiff(
    required_meta,
    colnames(obj@meta.data)
)

if (length(missing_meta) > 0) {

    warning(
        paste(
            "Missing metadata:",
            paste(missing_meta, collapse = ", ")
        )
    )

}


# ============================================================
# 4. Make sure clusters exist
# ============================================================

if (is.null(Idents(obj))) {
    stop("No identities found in Seurat object.")
}


cluster_ids <- levels(Idents(obj))

cat("\nClusters detected:\n")
print(cluster_ids)


# ============================================================
# 5. Canonical marker genes
#
# These are intentionally conservative marker sets.
# They are used for module scoring rather than hard filtering.
# ============================================================

marker_genes <- list(

    Microglia = c(
        "CX3CR1",
        "C1QA",
        "C1QB",
        "C1QC",
        "CSF1R",
        "AIF1",
        "P2RY12",
        "TMEM119",
        "CST3",
        "TYROBP"
    ),

    Astrocyte = c(
        "AQP4",
        "SLC1A2",
        "SLC1A3",
        "ALDH1L1",
        "GFAP",
        "GJA1",
        "SOX9",
        "Aldoc"
    ),

    Neuron = c(
        "RBFOX1",
        "SYT1",
        "SNAP25",
        "GRIA1",
        "GRIA2",
        "GRIN1",
        "GRIN2A",
        "GRIN2B",
        "CAMK2A",
        "MAP2"
    ),

    Oligodendrocyte = c(
        "MBP",
        "MOBP",
        "PLP1",
        "MAG",
        "MOG",
        "OPALIN",
        "CNP",
        "CLDN11",
        "MAL"
    ),

    OPC = c(
        "PDGFRA",
        "CSPG4",
        "VCAN",
        "OLIG1",
        "OLIG2",
        "SOX10",
        "BCAN"
    ),

    Endothelial = c(
        "CLDN5",
        "FLT1",
        "KDR",
        "VWF",
        "PECAM1",
        "EMCN",
        "RAMP2",
        "EPAS1"
    )
)


# ============================================================
# 6. Keep genes actually present in dataset
# ============================================================

feature_names <- rownames(obj)


marker_genes_present <- lapply(
    marker_genes,
    function(genes) {
        intersect(genes, feature_names)
    }
)


cat("\n")
cat("Marker genes found in dataset:\n")

for (ct in names(marker_genes_present)) {

    cat(
        sprintf(
            "%-18s %d / %d genes\n",
            ct,
            length(marker_genes_present[[ct]]),
            length(marker_genes[[ct]])
        )
    )

    if (length(marker_genes_present[[ct]]) == 0) {
        warning(
            paste0(
                "No marker genes found for cell type: ",
                ct
            )
        )
    }

}


# ============================================================
# 7. Add module scores
#
# AddModuleScore creates one score for each cell type.
# ============================================================

cat("\n")
cat("Calculating module scores...\n")


score_names <- c(
    "Microglia",
    "Astrocyte",
    "Neuron",
    "Oligodendrocyte",
    "OPC",
    "Endothelial"
)


score_columns <- character(0)


for (ct in score_names) {

    genes_use <- marker_genes_present[[ct]]

    if (length(genes_use) < 2) {

        warning(
            paste0(
                "Too few genes for module score: ",
                ct
            )
        )

        next
    }

    old_col <- paste0(ct, "1")

    obj <- AddModuleScore(
        object = obj,
        features = list(genes_use),
        name = paste0("Score_", ct, "_")
    )

    # AddModuleScore appends "1"
    created_col <- paste0(
        "Score_",
        ct,
        "_1"
    )

    if (created_col %in% colnames(obj@meta.data)) {

        score_columns <- c(
            score_columns,
            created_col
        )

    }

}


# ============================================================
# 8. Rename score columns to clean names
# ============================================================

score_rename <- c(
    "Score_Microglia_1" = "Score_Microglia",
    "Score_Astrocyte_1" = "Score_Astrocyte",
    "Score_Neuron_1" = "Score_Neuron",
    "Score_Oligodendrocyte_1" = "Score_Oligodendrocyte",
    "Score_OPC_1" = "Score_OPC",
    "Score_Endothelial_1" = "Score_Endothelial"
)


for (old_name in names(score_rename)) {

    new_name <- score_rename[[old_name]]

    if (old_name %in% colnames(obj@meta.data)) {

        obj@meta.data[[new_name]] <-
            obj@meta.data[[old_name]]

    }

}


score_columns_clean <- intersect(
    unname(score_rename),
    colnames(obj@meta.data)
)


# ============================================================
# 9. Cell-level preliminary annotation
#
# Assign the cell to the cell type with the highest module score.
# This is a preliminary automatic annotation.
# ============================================================

if (length(score_columns_clean) >= 2) {

    score_matrix <- obj@meta.data[
        ,
        score_columns_clean,
        drop = FALSE
    ]

    score_matrix_numeric <- as.matrix(
        score_matrix
    )

    max_idx <- max.col(
        score_matrix_numeric,
        ties.method = "first"
    )

    celltype_from_score <- colnames(
        score_matrix_numeric
    )[max_idx]

    celltype_from_score <- sub(
        "^Score_",
        "",
        celltype_from_score
    )

    obj$celltype_score <- celltype_from_score

} else {

    stop(
        "Not enough valid module scores for annotation."
    )

}


# ============================================================
# 10. Cluster-level annotation
#
# For each cluster:
#   calculate mean module score
#   assign cluster to highest scoring cell type
#
# This gives a preliminary cluster annotation.
# ============================================================

meta <- obj@meta.data

meta$cluster <- as.character(
    Idents(obj)
)


cluster_score_table <- meta %>%
    group_by(cluster) %>%
    summarise(
        across(
            all_of(score_columns_clean),
            mean,
            na.rm = TRUE
        ),
        .groups = "drop"
    )


cluster_score_matrix <- as.matrix(
    cluster_score_table[
        ,
        score_columns_clean,
        drop = FALSE
    ]
)

cluster_max_idx <- max.col(
    cluster_score_matrix,
    ties.method = "first"
)

cluster_annotation <- colnames(
    cluster_score_matrix
)[cluster_max_idx]

cluster_annotation <- sub(
    "^Score_",
    "",
    cluster_annotation
)


cluster_annotation_table <- data.frame(
    cluster = cluster_score_table$cluster,
    annotation = cluster_annotation,
    cluster_score_table[
        ,
        score_columns_clean,
        drop = FALSE
    ],
    stringsAsFactors = FALSE
)


# ============================================================
# 11. Add cluster annotation to each cell
# ============================================================

cluster_map <- setNames(
    cluster_annotation_table$annotation,
    cluster_annotation_table$cluster
)


obj$celltype <- unname(
    cluster_map[
        as.character(
            Idents(obj)
        )
    ]
)


# ============================================================
# 12. Print annotation result
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Cluster annotation\n")
cat("============================================================\n")

print(
    cluster_annotation_table
)


cat("\n")
cat("Cell type counts:\n")

print(
    table(
        obj$celltype
    )
)


# ============================================================
# 13. Find marker genes for each cluster
# ============================================================

cat("\n")
cat("Finding cluster marker genes...\n")


markers <- FindAllMarkers(
    object = obj,
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.25
)


# ============================================================
# 14. Add annotation to marker table
# ============================================================

if (nrow(markers) > 0) {

    markers$celltype <- unname(
        cluster_map[
            as.character(markers$cluster)
        ]
    )

}


# ============================================================
# 15. Save marker table
# ============================================================

write.csv(
    markers,
    marker_file,
    row.names = FALSE
)


# ============================================================
# 16. Save cluster annotation table
# ============================================================

write.csv(
    cluster_annotation_table,
    annot_file,
    row.names = FALSE
)


# ============================================================
# 17. Module score UMAP
#
# Plot the six module scores.
# ============================================================

pdf(
    umap_file,
    width = 12,
    height = 8
)


if (length(score_columns_clean) > 0) {

    plots <- FeaturePlot(
        object = obj,
        features = score_columns_clean,
        reduction = "umap",
        ncol = 3
    )

    print(plots)

}


# ============================================================
# 18. Cell type annotation UMAP
# ============================================================

p_celltype <- DimPlot(
    object = obj,
    reduction = "umap",
    group.by = "celltype",
    label = TRUE,
    repel = TRUE
) +
    ggtitle(
        "Cell type annotation"
    )


print(
    p_celltype
)


# ============================================================
# 19. Cluster UMAP
# ============================================================

p_cluster <- DimPlot(
    object = obj,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
) +
    ggtitle(
        "Seurat clusters"
    )


print(
    p_cluster
)


dev.off()


# ============================================================
# 20. Save cell-level module scores
# ============================================================

module_score_output <- obj@meta.data[
    ,
    c(
        "library",
        "condition",
        "seurat_clusters",
        "celltype",
        score_columns_clean
    ),
    drop = FALSE
]


write.csv(
    module_score_output,
    "module_scores.csv",
    row.names = TRUE
)


# ============================================================
# 21. Save annotated Seurat object
# ============================================================

saveRDS(
    obj,
    output_file
)


# ============================================================
# 22. Final summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Annotation finished\n")
cat("============================================================\n")

cat(
    "Annotated Seurat object:",
    output_file,
    "\n"
)

cat(
    "Marker table:",
    marker_file,
    "\n"
)

cat(
    "Annotation table:",
    annot_file,
    "\n"
)

cat(
    "UMAP:",
    umap_file,
    "\n"
)

cat(
    "Total cells:",
    ncol(obj),
    "\n"
)

cat("\nCell type distribution:\n")

print(
    prop.table(
        table(
            obj$celltype
        )
    )
)

cat("\n")
cat("Done.\n")