#!/usr/bin/env Rscript

# ============================================================
# GSE138852 snRNA-seq
# Manual cell-type annotation based on cluster marker genes
#
# Workflow:
#   merged_qc_filtered.rds
#          ↓
#   FindAllMarkers()
#          ↓
#   Top marker genes per cluster
#          ↓
#   DotPlot / Heatmap
#          ↓
#   Manual annotation by researcher
#          ↓
#   annotated_manual.rds
#
# IMPORTANT:
#   This script DOES NOT automatically assign cell types.
#   Cell-type labels are provided manually through
#   cluster_annotation.csv
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


input_file <- get_arg(
    "--input"
)

output_file <- get_arg(
    "--output",
    "annotated_manual.rds"
)

marker_file <- get_arg(
    "--markers",
    "cluster_markers.csv"
)

top_marker_file <- get_arg(
    "--top_markers",
    "top20_markers.csv"
)

dotplot_file <- get_arg(
    "--dotplot",
    "marker_dotplot.pdf"
)

heatmap_file <- get_arg(
    "--heatmap",
    "marker_heatmap.pdf"
)

annotation_file <- get_arg(
    "--annotation",
    "cluster_annotation.csv"
)


if (is.null(input_file)) {
    stop("Missing --input")
}


# ============================================================
# 2. Load Seurat object
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Manual cell-type annotation\n")
cat("============================================================\n")

cat("Input:", input_file, "\n")


obj <- readRDS(input_file)


cat("Cells:",ncol(obj),"\n")

cat("Features:",nrow(obj),"\n")

obj <- JoinLayers(obj)

# ============================================================
# 3. Check clustering
# ============================================================

if (!"seurat_clusters" %in% colnames(obj@meta.data)) {

    stop(
        "seurat_clusters not found in metadata."
    )

}


Idents(obj) <- "seurat_clusters"


clusters <- levels(
    Idents(obj)
)


cat("\n")
cat("Clusters:\n")
print(clusters)


# ============================================================
# 4. Find markers for all clusters
# ============================================================

cat("\n")
cat("Finding marker genes for all clusters...\n")


markers <- FindAllMarkers(
    object = obj,

    only.pos = TRUE,

    min.pct = 0.10,

    logfc.threshold = 0.25
)


if (nrow(markers) == 0) {

    stop(
        "No marker genes were detected."
    )

}


# ============================================================
# 5. Save complete marker table
# ============================================================

write.csv(
    markers,
    marker_file,
    row.names = FALSE
)


cat(
    "\nAll marker genes saved to:",
    marker_file,
    "\n"
)


# ============================================================
# 6. Obtain top 20 markers for each cluster
#
# Ranking:
#   1. avg_log2FC
#   2. pct.1
# ============================================================

top_markers <- markers %>%
    group_by(cluster) %>%
    arrange(
        desc(avg_log2FC),
        desc(pct.1)
    ) %>%
    slice_head(
        n = 20
    ) %>%
    ungroup()


write.csv(
    top_markers,
    top_marker_file,
    row.names = FALSE
)


cat(
    "Top marker genes saved to:",
    top_marker_file,
    "\n"
)


# ============================================================
# 7. Print top markers to terminal
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Top 20 marker genes per cluster\n")
cat("============================================================\n")


for (cl in clusters) {

    cat("\n")
    cat(
        "---------------- Cluster ",
        cl,
        " ----------------\n",
        sep = ""
    )

    genes <- top_markers %>%
        filter(
            cluster == cl
        ) %>%
        pull(
            gene
        )

    cat(
        paste(
            genes,
            collapse = ", "
        ),
        "\n"
    )
}


# ============================================================
# 8. Make a compact marker matrix
#
# Used for heatmap.
# ============================================================

heatmap_genes <- top_markers %>%
    group_by(cluster) %>%
    slice_head(
        n = 10
    ) %>%
    pull(gene) %>%
    unique()


heatmap_genes <- intersect(
    heatmap_genes,
    rownames(obj)
)


cat(
    "\nNumber of heatmap marker genes:",
    length(heatmap_genes),
    "\n"
)


# ============================================================
# 9. Marker DotPlot
#
# Important:
# We DO NOT use predefined cell-type markers here.
# We simply visualize cluster marker genes.
# ============================================================

pdf(
    dotplot_file,
    width = 16,
    height = 12
)


if (length(heatmap_genes) > 0) {

    p_dot <- DotPlot(
        object = obj,
        features = heatmap_genes,
        group.by = "seurat_clusters"
    ) +
        RotatedAxis() +
        ggtitle(
            "Cluster marker genes"
        ) +
        theme(
            axis.text.x = element_text(
                angle = 45,
                hjust = 1
            )
        )

    print(
        p_dot
    )

}


dev.off()


# ============================================================
# 10. Marker Heatmap
# ============================================================

pdf(
    heatmap_file,
    width = 14,
    height = 12
)


if (length(heatmap_genes) > 0) {
    # Ensure all marker genes are available in scale.data
    heatmap_genes <- intersect(
       heatmap_genes,
       rownames(obj))

    obj <- ScaleData(
        obj,
        features = heatmap_genes,
        verbose = FALSE)
        
    p_heat <- DoHeatmap(
        object = obj,
        features = heatmap_genes,
        group.by = "seurat_clusters"
    ) +
        ggtitle(
            "Top cluster marker genes"
        )

    print(
        p_heat
    )

}


dev.off()


# ============================================================
# 11. Create empty annotation table
#
# IMPORTANT:
# The labels below are intentionally EMPTY.
#
# You will manually fill:
#
# cluster,celltype
#
# Example:
#
# 0,Microglia
# 1,Oligodendrocyte
# 2,Astrocyte
#
# ============================================================


annotation_table <- data.frame(
    cluster = clusters,
    celltype = "",
    stringsAsFactors = FALSE
)


write.csv(
    annotation_table,
    annotation_file,
    row.names = FALSE
)


cat("\n")
cat(
    "Annotation template created:",
    annotation_file,
    "\n"
)


# ============================================================
# 12. Check whether annotation file has been manually filled
# ============================================================

annotation_data <- read.csv(
    annotation_file,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
)

# Treat NA and empty strings as unfilled annotations
annotation_filled <- !is.na(annotation_data$celltype) &
    trimws(annotation_data$celltype) != ""

if (!any(annotation_filled)) {

    cat("\n")
    cat("============================================================\n")
    cat("Manual annotation required\n")
    cat("============================================================\n")

    cat(
        "Please edit:",
        annotation_file,
        "\n"
    )

    cat(
        "Fill in the 'celltype' column for every cluster.\n"
    )

    cat(
        "Then rerun this script.\n"
    )

    cat(
        "\nNo automatic cell-type annotation was performed.\n"
    )

    # Still save the original object.
    saveRDS(
        obj,
        output_file
    )

} else {


    # ========================================================
    # 13. Validate annotation table
    # ========================================================

    if (
        any(
            annotation_data$celltype == ""
        )
    ) {

        stop(
            "Some clusters still have empty celltype labels."
        )

    }


    if (
        !all(
            as.character(
                annotation_data$cluster
            ) %in%
            clusters
        )
    ) {

        stop(
            "Annotation table contains unknown cluster IDs."
        )

    }


    # ========================================================
    # 14. Create cluster → cell type mapping
    # ========================================================

    cluster_map <- setNames(
        annotation_data$celltype,
        annotation_data$cluster
    )


    # ========================================================
    # 15. Add manual annotation to Seurat object
    # ========================================================

    obj$celltype_manual <- unname(
        cluster_map[
            as.character(
                Idents(obj)
            )
        ]
    )


    # ========================================================
    # 16. Print final annotation
    # ========================================================

    cat("\n")
    cat("============================================================\n")
    cat("Final manual annotation\n")
    cat("============================================================\n")


    print(
        annotation_data
    )


    cat("\n")
    cat("Cell numbers by cell type:\n")


    print(
        table(
            obj$celltype_manual
        )
    )


    # ========================================================
    # 17. Save annotated object
    # ========================================================

    saveRDS(
        obj,
        output_file
    )


    cat("\n")
    cat(
        "Annotated Seurat object saved:",
        output_file,
        "\n"
    )

}


# ============================================================
# 18. Finished
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Manual annotation workflow finished.\n")
cat("============================================================\n")