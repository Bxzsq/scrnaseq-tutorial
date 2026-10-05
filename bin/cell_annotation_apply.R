#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
})

# ============================================================
# 1. 参数解析
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
    idx <- match(flag, args)

    if (is.na(idx) || idx == length(args)) {
        stop(paste0("Missing argument: ", flag))
    }

    args[idx + 1]
}

input_file <- get_arg("--input")
annotation_file <- get_arg("--annotation")
output_file <- get_arg("--output")
annotation_output <- get_arg("--annotation_output")
umap_output <- get_arg("--umap_output")

cat("============================================================\n")
cat("Apply manual cell type annotation\n")
cat("============================================================\n")

cat("Input Seurat object :", input_file, "\n")
cat("Annotation file     :", annotation_file, "\n")
cat("Output Seurat RDS   :", output_file, "\n")
cat("\n")

# ============================================================
# 2. 读取 Seurat 对象
# ============================================================

obj <- readRDS(input_file)

cat("Cells:", ncol(obj), "\n")
cat("Clusters:\n")
print(table(Idents(obj)))

# ============================================================
# 3. 读取人工注释表
# ============================================================

annotation_data <- read.csv(
    annotation_file,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
)

cat("\nAnnotation table:\n")
print(annotation_data)

# ============================================================
# 4. 检查 annotation 文件
# ============================================================

required_columns <- c("cluster", "celltype")

if (!all(required_columns %in% colnames(annotation_data))) {
    stop(
        "Annotation file must contain columns: ",
        paste(required_columns, collapse = ", ")
    )
}

annotation_data$cluster <- as.character(annotation_data$cluster)
annotation_data$celltype <- trimws(annotation_data$celltype)

# 检查空注释
bad_annotation <- is.na(annotation_data$celltype) |
    annotation_data$celltype == ""

if (any(bad_annotation)) {

    cat("\nERROR: Some clusters do not have cell type annotations:\n")

    print(annotation_data$cluster[bad_annotation])

    stop(
        "Please fill in the 'celltype' column for every cluster."
    )
}

# 检查重复 cluster
if (anyDuplicated(annotation_data$cluster)) {

    dup_clusters <- unique(
        annotation_data$cluster[
            duplicated(annotation_data$cluster)
        ]
    )

    stop(
        "Duplicated cluster IDs found: ",
        paste(dup_clusters, collapse = ", ")
    )
}

# ============================================================
# 5. 检查 Seurat cluster 与 annotation 是否匹配
# ============================================================

obj_clusters <- levels(Idents(obj))

annotation_clusters <- annotation_data$cluster

missing_annotation <- setdiff(
    obj_clusters,
    annotation_clusters
)

extra_annotation <- setdiff(
    annotation_clusters,
    obj_clusters
)

if (length(missing_annotation) > 0) {

    stop(
        "Missing annotation for cluster(s): ",
        paste(missing_annotation, collapse = ", ")
    )
}

if (length(extra_annotation) > 0) {

    warning(
        "Annotation file contains cluster(s) not present in Seurat object: ",
        paste(extra_annotation, collapse = ", ")
    )
}

# ============================================================
# 6. 建立 cluster -> cell type 映射
# ============================================================

annotation_map <- setNames(
    annotation_data$celltype,
    annotation_data$cluster
)

cat("\nCluster annotation mapping:\n")

for (cl in names(annotation_map)) {

    cat(
        "Cluster ",
        cl,
        " -> ",
        annotation_map[[cl]],
        "\n",
        sep = ""
    )
}

# ============================================================
# 7. 将 cluster 映射为 cell type
# ============================================================

cluster_ids <- as.character(Idents(obj))

celltype_labels <- unname(
    annotation_map[cluster_ids]
)

if (any(is.na(celltype_labels))) {

    stop(
        "Some cells could not be mapped to a cell type."
    )
}

obj$celltype_manual <- celltype_labels

# ============================================================
# 8. 设置 cell type factor 顺序
# ============================================================

celltype_levels <- unique(annotation_data$celltype)

obj$celltype_manual <- factor(
    obj$celltype_manual,
    levels = celltype_levels
)

# ============================================================
# 9. 输出细胞类型统计
# ============================================================

cat("\n============================================================\n")
cat("Cell type distribution\n")
cat("============================================================\n")

celltype_table <- as.data.frame(
    table(obj$celltype_manual)
)

colnames(celltype_table) <- c(
    "celltype",
    "n_cells"
)

print(celltype_table)

write.csv(
    celltype_table,
    annotation_output,
    row.names = FALSE
)

# ============================================================
# 10. 生成 Cell Type UMAP
# ============================================================

if (!"umap" %in% Reductions(obj)) {

    warning(
        "UMAP reduction not found. celltype_umap.pdf will not be generated."
    )

} else {

    p <- DimPlot(
        obj,
        reduction = "umap",
        group.by = "celltype_manual",
        label = TRUE,
        repel = TRUE
    ) +
        ggtitle("Manual Cell Type Annotation") +
        theme_classic()

    ggsave(
        filename = umap_output,
        plot = p,
        width = 10,
        height = 8
    )
}

# ============================================================
# 11. 保存最终 Seurat 对象
# ============================================================

saveRDS(
    obj,
    output_file
)

cat("\n============================================================\n")
cat("Manual annotation completed\n")
cat("============================================================\n")

cat("Annotated Seurat object:", output_file, "\n")
cat("Cell type table         :", annotation_output, "\n")
cat("Cell type UMAP          :", umap_output, "\n")
cat("\n")

cat("Final cell types:\n")
print(table(obj$celltype_manual))