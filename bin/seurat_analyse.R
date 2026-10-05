# #!/usr/bin/env Rscript
# # =============================================================
# # seurat_analyse.R —— 单细胞下游分析（Seurat 版）
# #
# # 读取 STARsolo 输出的表达矩阵，合并所有样本后完成：
# # 细胞过滤 -> 线粒体比例质控 -> 标准化 -> 高变基因 ->
# # PCA/UMAP 降维 -> 聚类 -> marker 基因导出 -> 保存 rds。
# #
# # 用法（由 modules/seurat 调用，也可在 RStudio 里独立调试）：
# #   Rscript seurat_analyse.R a.Solo.out b.Solo.out --out merged.rds
# # =============================================================

# #!/usr/bin/env Rscript

# suppressPackageStartupMessages({
#     library(Seurat)
# })

# args <- commandArgs(trailingOnly = TRUE)

# # ============================================================
# # 参数解析
# # ============================================================

# get_arg <- function(flag, default = NULL) {

#     idx <- which(args == flag)

#     if (length(idx) == 0) {
#         return(default)
#     }

#     if (idx[1] == length(args)) {
#         stop("Missing value for ", flag)
#     }

#     args[idx[1] + 1]
# }


# samplesheet <- get_arg("--samplesheet")
# out_rds     <- get_arg("--out", "merged_raw.rds")
# out_qc      <- get_arg("--qc", "qc_summary.csv")

# # 去掉已经解析的参数
# remove_args <- c(
#     "--samplesheet",
#     samplesheet,
#     "--out",
#     out_rds,
#     "--qc",
#     out_qc
# )

# h5_files <- args[!args %in% remove_args]

# if (length(h5_files) == 0) {
#     stop("No Cell Ranger H5 files found.")
# }


# # ============================================================
# # 读取 samplesheet
# # ============================================================

# meta <- read.csv(
#     samplesheet,
#     stringsAsFactors = FALSE,
#     check.names = FALSE
# )

# required_cols <- c(
#     "library",
#     "condition",
#     "fastq_dir"
# )

# if (!all(required_cols %in% colnames(meta))) {
#     stop(
#         "samplesheet must contain columns: ",
#         paste(required_cols, collapse = ", ")
#     )
# }


# # ============================================================
# # 根据 library 获得 condition
# # ============================================================

# condition_map <- setNames(
#     meta$condition,
#     meta$library
# )


# # ============================================================
# # 读取 6 个 Cell Ranger H5
# # ============================================================

# objects <- list()

# for (h5 in h5_files) {

#     message("Reading: ", h5)

#     library_name <- sub(
#         "_filtered_feature_bc_matrix\\.h5$",
#         "",
#         basename(h5)
#     )

#     message("Library: ", library_name)

#     counts <- Read10X_h5(h5)

#     # Seurat Read10X_h5 可能返回 list
#     if (is.list(counts)) {

#         if ("Gene Expression" %in% names(counts)) {
#             counts <- counts[["Gene Expression"]]
#         } else {
#             counts <- counts[[1]]
#         }
#     }

#     obj <- CreateSeuratObject(
#         counts = counts,
#         project = library_name,
#         min.cells = 3,
#         min.features = 0
#     )

#     obj <- RenameCells(
#         obj,
#         add.cell.id = library_name
#     )

#     obj$library <- library_name

#     obj$condition <- unname(
#         condition_map[library_name]
#     )

#     objects[[library_name]] <- obj
# }


# # ============================================================
# # 检查是否获得 6 个 library
# # ============================================================

# message(
#     "Libraries detected: ",
#     paste(names(objects), collapse = ", ")
# )

# if (length(objects) != 6) {
#     warning(
#         "Expected 6 libraries, but detected ",
#         length(objects)
#     )
# }


# # ============================================================
# # 合并
# # ============================================================

# seu <- objects[[1]]

# if (length(objects) > 1) {

#     seu <- merge(
#         x = seu,
#         y = objects[-1]
#     )
# }


# # ============================================================
# # 计算 mitochondrial percentage
# # ============================================================

# seu[["percent.mt"]] <- PercentageFeatureSet(
#     seu,
#     pattern = "^MT-"
# )


# # ============================================================
# # 生成 QC summary
# # ============================================================

# qc_summary <- data.frame(
#     cell = colnames(seu),
#     library = seu$library,
#     condition = seu$condition,
#     nCount_RNA = seu$nCount_RNA,
#     nFeature_RNA = seu$nFeature_RNA,
#     percent_mt = seu$percent.mt,
#     stringsAsFactors = FALSE
# )

# write.csv(
#     qc_summary,
#     out_qc,
#     row.names = FALSE
# )


# # ============================================================
# # library-level summary
# # ============================================================

# library_summary <- aggregate(
#     cbind(
#         nCount_RNA,
#         nFeature_RNA,
#         percent_mt
#     ) ~ library + condition,
#     data = qc_summary,
#     FUN = median
# )

# write.csv(
#     library_summary,
#     "qc_library_median.csv",
#     row.names = FALSE
# )


# # ============================================================
# # 当前阶段暂不做 QC 过滤
# #
# # 原因：
# # 我们先看 6 个 library 的 QC 分布，
# # 再决定统一的：
# #
# # nFeature_RNA
# # nCount_RNA
# # percent.mt
# #
# # 阈值。
# # ============================================================


# # ============================================================
# # 保存 merged Seurat object
# # ============================================================

# saveRDS(
#     seu,
#     out_rds
# )


# message(
#     "Finished: ",
#     ncol(seu),
#     " cells x ",
#     nrow(seu),
#     " features"
# )

# message(
#     "Saved: ",
#     out_rds
# )


suppressPackageStartupMessages({
    library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {

    idx <- which(args == flag)

    if (length(idx) == 0) {
        return(default)
    }

    if (idx[1] == length(args)) {
        stop("Missing value for ", flag)
    }

    args[idx[1] + 1]
}

samplesheet <- get_arg("--samplesheet")
out_rds     <- get_arg("--out", "merged_qc_filtered.rds")
qc_before   <- get_arg("--qc_before", "qc_summary_before.csv")
qc_after    <- get_arg("--qc_after", "qc_summary_after.csv")
qc_library  <- get_arg("--qc_library", "qc_library_summary.csv")

# ============================================================
# 剩余参数全部视为 H5 文件
# ============================================================

remove_idx <- c()

for (flag in c(
    "--samplesheet",
    "--out",
    "--qc_before",
    "--qc_after",
    "--qc_library"
)) {
    idx <- which(args == flag)
    if (length(idx) > 0) {
        remove_idx <- c(remove_idx, idx, idx + 1)
    }
}

h5_files <- args[-sort(unique(remove_idx))]

if (length(h5_files) == 0) {
    stop("No Cell Ranger H5 files found.")
}

message("Number of H5 files: ", length(h5_files))

# ============================================================
# samplesheet
# ============================================================

meta <- read.csv(
    samplesheet,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_cols <- c(
    "library",
    "condition"
)

if (!all(required_cols %in% colnames(meta))) {
    stop(
        "samplesheet must contain: ",
        paste(required_cols, collapse = ", ")
    )
}

condition_map <- setNames(
    meta$condition,
    meta$library
)

# ============================================================
# 读取 Cell Ranger matrix
# ============================================================

objects <- list()

for (h5 in h5_files) {

    message("==============================================")
    message("Reading: ", h5)

    library_name <- sub(
        "_filtered_feature_bc_matrix\\.h5$",
        "",
        basename(h5)
    )

    message("Library: ", library_name)

    counts <- Read10X_h5(h5)

    if (is.list(counts)) {

        if ("Gene Expression" %in% names(counts)) {
            counts <- counts[["Gene Expression"]]
        } else {
            counts <- counts[[1]]
        }
    }

    obj <- CreateSeuratObject(
        counts = counts,
        project = library_name,
        min.cells = 3,
        min.features = 0
    )

    obj <- RenameCells(
        obj,
        add.cell.id = library_name
    )

    obj$library <- library_name
    obj$condition <- unname(
        condition_map[library_name]
    )

    objects[[library_name]] <- obj

    message(
        "Cells: ", ncol(obj),
        " | Features: ", nrow(obj)
    )
}

# ============================================================
# 检查 library 数量
# ============================================================

if (length(objects) == 0) {
    stop("No Seurat objects created.")
}

if (length(objects) != 6) {
    warning(
        "Expected 6 libraries, detected ",
        length(objects)
    )
}

# ============================================================
# Merge
# ============================================================

message("Merging libraries...")

seu <- objects[[1]]

if (length(objects) > 1) {
    seu <- merge(
        x = seu,
        y = objects[-1]
    )
}

message(
    "Merged object: ",
    ncol(seu),
    " cells x ",
    nrow(seu),
    " genes"
)

# ============================================================
# Mitochondrial percentage
# ============================================================

seu[["percent.mt"]] <- PercentageFeatureSet(
    seu,
    pattern = "^MT-"
)

# ============================================================
# QC before filtering
# ============================================================

qc_before_df <- data.frame(
    cell = colnames(seu),
    library = seu$library,
    condition = seu$condition,
    nCount_RNA = seu$nCount_RNA,
    nFeature_RNA = seu$nFeature_RNA,
    percent_mt = seu$percent.mt,
    stringsAsFactors = FALSE
)

write.csv(
    qc_before_df,
    qc_before,
    row.names = FALSE
)

# ============================================================
# Library-level QC summary before filtering
# ============================================================

library_summary_before <- aggregate(
    cbind(
        nCount_RNA,
        nFeature_RNA,
        percent_mt
    ) ~ library + condition,
    data = qc_before_df,
    FUN = median
)

write.csv(
    library_summary_before,
    qc_library,
    row.names = FALSE
)

# ============================================================
# QC plots
# ============================================================

pdf(
    "qc_metrics.pdf",
    width = 14,
    height = 7
)

print(
    VlnPlot(
        seu,
        features = c(
            "nFeature_RNA",
            "nCount_RNA",
            "percent.mt"
        ),
        group.by = "library",
        pt.size = 0,
        ncol = 3
    )
)

print(
    FeatureScatter(
        seu,
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA"
    )
)

dev.off()

# ============================================================
# QC filtering
# ============================================================

message("Applying QC filters:")

message("nFeature_RNA >= 200")
message("nFeature_RNA <= 3000")
message("percent.mt <= 5")

seu <- subset(
    seu,
    subset =
        nFeature_RNA >= 200 &
        nFeature_RNA <= 3000 &
        percent.mt <= 5
)

message(
    "After QC: ",
    ncol(seu),
    " cells"
)

# ============================================================
# QC after
# ============================================================

qc_after_df <- data.frame(
    cell = colnames(seu),
    library = seu$library,
    condition = seu$condition,
    nCount_RNA = seu$nCount_RNA,
    nFeature_RNA = seu$nFeature_RNA,
    percent_mt = seu$percent.mt,
    stringsAsFactors = FALSE
)

write.csv(
    qc_after_df,
    qc_after,
    row.names = FALSE
)

# ============================================================
# Normalize
# ============================================================

message("NormalizeData...")

seu <- NormalizeData(
    seu,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
)

# ============================================================
# Variable features
# ============================================================

message("FindVariableFeatures...")

seu <- FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
)

# ============================================================
# Scale
# ============================================================

message("ScaleData...")

seu <- ScaleData(
    seu,
    features = VariableFeatures(seu),
    verbose = FALSE
)

# ============================================================
# PCA
# ============================================================

message("RunPCA...")

seu <- RunPCA(
    seu,
    features = VariableFeatures(seu),
    npcs = 50,
    verbose = FALSE
)

# ============================================================
# Elbow plot
# ============================================================

pdf(
    "pca_elbow.pdf",
    width = 8,
    height = 6
)

print(
    ElbowPlot(
        seu,
        ndims = 50
    )
)

dev.off()

# ============================================================
# Neighbors
# ============================================================

message("FindNeighbors...")

seu <- FindNeighbors(
    seu,
    dims = 1:30,
    verbose = FALSE
)

# ============================================================
# Clustering
# ============================================================

message("FindClusters...")

seu <- FindClusters(
    seu,
    resolution = 0.5,
    verbose = FALSE
)

# ============================================================
# UMAP
# ============================================================

message("RunUMAP...")

seu <- RunUMAP(
    seu,
    dims = 1:30,
    verbose = FALSE
)

# ============================================================
# UMAP: library
# ============================================================

pdf(
    "umap_library.pdf",
    width = 8,
    height = 7
)

print(
    DimPlot(
        seu,
        group.by = "library"
    )
)

dev.off()

# ============================================================
# UMAP: condition
# ============================================================

pdf(
    "umap_condition.pdf",
    width = 8,
    height = 7
)

print(
    DimPlot(
        seu,
        group.by = "condition"
    )
)

dev.off()

# ============================================================
# UMAP: clusters
# ============================================================

pdf(
    "umap_clusters.pdf",
    width = 8,
    height = 7
)

print(
    DimPlot(
        seu,
        group.by = "seurat_clusters",
        label = TRUE
    )
)

dev.off()

# ============================================================
# Cluster markers
# ============================================================

message("FindAllMarkers...")

markers <- FindAllMarkers(
    seu,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    verbose = FALSE
)

write.csv(
    markers,
    "cluster_markers.csv",
    row.names = FALSE
)

# ============================================================
# 保存
# ============================================================

saveRDS(
    seu,
    out_rds
)

message("==============================================")
message("Finished!")
message(
    "Final object: ",
    ncol(seu),
    " cells x ",
    nrow(seu),
    " genes"
)
message("Saved: ", out_rds)
