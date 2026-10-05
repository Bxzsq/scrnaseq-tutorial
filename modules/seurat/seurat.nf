// // SEURAT —— 下游分析：合并样本、过滤、降维聚类、marker 基因
// // 输入是 collect() 之后的"全部样本矩阵目录列表"，因此只运行一次。
// // 实际分析逻辑在 bin/seurat_analyse.R ——
// // 放在项目 bin/ 目录下的可执行脚本会被 Nextflow 自动加入 PATH，
// // 这样你用熟悉的 R 写分析，Nextflow 只管调度与环境。

// process SEURAT {

//     tag "Seurat_merge"

//     container 'satijalab/seurat:5.5.1'

//     publishDir "${params.outdir}/seurat", mode: 'copy'

//     input:
//     path matrix_h5
//     path samplesheet
//     path analysis_script

//     output:
//     path "merged_raw.rds",
//          emit: rds

//     path "qc_summary.csv",
//          emit: qc

//     path "qc_library_median.csv"

//     script:
//     """
//     echo "============================================================"
//     echo "Seurat"
//     echo "Matrices:"
//     ls -lh *.h5
//     echo "============================================================"

//     Rscript ${analysis_script} \
//         --samplesheet ${samplesheet} \
//         --out merged_raw.rds \
//         --qc qc_summary.csv \
//         ${matrix_h5}
//     """
// }




process SEURAT {

    tag "Seurat_downstream"

    container 'satijalab/seurat:5.5.1'

    cpus 4
    memory '32 GB'

    publishDir "${params.outdir}/seurat", mode: 'copy'

    input:
    path matrix_h5
    path samplesheet
    path analysis_script

    output:
    path "merged_qc_filtered.rds",
         emit: filtered_rds

    path "qc_summary_before.csv",
         emit: qc_before

    path "qc_summary_after.csv",
         emit: qc_after

    path "qc_library_summary.csv",
         emit: qc_library

    path "qc_metrics.pdf",
         emit: qc_plot

    path "pca_elbow.pdf"

    path "umap_library.pdf"
    path "umap_condition.pdf"
    path "umap_clusters.pdf"

    path "cluster_markers.csv"

    script:
    """
    echo "============================================================"
    echo "Seurat downstream analysis"
    echo "============================================================"

    echo "Input H5 files:"
    ls -lh *.h5

    Rscript ${analysis_script} \
        --samplesheet ${samplesheet} \
        --out merged_qc_filtered.rds \
        --qc_before qc_summary_before.csv \
        --qc_after qc_summary_after.csv \
        --qc_library qc_library_summary.csv \
        ${matrix_h5}

    echo
    echo "Seurat finished."
    """
}