nextflow.enable.dsl = 2

include { FASTQC } from './modules/fastqc/fastqc.nf'
include { CELLRANGER } from './modules/cellranger/cellranger.nf'
include { SEURAT } from './modules/seurat/seurat.nf'
include { ANNOTATION } from './modules/annotation/annotation.nf'
include { ANNOTATION_APPLY } from './modules/annotation/annotation_apply.nf'
include { ANNOTATION_AUTO } from './modules/annotation/annotation_auto.nf'

workflow {


    /*1. 读取 samplesheet*/

    ch_samples = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.library.toString(),
                row.condition.toString(),
                file(row.fastq_dir.toString(), checkIfExists: true)
            )
        }


    /*
     * ============================================================
     * 2. FASTQ 模式
     *
     * FASTQ
     *   ↓
     * FASTQC
     *   ↓
     * CELLRANGER
     *   ↓
     * SEURAT
     *   ↓
     * ANNOTATION_PREPARE
     * ============================================================
     */

    if (params.input_mode == 'fastq') {

        FASTQC(ch_samples)

        CELLRANGER(ch_samples)

        ch_matrix = CELLRANGER.out.matrix
            .map { library, condition, h5 -> h5 }
            .collect()

        SEURAT(
            ch_matrix,
            file(params.samplesheet),
            file("$projectDir/bin/seurat_analyse.R")
        )

        ch_seurat = SEURAT.out.filtered_rds

        ANNOTATION(
            ch_seurat,
            file("$projectDir/bin/cell_annotation_manual.R")
        )
    }


    /*
     * ============================================================
     * 3. H5 模式
     *
     * 已经存在 Cell Ranger H5
     *
     * H5
     *   ↓
     * SEURAT
     *   ↓
     * ANNOTATION_PREPARE
     * ============================================================
     */

    else if (params.input_mode == 'h5') {

        ch_matrix = Channel
            .fromPath(
                "${params.h5_dir}/*_filtered_feature_bc_matrix.h5",
                checkIfExists: true
            )
            .collect()

        SEURAT(
            ch_matrix,
            file(params.samplesheet),
            file("$projectDir/bin/seurat_analyse.R")
        )

        ch_seurat = SEURAT.out.filtered_rds

        ANNOTATION(
            ch_seurat,
            file("$projectDir/bin/cell_annotation_manual.R")
        )
    }


    /*
     * ============================================================
     * 4. annotation_prepare 模式
     *
     * 只做：
     *
     * merged_qc_filtered.rds
     *        ↓
     * FindAllMarkers
     *        ↓
     * marker files
     *        ↓
     * cluster_annotation.csv
     *
     * 不应用人工标签
     * ============================================================
     */

    else if (params.input_mode == 'annotation_manual') {

        ch_seurat = Channel
            .fromPath(
                params.seurat_rds,
                checkIfExists: true
            )

        ANNOTATION(
            ch_seurat,
            file("$projectDir/bin/cell_annotation_manual.R")
        )
    }



    else if (params.input_mode == 'annotation_auto') {

    ch_matrix = Channel
        .fromPath(
            "${params.h5_dir}/*_filtered_feature_bc_matrix.h5",
            checkIfExists: true
        )
        .collect()

    ANNOTATION_AUTO(
        ch_matrix,
        file(params.samplesheet),
        file("$projectDir/bin/celltypist_annotation.py")
    )
}
    /*
     * ============================================================

     * merged_qc_filtered.rds
     *        +
     * cluster_annotation.csv
     *        ↓
     * ANNOTATION_APPLY
     *        ↓
     * annotated_manual.rds
     * ============================================================
     */

    else if (params.input_mode == 'annotation_apply') {

        ch_seurat = Channel
            .fromPath(
                params.seurat_rds,
                checkIfExists: true
            )

        ch_annotation = Channel
            .fromPath(
                params.annotation_file,
                checkIfExists: true
            )

        ANNOTATION_APPLY(
            ch_seurat,
            ch_annotation,
            file("$projectDir/bin/cell_annotation_apply.R")
        )
    }






    /*
     * ============================================================
     * 7. 非法模式
     * ============================================================
     */

    else {

        error """
        Invalid input_mode: ${params.input_mode}

        Valid values:

          --input_mode fastq
          --input_mode h5
          --input_mode annotation_manual
          --input_mode annotation_apply
          --input_mode annotation_auto
        """
    }
}

