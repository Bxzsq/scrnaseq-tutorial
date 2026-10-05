process CELLRANGER {

    tag "$library"
    container 'cellranger:10.1.0'
    cpus params.cellranger_cpus
    memory "${params.cellranger_mem_gb} GB"

    /*
     * Cell Ranger 很吃资源，一次只运行一个 library
     */
    maxForks 1

    publishDir "${params.outdir}/cellranger", mode: 'copy'

    input:
    tuple val(library), val(condition), path(fastq_dir, stageAs: 'fastqs')

    output:
    tuple val(library), val(condition),
          path("${library}_filtered_feature_bc_matrix.h5"),
          emit: matrix

    tuple val(library), val(condition),
          path("cr_${library}/outs/metrics_summary.csv"),
          emit: metrics

    path "cr_${library}/outs",
         emit: outs

    script:
    """
    echo "============================================================"
    echo "Cell Ranger"
    echo "Library   : ${library}"
    echo "Condition : ${condition}"
    echo "FASTQ dir : fastqs"
    echo "Reference : ${params.transcriptome}"
    echo "Output    : cr_${library}"
    echo "============================================================"

    cellranger count \
        --id="cr_${library}" \
        --transcriptome="${params.transcriptome}" \
        --fastqs="fastqs" \
        --sample="${library}" \
        --create-bam=${params.create_bam} \
        --localcores=${task.cpus} \
        --localmem=${params.cellranger_mem_gb}
      cp \
        "cr_${library}/outs/filtered_feature_bc_matrix.h5" \
        "${library}_filtered_feature_bc_matrix.h5"
    echo
    echo "Cell Ranger finished: ${library}"
    """
}