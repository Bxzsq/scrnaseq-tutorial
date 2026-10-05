// FASTQC —— 读段质量检查
// 每个样本运行一次（由输入通道的元素个数决定）

process FASTQC {

    tag "$library"

    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(library), val(condition), path(fastq_dir)

    output:
    path "fastqc/*_fastqc.html",
         emit: html

    path "fastqc/*_fastqc.zip",
         emit: zip

    script:
    """
    mkdir -p fastqc

    echo "Running FastQC for ${library}"

    fastqc \
        -t ${task.cpus} \
        --outdir fastqc \
        ${fastq_dir}/*.fastq.gz
    """
}
