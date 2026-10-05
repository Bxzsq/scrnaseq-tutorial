// MULTIQC —— 把 FastQC、STAR 等工具的日志汇总成一份网页报告

process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.21--pyhdfd78af_0'
    conda 'bioconda::multiqc=1.21'   // 无容器时的备选（-profile conda）
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files   // collect() 之后的全部质控文件

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc . -n multiqc_report.html
    """
}
