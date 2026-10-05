// STARSOLO —— 比对 + 定量，产出 基因x细胞 表达矩阵
// STARsolo 是 Cell Ranger 的开源替代，结果高度一致且快很多。
// 10x 数据约定：R1 = barcode+UMI，R2 = cDNA 插入序列。
// 对 10x 3' v3 / v3.1：CB whitelist 应为 3M-february-2018.txt，UMI 长度应为 12。

process STARSOLO {
    tag "$sample"
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_2'
    conda 'bioconda::star=2.7.11b'   // 无容器时的备选（-profile conda）
    cpus 8
    memory '32 GB'   // STAR 需要装载基因组索引，内存大头在这里
    publishDir "${params.outdir}/starsolo", mode: 'copy'

    input:
    tuple val(sample), path(read1), path(read2)
    path star_index

    output:
    // Solo.out/Gene/filtered/ 里就是 10x 格式的表达矩阵
    // （matrix.mtx + barcodes.tsv + features.tsv）
    tuple val(sample), path("${sample}.Solo.out"), emit: matrix
    path "${sample}.Log.final.out",                emit: log

    script:
    // 注意 readFilesIn 的参数顺序：cDNA 在前，barcode 在后
    """
    STAR --runThreadN ${task.cpus} \\
         --genomeDir ${star_index} \\
         --readFilesIn ${read2} ${read1} \\
         --readFilesCommand zcat \\
         --soloType CB_UMI_Simple \\
         --soloCBstart 1\\
         --soloCBlen 16\\
         --soloUMIstart 17\\
         --soloUMIlen 12\\
         --soloCBwhitelist ${params.whitelist} \\
         --outFileNamePrefix ${sample}.
    """
}
