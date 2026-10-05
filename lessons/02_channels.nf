// =============================================================
// 第 2 课：最常用的 channel 操作符
// 运行方式：nextflow run lessons/02_channels.nf
//
// 流程开发 80% 的工作是"数据整形"——把通道里的元素整理成
// 下游 process 需要的形状。这 4 个操作符最常用。
// =============================================================

nextflow.enable.dsl = 2

workflow {
    // ---- 1) of：手动列举元素 ----
    Channel.of(1, 2, 3)
        .view { v -> "of        -> ${v}" }

    // ---- 2) fromPath：按通配符匹配文件，产出 path 元素 ----
    // Channel.fromPath('data/*.fastq.gz')
    //     .view { f -> "fromPath  -> ${f}" }

    // ---- 3) splitCsv + map：解析 samplesheet（流程开发标配）----
    // 把 CSV 每行变成一个记录，再整形成 (样本名, R1, R2) 三元组
    // 注意：fromPath 的相对路径相对的是"你运行 nextflow 时所在的目录"，
    // 不是脚本所在目录——在别处运行报 No such file 就是这个问题
    Channel
        .fromPath('samplesheet.csv')
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
        .view { sample, r1, r2 -> "splitCsv  -> 样本 ${sample}: ${r1} + ${r2}" }

    // ---- 4) collect：把散开的元素收成一个列表 ----
    // 典型场景：上游逐样本并行，下游需要"全部样本到齐"再合并分析
    Channel.of('a', 'b', 'c')
        .collect()
        .view { list -> "collect   -> ${list}" }

    // 其他以后会常用到的：mix（合并通道）、join（按键配对）、
    // groupTuple（按键分组）、branch（按条件分流）
}
