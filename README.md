# 从零开始：用 Nextflow 搭建单细胞 RNA-seq 分析流程

> 面向人群：有自己的 Linux 服务器，平时在 RStudio 里用 Seurat 做单细胞分析，但没接触过 Nextflow。
> 目标：跑通一条 `FASTQ → 表达矩阵 → Seurat 聚类 → rds` 的完整流程，下游分析就是你熟悉的 R。

## 你将做出什么

```
samplesheet.csv
    ├─> FASTQC    读段质控（逐样本并行）
    ├─> STARSOLO  比对 + 定量，产出 基因×细胞矩阵（逐样本并行）
    ├─> SEURAT    合并样本：过滤/降维/聚类/marker（R 脚本，只运行一次）
    └─> MULTIQC   汇总所有质控日志为一份网页报告
```

核心思想：**分析逻辑用你擅长的 R 写（`bin/seurat_analyse.R`），Nextflow 只负责并行调度、环境和容错。**

## 第 0 步：准备运行环境（Linux 服务器，10 分钟）

```bash
# 1. Nextflow 需要 Java 17+，先检查
java -version

# 2. 安装 Nextflow（二选一）
curl -s https://get.nextflow.io | bash && mv nextflow ~/bin/   # 官方脚本
# 或
conda create -n nextflow -c bioconda -c conda-forge nextflow   # conda 方式

nextflow -version   # 验证
```

然后确认软件环境的运行方式，三种选一种（对应 `-profile`）：

| 方式 | 适用场景 | 检查命令 |
|------|----------|----------|
| **docker** | 你有 root 或在 docker 用户组 | `docker run hello-world` |
| **singularity** | 共享集群/HPC 最常见 | `singularity --version` |
| **conda** | 什么权限都没有的退路（模块里已写好 conda 指令，自动用 bioconda 装环境） | `conda --version` |

> 如果服务器上 RStudio 那套 R 环境已装好 Seurat，也可以不用任何容器——但建议至少用 conda 保持可重复性。

## 第 1 步：10 分钟理解两个核心概念

Nextflow 只有两块积木：

- **process**：一个任务 = 输入 + 输出 + 一段脚本（Bash/Python/R 都行）
- **channel**：连接 process 的异步队列，数据在里面流动

把项目传到服务器后，在项目根目录跑两个小课：

```bash
nextflow run lessons/01_hello.nf
nextflow run lessons/02_channels.nf
```

观察输出：`01` 里 SAY_HELLO 自动执行了 3 次——**通道里有几个元素，process 就运行几次**，这就是 Nextflow 的并行方式，你不需要写任何循环或并发代码。`02` 演示了最常用的 4 个 channel 操作符，其中 `splitCsv + map`（解析样本表）和 `collect`（先并行后汇总）马上就会在主流程里用到。

## 第 2 步：认识项目结构

```
scrnaseq-tutorial/
├── main.nf            # 编排：把模块连成流程（先读我）
├── nextflow.config    # 参数 + 运行环境 profile
├── samplesheet.csv    # 样本表：sample,fastq_1,fastq_2
├── modules/           # 每个工具一个可复用模块，注释详细
│   ├── fastqc/main.nf
│   ├── starsolo/main.nf
│   ├── seurat/main.nf
│   └── multiqc/main.nf
├── bin/
│   └── seurat_analyse.R    # 下游分析本体（你熟悉的 Seurat）
└── lessons/           # 第 1 步的两个小课
```

建议阅读顺序：`main.nf` → `modules/fastqc/main.nf`（最简单的模块）→ `modules/starsolo/main.nf` → `bin/seurat_analyse.R` → `nextflow.config`。

几个关键设计，看懂它们就看懂了 Nextflow 流程开发：

| 设计 | 位置 | 说明 |
|------|------|------|
| 样本表驱动 | `main.nf` 第 1 步 | `splitCsv` 解析 CSV，多样本自动并行 |
| 先并行后汇总 | `main.nf` 第 4 步 | `collect()` 等所有矩阵到齐再做合并分析 |
| 分析脚本外置 | `bin/` + `modules/seurat` | bin/ 下脚本自动加入 PATH，分析逻辑用 R 写 |
| 环境收进 profile | `nextflow.config` | 换环境只换 `-profile`，不改流程代码 |

## 第 3 步：准备数据与参考基因组

**测试数据**：10x Genomics 官网（10xgenomics.com/datasets）下载任意 PBMC 数据集的 FASTQ（如 "1k PBMCs from a Healthy Donor, v3 chemistry"），改好 `samplesheet.csv` 里的路径。

**STAR 索引**（只需构建一次）：

```bash
STAR --runMode genomeGenerate \
     --genomeDir /ref/star_index \
     --genomeFastaFiles /ref/genome.fa \
     --sjdbGTFfile /ref/genes.gtf \
     --sjdbOverhang 99 \
     --runThreadN 8
```

**barcode 白名单**：10x 3' v3 / v3.1 对应 `3M-february-2018.txt.gz`（给 STARsolo 用前先解压成 `.txt`）；`737K-august-2016.txt` 是旧版 v2，不适合这个 pbmc_1k_v3 数据集。另一个容易漏掉的参数是 **UMI 长度**：v3 / v3.1 应为 `12`，脚手架里已经在 STARsolo 模块中加好了 `--soloUMIlen ${params.solo_umi_len}`。把 `nextflow.config` 里的 `star_index`、`whitelist` 改成你的实际路径即可。

## 第 4 步：运行！

```bash
# 按你服务器的条件选一个 profile
nextflow run main.nf -profile docker
nextflow run main.nf -profile singularity        # 或集群: -profile singularity,slurm
nextflow run main.nf -profile conda              # 无 root 的退路

# 常用附加参数
nextflow run main.nf -profile singularity \
    --input my_samples.csv \     # 覆盖 config 里的参数
    -resume \                    # 失败后从断点续跑（最重要的参数！）
    -with-report report.html \   # 资源使用报告
    -with-timeline timeline.html # 任务时间线
```

结果在 `results/` 下：`fastqc/`、`starsolo/`（每样本矩阵）、`seurat/merged.rds` + `figures/` + `cluster_markers.csv`、`multiqc/multiqc_report.html`。

## 第 5 步：在 RStudio 里接着分析

流程跑到这里交付的是一个**已经质控、聚类好的 Seurat 对象**，正好接进你日常的交互式分析：

```r
seu <- readRDS("results/seurat/merged.rds")
DimPlot(seu, label = TRUE)                       # 检查聚类
VlnPlot(seu, features = c("nFeature_RNA", "percent.mt"))
# 然后按你自己的习惯做细胞注释、差异分析、画图……
```

这正是推荐的分工：**Nextflow 管"重计算、要并行、要可重复"的部分（比对、定量、批量质控），RStudio 管"需要边看边调"的部分（注释、可视化、探索）。** 调好了 `seurat_analyse.R` 里的过滤阈值，就 `-resume` 重跑一次，全团队拿到完全一致的结果。

如果你用的是 RStudio Server，流程跑在后台，浏览器里随时读结果，互不影响。

## 第 6 步：排错（迟早要用）

- 每个任务在 `work/` 下有独立目录，里面有 `.command.sh`（实际执行的脚本）、`.command.err`（报错信息）、`.command.out`
- `nextflow log` 查看历史运行；终端里 `[tag 名]` 对应模块里的 `tag "$sample"`
- R 脚本报错时，可以先在 RStudio 里用同一份矩阵逐行调 `seurat_analyse.R`，改好再 `-resume` 重跑

## 进阶路线

1. **想对比 Python 生态**：下游换成 scanpy 只需改一个模块——写个 `bin/analyse.py` + 换容器镜像，调度逻辑完全不动
2. **给模块写测试**：学习 nf-test
3. **复用社区模块**：nf-core/modules 有几百个写好的模块可直接 include
4. **对照生产级实现**：`nextflow run nf-core/scrnaseq -profile test,singularity` 跑一遍官方流程，对比学习

## 命令速查

```bash
nextflow run <流程> -profile <环境>   # 运行
-resume                               # 断点续跑
--参数 值                             # 覆盖 config 参数
nextflow log                          # 历史运行
nextflow clean -f                     # 清理 work 缓存
nextflow pull nf-core/scrnaseq        # 拉取社区流程
```
