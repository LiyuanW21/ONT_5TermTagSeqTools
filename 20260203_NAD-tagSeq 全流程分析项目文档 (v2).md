# NAD-tagSeq 全流程分析项目文档 (v2)

## 第一部分：项目背景与实验设计

### 1.1 项目目标

本项目旨在利用 **NAD-tagSeq** 技术，结合牛津纳米孔（Oxford Nanopore）长读长测序，对特定样本中的 **NAD 加帽 RNA (NAD-capped RNA)** 进行全长转录组水平的鉴定和定量。通过对 NAD 标记（tagged）和非标记（non-tagged）的 RNA reads 进行对比分析，本项目致力于揭示 NAD-capping 在转录本水平的调控规律，包括其在不同基因和异构体上的分布、丰度比例，以及与特定序列基序（motif）的潜在关联。

### 1.2 实验技术原理

本项目采用的 NAD-tagSeq 方法，其核心流程基于对 NAD-capped RNA 的特异性化学标记。具体技术细节如下：

- **化学标记**: 通过 **DBCO-PEG4-Azide** 与 **末端炔基修饰的 RNA 标签** 之间发生的 **菌株促进的叠氮-炔环加成反应 (SPAAC)**，将一个 39bp 的已知序列 RNA 标签（NAD tag）共价连接到 RNA 5' 端的 NAD 帽子上。
- **测序平台**: 使用 **Oxford Nanopore PromethION 48** 平台进行直接 cDNA 测序。
- **建库策略**: 采用 **SQK-PCS109 (cDNA-PCR Sequencing Kit)** 结合 **SQK-PBK004 (PCR Barcoding Kit)** 进行文库构建，实现了全长转录本的捕获和多样本混合测序。

由此产生的测序 reads 具有明确的结构特征，即在 5' 端携带 NAD tag 的序列为 NAD-capped RNA，反之则为 non-capped RNA。其典型的 read 结构如下：

```
5' --- [SSP Adapter] --- [NAD Tag] --- [cDNA/转录本序列] --- [PolyA] --- [VNP Adapter反向互补] --- 3'
```

### 1.3 样本与原始数据
    ###样本为致病疫霉菌株1306三个菌丝阶段重复以及一个侵染马铃薯desiree 2天的样本。主要分析1306菌丝样本，顺带看看侵染样本。

| 样本名称         | 原始数据文件路径 (FASTQ)                |
| ---------------- | --------------------------------------- |
| `1306_1_pass`    | `/home/liyuan/NAD/20230522DATA/1306_1_pass.fq` |
| `1306_2_pass`    | `/home/liyuan/NAD/20230522DATA/1306_2_pass.fq` |
| `1306_3_pass`    | `/home/liyuan/NAD/20230522DATA/1306_3_pass.fq` |
| `1306_2dpi_1_pass` | `/home/liyuan/NAD/20230522DATA/1306_2dpi_1_pass.fq` |

---
参考基因组：/home/liyuan/001_analysis/001_NAD_21Project/Reference/1306reference/Pinf1306_UCR_1_genomic.fna
基因gff/gtf注释：/home/liyuan/001_analysis/001_NAD_21Project/Reference/1306reference/1306.gff  /home/liyuan/001_analysis/001_NAD_21Project/Reference/1306reference/1306.gtf
mRNA bed文件：/home/liyuan/001_analysis/001_NAD_21Project/Reference/1306reference/1306_mRNA_gene.bed

## 第二部分：已完成的数据预处理流程

在进入核心生物信息学分析之前，原始测序数据已经过一系列标准化的预处理步骤。整个流程由一系列自动化脚本完成，关键步骤、参数和文件路径总结如下。

### 2.1 步骤一：全长 cDNA reads 筛选 (`001_run_pychopper_pipeline.sh`)

- **目的**: 从原始 FASTQ 文件中筛选出包含完整 5' 和 3' 端接头的全长 cDNA reads。
- **核心工具**: `pychopper`
- **策略**: 采用两阶段筛选策略（`edlib` 快速筛选 + `phmm` 精确复筛）。
- **输入目录**: `/home/liyuan/NAD/20230522DATA/`
- **输出目录**: `/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/doublefilter_pychopper/`
- **关键输出**: 在输出目录中，生成 `${SAMPLE}_final_full_length.fq` 文件。

### 2.2 步骤二：NAD-tag 识别参数优化 (`003_cutadapt_grid_tuning.py`)

- **目的**: 为 `cutadapt` 工具寻找最优的参数组合，以最准确地区分带有 NAD tag 和不带 NAD tag 的 reads。
- **策略**: 通过网格搜索，测试了不同的错误率 (`-e`) 和最小重叠长度 (`-O`) 组合。
- **最终参数**: 确定 **`-e 0.2`** 和 **`-O 12`** 为后续分析的最佳参数。

### 2.3 步骤三：NAD-tag 识别与分离 (`002_run_cutadapt_tag_pipeline.sh`)

- **目的**: 应用优化后的参数，对全长 reads 进行切割，分离出 NAD-tagged 和 non-tagged reads。
- **核心工具**: `cutadapt`
- **输入目录**: `/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/doublefilter_pychopper/`
- **输出目录**: `/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/cutadapt_tag_e_0.2_O_12_out/`
- **关键输出**: 在输出目录中，生成 `${sample}_tagged.fq` 和 `${sample}_nontagged.fq` 文件。

### 2.4 步骤四：全基因组比对 (`004_run_minimap2_genome_mapping.sh`)

- **目的**: 将分离后的 tagged 和 nontagged reads 比对到参考基因组上。
- **核心工具**: `minimap2`
- **策略**: 采用 `map-ont` 预设参数进行比对。
- **参考基因组**: `/home/liyuan/001_analysis/001_NAD_21Project/Reference/1306reference/Pinf1306_UCR_1_genomic.fna`
- **输入目录**: `/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/cutadapt_tag_e_0.2_O_12_out/`
- **输出目录**: `/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/cutadapt_tag_e_0.2_O_12_out/mapping_1306_genome/`
- **关键输出**: 在输出目录中，生成 `${sample}.sorted.bam` 和 `${sample}.sorted.bam.bai` 文件。

---

## 第三部分：关键文件与目录路径汇总

为了方便 Agent 理解和执行，此处集中列出所有关键路径。

| 类别 | 路径 |
| :--- | :--- |
| **参考基因组 (FNA)** | `/home/liyuan/001_analysis/001_NAD_21Project/Reference/1306reference/Pinf1306_UCR_1_genomic.fna` |
| **基因组注释 (GFF/GTF)** | `/home/liyuan/001_analysis/001_NAD_21Project/Reference/1306reference/1306.gff`|
| **原始数据 (FASTQ)** | `/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/data/seqs/20230523_1306_3reps_2dpi_1rep` |
| **全长 Reads (Pychopper)** | `/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/doublefilter_pychopper/` |
| **Tagged/Nontagged Reads (Cutadapt)** | `/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/cutadapt_tag_e_0.2_O_12_out/` |
| **比对结果 (Minimap2 BAM)** | `/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/cutadapt_tag_e_0.2_O_12_out/mapping_1306_genome/` |

---

## 第四部分：后续深度分析计划

基于上述预处理得到的高质量比对文件，我们计划进行以下五步深度分析。

### 4.1 步骤一：比对后质量评估

- **输入**: `${sample}.sorted.bam` 文件，位于 Minimap2 输出目录。
- **操作**: 
    1. 使用 `Qualimap` 对 BAM 文件进行质控。
        ##安装：
        conda install -c conda-forge -c bioconda qualimap -y --solver libmamba

### 4.2 步骤二：基于对称原则的“完整 Read”筛选
- **输入**: `${sample}.sorted.bam` 文件和基因组注释文件。
- **操作**: 
    1. 定义“完整 Read”标准保证reads 5' 端序列大概率完整，确保tagReads和notagreads之间的对称性。
    2. 编写 Python 脚本，对 tagged 和 non-tagged 两组 BAM 文件执行4.2.1筛选逻辑。
    3. 最终参数记录: MAPQ >= 0, CDS overlap >= 50 bp, NAD-TSS window = 50 bp, 仅保留 primary 比对 (排除 secondary/supplementary)。
- **输出**: 
    1. 筛选后的 reads 子集（BAM 格式）以及reads筛选后的统计表。

    #### 4.2.1 筛选原则
    1. 比对可信度 (Mapping Confidence): 排除 secondary/supplementary 比对记录；保留所有 primary 比对 (MAPQ >= 0)。
    3. CDS 起始位点定义 (CDS Start Definition): 以基因组注释的 CDS 起始密码子 (Start Codon) 坐标作为参考锚点（若无 CDS 注释，则该基因标记为不适用）。注意 CDS 为 mRNA 的编码区，不包含 UTR。
    4. 5' 端完整性校验 (CDS Overlap): 所有保留的 reads，其 5' 端必须覆盖 CDS 起始位点；最小重叠阈值为向 CDS 内部延伸 >= 50 bp。
    5. 对称原则过滤 (Symmetric 5' Filtering): 为贯彻 TagSeqTools 的对称原则，消除截断 reads 对定量带来的偏差，采用以下分级策略：
        - 首选方案（数据驱动 NAD-TSS）: 利用 tagged reads 在基因组上的堆叠峰值定义每个基因的经验性 NAD-TSS。对于 untagged reads（分母）:
            - 保留: 起始位置位于 NAD-TSS 上游 (upstream) 的 reads（允许更长的 5' UTR），或位于 NAD-TSS 下游 50 bp 窗口内的 reads。
            - 丢弃: 起始位置位于 NAD-TSS 下游 >50 bp 的 reads（视为 5' 端严重截断）。
        - 备选方案（无 tag 时的 CDS 补救）: 仅当某基因未检测到 tagged reads 时，参考注释的 CDS 起始位点。要求 untagged reads 必须完整覆盖 CDS 5' 端，否则予以丢弃。

    **reads 筛选结果汇总（基于整合统计表）**
    - 来源文件：`/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/003_cutadapt_tag_e_0.2_O_12_out/mapping_1306_genome/001_complete_reads_no_strand_mapq0_cds50/mapping_rate_with_complete_reads.tsv`

    | sample | Total_reads | group | total_reads | primary_mapped_reads | primary_mapping_rate(%) | primary_map | cds_overlap | nad_tss | kept |
    | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
    | 1306_1 | 4267408 | tagged | 55334 | 51256 | 92.63 | 51256 | 27031 | 0 | 27031 |
    | 1306_1 | 4267408 | nontagged | 4212074 | 3810009 | 90.45 | 3810009 | 129281 | 108333 | 108333 |
    | 1306_2 | 4101913 | tagged | 60848 | 55579 | 91.34 | 55579 | 27146 | 0 | 27146 |
    | 1306_2 | 4101913 | nontagged | 4041065 | 3596918 | 89.01 | 3596918 | 130586 | 104834 | 104834 |
    | 1306_3 | 3594173 | tagged | 128344 | 122383 | 95.36 | 122383 | 74756 | 0 | 74756 |
    | 1306_3 | 3594173 | nontagged | 3465829 | 3321433 | 95.83 | 3321433 | 228941 | 183723 | 183723 |
    | 1306_2dpi_1 | 17710391 | tagged | 709046 | 7021 | 0.99 | 7021 | 2736 | 0 | 2736 |
    | 1306_2dpi_1 | 17710391 | nontagged | 17001345 | 216404 | 1.27 | 216404 | 10851 | 9953 | 9953 |


    #### 4.2.2 CDS 覆盖度分析
    -**输入**: `${sample}.sorted.bam` 文件和基因组注释文件。
    - **操作**: 
        1. 结合基因组注释文件 (GFF/GTF)，使用 deepTools 绘制“元基因（metagene）”覆盖度曲线。
            脚本: `005_metagene_deeptools.sh`（先用 `bamCoverage` 生成 BPM 归一化 BigWig，再用 `computeMatrix`/`plotProfile` 画图）。
        2. 合并同一样本重复 BigWig 并取平均: `005_merge_bigwig.sh`（输出到 `bamcoverage_bpm_merged/`）。
        3. 对合并后的 BigWig 计算矩阵（scale-regions）并输出：
            - raw 合并样本矩阵: `bamcoverage_bpm_merged/raw/` -> `metagene_matrix_merged/raw/`
            - complete 合并样本矩阵: `bamcoverage_bpm_merged/complete/` -> `metagene_matrix_merged/complete/`
        4. 使用单独 R 脚本作图（便于调整样式）: `006_plot_metagene_profile.R`
            - 输入: `metagene_matrix_merged/raw/metagene_profile.tsv` 与 `metagene_matrix_merged/complete/metagene_profile.tsv`
            - 输出: `metagene_matrix_merged/plots/metagene_profile_merged.png` 与 `metagene_matrix_merged/plots/metagene_profile_merged.eps`
            -要求：左右两个面板，共用图例，画1306-1，1306-2,1306-3，Tagged和Nontagged两组分别用相似色系。

### 4.3 步骤三：基因counts统计与 NAD Ratio 分析

- **输入**: 筛选后的“完整 Read”子集。
- **操作**: 
    1. 计量每个基因有多少 NAD read 和多少 non-NAD reads。**NAD Ratio = (Tagged Reads) / (Total Reads)**
        - 方法：使用 `bedtools coverage -counts` 对 complete reads BAM 进行基因级计数（no strand）。
        - 计数注释：/home/liyuan/001_analysis/001_NAD_21Project/Reference/1306reference/1306_mRNA_gene.bed
        - 脚本：007_gene_counts_nad_ratio.py
        - Criteria of NAD capping genes：1. NAD reads ≥2，Enriched in at least two copies.
        - 说明：2‰（tagged/total > 0.002）在当前数据上无额外筛选作用，已移除。
        - 输出目录：/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/007_gene_counts_nad_ratio/
        - 汇总表（横向对比）：gene_counts_all_samples.tsv（同一基因一行，各样本横向展开）
        - 样本级 NAD capping 列表：nad_capping_by_sample/{sample}.nad_capping.tsv
        -1.1 新增统计 nonNAD-capped genes:
         -在1306-1, 2, 3三个样本中，至少有2个replicates的tag+nontag reads总和 ≥ 2 的基因共有 8265 个。其中3091个NAD-tagged gene,所以Non-tagged gene为5174个。
        
    
    2. 计算3个重复的NAD ratio平均值，列出 top10。
        - 排序规则：先按 mean NAD ratio 降序，再按 mean NAD reads（tagged）降序。
        - 输出：nad_capping_genes_replicates.tsv（包含每个样本与均值）与 nad_ratio_top10_replicates.tsv

    3. 写 R 脚本绘制筛选 NAD capping genes 的维恩图。
        - 脚本：007_plot_nad_capping_venn.R
        - 输出：004_gene_counts_nad_ratio/plots/nad_capping_genes_venn.svg
        - 说明：SVG 可在 Adobe Illustrator 里直接编辑。
    
    4. 对3091个NAD capping genes进行GO富集分析。
        -新建文件夹005_GO_anno,放1306基因的GO注释和蛋白注释。
            #GO注释主要参考20240125用blast2GO注释的信息。
                #有两个版本：1.原始noslim的版本: 1306_noslim_OOmycota. 2.Pathoslim: 1306_Pathoslim_OOmycota.
            #蛋白注释参考20260207_EGGNOG完成的注释: proteins.eggnog.emapper.annotations.1306.txt
        -PInf<->KAI ID对应关系来自GFF，文件：005_GO_anno/pinf_to_kai_from_gff.tsv（统一使用PInf名称）。
        -对3091个nad_capping_genes进行GO富集分析，两个GO版本都做一遍。
        -脚本：/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/005_GO_anno/008_run_go_enrichment_topgo.R
        -输出目录：/home/liyuan/001_analysis/001_NAD_21Project/001_NAD_seq_analysis/analysis/001_20230522data1306/005_GO_anno/go_enrichment/
            - noslim: go_enrichment_noslim_{BP,MF,CC}.tsv + summary
            - pathoslim: go_enrichment_pathoslim_{BP,MF,CC}.tsv + summary
            - genericslim: go_enrichment_genericslim_{BP,MF,CC}.tsv + summary
            - Top20图：go_enrichment/plots/*_top20.png
        -背景基因数（GO注释覆盖后的PInf基因数）:
            - noslim: 12973
            - pathoslim: 19981 ❌️
            - genericslim: 13919 ✅️
        -在3091个NAD capping genes中的GO覆盖数:
            - noslim: 2494
            - pathoslim: 3091 ❌️
            - genericslim: 2522 ✅️
        -补充：3个replicates标准（1640基因）已做noslim与genericslim GO分析，结果未见比2个replicates标准更集中于致病相关基因。
            
    5. Density plot 或 violin plot 绘制NAD ratio分布
    6. 分析NADcappedgene和nonNADcappedgene之间reads长度是否一样？
        

### 4.4 步骤四：高级分析——Motif 探索与异构体分析
    
    ####4.4.1 Motif鉴定

    - **输入**: 定量结果、BAM 文件、基因组序列。

    - **操作**: 
        1. 精确鉴定 NAD capping gene RNA 的 TSS/帽子下游100bp，并使用 **MEME Suite** 进行 motif 富集分析。
        

    #####4.4.2 异构体水平的定量
        1. 进行差异异构体使用 (DIU) 分析。

### 4.5 步骤五：结果可视化与报告整合

- **操作**: 将所有分析结果进行系统性整合，生成火山图、Sashimi plot、Motif logo 等图表，并更新至本报告。

## 第五部分：参考文献

[1] Huang, Z., Cai, Z., Yang, Z., & Xia, Y. (2020). TagSeqTools: a flexible and comprehensive analysis pipeline for NAD tagSeq data. *bioRxiv*. [https://doi.org/10.1101/2020.03.09.982934](https://doi.org/10.1101/2020.03.09.982934)
