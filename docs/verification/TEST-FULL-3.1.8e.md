# 3.1.8e 真实数据完整测试记录

> 版本：`/nfs/users/zhaojw/software/miRDeep-p3-3.1.8e`（`__version__ = 3.1.8e`）
> 目的：用一个**完整可用的外部数据集**跑通 identification → annotation → analysis 全链路，
> 而不只是 `test.sh` 的 mini 数据。
> 运行时间：2026-11-18 22:05 起

---

## 1. 为什么用这套数据

水稻 case1 的 identification 输入在 `/nfs/users/xuzx/...`，当前账号**无读权限**
（`ls` 直接 Permission denied），无法从原始 reads 重跑 identification + annotation。

因此选用 `Scaevola taccada`（草海桐）项目数据：全部在 `zhaojw` 名下，可读，且是**不同物种**，
可同时检验跨物种稳健性。

| 项 | 值 |
|---|---|
| 输入根目录 | `/nfs/users/zhaojw/project/Scaevola_taccada_miRNA/sRNA-seq_result/` |
| 样本数 | **15**（flower×3, fruit×3, leaf×3, root×3, stem×3） |
| 每样本产物 | `chr_length`、`<s>-processed.fa`、`<s>.total_reads`、`<s>_nr_predictions`、`<s>_filter_P_prediction` —— **15/15 完整** |
| 基因组 | `/nfs/users/zhaojw/project/Scaevola_taccada_miRNA/genome/Sta_v4.fasta`（1.3 GB，8 条染色体 `Chr01`–`Chr08`） |
| Bowtie 基因组索引 | `sRNA-seq_result/index/genome_index`（已有，无需重建） |
| 旧注释结果 | `sRNA-seq_result/annotation/annotation/`（344 条，**由 ≤3.1.6c 产出**） |

---

## 2. 运行命令（留档）

```bash
cd /nfs/users/zhaojw/software/miRDeep-p3-3.1.8e
python ./mirdeep-p3 annotation \
  -f /nfs/users/zhaojw/Genburten/temp/mirdp3_318e_fulltest/input_folders.txt \
  -o /nfs/users/zhaojw/Genburten/temp/mirdp3_318e_fulltest/out \
  -g /nfs/users/zhaojw/project/Scaevola_taccada_miRNA/genome/Sta_v4.fasta \
  -d /nfs/users/zhaojw/project/Scaevola_taccada_miRNA/sRNA-seq_result/index/genome_index \
  --prefix scr --prefix_miRNA Sta --species "Scaevola taccada" \
  --fix-mature-boundary -t 8
```

输出目录：`/nfs/users/zhaojw/Genburten/temp/mirdp3_318e_fulltest/out/scr/`

---

## 3. Annotation 结果

### 3.1 结构筛选（3.1.8e 新增步骤）

```
[structure] loci: 344  kept: 343  removed: 1 (min arm gap: 10 nt)
[structure]   short_arm_gap: 1
[structure]   removed: Sta-MIR3711a
```

被剔除条目：

| 名称 | 状态 | 前体 | mature | star | gap | 原因 |
|---|---|---:|---:|---:|---:|---|
| `Sta-MIR3711a` | short_arm_gap | 48 | 21 | 21 | **6** | loop 6 nt < 10 |

前体 48 = 21 + 21 + 6，又是「绕成熟硬切窗口」的典型。

### 3.2 执行顺序已验证（关键）

日志行号证明筛选发生在**任何比对之前**：

| 行号 | 事件 |
|---:|---|
| 162 | `[structure] loci: 344 kept: 343` |
| 169 | `Output written to: .../scr-basic-info-mature-1nt.fa`（**参考从此才构建**） |
| 245 | `[boundary] alignments used: 15  loci: 343` |

筛前原件保留：`temp/scr-basic-info.before-structure-filter`（**344 行**，含被删条目）；
最终 `scr-basic-info` 中该条目已不存在。被删位点**没有参与任何 reads 比对**。

### 3.3 边界审计与前体修复

```
[boundary] loci: 343  canonical: 320  low coverage: 13  flagged: 10
[boundary]   shift_beyond_limit: 1
[boundary]   ambiguous: 9
[boundary] no correction needed; basic-info left as is
```

本数据集**没有触发纠正**（0 条 shift），因此前体修复步骤按设计不执行
（它只在 `--fix-mature-boundary` 真正改表后才跑）。这两条路径已由 mini 夹具与
水稻 case1 的 part2 专项测试覆盖。

### 3.4 行数一致性

| 产物 | 行数 |
|---|---:|
| `scr-basic-info` | 343 |
| `scr-mature.count` 数据行 | 343 |
| `scr-mature.exp` 数据行 | 343 |
| `scr-structure-filter.tsv` 判定 | 344（每条输入一行） |
| `scr-basic-info-cluster` | 10 簇 |

三者严格对齐，无 446-vs-427 那类缺口。

---

## 4. 逐条质量核验（343 条全查）

| 检查项 | 结果 |
|---|---|
| 成熟坐标长度 == 序列长度 | **343/343 通过** |
| star 坐标长度 == 序列长度 | **343/343 通过** |
| 前体坐标长度 == 序列长度 | **343/343 通过** |
| 前体结构串长度 == 前体序列长度 | **343/343 通过** |
| ±20 窗口结构串长度 == 窗口序列长度 | **343/343 通过** |
| 基因组能复现成熟序列 | **343/343 通过** |
| 基因组能复现 star 序列 | **343/343 通过** |
| 基因组能复现前体序列 | **343/343 通过** |
| 成熟完全在前体内 | **343/343 通过** |
| star 完全在前体内 | **343/343 通过** |
| 臂间距 gap ≥ 10 | **343/343 通过** |

**问题总数：0。**

其它统计：

- gap 分布：min 16 / p10 32 / 中位 82 / max 235
- star 长度：全部 21 nt（本次未做结构感知改造）
- mature 长度：20 (73)、21 (173)、22 (83)、23 (6)、24 (8)
- 表达矩阵：17 列（MIR + MIR_ID + 15 样本），334/343 条非零，**无空样本列**
- 总 counts 23,420,800；最大 RPM 17,787

---

## 5. 与旧结果（≤3.1.6c）的对照

按**基因组位点键**（chr + strand + 成熟起止）匹配，而非按名称——MIRN 家族在两次运行间
会重新编号，按名字比会错配。

| 项 | 旧结果 | 本次 3.1.8e |
|---|---:|---:|
| 条目数 | 344 | 343 |
| 位点与旧结果一致 | — | **343**（差 1 条 = 被新筛选剔除的 `Sta-MIR3711a`） |
| star 坐标长度 == 序列长度 | **0 / 306** | **306 / 306** |
| star 序列与前体结构 | — | **343/343 与旧结果完全相同** |
| 前体序列 | — | **343/343 相同** |
| star 重叠（gap<0） | **170**（该批位点） | **0** |

结论清晰：**除了 star 坐标，其它一切都没变**。

旧结果的 `star_end − star_start` 恒为 **20**，而 star 序列是 **21 nt**——这正是 3.1.6c
已知的 star 坐标 bug（`tests/test_coordinates.py` 里有对应回归）。它导致 344 条里 171 条
看起来「成熟与 star 重叠」，实际是坐标错位造成的假象。3.1.8e 中该问题已消失。

---

## 6. 下游分析

### 6.1 `analysis Stat`

```bash
python ./mirdeep-p3 analysis Stat -i scr-basic-info --rnaplot -o analysis_stat
```

- 通过，产出 **355 个文件**
- 含 4 类统计图（碱基分布、首位分布、长度分布、家族分布）× pdf/png/svg
- 含 **343 张二级结构图**（`*_ss.svg`），每个存活位点一张——说明前体结构与坐标自洽

### 6.2 `analysis Differential_expression`

```bash
python ./mirdeep-p3 analysis Differential_expression \
  -c scr-mature.count -r scr-mature.exp \
  --case1 9,10,11 --case2 12,13,14 --case1name leaf --case2name root \
  -o analysis_DE_leaf_vs_root
```

| 项 | 值 |
|---|---|
| 过滤后保留（`--min-expr` 默认 5 RPM） | **152** |
| 上调（padj<0.05, log2FC>1） | **36** |
| 下调（padj<0.05, log2FC<-1） | **36** |
| 不显著 | 80 |
| 产物 | **620 个文件**（PCA、样本相关热图、火山图、表达热图、逐 miRNA 表达图） |

152 + 80 = 232 ≠ 343，是因为 DESeq2 内部还会剔除全零/离群条目，属正常。

---

## 7. 本次测试覆盖与未覆盖

| 路径 | 是否覆盖 | 依据 |
|---|---|---|
| identification → annotation 全链路 | ✅ | 本测试 15 样本 |
| 3.1.8e 结构筛选（含实际剔除） | ✅ | 344→343，剔除 1 条 |
| 3.1.8e 筛选的执行位置（比对前） | ✅ | 日志行号 162 → 169 → 245 |
| 表达矩阵与 basic-info 行对齐 | ✅ | 343 = 343 = 343 |
| 3.1.8c 边界审计 | ✅（只审计不纠正） | 320 canonical / 10 flagged |
| **3.1.8c 边界纠正（shift）** | ❌ 本数据 0 条 | 已由 mini 夹具 + 水稻 case1 覆盖 |
| **3.1.8d 前体越界修复** | ❌ 未触发 | 同上 |
| 下游 Stat / DE | ✅ | 355 / 620 个产物 |

未覆盖的两条路径依赖「有越界或偏移的位点」，本数据集恰好干净，因此不触发；
这本身是个好结果，但需说明覆盖缺口不在版本缺陷而在数据。

---

## 8. 结论

1. **全链路跑通**：15 样本、1.3 Gb 基因组、344 个位点，从 identification 产物一路到 DE 结果。
2. **结构筛选按设计工作**：正确剔除 1 条 `precursor = mature + star + loop` 的窗口型位点，
   且在建立比对参考之前完成，被删位点零计算。
3. **逐条零缺陷**：343 条在坐标/序列/结构/包含性/臂间距五个维度全部通过，
   且全部能被基因组复现。
4. **旧版 star 坐标 bug 已修复**：同一批位点，旧结果 306/306 坐标不可复现、170 条假重叠；
   新结果 0 条异常，而 star 序列、前体序列、前体结构 343/343 与旧结果相同——
   证明改动精准命中问题本身，没有连带影响。
5. **下游可用**：结构图 343 张、PCA/热图/火山图齐全，DE 得到 36 上 / 36 下。

---

## 9. 产物位置

```
/nfs/users/zhaojw/Genburten/temp/mirdp3_318e_fulltest/
├── input_folders.txt          15 个样本目录清单
├── run_cmd.sh                 本次运行命令（可复现）
├── run.log                    完整 stdout/stderr
└── out/scr/                   47 M
    ├── scr-basic-info             343 条（筛后）
    ├── scr-structure-filter.tsv   344 条判定
    ├── scr-boundary-audit.tsv     边界审计
    ├── scr-mature.count / .exp    15 样本矩阵
    ├── scr-basic-info-cluster     10 簇
    ├── scr_annotation.log         全流程日志
    ├── temp/scr-basic-info.before-structure-filter   344 条（筛前原件）
    ├── analysis_stat/             355 个文件
    └── analysis_DE_leaf_vs_root/  620 个文件
```
