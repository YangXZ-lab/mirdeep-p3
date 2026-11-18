# 3.1.8c 验证报告 — 成熟序列 5' 边界审计与校正

> 基线：`/nfs/users/zhaojw/software/miRDeep-p3-3.1.7a`（实际版本 3.1.8b）
> 工作目录：`/nfs/users/zhaojw/software/miRDeep-p3-3.1.8c`
> 日期：2026-11-18
> 状态：本地完成，**未提交 GitHub**

---

## 1. 要解决的问题

3.1.8b（及之前版本）的成熟序列区间来自预测，再经 `utils/coords.py` 按臂做 ±1 调整。
该调整可能让区间相对支撑它的 reads 在 5' 端错位 1 个碱基。**错位时不报错**：
表达矩阵把它当作 1 nt offset 吸收掉，表和矩阵都"看起来正常"。

四物种实测（Arabidopsis 135 / Gnetum 513 / rice 446 / Scaevola 344，共 1420 条可定量）：

| 物种 | 可定量 | 5' 主导偏移非 canonical | 占比 |
|---|---:|---:|---:|
| Arabidopsis thaliana | 135 | 1 | 0.7% |
| Gnetum montanum | 513 | 11 | 2.1% |
| Oryza sativa | 428 | 22 | 5.1% |
| Scaevola taccada | 344 | 13 | 3.8% |
| **合计** | **1420** | **47** | **3.3%** |

其中与 miRBase 交叉验证确认为**注释错误**的有 3 条：`Nip-MIR166a`、`Nip-MIR166b`、
`Gmo-MIR159a`，全部为注释 5' 端多 1 nt。`Nip-MIR166a/b` 是该家族最高丰度成员
（各约 160 万 counts），错误并非数值上可忽略。

---

## 2. 两个关键发现（决定了实现方式）

### 2.1 basic-info 的区间是「1-based 起点 + 排他终点」

`Nip-MIR166a` 成熟为 `3487876..3487897`，即 `3487897 − 3487876 = 21`，
**等于成熟序列长度**。这不是 1-based 闭区间（那样会是 22）。前体、region、star 同样如此。

实现初期按闭区间计算，导致整体偏 1，修正后的序列取到了前体 5' 端——
被 `Nip-MIR166a` 的断点直接暴露。现在的 `mature_rel_bounds()` / `rel_to_genomic()`
按此约定实现，并在读取前体坐标时用 `end − start == len(pre_seq)` 做一致性校验，
不匹配就停止而不是"猜"。

### 2.2 修法是「平移」而非「裁 5'」

只看 offset 无法决定怎么改：

| 记录 | 注释成熟 | 主导 reads | 正确动作 |
|---|---|---|---|
| `Nip-MIR166a` | `TTCGGACCAGGCTTCATTCCC` (21) | `TCGGACCAGGCTTCATTCCCC` (21) | **平移**（5' 去 1、3' 补 1）|
| `Gmo-MIR159a` | `GTTGGATTGAAGGGAGCTCCA` (21) | `TTGGATTGAAGGGAGCTCCA` (20) | **裁 5'** |

两者长度关系不同。因此实现同时构造 `trim`（只动 5'）与 `shift`（整体平移）两个候选，
**只接受能精确复现主导 read 序列的那个**；都不匹配则不动。

### 2.3 只裁不补会更糟（反例）

3.1.8c 第一版只做 5' 裁切，`Nip-MIR166a` 得到 `TCGGACCAGGCTTCATTCCC`（20 nt），
比修正前更偏离 miRBase。加入 `shift` 候选后得到 `TCGGACCAGGCTTCATTCCCC`，
与 miRBase `osa-miR166a-3p` 一致。

---

## 3. 改动清单

| 文件 | 类型 | 内容 |
|---|---|---|
| `src/utils/mature_boundary.py` | 新增 | 纯函数：区间换算、主导 offset 选择（tie 偏向 canonical）、候选区间生成与裁决、安全护栏 |
| `src/audit_mature_boundary.py` | 新增 | CLI：读 `.mature.aln` → 出审计表；`--apply` 写修正表 |
| `tests/test_mature_boundary.py` | 新增 | stdlib 单测，含 `Nip-MIR166a` 真实记录回归 |
| `src/commands/annotation.py` | 修改 | 新增 3 个参数；part2 改为 对齐 → 审计(+可选修正) → 聚类 → 定量 |
| `src/version.py` | 修改 | 3.1.8b → 3.1.8c |
| `CHANGELOG.md` | 修改 | 新增 `[Unreleased] - 3.1.8c` |
| `test.sh` | 修改 | 第 11 节接入新单测 |
| `README.md` | 修改 | `Current versions` 版本串 → v3.1.8c |
| `docs/DATA-INDEX.md` | 修改 | 软件版本行 → `3.1.8c` |

### 新增参数

| 参数 | 默认 | 作用 |
|---|---|---|
| `--fix-mature-boundary` | off | 打开才真正改写 basic-info |
| `--boundary-min-fraction` | 0.60 | 主导 offset 的最低占比 |
| `--boundary-min-weight` | 10.0 | 最低 read 权重 |

`audit_mature_boundary.py` 另有 `--max-shift`(1)、`--max-3p-gap`(2)、`--allow-extend`(off)。

### 保守边界（重要）

- **只接受 ±1 nt**。≥2 nt 归为 `shift_beyond_limit`，仅报告。理由：那个距离上
  主导 offset 更可能是 5' 降解而非边界错判。
- **不自动向上游延伸**。offset 0（reads 起点早于注释成熟）归为 `shift_upstream`，
  仅报告。1 nt flank 的存在就是为了吸收这类 reads；据此反推"注释少算 1 nt"
  是比"裁掉 1 个 reads 不支持的碱基"强得多的断言。需要时可 `--allow-extend`。

---

## 4. 验证

### 4.1 单元测试

```
python3 tests/test_mature_boundary.py
test_mature_boundary: all checks passed
```

覆盖：双链相对/基因组往返、offset→shift 换算、tie 决断（50/50 不纠偏）、
候选裁决（trim / shift / 都不匹配 / 越出前体）、臂归属、以及 `Nip-MIR166a` 真实回归。

### 4.2 端到端（`test.sh`，mini 数据集）

```
[test] Mature-boundary unit tests passed.
[test] Version block in README matches src/version.py (3.1.8c) and data-index.env (...)
[test]  ALL TESTS PASSED
```

**默认行为与 3.1.8b 逐字节一致**（md5 对比 3.1.8b 发布记录）：

| 产物 | 3.1.8b | 3.1.8c 默认 |
|---|---|---|
| `mirdp3-test-basic-info` | `f9223f3c2b120d62ec24fa60095aca09` | 相同 |
| `mirdp3-test-basic-info-cluster` | `db799e4e78ea901bcb75fb2d0f160dbd` | 相同 |

新增产物 `mirdp3-test-boundary-audit.tsv`，mini 数据 19 条全部 canonical。

### 4.3 真实数据（case1，21 样本，446 条）

```
[boundary] loci: 446  canonical: 405  low coverage: 18  flagged: 23
[boundary]   shift: 2
[boundary]   shift_upstream: 1
[boundary]   shift_beyond_limit: 12
[boundary]   ambiguous: 7
[boundary]   out_of_range: 1
[boundary] corrected 2 mature intervals
```

| miRNA | 旧 | 新 | 说明 |
|---|---|---|---|
| `Nip-MIR166a` | `3487876..3487897` `TTCGGACCAGGCTTCATTCCC` | `3487877..3487898` `TCGGACCAGGCTTCATTCCCC` | = miRBase `osa-miR166a-3p` |
| `Nip-MIR166b` | `30392798..30392819` `TTCGGACCAGGCTTCATTCCC` | `30392797..30392818` `TCGGACCAGGCTTCATTCCCC` | 同上（负链，坐标同步下移）|

**改动面核验**：446 行中仅 2 行变化，且全部只落在第 16/17/18 列
（成熟 start / end / seq）；其余 444 行逐字节不变；原始 `case1-basic-info`
的 md5 未变（`034b74ad7edb59dbc891274cb9641d2f`）。

**被拒绝的例子**（体现保守性）：`Nip-MIR11339a` offset=2、占比 75%，
但两个候选都无法复现主导 read，判为 `out_of_range` 不动——
它是真实的 5' 截短分子，不是边界错判。

### 4.4 修正分支端到端（构造输入）

mini 数据全部 canonical，默认路径不会触发修正。因此构造了一份"reads 整体后移 1 nt"
的 identification 输出（每条 read 去掉 5' 首碱基，头部与 `_xN` 计数保持），
再跑完整 annotation：

```
[boundary] loci: 19  canonical: 0  low coverage: 0  flagged: 19
[boundary]   shift: 19
[boundary] corrected 19 mature intervals
[boundary] loci: 19  canonical: 19  low coverage: 0  flagged: 0   ← 修后验证
```

- 19 条成熟序列全部按 reads 裁掉首碱基（`CCCCAAATGTAGACAAAGCA` → `CCCAAATGTAGACAAAGCA` 等）
- `temp/fixtest-basic-info.before-boundary-fix` 保留了修正前原文
- 修后审计 19/19 canonical，说明重建参考 → 重新比对 → 复核这条链闭合

---

## 5. 使用方式

```bash
export PATH=/nfs/users/zhaojw/miniconda3/envs/mirdp3_conda/bin:$PATH
cd /nfs/users/zhaojw/software/miRDeep-p3-3.1.8c

# 默认：只审计，不改数据（输出 <prefix>-boundary-audit.tsv）
./mirdeep-p3 annotation -i <identification输出> -g <genome> -d <genome_index> \
    --prefix_miRNA <prefix> --prefix <outname> -t 16

# 启用修正
./mirdeep-p3 annotation ... --fix-mature-boundary

# 单独复核已有结果
python src/audit_mature_boundary.py -i <prefix>-basic-info \
    -aln-dir <prefix>/temp -o audit.tsv [--apply -oc fixed-basic-info]
```

---

## 6. 过程记录（供后续参考）

1. **区间约定踩坑**：先按 1-based 闭区间实现，`Nip-MIR166a` 修正后序列取到前体 5' 端。
   用真实记录做回归用例后暴露并修正。
2. **`test.sh` 补丁吃掉 `fi`**：首次锚点选在 `log "...O-10"` 行，导致原 `if` 的 `fi`
   被新块占用，`bash -n` 报 `unexpected end of file`。改为锚定完整 `if…fi` 块。
3. **README 版本守卫生效**：首次端到端在第 11 节以 exit 7 失败
   （`README 'Current versions' says v3.1.8b but src/version.py says 3.1.8c`）。
   这是 3.1.8b 新增的文档守卫，说明它不是装饰。同步 README 与 `docs/DATA-INDEX.md` 后通过。
4. **不做的事**：未改 `utils/coords.py` 的臂相关 ±1 逻辑。那是全量坐标调整，
   影响所有物种的 basic-info / star / 靶标 / 聚类；本次是在注释之后做**有证据的局部校正**，
   改动面可控、可回滚。

---

## 7. 未做 / 后续

- 未提交 Git、未打 tag、未改动 conda/docker/bioconda（按赵哥要求）。
- `data/index` 仍是指向旧索引的软链接（`software/github/miRDeep-P3/data/index`），
  打包前需处理。
- 建议后续把本模块的判定逻辑并入 `utils/coords.py` 的长期方案讨论：
  根治点在于臂调整时就用 reads 定位，而不是事后校正。

---

## 8. 改动范围核查（2026-11-18，回答"是否会波及其他 miRNA"）

结论：**只有被修正的那几条会变，其余一律不动**。下面每一项都是实测，不是推断。

### 8.1 basic-info 逐格比对（case1，446 行 × 22 列 = 9,812 个单元格）

| 项目 | 结果 |
|---|---|
| 行数 | 446 → 446，行序不变 |
| 变化单元格 | **6 个（0.061%）** |
| 变化行 | **2 行**：`Nip-MIR166a`、`Nip-MIR166b` |
| 变化列 | **仅第 16/17/18 列**（成熟 start / end / seq） |
| 身份列（名称、MIR_ID） | 未改动 |
| 其余 444 行 | **逐格完全相同** |

### 8.2 star 列（第 19–22 列）不需要同步改

star 是前体另一条臂，判定只依赖成熟在前体中的相对位置，容差 ±2 nt。成熟 5' 移动 1 nt
不改变臂归属；实测修正后 star 与成熟区间**无重叠**（MIR166a 相距 57 nt，MIR166b 22 nt），
star 的坐标与序列保持原值是正确的。

### 8.3 下游文件

| 下游产物 | 结果 |
|---|---|
| `*-basic-info-cluster` | 两次运行输出 **0 行差异**（17 个簇完全相同） |
| `*-basic-info-mature-1nt.fa`（定量参考） | 446 条，顺序一致，**仅 2 条内容不同**，长度不变（23 → 23） |
| 表达矩阵（对两个参考分别建索引、比对、定量） | 446 行中 **仅 2 行数值变化**，其余 444 行**逐值相同** |

定量层之所以必然隔离：`bowtie -a -v 0` 对每条 read 独立报告所有参考命中，某条参考
序列改变不会改变其他参考的命中集合，`expression_matrix.py` 又按 (name, MIR_ID) 分别累加。

### 8.4 本次核查中发现并修掉的两个问题

**① 审计脚本把"±20 nt 窗口"当成了前体（22 列布局）**

22 列表里第 7–9 列是前体 ±20 nt 窗口（`Precursor_20nt_*`），第 11–13 列才是真实前体。
原实现（以及 24 列布局）读的是真实前体，22 列却读成窗口，导致：

- `arm` 标签用了一个比发夹大 20 nt 的尺子量，结论偏乐观；
- `resolve_candidate()` 本就要求候选落在前体内，但因为框选错了，越界检查形同虚设。

已改为**双窗口**：候选解析用 ±20 窗口（一定有腾挪空间），臂边界与包含性校验用真实前体。
两种列布局（22 / 24）现在对同一批 loci 给出**完全相同的判定与坐标**（446 条逐一比对，0 差异）。

**② 修正后成熟会超出注释前体 1 nt（真实存在）**

`Nip-MIR166a/b` 的前体当初是贴着（错误的）成熟 3' 端切出来的，成熟后移 1 nt 后末端就
落到前体之外：

```
Nip-MIR166a  mature 3487877..3487898   pre 3487799..3487897
Nip-MIR166b  mature 30392797..30392818  pre 30392798..30392861
```

原始表里没有任何一条这个现象（0/446），修正后出现 2 条。判定为**前体本身短 1 nt**——
突出那一位在 ±20 窗口里能取到，属基因组编码、非非模板加成。

处理方式：**照常修正**（reads 是证据），但在运行日志里显式报出，不静默：

```
[boundary]   corrected mature overhangs the annotated precursor: 2 (Nip-MIR166a, Nip-MIR166b);
             the precursor end is short by the same amount
```

不直接改前体的原因：改前体边界就得同步重算 `Precursor_struc`（RNAfold 结构串），
否则序列与结构长度对不上，会把一个 1 nt 的坐标问题换成一个结构列损坏问题。
下游实测都不受影响：`RNA_plot.py` 画图用的是 ±20 窗口（第 9/10 列），
`target_finder` 用成熟序列本身，定量也用窗口。若将来要彻底闭合，正确做法是
"扩前体 + 重折叠 + 重算 star"，属独立改动。

### 8.5 复跑与一致性

- 用修好的审计脚本重跑 case1，修正结果与修前**逐字节一致**（`cmp` 通过）。
- 单元测试 `tests/test_mature_boundary.py` 全通过。
- `test.sh` 端到端重跑通过（见第 4 节）。
- 顺带把新增文件的 CRLF 换行统一为 LF，与仓库其余文件一致。
