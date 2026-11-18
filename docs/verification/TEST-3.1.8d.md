# 3.1.8d 验证报告 — 前体越界修复与重新折叠

> 基线：3.1.8c（`/nfs/users/zhaojw/software/miRDeep-p3-3.1.8c`）
> 交付：`/nfs/users/zhaojw/software/miRDeep-p3-3.1.8d`
> **未提交 GitHub**，未改 conda / docker / bioconda。

---

## 1. 要解决的问题

3.1.8c 能把成熟 5' 边界按 reads 证据纠正过来，但纠正会带来一个副作用：
**成熟可能落到注释前体之外 1 nt**。

原因是前体当初是贴着（错误的）成熟 3' 端切出来的。成熟往后挪 1 nt，末端就越界了。
原始数据 446 条里 0 条越界；纠正后有 2 条：

```
Nip-MIR166a  mature 3487877..3487898   pre 3487799..3487897   → 末端差 1
Nip-MIR166b  mature 30392797..30392818  pre 30392798..30392861 → 起点差 1
```

后果不是崩溃，而是**内部不一致**：成熟不再是前体序列的子串，存下来的点括号结构
描述的是一个装不住这条成熟的窗口，star 也是相对一个偏 1 的坐标系放的。

---

## 2. 做法

新增一步，只在 `--fix-mature-boundary` 真正改写了表之后运行：

1. **扩前体**：在短的那一侧按需求补碱基，碱基从基因组读。
2. **验证"只增不移"**：旧前体必须仍出现在坐标指示的位置上。若不是，说明基因组版本不符，
   该条不动（`sequence_mismatch`）。
3. **重新折叠**：对扩后的前体和它的 ±20 nt 窗口各跑一次 RNAfold。
4. **重算 star**：用新前体调 `compute_star`，并报告是否移动。

成熟序列完全不动，所以 bowtie 参考与比对结果继续有效，**不需要第三轮比对**。

### 写回哪些列

| 列 | 内容 |
|---|---|
| 7–10 | ±20 nt 窗口（起止 / 序列 / 结构） |
| 11–14 | 前体（起止 / 序列 / 结构） |
| 20–22 | star 坐标 / 序列（**仅当它移动时**） |

成熟三列（16/17/18）、名称、家族、MIR_ID 一律不动。

---

## 3. 改动清单（4 新 + 6 改）

| 文件 | 类型 | 内容 |
|---|---|---|
| `src/utils/genome_fasta.py` | 新增 | 自建字节偏移索引的 FASTA 随机访问，索引缓存写在 fasta 旁边，写不了就落到临时目录 |
| `src/utils/precursor_repair.py` | 新增 | 扩增决策、转录向偏移校验、护栏、star 重算，纯 stdlib |
| `src/repair_precursor.py` | 新增 | CLI，可独立对任意 22 列 basic-info 使用 |
| `tests/test_precursor_repair.py` | 新增 | 单测：基因组读取（多行、双链）、上下侧扩增、拒绝路径、写回列范围 |
| `src/commands/annotation.py` | 改 | part2 在边界纠正之后、聚类之前调用修复 |
| `src/utils/mature_boundary.py` | 改 | 文档：basic-info 区间由 "1-based" 更正为 **0-based** |
| `src/audit_mature_boundary.py` | 改 | 同上 |
| `src/version.py` | 改 | 3.1.8d |
| `CHANGELOG.md` / `README.md` / `docs/DATA-INDEX.md` / `test.sh` | 改 | 版本串 + 单测接入 |

---

## 4. 为什么自带 FASTA 读取

`samtools faidx` 与 `bedtools getfasta` 都想在 FASTA 旁边写 `.fai`，而基因组常在只读共享目录
（本例就在 `xuzx` 名下）；pyfaidx 又是新依赖。所以自己扫一遍建索引：360 MB 基因组一次扫描
约 2 秒，之后每次取序列只读几百字节。

**注意坐标**：`fetch()` 用 **0-based 半开区间**，与 basic-info 一致。

---

## 5. 顺带纠正的一个文档错误

3.1.8c 的注释写 basic-info 是 "1-based start, exclusive end"。**实测是 0-based。**

用基因组检验三种约定：

| 约定 | 前体序列命中 | 成熟序列命中 |
|---|---|---|
| 1-based 起点 + 排他终点 | 0/446 | 0/446 |
| **0-based 起点 + 排他终点** | **446/446** | **446/446** |
| 1-based 双端包含 | 0/446 | 0/446 |

mini 数据集（拟南芥）同样只认 0-based（19/19）。

3.1.8c 的计算没出错，因为它只用两个坐标的**差值**，前后惯例一致；但**任何拿这些坐标去切
基因组的代码**（包括本次新增的修复步骤）必须按 0-based，否则整体偏 1 nt。注释已改正。

---

## 6. 验证

### 6.1 单元测试

`tests/test_precursor_repair.py`：all checks passed（本地 + 服务器各跑一次，均已接入 `test.sh`）。

覆盖：
- 基因组随机访问：整条 / 片段 / 跨行 / 反向互补 / 末位 / 越界 / 未知染色体 / 空区间
- 扩增决策：已包含 / 单侧越界（高低各一）/ 双侧越界
- 转录向偏移：`+` / `-` 各两种
- 拒绝路径：超限（`extension_too_long`）、基因组不符（`sequence_mismatch`）、折叠失败（`fold_failed`）
- 写回：只有窗口与 star 列被改，成熟列与身份列不动（逐列断言）

### 6.2 端到端（`test.sh`，mini 数据集）

**ALL TESTS PASSED**，含新增的 `Precursor-repair unit tests passed.`

默认运行（不启用纠正）与 3.1.8b 基线**逐字节一致**：

| 产物 | md5 | 基线 |
|---|---|---|
| `mirdp3-test-basic-info` | `f9223f3c2b120d62ec24fa60095aca09` | 同 |
| `mirdp3-test-basic-info-cluster` | `db799e4e78ea901bcb75fb2d0f160dbd` | 同 |

且**不生成** `*-precursor-repair.tsv`（未触发即不动）。

### 6.3 真实数据（case1，21 样本，446 条）

| 步骤 | 结果 |
|---|---|
| 边界纠正 | `shift: 2`，并报 `overhangs the annotated precursor: 2 (Nip-MIR166a, Nip-MIR166b)` |
| 前体修复 | `mature inside precursor: 444 / needing repair: 2 / extended: 2` |
| 最终表自检 | 成熟越出前体 **0** 条；成熟非前体子串 **0** 条；结构串长度不符 **0** 条；前体区间与序列长度不符 **0** 条 |
| 改动范围 | 原始 → 最终只变 2 行，变化列 = {6,7,8,9,10,11,12,13,15,16,17}（后者来自前一步的边界纠正） |
| 修复步骤本身 | 修正后 → 修复后只变 2 行，变化列 = {6,7,8,9,10,11,12,13}，**成熟四列与修正后完全一致（0 行差异）** |
| 运行耗时 | 报告模式 2.1 s（含索引构建） |

### 6.4 修复分支端到端（构造输入）

要触发修复，需要"成熟整体后移"的 reads，而旧夹具（掐掉 read 首碱基）产生的是 trim，不会越界。
所以新写了一个夹具：先用 bowtie 定位每条 read，再取**下移 1 nt、长度不变**的基因组窗口
（3235 条），这才是 MIR166a 的真实情形。

```
[boundary] shift: 18, ambiguous: 1
[boundary] corrected mature overhangs the annotated precursor: 14 (...)
[precursor] loci: 19  mature inside precursor: 5  needing repair: 14
[precursor]   extended: 14
[boundary] 复检: canonical: 18, flagged: 1 (ambiguous)
```

核验（修复前 vs 修复后）：**问题计数 无**；改动 14 行、5 行逐格不动；变化列全在允许集合内；
14 条序列均与基因组一致、结构串长度自洽、成熟全部进入前体；**star 全部未移动**。
下游 19 = 19 = 19 行，聚类正常。

---

## 7. star 为什么不动（这不是巧合）

两条都是**在转录 3' 端**差 1 nt：

- `Nip-MIR166a`（`+`）：前体末端扩 1
- `Nip-MIR166b`（`-`）：前体起点减 1 = 转录 3' 端扩 1

而 star 位于**对侧臂（转录 5' 端）**，扩 3' 端不影响它。14 条端到端夹具里同样全部 `star_moved=no`，
与 `compute_star` 的重算结果一致。

若将来遇到 5' 端扩增（`shift_upstream` 被 `--allow-extend` 放行时），star 会真的移动，
那时 20–22 列会被一并更新 —— 这条路径已实现，但默认不开。

---

## 8. 护栏

| 护栏 | 默认 | 行为 |
|---|---|---|
| `--max-extension` | 5 nt / 侧 | 超过则 `extension_too_long`，不动 |
| 基因组不符 | — | `sequence_mismatch`，不动 |
| RNAfold 失败 | — | `fold_failed`，不动 |
| 成熟不在前体臂上 | — | 扩增照做，star 保留原值并报 `star_unresolved` |
| 前体越出染色体 | — | `out_of_genome`，不动 |

---

## 9. 使用方式

### 正常流程（推荐）

```bash
./mirdeep-p3 annotation ... --fix-mature-boundary
```

修复自动跟随纠正执行，产出：

```
<prefix>-boundary-audit.tsv                 审计（纠正前）
<prefix>-boundary-audit.after-fix.tsv       复检（纠正后）
<prefix>-precursor-repair.tsv               前体修复报告
temp/<prefix>-basic-info.before-boundary-fix      纠正前原件
temp/<prefix>-basic-info.before-precursor-repair  修复前原件
```

### 单独使用

```bash
python src/repair_precursor.py \
    -i sample-basic-info -g genome.fasta \
    -o sample-precursor-repair.tsv --apply -oc sample-basic-info.repaired
```

先不加 `--apply` 跑一遍看报告，确认 `needing repair` 与预期一致再写入。

---

## 10. 未做 / 后续

- 未提交 Git、未打 tag、未改 conda / docker / bioconda（按赵哥要求）。
- `data/index` 仍是指向旧索引的软链接，打包前需处理。
- 端到端只覆盖 22 列布局（`repair_precursor.py` 对 24 列会明确报错而不是猜）。
- `Nip-MIR12457a` 存在 star 与成熟重叠 1 nt 的**既有**现象（原始表就有，与本次无关），未处理。
