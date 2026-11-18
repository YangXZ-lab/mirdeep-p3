# 3.1.8e 验证报告 — 成熟/star 臂间距筛选

> 基线：3.1.8d（`/nfs/users/zhaojw/software/miRDeep-p3-3.1.8d`）
> 交付：`/nfs/users/zhaojw/software/miRDeep-p3-3.1.8e`
> **未提交 GitHub**，未改 conda / docker / bioconda。

---

## 1. 问题

`compute_star()` 从前体的另一端取 `star_len` 个碱基，只校验成熟贴在某一端，
**从不检查两条臂之间是否还留得下 loop**。当前体不是折叠出来的发夹、而是绕成熟
硬切出来的窗口时，两条臂就会互相压上去：

```
Nip-MIR12457a   前体 3742249..3742289   40 nt
                成熟 3742249..3742269   20 nt
                star 3742268..3742289   21 nt      20 + 21 = 41 > 40
```

star 起点落在成熟内部，存下来的点括号描述的是两条臂都对不上的结构。
下游不会报错——这是最麻烦的地方。

---

## 2. 判据

**净间距 gap**（即 loop 长度）：

| 位置关系 | gap |
|---|---|
| star 在 mature 之后 | `star_start − mature_end` |
| star 在 mature 之前 | `mature_start − star_end` |
| 重叠 | 负值（重叠碱基数） |

剔除条件（任一命中）：

| 状态 | 含义 |
|---|---|
| `arm_overlap` | gap < 0，两臂互相覆盖 |
| `short_arm_gap` | 0 ≤ gap < `--min-arm-gap`（默认 10） |
| `missing_precursor_structure` | 前体结构缺失或长度与序列不符 |
| `missing_stem_structure` | ±20 窗口结构缺失或长度不符 |
| `missing_star` | star 序列为空 |
| `bad_coordinates` | 坐标非整数或区间空/倒置 |

阈值 10 nt 是**宽松**的：真实植物发夹 loop 通常 15–100+ nt。

**判据的实际含义**要写清楚：它检验的是"前体有没有空间装下两条臂加一个 loop"，
**不是**"star 是否真实表达"。star 在多数植物里是推算而非实测，这是结构合理性检查。
`--min-arm-gap 0` 可退化为"仅要求不重叠"。

---

## 3. 位置：part 1 末尾，任何比对之前

```
part1: ... 4.9 结构提取 → 4.10 basic-info → **[4.11 结构筛选]** → 返回
part2: 比对 → 边界审计 → 前体修复 → 聚类 → 定量
```

这样 part2 的**所有**计算都只针对存活条目：不为将被剔除的位点做比对、不做边界审计、
不折叠发夹、不进表达矩阵。这是赵哥特别指出的要点。

`temp/<prefix>-basic-info.before-structure-filter` 保留筛前原件，可追溯。

---

## 4. 改动清单（3 新 + 6 改）

| 文件 | 类型 | 内容 |
|---|---|---|
| `src/utils/structure_filter.py` | 新增 | 判据本体，纯 stdlib；`arm_gap` / `assess` / `filter_rows` / `report_line` / `summarize` |
| `src/filter_structure.py` | 新增 | CLI，可对任意 22 列 basic-info 独立使用 |
| `tests/test_structure_filter.py` | 新增 | 单测：四种几何情形、阈值边界、各拒绝分支、报告列数 |
| `src/commands/annotation.py` | 改 | part1 末尾调用筛选；新增 `--min-arm-gap`（含默认值与负值校验） |
| `src/version.py` | 改 | 3.1.8e |
| `CHANGELOG.md` / `README.md` / `docs/DATA-INDEX.md` / `test.sh` | 改 | 版本串 + 单测接入 |

---

## 5. 验证

### 5.1 单元测试

本地与服务器各跑一次，均 `all checks passed`。覆盖：

- `arm_gap`：star 在后 / star 在前 / 紧邻（0）/ 重叠 1 / 重叠 5，含 `Nip-MIR12457a` 真实几何
- 阈值边界：gap 正好 10 保留、gap 9 剔除；`min_gap=0` 时 gap 9 保留；`min_gap=15` 时剔除
- 各拒绝分支：无 star、无前体结构、前体结构长度不符、无窗口结构、窗口结构长度不符、空区间、非整数
- `filter_rows`：保留数、报告覆盖每一输入行、顺序保持、汇总计数
- 报告行：13 列、名称/状态/gap 正确
- 保留行**逐字节不变**

### 5.2 端到端（`test.sh`，mini 数据集）

**ALL TESTS PASSED**，含新增 `Structure-filter unit tests passed.`

mini 数据 19 条全部通过筛选：

```
[structure] loci: 19  kept: 19  removed: 0 (min arm gap: 10 nt)
```

产物 `mirdp3-test-basic-info` md5 = `f9223f3c2b120d62ec24fa60095aca09`，与 3.1.8b
基线**一致**；`mature.count` / `mature.exp` 均 19 行，与 basic-info 对齐。
说明筛选在小数据上零影响，不误杀。

### 5.3 真实数据（case1，446 条）

| 项目 | 结果 |
|---|---|
| 筛选 | `446 → 432`，剔除 14 条（`arm_overlap` 1、`short_arm_gap` 13） |
| 保留行完整性 | 432 行与原始**逐字节一致**（筛选只删不改） |
| 剩余集复检 | 最小 gap = 13，`gap < 10` 的条目 **0** |
| 运行耗时 | 秒级 |

剔除名单（14 条）：

| 名称 | 状态 | gap | 前体 |
|---|---|---:|---:|
| `Nip-MIR12457a` | arm_overlap | −1 | 40 |
| `Nip-MIR2120a/b/d` | short_arm_gap | 2 | 45 |
| `Nip-MIR12389a` | short_arm_gap | 8 | 50 |
| `Nip-MIR12753a–e` | short_arm_gap | 8 | 50 |
| `Nip-MIR12325a–d` | short_arm_gap | 9 | 50 |

**共同特征**：前体长度恰好等于 `mature + star + loop`（如 `MIR12389a`：50 = 21 + 21 + 8），
远低于全库前体长度 p5 分位（64 nt）。前体不是折叠产物，是绕成熟硬切的窗口。

### 5.4 为什么这 14 条只能删、不能扩

3.1.8d 有扩前体的能力，本来可以试试救回。实测否掉：**14 条的成熟全部落在 ±20 nt
窗口的中部，一条都不贴臂**（`MIR12457a` 40/80、`MIR2120a` 43/85、`MIR12389a` 20/90、
`MIR12753a` 49/90）。而 `compute_star` 要求成熟贴在某一端（±2 nt），所以把前体扩到
整个窗口后它会直接返回 `None`——这些位点根本没有"以该成熟为臂"的发夹。
扩前体无效，删除是唯一正确处理。

---

## 6. 与 3.1.8d 的兼容性

| 场景 | 行为 |
|---|---|
| mini 测试数据 | 19/19 通过，产物与 3.1.8b 基线**逐字节一致** |
| case1 真实数据 | 446 → 432 行 |
| 默认值 | `--min-arm-gap 10`；`0` 退化为"仅剔除重叠" |
| 保留行 | 永不修改，逐字节原样写出 |

**注意**：这是 3.1.8x 系列里**第一个默认会改变行数**的版本。若下游脚本假定
basic-info 行数等于 identification 的 `-all-mod_fp_prediction` 条数，需注意差值来源，
并读 `<prefix>-structure-filter.tsv` 确认。

---

## 7. 使用方式

### 正常流程（默认已启用）

```bash
./mirdeep-p3 annotation ...          # 默认 --min-arm-gap 10
./mirdeep-p3 annotation ... --min-arm-gap 0    # 只剔除重叠
./mirdeep-p3 annotation ... --min-arm-gap 15   # 更严格
```

产出：

```
<prefix>-structure-filter.tsv                     每条输入位点的判定（含剔除原因）
<prefix>-basic-info                               筛后表
temp/<prefix>-basic-info.before-structure-filter  筛前原件
```

### 单独使用

```bash
python src/filter_structure.py -i sample-basic-info \
    -o sample-structure-filter.tsv --apply -oc sample-basic-info.filtered
```

先不加 `--apply` 看报告，确认剔除数量与名单符合预期再写入。

---

## 8. 未做 / 后续

- 未提交 Git、未打 tag、未改 conda / docker / bioconda（按赵哥要求）。
- `data/index` 仍是指向旧索引的软链接，打包前需处理。
- 端到端只覆盖 22 列布局（`filter_structure.py` 对 24 列明确报错而不是猜）。
- 判据只做结构合理性，不评估 star 表达证据。若要进一步收紧，可考虑要求 star 区域
  有 reads 支持，但那会把大量正常条目误杀（star 常不表达），不建议。
- 3.1.8d 修好的 `Nip-MIR166a/b` 前体越界问题在 3.1.8e 中保持有效，两者判据不冲突
  （臂间距 vs 包含性）。
