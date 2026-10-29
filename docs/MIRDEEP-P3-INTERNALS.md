# miRDeep-P3 代码知识库（内部实现笔记）

> 用途：让 AI 助手（或任何接手的人）在**不重新通读全部源码**的前提下，准确回答"某功能在哪实现、某阈值在哪改、某输出列是什么含义"。
> 版本基准：`src/version.py: __version__ = "3.1.8a"`，开发树 `/nfs/users/zhaojw/software/miRDeep-p3-3.1.7a`（git，base = `fb0b754` / 3.1.6c）。
> 3.1.7a/b/c 与 3.1.8a 的增补点在各节内就地注明（如 §4.4.1、§6、§14）；其余结论来自 3.1.6c 源码静态通读，仍适用。
> 维护约定：改代码后请同步更新本文对应小节，尤其是第 4/5/7/8/10 节的流程与常量表。
> 待办与问题台账（编号 I-01…I-24 / O-01…O-10；其中 O-01/O-02/O-03/O-08 已于 3.1.8a 落地）见 `docs/MIRDEEP-P3-ISSUES.md`。

---

## 0. 速查（先看这一节）

| 项目 | 值 |
|---|---|
| GitHub | https://github.com/YangXZ-lab/mirdeep-p3 |
| GitHub 主干副本 | `/nfs/users/zhaojw/software/github/git/mirdeep-p3` |
| 当前开发树（3.1.7c） | `/nfs/users/zhaojw/software/miRDeep-p3-3.1.7a`（git 仓库，单个本地 commit + 后续改动） |
| 入口 | `./mirdeep-p3`（无扩展名，Python3，可执行） |
| 版本唯一真源 | `src/version.py` 的 `__version__` |
| 数据索引版本真源 | `data-index.env` 的 `DATA_INDEX_RELEASE`（当前 `mirdeep-p3-v3.1.4c-full`） |
| 三步 | `identification`(1) → `annotation`(2) → `analysis`(3) |
| 最小跑通 | `./test.sh -t 8`（用 tests/data/mini.fq + mini_genome.fasta） |
| 输出根 | `identification/`、`annotation/`、各 analysis 子命令自己的 `-o` |

**同机存在三份副本，改代码前务必确认运行命令指向哪一份**：

1. `/nfs/users/zhaojw/software/github/git/mirdeep-p3` —— 最新（本文基准）
2. `/nfs/users/zhaojw/software/github/miRDeep-P3` —— 较新副本（入口名无扩展名）
3. `/nfs/users/zhaojw/software/github/mirdp3` —— 较旧（入口 `mirdeep-p3.py`）

三份不是软链接，`PROJECT_ROOT = Path(__file__).resolve().parent` 各取自己的根，改一份另一份不生效。

---

## 1. 软件定位与设计哲学

- **来源**：miRDeep-P2（Perl 原版）的 **Python 重写 + 植物化改造 + 下游分析套件**。作者 Jiawen Zhao，版权 2026。
- **保留**：miRDeep 的核心统计打分（Gumbel 分布拟合 MFE、频次 log-odds、star 证据、保守性）。
- **改造**：把 miRDeep-P2 里"植物"专用的硬判据抽出来放进 `mod_rm_redundant_meet_plant.py`；把毛虫/果蝇才需要的 `--drosha`、`--randfold` 降级为可选；负链成熟体在核心算法里不支持（见 §7.3）。
- **下游**：annotation 做 PmiREN 家族注释与命名，analysis 做统计/靶标/TFBS/差异表达/富集，形成"从 reads 到文章图表"的闭环。
- **工程约定（3.1.6c 起严格执行）**：
  - 退出码语义：`0=成功`、`1=失败`、`2=PARTIAL_SUCCESS`、`3=NO_TERMS`。
  - "空结果"**不是错误**：写空文件 + exit 0。
  - 单样本失败不拖垮整批：按样本隔离 + 末尾 `[summary] ok/empty/failed` 汇总。
  - 所有拼进 shell 的路径必须过 `src/utils/shellquote.py: shq()`。

---

## 2. 仓库结构

```
mirdeep-p3/                     入口（argparse 分发，17979 B）
├── README.md                   用户手册（25958 B）
├── RELEASE.md                  发版 runbook（13966 B）
├── CHANGELOG.md                变更记录（6435 B）
├── CONTRIBUTING.md
├── test.sh                     本地验证（5541 B）
├── data-index.env              DATA_INDEX_RELEASE=...（勿删对应 Release）
├── environment.yml             conda 完整锁定（391 行）
├── mirdp3_environment.yml      人工维护的宽松约束（推荐建环境用）
├── Dockerfile
├── conda.recipe/{build.sh,meta.yaml}
├── .github/workflows/{ci.yml,docker-release.yml}
├── src/
│   ├── commands/{identification.py (49957 B), annotation.py (29932 B), analysis/*}
│   ├── mod_miRDP.py                    核心打分引擎（47794 B）
│   ├── mod_rm_redundant_meet_plant.py  植物判据 + 坐标映射（29511 B）
│   ├── excise_candidate.py             切前体候选（21829 B）
│   ├── primary-basic-info.py           basic-info 24 列生成（19437 B）
│   ├── filter_alignments.py            read 级硬门限（16382 B）
│   ├── convert_bowtie_to_blast.py      bowtie→BLAST 格式（14415 B）
│   ├── mirdp3_all_pred.py              多样本合并 RPM（11862 B）
│   ├── preprocess_files.py / preprocess_reads.py / remove_nr_reads.py
│   ├── extract_mature.py / extract_mature_flank.py / extract_struc.py
│   ├── expression_matrix.py / scan_cluster.py / rearrange_result.py
│   ├── predictions_annotate.py / predictions_annotate_filter.py
│   ├── all_pred_to_bed.py / all_pred_to_bed_sort.py / convert_to_bed.py
│   ├── check_python_deps.py
│   ├── version.py
│   ├── utils/{config,dependencies,help_display,shellquote}.py
│   └── MiRNATarget/{parse_mirna_targets.py, parse_ssearch.py, README.md}
├── scripts/                    R 为主的下游脚本 + 构建用 Python
├── bin/                        PmiREN core 构建 / 注释 / 去重工具箱
├── data/
│   ├── PmiREN-20260810-isoform.fa    core 数据集（94038 B）
│   ├── ko00001.json                  KEGG 层级（5.6 MB）
│   ├── motif/                        80+ 物种 × 3 文件（TF_list / meme / bg）+ sp.map
│   └── index/                        ★ 需从 data-index Release 下载（当前缺失）
├── examples/PmiREN2.0_basic_info_MIR.fasta / ..._unqiue.fasta（注意 "unqiue" 是原样拼写）
├── tests/data/{mini.fq 49MB, mini_genome.fasta 31MB}
└── docs/{core-dataset.md, file_process.md, reannotation.md}
```

> `data/index/` 当前**不存在**：任何 identification/annotation 首次运行都会尝试下载/构建索引。离线环境需先按 README 放好。

---

## 3. 入口层 CLI（`./mirdeep-p3`）

### 3.1 启动顺序
1. `PROJECT_ROOT = Path(__file__).resolve().parent`；`SRC_DIR = PROJECT_ROOT/"src"` 插到 `sys.path[0]`。
2. `from version import __version__`；`from utils.config import load_config`；`from commands import identification, annotation, analysis`。
3. 打印 BANNER + GitHub 链接。
4. `preprocess_argv()`：把老写法 `mirdeep-p3 step 1 ...` 改写成 `--step 1`（向后兼容）。
5. `main()` 按子命令分发：`identification|1`、`annotation|2`、`analysis|3`。

### 3.2 帮助系统
- 子命令都 `add_help=False`，自己处理 `-h`。
- `handle_help_request()`：
  - `--step` 同时含 1 和 2 且 `-h` → `show_combined_identification_annotation_help()`（内置大段手写帮助，逐参数标注 `[ID]` / `[AN]` / `[BOTH]` 归属）。
  - 单 step `-h` → 该 step 的 `parser.print_help()`。
  - 多 step `-h` → `utils/help_display.print_combined_step_help()`（**仍是占位实现**，文案为 "(Full combined help will be implemented later.)"）。
- `--step 3` 单独使用时只调 `analysis.main(["-h"])` 打印帮助，不执行任何分析。

### 3.3 配置合并
`merge_config_with_args()`：`--config` 读 key:value（格式见 §11.1），**命令行优先**；配置文件缺 `:` 直接抛 `ValueError("Invalid config line {n}: missing ':' separator")`。

### 3.4 一次跑完 1+2
`execute_combined_identification_annotation()`：
- 为两个 step 各建独立 parser，`parse_known_args`。
- 输出根：`<out>/identification` 与 `<out>/annotation`。
- identification 跑完后读 `<id_output>/mirdp3-identification-*.pipe`，**取每行第 2 列**得到各样本目录，逗号拼成 annotation 的 `-i`。
- annotation 缺失的参数用 identification 的 genome / index / replicate / threads / progress 补齐；要求显式给 `--prefix_miRNA` 与 `--species`。
- ⚠️ `main()` 里用 `del cli_args[step_idx - 1]` 连删两次，索引依赖"脚本名已从 argv 去掉"这一假设，边界情况（`--step` 出现在首位）易错位。改动此处务必加测试。

---

## 4. Step 1 · identification 全流程

文件：`src/commands/identification.py`（49957 B，全流程唯一编排者）。

### 4.1 默认参数（改默认值来这里）
```python
DEFAULT_MIN_LEN       = 18     # trim_galore --length
DEFAULT_MAX_LEN       = 26     # trim_galore --max_length
DEFAULT_RPM_THRESHOLD = 5      # 去 ncRNA 后的 RPM 门槛
DEFAULT_MAX_MAPPINGS  = 15     # 每 read 最多基因组比对位点
DEFAULT_PRE_LENGTH    = 300    # 前体窗口长度
DEFAULT_THREADS       = 1
DEFAULT_PROGRESS      = 1      # 并发样本数
```

### 4.2 输入预检（3.1.6c 新增，先验后跑）
- `INPUT_VALIDATION_MAX_RECORDS=4000`、`MAX_LINES=200000`、`FOOTER_BYTES=65536`。
- `validate_input_file()` → `{seq_type, records, eof, gzip, gzip_full, size}`。
- `_scan_fastq()`：校验 `@`/`+`、长度一致、非核苷酸比例 >50% 报错；`_scan_fasta()` 同理。
- `_check_plain_fastq_footer()` 读尾部 64KB 防**截断文件**；`_verify_gzip_stream()` 全量解压校验（`--verify-gzip` 开关）。
- `_NUC_CHARS = set("ACGTUNacgtun-.*")`；错误类型 `InputValidationError`，消息格式 `"<原因> -- <path>"`。
- 全部输入都坏 → `sys.exit`；部分坏 → 只跳过并打 `[skip]`。genome 同样预检。

### 4.3 样本命名
- `decompress_file()`：进程内解压到 `temp/<stem>`，先写 `.part` 再 `os.replace`（防半文件被复用）。
- `_uniquify_prefixes()`：重名加 `-1/-2` 后缀，再兜底去重。
- `resolve_prefixes()` → `(prefixes, group_prefixes)`；`--prefix` 支持三种给法：给 1 个（按组展开 `<p>-<g>-<n>`）、给"组数"个（每组一个）、给"输入数"个（每输入一个）。

### 4.4 单样本流水线（声明式 stage，由 `src/commands/stage_runner.py` 驱动）

3.1.8a 把原来 15 个内联动作收成 **12 个 `Stage`**：每个 stage 声明 `sid`（沿用本文小节号）、`name`、`artifacts(ctx)`（产物路径，resume 判定用）、`run(ctx)`。`run_stages()` 顺序执行、计时、写 manifest；stage 内 `raise StageEmpty(reason)` 表示"空结果"（不是错误）。

| sid | 动作 | 关键命令 / 产物 |
|---|---|---|
| 4.4-4.5 | 解压 + 修剪 + 格式整理 | `temp/<stem>`；FASTA 直接 copy → `<prefix>.processed.fa`；FASTQ 且 `reads_clean` → `trim_galore --small_rna --length {min} --max_length {max} --dont_gzip --suppress_warn -j 8 --basename {prefix} -o {out}` → `<prefix>_trimmed.fq`；再 `src/preprocess_files.py -i <trimmed> -o <temp>` |
| 4.6 | ncRNA 比对 | `bowtie -v 0 data/index/rfam_index`；`bowtie -v 1 data/index/mature_index` |
| 4.7 | 去 ncRNA | `src/preprocess_reads.py <processed.fa> <rfam.aln> <mature.aln> <rpm_threshold> <prefix>.fa <prefix>.processed.fa <prefix>.total_reads` |
| 4.8 | 基因组比对 | `bowtie -a -v 0 [--large-index] <genome_index>` → `<prefix>-processed.aln`（`*.ebwtl` 判 large） |
| 4.9 | 转 BLAST + 过滤 | `src/convert_bowtie_to_blast.py <aln> <prefix>.fa <genome> -o <prefix>-processed.bst` → `src/filter_alignments.py -c {max_mappings}` → `*-processed-filter.bst`（空 → StageEmpty） |
| 4.10 | 切前体 + 折叠 | `src/excise_candidate.py <genome> <filter.bst> -l {pre_length} -o <prefix>_precursors.fa` → `RNAfold --noPS -j{threads}` → `*_precursors.struc` |
| 4.11 | 全量 read 也走一遍 | 对 `<prefix>.fa`（含 ncRNA 的）重复 bowtie→bst→filter→`filter_alignments.py -b <prefix>.fa` → `<prefix>_filtered.fa` |
| 4.12 | 前体索引 + 计数 | `bowtie-build` 建 precursor 索引 → bowtie → bst → `sort +3 -25` → `<prefix>_signatures` |
| 4.13 | 核心预测 | `src/mod_miRDP.py <signatures> <precursors.struc> -o <prefix>_predictions`（空 → StageEmpty） |
| 4.14 | 植物化过滤 | `samtools faidx <genome> --fai-idx temp/<genome>.fai` → `<prefix>/chr_length`；`src/mod_rm_redundant_meet_plant.py <chr_length> <precursors.fa> <predictions> <total_reads> -n <prefix>_nr_predictions -f <prefix>_filter_P_prediction`（非 0 退出 → RuntimeError；空 → StageEmpty） |
| 4.15 | 转 BED | `src/convert_to_bed.py -i <nr_pred> -o <prefix>_nr_predictions.bed` |
| 4.16 | 收尾 | 复制 `<prefix>.processed.fa` 到输出根，并尝试建兼容软链 `<prefix>-processed.fa` |

> ⚠️ 第 4.13 步（`_stage_prediction()`）**不传** `-s/-x/-y/-z`。即：阈值恒为 **1.0**、保守性**不打分**（known_nuclei 为空）、不跑 randfold、不跑 drosha。这是"植物精简模式"，与 Perl miRDeep 默认行为不同，是最容易被误认为 bug 的地方（改进机会见台账 I-16）。

### 4.4.1 resume 与 run manifest（3.1.8a）

- **每样本 manifest**：`<out>/<prefix>/<prefix>.run.json`，原子写（`.part` → `os.replace`）。
- **内容**：`schema`(1)、`version`、`signature`、`input{path,size,mtime,md5,seq_type}`、`genome{path,size,mtime}`、`params`、`tools{bowtie,bowtie-build,RNAfold,samtools,trim_galore}`、`env{python,platform}`、`stages{sid:{name,status,artifacts,artifacts_info,message,seconds}}`。
- **批次 manifest**：`<out>/mirdp3-identification-<mmddyy-HHMM>.manifest.json`，每样本一条（prefix / status / message / manifest 路径），status 取 `ok|empty|failed|partial`。
- **`--resume` 命中条件**（三者同时满足）：① manifest 记的 status ∈ {ok, skipped}；② 记的 artifacts 列表 == 当前算出的列表（防 prefix/参数漂移）；③ 磁盘上产物齐全且非空（`allow_empty` 的 stage 除外）。
- **签名** = md5(version + 参数 + genome size/mtime + 输入 size/mtime + genome_index)。不一致打印 `[resume] <prefix>: input or parameters changed, not resuming` 并整轮重跑；**输入 md5 只在建 manifest 时算一次**，续跑判定只用 size/mtime，所以大文件也是秒级决策。
- **与 `--clean` 冲突**：resume 时强制 `clean=False` 并打印 `[resume] --clean ignored: temp/ is needed to resume.`（resume 依赖 temp 里的中间产物）。
- **已知小瑕疵（台账 I-24）**：批次 manifest 文件名精度到分钟，同一分钟内重跑会覆盖上一条（每样本 manifest 不受影响）。


### 4.5 状态机与汇总
- `PIPELINE_OK="ok"` / `PIPELINE_EMPTY="empty"` / `PIPELINE_FAILED="failed"`。
- worker 内 `try/except` + traceback 写日志，单样本失败不影响其他样本。
- `_write_empty_results()` 触碰三个空文件：`<prefix>_nr_predictions`、`.bed`、`_filter_P_prediction`，并记 `"[tag] {prefix}: {reason}"`。
- 结尾打印 `[summary] ok/empty/failed`；**全部失败才 exit 非 0**。

### 4.6 索引与并行
- 默认建在 `<out>/index/genome_index`；`is_bowtie_index_complete()` 检查 `.1/.2/.3/.4/.rev.1/.rev.2` + `.ebwt`/`.ebwtl`；`-d` 可跳过建索引。日志 `index_build.log`。
- 并行：`multiprocessing.Pool(min(args.progress, len(runnable)))`，args 对象跨进程传。
- `--clean` 只删各样本 `temp/`；索引由 `--clean-index` 删（`-d` 传入的索引不动）。3.1.7c 起帮助文案已同步。

### 4.7 产物清单
| 文件 | 含义 |
|---|---|
| `<out>/mirdp3-identification-<mmddyy-HHMM>.pipe` | 清单，每行 tab 分隔：`<input>\t<out_dir>\t<group_prefix>` |
| `<prefix>/<prefix>.fq_trimming_report.txt` | trim_galore 报告 |
| `<prefix>_identification.log` | 每样本全流程日志（含各工具 stdout） |
| `<prefix>.total_reads` | 该样本有效 read 总数（RPM 分母） |
| `<prefix>_filter_P_prediction` | **最终预测**（8 列，见 §10.2） |
| `<prefix>_nr_predictions` / `.bed` | 去冗余 8 列 / BED6 前体区间 |
| `<prefix>/chr_length` | 染色体长度（下游坐标校验用） |
| `<prefix>/<prefix>.run.json` | 每样本 run manifest（版本/参数/工具版本/输入指纹/各 stage 状态与耗时） |
| `<out>/mirdp3-identification-<mmddyy-HHMM>.manifest.json` | 批次 manifest（每样本一条：status + manifest 路径） |

---

## 5. Step 2 · annotation 全流程

文件：`src/commands/annotation.py`（29932 B）。

### 5.1 参数要点
`-i` 输入目录（identification 产物）、`-f` 列表、`-o`、`--prefix`、`--prefix_miRNA`（**必填**）、`--species`（**必填**）、`-g` genome、`--pmiren`（默认 `PmiREN-20260810-isoform.fa`）、`--pmiren-index`（默认 `index/isoform-in`，相对 `data/`）、`-d`、`-r`、`-t`、`-p`、`-s/--same`、`--consistency`、`--common`、`--clean`，以及 `--bowtie/--blastn/--bowtie-build/--makeblastdb/--RNAfold/--samtools/--bedtools/--seqkit` 自定义可执行文件。

- `validate_input_folders()` 要求目录内存在 `*_filter_P_prediction`；`detect_file_prefix()` 从 stem 去后缀得前缀。
- `resolve_annotation_prefixes()` → `(folder_prefixes, file_prefixes, output_dirs)`；`--same` 时文件夹固定 `case1..caseN` 且 `--prefix` 必须每输入一个；未给 `--prefix` 也自动 `case1..caseN`。

### 5.2 Part 1（保守性注释 + basic-info）
1. 写 `read_total_file`：`<fprefix>\t<total_reads>`。
2. 每样本：`src/mirdp3_all_pred.py -i <fprefix>_filter_P_prediction -c read_total_file -o <fprefix>-all-mod_fp_prediction -m <fprefix>_count_to_rpm`。
3. 多样本按序 `cat` 合并 → `<group>-all-mod_fp_prediction`。
4. `all_pred_to_bed.py` → `all_pred_to_bed_sort.py` → `bedtools merge -d 0 -c 4,5,6 -o collapse,distinct,distinct` → `<group>-merge_output.bed`。
5. `remove_nr_reads.py -m merge -i <all_mod> -o <group>-all-mod_fp_prediction-nr`。
6. `extract_mature.py -i <nr> -genome <g> -o <nr>-mature.bed -fo <nr>-mature.fasta`。
7. `bin/anno_miRNA.py -i <mature.fasta> -p data/<pmiren> -d data/<pmiren-index> -o temp -t <threads> --type MIRN` → 产物在 `temp/temp/anno.map`（注意**双层 temp**，这是历史遗留路径）。
8. `predictions_annotate.py -i <nr> -c <anno.map> -o <nr>-annotated`。
9. `predictions_annotate_filter.py -o <group>-filtered-all-nr-anno`（去掉非保守）。
10. `extract_struc.py -i <filtered> -genome <g> --threads <t> -struc <s> -struc_20nt <s20> -fasta_20nt <f20>`。
11. `primary-basic-info.py -f <filtered> -struc_20nt <s20> -struc <s> -chr_length <第一个输入目录>/chr_length -species <sp> -prefix_miRNA <p> -m <anno.map> -o <group>-primary-basic-info`。
12. `rearrange_result.py` → `<group>-second-basic-info` → 复制为 `<out>/<group>-basic-info`。

### 5.3 Part 2（表达矩阵 + 簇）
1. 从 basic-info 取**第 1 / 15 / 18 列**（0-based 索引 0/14/17，判断 `len(parts) >= 18`）拼 mature.fa：`>{parts[0]}\n{parts[14]}\n{parts[17]}`。
2. `scan_cluster.py -i <basic-info> -species <sp> -prefix <p> -o <group>-basic-info-cluster`。
3. `extract_mature_flank.py` 取成熟序列 ±1 nt → `bowtie-build` 建 `<index>/<group>-mature-1nt`。
4. 每样本 `bowtie -a -v 0` 比对 `<fprefix>.processed.fa`（兼容 legacy `<fprefix>-processed.fa`，再退化为 `rglob("*.processed.fa")`）。
5. `expression_matrix.py -aln a,b,c -c read_total_file -species <sp> -o <group>-mature.exp -oc <group>-mature.count`。

### 5.4 一致性分支（★ 3.1.7c 起可用）
`--consistency` / `--common` 调用 `src/miRNA_consistency.py`（544 行；3.1.7c 补回，此前缺失，触发即 `FileNotFoundError`）。三种模式：
1. **参考重命名**：`-i <basic-info> -b <参考 basic-info> -ob <out>`（有 `-b`、无 `-c`）；
2. **跨组统一命名**：`-i a,b --common -b <参考> -o c1,c2`；
3. **纯合并**：无 `-b`，`-o` 给多个输出文件。
容差：`--mature-distance`（默认 1）、`--pre-distance`（默认 2）。

### 5.5 产物
`<out>/mirdp3-annotation-<time>.pipe`、`<group>-basic-info`（22 列，见 §10.1）、`<group>_annotation.log`（Part1 与 Part2 追加同一文件）、`<group>-basic-info-cluster`、`<group>-mature.exp`、`<group>-mature.count`；`--clean` 删每组 `temp/` 与 `index/`。

---

## 6. Step 3 · analysis 六个子命令

分发器：`src/commands/analysis/__init__.py` 的 `main(argv, project_root)`，subparsers 名：`TFBS`、`Target_finder`、`Differential_expression`（别名 `DE`）、`Functional_analysis`（别名 `FA`）、`Stat`、`Onestep`；存进 `args.analysis_command`，并注入 `args.project_root = project_root`。

> 外部工具（★ 3.1.7c）：`src/utils/toolset.py` 统一解析。`add_tool_arguments(parser)` 给每个子命令加 `--rscript / --fimo / --bedtools / --ssearch36 / --seqkit / --emapper / --python`；`resolve_tools(args, required=[...], optional=[...], hint=...)` 先打印 `External dependency check:` 与每个工具的 `✔ (/path)`，必需工具不可用即 `sys.exit` 并提示该用哪个参数；`python` 默认 `sys.executable`，只有显式传入才校验。`Onestep` 用 `forward_tool_args(args)` 把用户给的路径透传给五步。

### 6.1 Stat
- 参数：`-i` basic-info（必填）、`-o`（默认 `mirdeep-stat-<mmddyy-HHMM>`）、`--rnaplot`、`--list`。
- 命令：`Rscript scripts/basic_stat.R -i <basic-info> -o <out>`；可选 `python scripts/RNA_plot.py --list`。
- 工具：`Rscript`（必需，`--rscript`）；`python` 默认 `sys.executable`。

### 6.2 TFBS
- 参数：`-i` basic-info **或** `-b` BED（互斥、必给其一）、`-o`、`--fai`（必填）、`--list`（列 `data/motif/sp.map` 第 1 列）、`-s/--species`（默认 `Arabidopsis_thaliana`）、`-u/--upstream`（默认 2000）、`-e/--evalue`（默认 1e-6）、`-g/--genome`（必填）、`-p/--picture`。
- 依赖：`bedtools` + `fimo`（必需，`--bedtools/--fimo`）；`-p/--picture` 另需 `Rscript`（`--rscript`），缺失时**警告并跳过作图**，不中断。
- 流程：basic-info → `awk` 取 `$5,$7,$8,$1,$2,$6` 成 BED（→ 正好是 chrom/start/end/name/score/strand 的 BED6）→ `scripts/get_upstream.py -u <up> -o <bed> <fai>` → `bedtools getfasta -s -nameOnly` → 从 `sp.map` 按物种取 `(tf_list, meme, bg)` → `fimo --max-strand --bgfile <bg> --thresh <e> -oc temp <meme> <upstream.fa>` → 两段 awk 关联 `fimo.tsv` 与 TF_list（`NR==FNR` 建 TF→motif 映射）、裁成 8 列（`$1,$NF,$3,$4-$5,$6,$8,$9,$10`）、`gsub(/\([+-]\)/,"",$3)` 去掉区段里的 `(+/-)` → `<out>/tfbs.tsv` → 可选 `Rscript scripts/tfbs_report.R` → `rmtree(temp)`。

### 6.3 Target_finder
- 参数：`-i` FASTA **或** `-b` basic-info（互斥）、`-c/--cds`（必填）、`-o`、`-e/--evalue`（默认 2.5）、`--engine {psrna,ssearch}`（★ 默认 `psrna`，3.1.7b 起）、`--allow-gap/--no-allow-gap`（★ psrna 默认 on）、`--psrna-detail`、`--max-hsp`（ssearch，默认 200）、`--GUs`（ssearch，默认 1）、`--keep-temp`（★ 3.1.7c）、`-t`。
- 依赖：`--engine psrna` **只需 numpy**（无需外部比对器）；`--engine ssearch` 需 `ssearch36` + `seqkit`（`--ssearch36/--seqkit`），缺失时报错并提示改用 psrna。
- 流程：从 basic-info `awk` 取 `$1/$18` 生成 mature.fa → `seqkit seq -r -p` 得反向互补 → 两次 `ssearch36`，选项 `-f -8 -g -3 -E 10000 -T {threads} -b 200 -r +4` 与 `-3 -n -U -W 10 -N 20000` → `scripts/parse_ssearch.py` 各转 tsv → cat 合并 → `src/MiRNATarget/parse_mirna_targets.py --E_cutoff --GUs_cutoff` → `seqkit fx2tab --length --name --header-line` 取长度 → awk 关联并过滤 `$(NF-6) >= $(NF-1)`（成熟长度 ≤ 比对长度）→ 末段 awk 重排 `print $1,$2,$NF,$8,($8+$NF-1),$(NF-2)` → `<out>/target_finder.tsv`（6 列）。
- temp：★ 3.1.7c 起跑完自动删 `<out>/temp`，`--keep-temp` 保留。
- 输出排序：★ 3.1.7b 起两条引擎都按 (mature.fa 顺序, target, start) 稳定排序，同参数两次运行 byte-identical。
- psrna 引擎实现：`src/MiRNATarget/psrna_scan.py`（369 行，numpy only）；规则与验证见 `docs/PSRNA-ENGINE.md`。

### 6.4 Differential_expression（DE）
- 参数：`-c` count（必填）、`-r` rpm（必填）、`--case1/--case2`（**列号**，必填；等长、不重叠、**不得含 1**；若表头第 2 列 == `MIR_ID` 则 2 也不允许）、`--case1name/--case2name`（默认 case1/case2）、`-o`、`--miRNA`（与 `-f` 互斥、与 `--DEOnly` 互斥）、`-f`、`--DEOnly`、`--min-expr`（默认 5.0）。
- 流程：`scripts/extract_col.py` 抽列 → `filter_low_expression()`（任一列 ≥5 才保留，**原地重写** exp/count）→ `Rscript scripts/pca_miRNA.R -i -r <replicates> -o` → `correlation_miRNA.R` → `miRNA_DEG.R -i <count> -r <replicates> --case1 <name1> --case2 <name2>` 产 `<case1name>_<case2name>_res.txt` → `miRNA_dotplot.R`（DEOnly 时先 `awk '$NF!="NOT"{print $1}'` 出列表）→ `miRNA_expression.R -i <exp> --deg <res.txt>` → bash 拼 `<out>/final_miRNA_expression.txt`（`printf` 表头 + `join -t $'\t'` 两个已排序文件，`executable='/bin/bash'` 因为用了进程替换 `<()`）→ `rmtree(temp)`。

- 工具：`Rscript`（必需，`--rscript`）；列抽取与结果拼接用 `python`（默认 `sys.executable`）。

### 6.5 Functional_analysis（FA）
- 两条互斥路线：`-p/--protein`（从零建 OrgDb）**或** `--orgdb` + `-f`；`-g/--gene` **或** `--target`（互斥必给其一）；`--chord` 只能配 `--target`；另有 `--EGGNOG_DATA_DIR`、`--kojson`、`-t`、`-o`。
- 会**先 `rmtree` 输出目录**（注意别指错 -o）。
- `os.environ.setdefault("MIRDEEP_R_LIB", <out>/Rlib)`（保留父进程已导出的值）。作用：`enrich_analysis.R` 会 `install.packages(orgdb, repos=NULL, type="source")` 且**不传 `lib=`**，并发会争写同一个 conda R 库 → 必须私有库 + 清 `00LOCK-*`。
- 建库路线：`sed '/^>/! s/\.//g'` 清洗蛋白 → `emapper.py --cpu <N> -m diamond --override --dbmem -d euk --tax_scope Viridiplantae [--data_dir <dir>]` → `<out>.emapper.annotations` → sed/awk 取 `$1,$9,$10,$12` 生成 `Go.eggnog.emapper.annotations` → `Rscript scripts/build_orgdb.R -i <go> --kojson data/ko00001.json -o <out>` → `org.Morg.eg.db`。
- 富集：gene 路线 `Rscript scripts/enrich_analysis.R -i <orgdb> -f <dir> -g <gene> --goterm 10 -o <out>`；target 路线 `python scripts/miRNA_enrich_analysis.py -i <target> --orgdb <db> -f <dir> -o <out>`，可选 `Rscript scripts/miRNA_chord_type2.R`。
- `run_enrichment_step()` 返回码分级：`0 ok` / `2 partial(PARTIAL_SUCCESS)` / `3 no_terms(NO_TERMS，空结果不算错误)` / 其他 `raise RuntimeError(f"{label} failed (exit {rc}); see the output above")`；partial 时 `sys.exit(2)`。

- 工具：`Rscript`（必需，`--rscript`）；建库路线还需 `emapper.py`（`--emapper`）。

### 6.6 Onestep
- 必填：`-b` basic、`-c` count、`-r` rpm、`--fai`、`-g` genome、`-t` transcript（**注意：这个 `-t` 是给 Target_finder 当 CDS 用的**，命名有历史包袱）、`-p` protein、`-o`、`--case1`、`--case2`。
- 可选：`--case1name/--case2name/--threads/--rnaplot/--tfbsplot/--DEOnly/--chord/--EGGNOG_DATA_DIR/--kojson`。
- 执行顺序：校验列号（禁 1、2）→ **预校验** `EGGNOG_DATA_DIR` 与 `ko00001.json`（失败早退，避免白跑五步）→ `rmtree` 输出目录 → 建 Rlib 并 export → Stat → TFBS → Target_finder → DE（`--DEOnly` 时用 `awk NR==FNR` 过滤 `target_finder.tsv` 得 `target_finder_DEG.tsv`）→ FA；FA 返回 2/3 分别标 partial/note，其它非 0 抛 `RuntimeError` → 清 temp → partial 时 `sys.exit(2)`。

- 工具：接受全部 `--<tool>` 参数并透传给五步；`--target_engine`（默认 `psrna`）透传给 Target_finder。

---

## 7. 核心算法 · `src/mod_miRDP.py`

Perl miRDeep 的逐行移植。类：`RNAfoldResult`、`Alignment`、`HairpinComponent`、`HairpinAnalysis`、`miRDPPredictor`。

### 7.1 BLAST 10 字段解析
`parts[0]=query`、`1=query_length`、`2=query beg..end`、`3=subject`、`4=subject_length`、`5=subject beg..end`、`6=e_value`、`7=identity`、`8=bitscore`、`9=other_info`。
- 链：`other_info` 含 `Minus` → 负链。
- 频次：query 名尾部 `_x(\d+)` → `frequency`（默认 1）。**这是上游 `sort` + `convert_bowtie_to_blast` 协同形成的约定**。

### 7.2 打分常量（改分数分布就改这里）
```python
nucleus_length = 7
score_star       =  3.9 ; score_star_not    = -1.3
score_nucleus    =  3.0 ; score_nucleus_not = -0.6
score_randfold   =  1.6 ; score_randfold_not = -2.2
score_intercept  =  0.3
scores_stem = [-3.1,-2.3,-2.2,-1.6,-1.5,0.1,0.6,0.8,0.9,0.9,0]  # 索引 = stem_bp
```
MFE 分布（Gumbel）：`test: scale 5.5 / location 32`、`background: scale 4.8 / location 23`。
频次分布：`parameter_test=0.999`、`parameter_control=0.6` → `log_odds = log(0.999/0.6)*freq + log(0.001/0.4)`。

`calculate_score()` = MFE 分 + 频次分 + star 分(±) + 保守分(仅当 `known_nuclei` 非空) + stem 分(仅 `--drosha`) + `0.3` + randfold 分(仅 `--randfold`)。

> `_score_mfe()`：`mfe >= 0` 直接给 `-10`；否则 `log(P_test/P_bg)`，`background==0` 时返回 5.4。
> ⚠️ `_cdf_gumbel` 写成 `e**(-(e**(-(x-location)/scale)))`，与标准 Gumbel CDF 的指数符号不同。**这是 Perl 原版照抄**，改动等于改整个分数分布，不要"顺手修正"。

### 7.3 结构解析与组件定位
- `parse_rnafold_output()`：序列转大写且 `U→T`；结构行 `re.sub(r'[^\.\(\)]','',line)` 只留 `.`/`(`/`)`；MFE 用 `re.search(r'\(([-\d\.]+)\)\s*$')` 从行尾抓。
- `parse_known_mirnas()`：`seed = seq[1:8]`（跳过第 1 位、取 7 nt），U 化后入 `known_nuclei`。
- `analyze_hairpin()`：先精确匹配 RNAfold ID，失败则去 `_\d+$`，再试 `_1.._5`；无 alignment 或无 arm → `None`。**mature = 按 `(frequency, query)` 降序排序后的第一条** alignment 的 subject 区间。
- `_determine_arm()`：mature 结构只含 `.`/`(` 且有 `(` → `"first"`（5' 臂）；只含 `.`/`)` 且有 `)` → `"second"`（3' 臂）；否则失败。
  **mature 在负链时直接返回 `""` → 负链 miRNA 在 identification 阶段不会产出**（植物负链 miRNA 极少，属有意取舍）。
- `_find_star_sequence()`：以成熟 3' 端留 2 nt overhang（`mature_end-2`）为锚，在 `base_pairs` 找配对位作 `star_beg`；`star_end = 配对位 + offset + 2`；裁剪到 `pri_end`。
- `_find_loop_sequence()` / `_find_flanking_sequences()`：按 arm 组合推边界；`flank_first` 从 1 到 `mature.beg-1` 或 `star.beg-1`；`flank_second` 从 `mature.end+1`/`star.end+1` 到 `pri_end`。
- `_score_stems_drosha()`：取 `flank_first` 末 10 nt 与 `flank_second` 前 10 nt，`stem_bp = min(两段括号计数)`。

### 7.4 两级初始过滤
**结构级** `_pass_filtering_structure()`：
1. 组件齐全（无 "no mature / no star / no loop / no flanks"）。
2. `_check_no_bifurcations()`：loop ≥20 nt 时检查 `mature + loop[:5]` 与 `loop[-5:] + star` 是否为纯 `.`+`(` 或 `.`+`)`；loop 不足 20 只查 mature/star 本身。
3. duplex 碱基对 `< 14` → `"too few pairings in duplex"`。
4. `|mature_len - star_len| >= 6` → `"too big difference between mature and star length"`。

**signature 级** `_pass_filtering_signature()`：
- 按 `subject_beg` 排序遍历 alignment，`_classify_read()` 判定 mature/star/loop（`fuzz_beg=2`、`fuzz_end=5` 的 contained 判定）；**负链 read 直接算 inconsistent**。
- `inconsistent/total > 0.2` → `"inconsistent\t{n}\nconsistent\t{n}"` 淘汰。
- star 证据：star 类 read 且 `_check_star_overhang()`（`beg - star.beg ∈ [-1,0,1]`）的频次 > fuzzy 频次 → `star_read_present=True`。

### 7.5 入口参数（★ 调整打分的地方）
```
positional: signature_file  structure_file
-o/--output (required)      -s/--known-mirnas      -t/--filtered-output
-u/--limited-output          -v/--threshold (1.0)   -x/--sensitive
-y/--randfold                -z/--drosha            --debug
```
- `process_all()`：`filtered_output` 时被淘汰的也写；否则 `score >= threshold` 才保留。
- `write_output()`：每条写 filter_message（如有）、`score_message`（`score_mfe\t.., score_freq\t.., score_star\t.., score_nucleus\t.., score\t..`）、各组件 beg/end/seq/struct、`pri_*`、原始 alignment 行，**空行分隔**。

---

## 8. 植物化后处理 · `src/mod_rm_redundant_meet_plant.py`

### 8.1 输入解析
- `parse_miRDeep_result()`：按空行切块；含 `_x\d+\s+\d+\s+\d+\.\.\d+\s+` 的行归 `alignment_lines`，其余 `key\tvalue`（或双空格）归 `pred` dict；块里有 `pri_id` 或 `mature_query` 就算一条预测。

### 8.2 坐标映射
`calculate_genomic_positions()`：
- precursor 名匹配顺序：精确 → `startswith(pri_id)` → 去掉最后一段下划线。
- 正链 `mature_bg_abs = excise_beg + mature_beg_rel - 1`；负链 `mature_bg_abs = excise_end - mature_end_rel`。
- `mature_arm == 'first'`：mature 区间整体 `-1` 且 `end+1`（补偿 3' overhang），`pre_bg..pre_end` 取 `mature_bg..star_end`；`'second'` 时 `pre` 取 `star_bg-1..mature_end`。
- 任一坐标 `<= 0` 或 chrom 不在 `chr_length` → 丢弃该条；**全部丢弃 → 写空文件 + exit 0**（3.1.6c 之前是报错）。

### 8.3 read 分布与 RPM
`analyze_read_distribution()`：从 `alignment_lines` 解析 `_xN` 频次；位置区间 `re.match(r'(\d+)\.\.(\d+)$')`；落在 star 区间 ±1 → `star_count`；落在 mature 区间 ±1 → `mature_count`；其余取最大值为 `biggest_other`；`mature_rpm = mature_count/total_reads*1e6`。

### 8.4 植物判据（★ 硬编码，全在这里）
`apply_plant_criteria()`：
1. `mature_struct` / `star_struct` **去掉末 2 nt**（3' overhang 区）后统计 `.` 数作为错配数。
2. mature 错配 `<=5`、star 错配 `<=5`、且 `|两者差| <= 3`。
3. `total_count != 0`。
4. `(mature_count + star_count) >= 0.75 * total_count`。
5. 成熟长度必须 **20–24 nt**。
6. 若长度是 **23 或 24**：额外要求 mature/star 错配都 `!= 0`、`mature_rpm >= 20`、必须有 star read。

### 8.5 输出
- `write_predictions()` 8 列：`chrom, strand, mature_query, pri_id, "mature_bg..mature_end", "pre_bg..pre_end", mature_seq, pre_seq`。
- `write_bed_file()`：BED6，`start = pre_bg-1`、`end = pre_end`、`name = {pri_id}_{arm}`（正链 first→`'5'`、second→`'3'`；负链反过来），score 固定 0。
- CLI：`chr_length precursor miRDeep_result total_reads` + `-n/--nr-output`（必填）、`-f/--filtered-output`（必填）、`-b/--bed-output`、`--debug`、`--stats`。

---

## 9. 坐标与链的约定（最容易出 bug 的地方）

1. **bowtie → 1-based**：`convert_bowtie_to_blast.calculate_alignment_coordinates()` 中**正负链一致**：`start = position+1`、`end = position+read_length`，负链只把 `strand_info` 标成 `"Plus/Minus"`，**坐标统一报在正链上**。
2. **负链切前体**：`excise_candidate.excise_candidates()` 先 `revcom` 整条染色体，再把正链坐标换成 `L-end+1 .. L-beg+1` 在反向序列上切，切完再换回（`is_negative_strand=True, original_seq_length`）。
3. **相对 → 绝对**：`mod_rm_redundant_meet_plant.calculate_genomic_positions()`，正负链公式不同（§8.2），且 `first` arm 有 `-1/+1` 的 overhang 补偿。**动这里必须回归 star 坐标**（3.1.6c 修的就是负链 star 被放到前体错误一端，测试集 171/171 条 star 与 mature 重叠）。
4. **BED 是从"pre_bg/pre_end"推的**（`start = pre_bg-1`），而 basic-info 的 `ext_start/ext_end` 是前体 ±20 nt 的扩展窗口 —— 两者不是一回事，别混用。
5. **`chr_length` 文件**由 `samtools faidx` 生成，annotation 阶段**取第一个输入目录的** `chr_length`（多组不同物种/基因组混跑时需注意）。

---

## 10. 关键数据格式

### 10.1 basic-info（最终 22 列；中间文件 24 列）
**两个阶段，别混**：`src/primary-basic-info.py` 先写 24 列（`out_fields`），再由 `src/rearrange_result.py: reorder_line()` 重排并丢掉 `confidence`/`source`：

```python
order = [0, 1, 2, 3, 4, 7, 5, 6, 8, 9, 12, 13, 10, 11, 14, 15, 16, 17, 18, 20, 21, 19]
```

`<group>-second-basic-info`（= 最终 `<group>-basic-info`）是重排后的 22 列文件，**下游一律按这 22 列取**（已用 `tests/test_output/annotation/mirdp3-test/mirdp3-test-basic-info` 逐列核对）：

| 列 | 字段 | 中间文件 0-based 索引 |
|---|---|---|
| 1 | `miR_id` | 0 |
| 2 | `MIR_ID` | 1 |
| 3 | `family_name` | 2 |
| 4 | `species` | 3 |
| 5 | `chrom` | 4 |
| 6 | `strand` | 7 |
| 7 | `ext_start` | 5 |
| 8 | `ext_end` | 6 |
| 9 | `stem_20` seq | 8 |
| 10 | `stem_20` struc | 9 |
| 11 | `pre_start` | 12 |
| 12 | `pre_end` | 13 |
| 13 | `pre_seq` | 10 |
| 14 | `pre_full` struc | 11 |
| 15 | `mature_mir_id` | 14 |
| 16 | `mature_start` | 15 |
| 17 | `mature_end` | 16 |
| 18 | `mature_seq` ★ | 17 |
| 19 | `star_id` | 18 |
| 20 | `star_start` | 20 |
| 21 | `star_end` | 21 |
| 22 | `star_seq` | 19 |

（`confidence`=索引 22、`source`=23，重排时丢弃，所以最终只有 22 列。）

**下游取列语义**（按最终 22 列，不要凭感觉改列号）：
- TFBS：`$5,$7,$8,$1,$2,$6` = chrom/ext_start/ext_end/name/MIR_ID(score)/strand → BED6。
- Target_finder：`$1`（名）+ `$18`（成熟序列）。
- annotation Part2 拼 mature.fa：`parts[0]`、`parts[14]`、`parts[17]`（0-based）= 第 1/15/18 列 → 名/成熟名/成熟序列。
- `bin/basic2bed.py` 抽 `5,7,8,1,4,6`（score 列取 species 而非 MIR_ID），与 tfbs.py 的内联 awk **只差 score 列**，属两套并存实现。

### 10.2 预测文件 `*_filter_P_prediction`（8 列，tab）
`chrom / strand / mature_query(name_xN) / pri_id / mature_bg..mature_end / pre_bg..pre_end / mature_seq / pre_seq`

### 10.3 `.pipe` 清单
`<input>\t<out_dir>\t<group_prefix>`，每行一个样本；annotation 的 `-i` 就是从这里第 2 列拼出来的。

### 10.4 `anno.map`
两/三列映射（`原名` → `新家族名`[`值`]），由 `bin/anno_miRNA.py` 生成，AmiGO 无关；注释与命名都用它。

### 10.5 count / exp 矩阵
首列 miRNA 名（或 `MIR_ID`），其后每样本一列；DE 的 `--case1/--case2` 是**列号**（1-based，禁 1；表头第 2 列是 `MIR_ID` 时禁 2）。

---

## 11. 辅助脚本清单（scripts/ 与 bin/）

### 11.1 `src/utils/`
- `config.py`：`load_config(path) -> Dict[str,str]`，格式 `key: value`，`#` 注释，缺 `:` 抛 `ValueError("Invalid config line {n}: missing ':' separator")`。
- `dependencies.py`：`check_external_tools(Dict[str, Optional[str]]) -> (all_ok, missing)`，先 `shutil.which` 或校验可执行位，再依次试 `--version` / `-version` / 空参，`timeout=20`，逐行打 ✔/✘。
- `shellquote.py`：`shq(value) = shlex.quote(str(value))`；只对路径/值用，别包 shell 语法。
- `help_display.py`：`print_combined_step_help(step_str)` **仍为占位**。

### 11.2 `scripts/`（按用途）
| 脚本 | 职责 |
|---|---|
| `basic_stat.R` | 基础统计（长度分布、基因组分布等） |
| `classify_scores.R` | core 评估分类 category 1 / 2.1 / 2.2 / 2.3 |
| `correlation_miRNA.R` / `pca_miRNA.R` | 样本相关 / PCA |
| `enrich_analysis.R` | 富集主脚本（返回 0/1/2/3，见 §6.5） |
| `miRNA_DEG.R` / `miRNA_dotplot.R` / `miRNA_expression.R` | 差异表达 / 点图 / 表达条形图 |
| `miRNA_chord_type2.R` | 和弦图（target 路线专用） |
| `tfbs_report.R` | TFBS 汇总图 |
| `build_orgdb.R` | 用 emapper 注释 + `ko00001.json` 建 `org.Morg.eg.db` |
| `MIR_cluster.py` | 两列关系做连通分量聚类，输出 miRNA / CL1.. / 簇内连通率 |
| `assign_family_by_score.py` | 按最佳比对分数定家族；MIR156/157、MIR165/166、MIR170/171 三对按序列特征 |
| `blast-score-filter.py` | 每 query 最高分 + 阈值过滤（`-t 55 -d 20`） |
| `blast-score-stat.py` | 每 query 只留最高分，输出 query/家族/分数 |
| `bowtie_score_filter.py` | 三项分数排序取最高，全平局全留 |
| `bowtie_scoring_blastn.py` | BLASTn-short 打分：左 overhang 单边 gap=1.0、内部错配=1.0、右 overhang 单边 gap=0.5、参考第 10/11 位差异 +1、长度差 ×1.0；`总分 = 100 - 10*(mismatch+position+length)` |
| `dedup_seq.py` | 全局序列去重（首个→`--uni`，后续→`--dup`，忽略 header、大小写不敏感） |
| `extract_col.py` | 按 case1/case2 列号抽列（第 1 列 + case1 全部 + case2 全部） |
| `get_upstream.py` | 读 `.fai`，按 BED + 上游长度取启动子区 |
| `miRNA_enrich_analysis.py` | 逐 miRNA 建子目录跑 `enrich_analysis.R`，聚合返回码 |
| `parse_ssearch.py` | ssearch36 默认输出 → tsv；docstring 记录了复现 psRNATarget 的参数 `-f -8 -g -3 -E 10000 -T 8 -b 200 -r +4/-3 -n -U -W 10 -N 20000` |
| `RNA_plot.py` | TSV（序列+结构）→ `RNAplot` 出 SVG |
| `split_fasta_gh_v2.py` | 按家族大小拆 host/guest（≤5 取 1、6–10 取 2、>10 取 3） |
| `reanno_fasta_by_map_v2.py` / `v3.py` | 按两/三列 mapping 重命名 FASTA |
| `reanno_remap_v3.py` / `v5.py` | 合并 cluster/singleton/alignment 为三列 miRNA→家族映射；v5 优先复用内在家族号，v3 支持 `-f` 已知家族与冲突消解 |

### 11.3 `bin/`（PmiREN core 构建工具箱）
| 脚本 | 职责 |
|---|---|
| `anno_miRNA.py` | PmiREN 注释主流程：建 BLAST 库（或用现成索引）→ 保留 ≥13 bp → 打分/过滤/定家族 → 未命中做自比对聚类 → 合并映射并重命名，产物 `anno.fasta` + `anno.map` |
| `basic2bed.py` | 抽 `5,7,8,1,4,6` 列（basic-info → BED） |
| `build_orgdb.R` | 与 `scripts/build_orgdb.R` 同内容（6435 B） |
| `dedup_family.py` | 按家族（只取数字部分，`miR156a`/`miR156d*` 合并为 `miR156`）内部去重并标 `isoform1/2…` |
| `mirdp3_core_build.py` | PmiREN-core V2 构建：比对→打分→聚类→统一重命名→拆 host/guest→与参考库合并 |
| `remove_matched_seq.py` | 从 `-i` 删掉与 `-r` 完全同长同序列的记录 |

### 11.4 core 数据集文档口径（`docs/core-dataset.md`）
```bash
python scripts/mirdp3_core_build.py -p data/PmiREN-20260810-isoform.fa \
  -i examples/PmiREN2.0_basic_info_MIR_unqiue.fasta -s 12178 -t 14 -o test/PmiREN-v2
# 评估：makeblastdb → blastn -task blastn-short（outfmt 含 qseq sseq）
#       → awk '$4>=13' → bin/bowtie_scoring_blastn.py → scripts/dedup_score_pair.py
#       → Rscript scripts/classify_scores.R
```
> `scripts/dedup_score_pair.py` 与文档里的 `bin/extract_columns.py` 当前**不在仓库**（文档滞后，见 §14-3）。

---

## 12. 依赖与环境

### 12.1 推荐建环境
```bash
conda env create -f mirdp3_environment.yml -n mirdp3
```
宽松约束要点：`python=3.12.*`、**`bowtie=1.3.1`（注释明确 "do NOT use bowtie2 here"）**、`bowtie2>=2.5`（仅个别环节）、`samtools>=1.17`、`bedtools>=2.30`、`seqkit>=2.4`、**`viennarna=2.7.*`（注释 "RNAfold; output format sensitive"）**、`trim-galore>=0.6.7`、`blast>=2.13`、`hmmer>=3.3`、`fasta3>=36.3.8`（ssearch）、`mmseqs2>=14`、`diamond>=2.0`、`eggnog-mapper>=2.1`、`meme>=5.4`（FIMO）、`r-base>=4.2`。

### 12.2 已锁定的实测版本（`environment.yml`）
`python=3.12.13`、`r-base=4.4.1`、`bowtie=1.3.1`、`bowtie2=2.5.5`、`viennarna=2.7.2`、`trim-galore=2.2.0`、`seqkit=2.13.0`、`meme=5.5.9`、`diamond=2.2.0`、`fasta3=36.3.8`、`eggnog-mapper=2.1.12`、`hmmer=3.4`、`bedtools=2.31.1`、`mmseqs2=18.8cc5c`、`samtools=1.23.1`。

### 12.3 三种安装
1. `conda install -c jaguares mirdeep-p3`（自带 `data/index`）。
2. 源码：`conda env create -f mirdp3_environment.yml -n mirdp3` + 下载 `data-index.tar.gz` 解到 `data/`。
3. Docker：`merc3dez/mirdeep-p3:3.1.6c-full`（另有阿里云镜像）。

### 12.4 数据索引
- 版本号唯一真源：`data-index.env` 的 `DATA_INDEX_RELEASE=mirdeep-p3-v3.1.4c-full`。
- 内容：`rfam_index`、`mature_index`、`isoform-in`（PmiREN blast 索引）。
- 该文件里明确警告：**不要删对应 GitHub Release/tag**，因为 `ci.yml`、`docker-release.yml`、`test.sh`、`README.md`、`CONTRIBUTING.md`、`conda.recipe/meta.yaml`、bioconda recipe 都读它。

---

## 13. 测试与验证

```bash
./test.sh [-t threads=8] [-q] [-k]     # -k 保留 tests/test_output
```
- 退出码：`0` 全通过、`1` 依赖缺失、`2` 数据下载失败、`3` identification 失败、`4` annotation 失败。
- 依赖检查：`bowtie`、`RNAfold`、`samtools`、`bedtools`、`seqkit`、`blastn`。
- 数据：`tests/data/mini.fq` + `mini_genome.fasta`。
- 校验点：`<out>/mini/mini-processed.fa`、`mini_filter_P_prediction`；annotation 侧 `annotation/mirdp3-test/mirdp3-test-basic-info`。
- ★ 3.1.7b 起还覆盖 Target_finder：默认引擎断言 120 行且 NF==6、连跑两次 `cmp` 一致、`--no-allow-gap` 与默认一致、有 `ssearch36` 时 `--engine ssearch` ≥100 行；失败退出码 5。
- ★ 3.1.7c 起 `test.sh` 依赖检查前会打印 `check_python_deps.py` 的 core/optional 状态（core 无三方依赖，numpy 供 psrna 引擎）。
- ⚠️ test.sh 用 `/tmp` 中转（`DATA_INDEX_RELEASE` 拼的 URL 下载到 /tmp），本机执行需改到工作目录。
- annotation 调用示例（README 口径）：
  ```bash
  ./mirdeep-p3 annotation -i <ID_OUT>/mini -d <ID_OUT>/index/genome_index \
    --prefix_miRNA Ath --prefix mirdp3-test --species "Arabidopsis thaliana" -o <OUT>
  ```

---

## 14. 已知问题 / 技术债（改之前先看）

| # | 问题 | 位置 | 现状 |
|---|---|---|---|
| 1 | `--consistency` / `--common` 调用不存在的 `src/miRNA_consistency.py` | `src/commands/annotation.py:582,590,600` | ★ 3.1.7c 已修（脚本补回，两种模式实测通过） |
| 2 | `--clean` 帮助说删 "temp, index"，identification 实际只删 temp | `identification.py` | ★ 3.1.7c 已修（文案改正 + 新增 `--clean-index`） |
| 3 | 文档引用的 `bin/extract_columns.py`、`scripts/dedup_score_pair.py` 缺失 | `docs/*.md` | 部分修：3.1.7c 补 `dedup_score_pair.py`；`extract_columns.py` 仍缺（I-09） |
| 4 | analysis 子命令硬编码 `Rscript/python/fimo/ssearch36/seqkit`，不检查、不可配 | `src/commands/analysis/*` | ★ 3.1.7c 已修（`src/utils/toolset.py`） |
| 5 | 仓库提交了 `__pycache__/*.pyc`（5 个 `__pycache__` 目录、多个 Python 版本） | 多处 | 仍存在（`__pycache__` 未入 `.gitignore`）；本地提交已排除 |
| 6 | `check_python_deps.py` 只列标准库 | `src/check_python_deps.py` | ★ 3.1.7c 已修（numpy 标为 optional） |
| 7 | `--step` 的 `del cli_args[step_idx-1]` 索引脆弱 | `mirdeep-p3` `main()` | ★ 3.1.7c 已修（按索引过滤） |
| 8 | `--step 3` 单用只打印帮助 | `mirdeep-p3` | 仍存在（I-06） |
| 9 | Target_finder temp 未清理（`rmtree` 被注释） | `src/commands/analysis/target_finder.py` | ★ 3.1.7c 已修（默认清理 + `--keep-temp`） |
| 10 | README 的 Citation/Contact/Acknowledgements 仍是占位符 | `README.md` | 仍存在（I-11） |
| 11 | 本文基线是"未跑过真实数据"的静态阅读结论；`data/index` 缺失 | — | 已跑通：`./test.sh -t 12 -k` ALL TESTS PASSED（3.1.7c） |
| 12 | basic-info 先写 24 列再重排成 22 列，重排逻辑只在 `rearrange_result.py` 的 docstring 里 | `src/rearrange_result.py` | 待记录（I-22）；列号含义见 §10.1 |

---

## 15. 常见改动的落点（改哪一行）

| 想改什么 | 改哪里 |
|---|---|
| 成熟体长度范围（18–26） | `identification.py` 的 `DEFAULT_MIN_LEN/DEFAULT_MAX_LEN`（同时影响 trim_galore `--length/--max_length`） |
| 植物成熟体 20–24 nt 的硬判据 | `mod_rm_redundant_meet_plant.py: apply_plant_criteria()` |
| 错配数、star 比例、RPM 门槛 | 同上（0.75、≤5、\|差\|≤3、rpm≥20） |
| miRDeep 打分常量 / Gumbel 参数 | `mod_miRDP.py` 的 `miRDPPredictor` 类常量区 |
| 预测阈值（默认 1.0） | `mod_miRDP.py -v`；但目前 identification **不透传**，需在 `identification.py:818` 加参数 |
| 保守性打分（known miRNAs） | 同上，加 `-s <mature.fa>` |
| read 级过滤强度（全长精确匹配、位点数上限） | `filter_alignments.py`（`identity==1.0`、`trunc_*`、`-c`） |
| 前体窗口长度（300） | `DEFAULT_PRE_LENGTH`（> `excise_candidate` 的 `-l`） |
| rpm_threshold（去 ncRNA） | `DEFAULT_RPM_THRESHOLD=5` |
| 新增 identification 参数 | `build_parser()` 加 `add_argument` → 在 `run_pipeline_for_input()` 对应 stage 拼命令（记得 `shq()`） |
| 新增 analysis 子命令 | `analysis/__init__.py` 的 `add_subparsers` 加一个 → 新建 `analysis/<name>.py` 暴露 `main(argv, project_root)` → 在 `if/elif` 分发 |
| 改 basic-info 列顺序 | `primary-basic-info.py` `out_fields` **+ 同步所有下游取列**（§10.1 清单） |
| 加数据索引文件 | 改 `data-index.env` 指向新 Release，并同步 `test.sh`/CI/conda recipe |
| 外部工具默认路径 / 新增可覆盖工具 | `src/utils/toolset.py` 的 `TOOL_SPECS`（key → 默认可执行名、CLI flag、用途）；子命令里 `resolve_tools(args, required=[...], optional=[...])` |
| 新增 analysis 参数别名（如 `--keep-temp`） | 对应子命令 `build_parser()`；`Onestep` 记得在 `forward_tool_args` 或各自 cmd 里透传 |
| Target_finder 引擎行为 | `src/commands/analysis/target_finder.py` 的 `DEFAULT_ENGINE`、`_run_psrna_engine` / `_run_ssearch_engine`；打分规则在 `src/MiRNATarget/psrna_scan.py` |

---

## 16. 与 miRDeep-P2 的关系（背景，避免误判"bug"）

- miRDeep-P2 是 Perl；本仓库把它整体 Python 化，因此 `mod_miRDP.py` 里许多"看起来奇怪"的写法（Gumbel CDF 的指数符号、`score_star_not` 等魔法数、`nucleus_length=7`）都是**原版移植**，不是笔误。
- 植物化差异主要在：负链 mature 不支持、加 `apply_plant_criteria()` 硬判据、默认关掉 randfold/drosha/conservation。
- 因此："某个 miRNA 没被检出"时，排查顺序应是：① 是否负链；② 是否被植物判据（长度/错配/star 比例）挡掉；③ 是否被 read 级全长精确匹配挡掉；④ 打分为负；⑤ 输入本身被 trim/RPM 过滤。

---

*文档生成时间：2026-11-15；基准 commit `fb0b754`（v3.1.6c）。*
