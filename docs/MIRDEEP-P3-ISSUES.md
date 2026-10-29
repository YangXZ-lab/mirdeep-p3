# miRDeep-P3 问题与改进台账（ISSUES & BACKLOG）

> **这份文件只记录，不动代码。** 等赵哥点名再做。
> 基准：commit `fb0b754` / `v3.1.6c`，记录日期 2026-11-15。
> 配套文档：`docs/MIRDEEP-P3-INTERNALS.md`（实现细节知识库，575 行）。
>
> **用法**：每条都有编号（如 `I-03`）。要动哪条直接说编号，例如"改 I-03"或"把 A 组全修了"。
> 每条都写了 **证据 / 影响 / 建议改法 / 验证方法**，动手时不必重新调研。
> 改完请在 §5 变更日志追加一行，并把该条状态改成"已修（commit）"。

---

## 1. 索引表

严重度：**P0** = 功能直接不可用；**P1** = 结果/体验受损但不致命；**P2** = 工程卫生/一致性；**P3** = 设计改进机会（不修也能用）。

| 编号 | 组 | 严重度 | 位置 | 一句话问题 | 状态 |
|---|---|---|---|---|---|
| I-01 | A 缺陷 | P0 | `src/commands/annotation.py:582,590,600` | `--consistency` / `--common` 调用的 `src/miRNA_consistency.py` 不存在 | 已修 (3.1.7c) |
| I-02 | A 缺陷 | P1 | `src/commands/analysis/target_finder.py:146` | temp 目录不清理（`rmtree` 被注释掉） | 已修 (3.1.7c) |
| I-03 | A 缺陷 | P1 | `src/commands/analysis/*.py` | 六个分析子命令硬编码 `Rscript/python/fimo/ssearch36/seqkit`，不检查、不可配 | 已修 (3.1.7c) |
| I-04 | A 缺陷 | P1 | `src/check_python_deps.py` | 依赖检查只列标准库，实际什么也没检查 | 已修 (3.1.7c) |
| I-05 | A 缺陷 | P2 | `mirdeep-p3` `main()` | `del cli_args[step_idx - 1]` 连删两次，索引依赖"脚本名已被移除"的假设 | 已修 (3.1.7c) |
| I-06 | A 缺陷 | P2 | `mirdeep-p3` | `--step 3` 单独使用时只打印帮助，与 `--step 1/2` 行为不一致 | 待办 |
| I-07 | A 缺陷 | P2 | `identification.py` `--clean` | 帮助文案说删 "temp, index"，实际只删 temp | 已修 (3.1.7c) |
| I-08 | A 缺陷 | P2 | `src/utils/help_display.py` | `print_combined_step_help()` 是占位实现 | 待办 |
| I-09 | B 文档 | P1 | `docs/file_process.md` | 引用 `bin/extract_columns.py`，仓库里只有 `bin/basic2bed.py` | 待办 |
| I-10 | B 文档 | P1 | `docs/core-dataset.md` | 流程里用到 `scripts/dedup_score_pair.py`，仓库缺失 | 已修 (3.1.7c) |
| I-11 | B 文档 | P2 | `README.md` | Citation / Contact / Acknowledgements 全是占位符 | 待办 |
| I-12 | C 工程 | P2 | 全仓库 | 提交了 5 个 `__pycache__` 目录的 `.pyc`，无 `.gitignore` 约束 | 待办 |
| I-13 | C 工程 | P2 | git 索引 | 入口脚本本地被 `chmod +x`（100644→100755）但未提交 | 已修 (3.1.7c) |
| I-14 | C 工程 | P2 | `test.sh` | 用 `/tmp` 中转下载，与"禁用 /tmp"的环境约定冲突 | 待办 |
| I-15 | C 环境 | P1 | `data/index/` | 本地索引缺失，首次运行会触发下载或现场建索引 | 待办 |
| I-16 | D 设计 | P1 | `identification.py:818-821` | 调 `mod_miRDP.py` **不传** `-v/-s/-x/-y/-z`，用户无法调灵敏度、保守性永不启用 | 待记录 |
| I-17 | D 设计 | P2 | `mod_rm_redundant_meet_plant.py` | 植物判据全部硬编码，无法按物种放宽 | 待记录 |
| I-18 | D 设计 | P2 | 多处 | basic-info 列号以裸 awk 索引散落在 5 处以上，无单一 schema 来源 | 待记录 |
| I-19 | D 设计 | P3 | `annotation.py` | 注释产物路径 `temp/temp/anno.map`（双层 temp） | 待记录 |
| I-20 | D 设计 | P3 | `bin/basic2bed.py` vs `analysis/tfbs.py` | 两套 basic-info→BED，第 5 列分别取 4 与 2 | 待记录 |
| I-21 | D 设计 | P3 | `mod_miRDP.py` `_determine_arm()` | 负链 mature 直接返回空 → 负链 miRNA 结构上不可能被检出 | 待记录 |
| I-22 | D 设计 | P3 | `src/rearrange_result.py` | basic-info 先写 24 列再重排成 22 列，重排顺序只在 docstring 里，易把两套列号搞混 | 待记录 |
| I-23 | A 缺陷 | P1 | `src/commands/analysis/__init__.py:76-86` | 子命令别名 `DE` / `FA` 落到 else 分支，只打印帮助并 exit 1 | 已修 (3.1.8a) |
| I-24 | C 工程 | P3 | `identification.py` 批次 manifest | 批次 manifest 文件名精度到分钟，同一分钟内重跑互相覆盖 | 待办 |

---

## 2. 逐条详情

> **3.1.7c 已修**：I-01、I-02、I-03、I-04、I-05、I-07、I-10、I-13。逐条证据见 §5 变更日志；未列出的条目状态不变。

### I-01（P0）`--consistency` / `--common` 指向不存在的脚本
- **证据**：`src/commands/annotation.py:582,590,600` 三处 `subprocess.run([sys.executable, str(project_root / "src" / "miRNA_consistency.py"), ...])`；`ls src/` 无此文件。
- **影响**：这两个参数一旦使用必然 `FileNotFoundError`，且报错指向一个"看起来应该存在"的路径，容易误判成环境问题。README 还专门描述了它们的用途。
- **待确认**：是改名（并入 `bin/anno_miRNA.py` 或 `scripts/reanno_remap_v*.py`？）、漏提交、还是功能下线但参数没删。
- **建议改法**（按确认结果三选一）：
  1. 找到真实文件 → 改路径 + 补端到端跑通；
  2. 漏提交 → 从另外两份副本 `../miRDeep-P3/`、`../mirdp3/` 里找回来；
  3. 功能下线 → 参数保留但改成 `parser.error("--consistency is not available in this build")`，README 标注。
- **验证方法**：带 `--consistency` 跑一次完整注释（需两组、其中一组有历史 basic-info）。

### I-02（P1）Target_finder 不清理 temp
- **证据**：`src/commands/analysis/target_finder.py:146` 为 `# shutil.rmtree(temp_dir)`。
- **影响**：每次跑都留全 CDS 比对中间产物（可能几百 MB）。
- **建议改法**：`shutil.rmtree(temp_dir, ignore_errors=True)`，与 DE/Stat/TFBS 对齐；如担心排查，加 `--keep-temp`。
- **验证方法**：跑一次后 `ls <out>` 只剩最终产物。

### I-03（P1）分析子命令的外部工具不可配、不检查
- **证据**：`tfbs.py` 直接拼 `fimo`/`bedtools`；`target_finder.py` 拼 `ssearch36`/`seqkit`；`stat.py`/`differential_expression.py`/`functional_analysis.py` 拼 `Rscript` 与 `python`。全仓库 `check_external_tools` 只在 identification/annotation 里用。
- **影响**：环境没 activate / PATH 不同 / 想指定特定版本时，报错是裸 "command not found"，且跨版本复现困难。
- **建议改法**：复用 `src/utils/dependencies.py: check_external_tools()` 做前置检查；把 identification 的 `--bowtie/--RNAfold/--samtools` 模式推广为 `--rscript/--fimo/--ssearch36/--seqkit/--bedtools/--emapper`。
- **验证方法**：清空 PATH 跑一次，应提示"缺少 X，请用 --X 指定"而非 traceback。

### I-04（P1）`check_python_deps.py` 形同虚设
- **证据**：`REQUIRED_PACKAGES` 全是标准库（sys/os/gzip/argparse/typing/collections/re/pathlib/importlib/subprocess/tempfile/traceback）。
- **影响**：让人误以为依赖已被检查。
- **建议改法**：当前确实零第三方依赖 —— 要么删掉模块与调用点（最诚实），要么真检查并声明"只用标准库"。
- **验证方法**：输出与实际依赖一致。

### I-05（P2）`--step` 参数删除的索引脆弱
- **证据**：`main()` 里两次 `del cli_args[step_idx - 1]`，只有在"`sys.argv` 已去掉脚本名"假设成立时索引才对。
- **影响**：`--step` 在 argv 首位、或与 `--config` 混排时可能删错元素，表现为参数消失或 "unrecognized arguments"，难复现。
- **建议改法**：改成构造新列表（跳过两个下标）或 `parse_known_args` + 显式剥离。
- **验证方法**：单测三种写法：`--step 1 -i x`、`identification --step 1 -i x`、`--config c.txt --step 1 -i x`。

### I-06（P2）`--step 3` 只打印帮助
- **证据**：`main()` 中 `--step` 单选到 "3"/"analysis" 时只调 `analysis.main(["-h"])`。
- **影响**：用户以为在跑分析，实际什么都没跑，且退出码 0。
- **建议改法**：改为 `parser.error("--step 3 requires an analysis subcommand; see Onestep")`。
- **验证方法**：`./mirdeep-p3 --step 3` 不静默 exit 0。

### I-07（P2）`--clean` 文案与实现不符
- **证据**：帮助写删 "temp, index"；实现只删各样本 `temp/`。
- **影响**：以为清了索引，实际 `<out>/index` 还在（bowtie 索引可能几 GB）。
- **建议改法**：改文案为"只删 temp"，或加 `--clean-index` 单独控制。
- **验证方法**：`--clean` 后 `ls <out>/index` 与文档一致。

### I-08（P2）合并帮助仍是占位
- **证据**：`print_combined_step_help()` 输出 "(Full combined help will be implemented later.)"。
- **影响**：`-h --step 1,3` 这类组合拿不到有效帮助（1+2 已有手写版 `show_combined_identification_annotation_help()`）。
- **建议改法**：把 1+2 手写帮助泛化为"任意 step 组合"，复用 `[ID]/[AN]/[BOTH]` 标注。
- **验证方法**：`./mirdeep-p3 -h --step 1,3` 输出含 Stat/TFBS 参数。

### I-09 / I-10（P1）文档引用不存在的脚本
- **证据**：`docs/file_process.md` 写 `bin/extract_columns.py`（实际是 `bin/basic2bed.py`，抽 5,7,8,1,4,6 列）；`docs/core-dataset.md` 流程含 `scripts/dedup_score_pair.py`（仓库无）。
- **影响**：照文档复现 core 构建会在中间失败。
- **建议改法**：改名对齐（I-09）；找回脚本或改写流程（I-10）。若 `dedup_score_pair` 确实丢了，可用 `blast-score-stat.py` + `dedup_seq.py` 近似替代，但需重验分类结果与 V2 一致。
- **验证方法**：严格照文档从头跑一遍 core 构建。

### I-11（P2）README 占位符
- **证据**：Citation / Contact / Acknowledgements 仍为 `[Project Name]`、`[Authors]`、`[email]`、`[Funding]`。
- **影响**：发版/投稿时引用信息缺失。
- **建议改法**：补作者、单位、通讯邮箱、基金号；Citation 段给 Software 引用格式与 Zenodo DOI（若已归档）。
- **验证方法**：人工核对。

### I-12（P2）`__pycache__` 入库
- **证据**：`find . -name "__pycache__" -type d | wc -l` → 5，`.pyc` 有多种 Python 版本后缀。
- **影响**：仓库脏、diff 噪音、跨版本误 import 风险。
- **建议改法**：`git rm -r --cached` + 加 `.gitignore`（`__pycache__/`、`*.pyc`）。
- **验证方法**：`git status` 干净，新克隆运行不产生被跟踪文件。

### I-13（P2）入口脚本执行位未提交
- **证据**：`git status` 显示 `M mirdeep-p3`，`git diff` 只有 mode 100644 → 100755。
- **影响**：新克隆者必须 `chmod +x`，README 里所有 `./mirdeep-p3` 都会 Permission denied。
- **建议改法**：`git update-index --chmod=+x mirdeep-p3` 提交 mode 变更。
- **验证方法**：`git ls-files -s mirdeep-p3` 显示 `100755`。

### I-14（P2）test.sh 用 /tmp
- **证据**：`test.sh` 下载 `data-index.tar.gz` 与解压中转走 `/tmp`。
- **影响**：与本机"禁止 /tmp"约定冲突；/tmp 空间或清理策略会致失败。
- **建议改法**：引入 `TMPDIR`（未设置时回落 /tmp）并允许 `-T <dir>` 覆盖；SLURM 里统一 export。
- **验证方法**：`TMPDIR=/nfs/users/zhaojw/Genburten/temp ./test.sh` 不往 /tmp 写。

### I-15（P1）data/index 缺失
- **证据**：`ls data/index` → No such file or directory；`data-index.env` 指向 Release `mirdeep-p3-v3.1.4c-full`。
- **影响**：本地任何 identification/annotation 首次运行都会卡在索引。
- **建议改法**：无需改代码 —— 按 README 下载 `data-index.tar.gz` 解到 `data/`，或直接用 conda/Docker 版；离线场景把解压后的 `data/index` 放共享路径，运行时用 `-d` 指过去。
- **验证方法**：`ls data/index` 含 `rfam_index*`、`mature_index*`、`isoform-in*`。

### I-16（P1）identification 不透传核心算法参数
- **证据**：`src/commands/identification.py:818-821` 的调用只有 `{signatures} {precursors.struc} -o {predictions}`。
- **影响**：① 阈值恒 1.0，无法用"降阈值提召回"这条最常规调参路径；② 保守性打分永远关闭（`known_nuclei` 空），PmiREN 已知 miRNA 的先验完全没用上；③ 无法启用 randfold/drosha 做对照实验。
- **建议改法**：`build_parser()` 增加 `--threshold/-v`（默认 1.0）、`--known-mirnas`、`--randfold`、`--drosha`、`--sensitive` 并透传到 stage 11；默认值不变以保证向后兼容。
- **风险**：开启 `-s` 会改变所有预测分数 → 与历史结果不可比；建议把"是否启用保守性打分"写进日志头。
- **验证方法**：同一 mini 数据跑 1.0 vs 0.5，条数单调增加；加 `-s` 后分数分布整体上移。

### I-17（P2）植物判据不可配
- **证据**：`apply_plant_criteria()` 内联常量：错配 ≤5/≤5、差 ≤3、`(mature+star) >= 0.75*total`、长度 20–24、23–24 nt 追加 `rpm >= 20` 且必须有 star read。
- **影响**：部分物种（成熟体 19 或 25 nt）被系统性漏检，且日志里没有分项淘汰统计，用户无从判断是哪条挡的。
- **建议改法**：抽成 `DEFAULT_PLANT_CRITERIA` dict + CLI 参数（`--plant-min-len/--plant-max-len/--mismatch-max/--star-frac`），日志输出各判据淘汰条数。
- **验证方法**：`--plant-min-len 18` 后原来因长度被挡的条目应出现；各判据淘汰数之和等于输入条数。

### I-18（P2）basic-info 列号散落
- **证据**：`tfbs.py:104` 用 `$5,$7,$8,$1,$2,$6`；`target_finder.py:67` 用 `$1,$18`；`annotation.py:372` 用 `parts[0]/[14]/[17]`（0-based）；`bin/basic2bed.py` 用 `5,7,8,1,4,6`。
- **影响**：改列顺序时漏改一处会产生**静默错误结果**（比报错更危险）。
- **建议改法**：建 `src/basic_info_schema.py` 定义列名与索引常量，Python 侧统一引用；shell 侧 awk 用 Python 生成的 `-v` 变量传列号。
- **验证方法**：改一次 schema 顺序，全链路仍跑通（同时也是回归测试）。

### I-19（P3）双层 temp 路径
- **证据**：`annotation.py` 传 `-o temp`，随后从 `temp/temp/anno.map` 取产物。
- **影响**：`bin/anno_miRNA.py` 内部再建一层 temp，路径耦合，换 `-o` 易找不到文件。
- **建议改法**：显式 `anno_dir` 变量；或让 `anno_miRNA.py` 支持 `--map-out` 直接指定映射路径。
- **验证方法**：`--clean` 后无残留；换输出目录仍能找到 `anno.map`。

### I-20（P3）两套 basic-info→BED
- **证据**：`bin/basic2bed.py` 抽 `5,7,8,1,4,6`（score=species）；`analysis/tfbs.py:104` 内联 awk 抽 `$5,$7,$8,$1,$2,$6`（score=MIR_ID）。
- **影响**：两路产出第 5 列语义不同，混用会让下游按 score 排序/过滤行为不一致。
- **建议改法**：统一为 `tfbs.py` 版本（MIR_ID 作 score 更有意义），`basic2bed.py` 改为调用同一实现或废弃。
- **验证方法**：两路对同一 basic-info 产出逐行相同的 BED。

### I-21（P3）负链 mature 不支持
- **证据**：`mod_miRDP.py: _determine_arm()` 中 mature 为负链时返回 `""`。
- **影响**：负链编码 miRNA 结构上不可能被检出（植物罕见但非不存在；也妨碍用动物数据做对照）。
- **建议改法**：需把 mature 区间镜像到前体坐标系（参考 `mod_rm_redundant_meet_plant.py` 的负链换算）并补回归用例。**风险高，非必要不做。**
- **验证方法**：构造已知负链 miRNA 前体 + reads，应被检出且 star/mature 不重叠。

### I-22（P3）basic-info 的 24 列 → 22 列重排未文档化
- **证据**：`src/primary-basic-info.py:313,412` 的 `out_fields` 写 24 列；`src/rearrange_result.py: reorder_line()` 用 `order = [0,1,2,3,4,7,5,6,8,9,12,13,10,11,14,15,16,17,18,20,21,19]` 重排，并丢弃 `confidence`/`source`，最终文件 22 列。
- **影响**：对着 `primary-basic-info.py` 数第 6 列会得到 `ext_start`，而真实文件第 6 列是 `strand`；改下游 awk 列号时极易错位（TFBS 的 `$5,$7,$8,$1,$2,$6` 只有按 22 列才成立）。
- **建议改法**：把 `order` 提升为共享常量（如 `src/utils/basic_info_schema.py`），并在 `docs/MIRDEEP-P3-INTERNALS.md §10.1` 保留双表；下游脚本改为引用常量或列名访问。
- **验证方法**：用 `tests/test_output/annotation/mirdp3-test/mirdp3-test-basic-info` 做列名断言。

### I-23（P1）analysis 子命令别名 DE / FA 完全不可用 —— 已修 (3.1.8a)
- **证据**：`src/commands/analysis/__init__.py:76-86` 原来逐个比较 `args.analysis_command == "Differential_expression"`；但 argparse 经别名调用时 dest 值是 `"DE"`，于是落到 `else: parser.print_help(); sys.exit(1)`。实测 `mirdeep-p3 analysis DE ...` 只打印帮助、返回码 1、stderr 为空。
- **影响**：README 把 `DE` 写成常用简写，照文档敲命令的人拿不到结果，也看不到报错。
- **修法**：parse 后先归一化 `command = {"DE": "Differential_expression", "FA": "Functional_analysis"}.get(args.analysis_command, args.analysis_command)`，六个分支改比 `command`。
- **验证方法**：`test.sh` 第 10 段故意用 `analysis DE` 调用，通过即覆盖。

### I-24（P3）批次 manifest 文件名精度到分钟 —— 待办
- **证据**：`identification.py: run()` 用 `datetime.now().strftime("%m%d%Y-%H%M")` 生成时间戳（与 pipe 文件同源），批次 manifest 复用之。
- **影响**：同一分钟内重跑会覆盖上一条批次记录（每样本 `*.run.json` 不受影响，信息不丢）。
- **建议改法**：时间戳加秒（`%m%d%Y-%H%M%S`），或存在同名时追加 `-2`。
- **验证方法**：连续两次 `--resume` 后 `ls <out>/*.manifest.json` 应有两条。

---

## 3. 设计改进机会（不影响正确性）

| 编号 | 想法 | 收益 | 代价 |
|---|---|---|---|
| O-01 | identification 的 15 个 stage 抽成声明式 Stage 列表（名称/命令/产物/空判定） | 加 stage、改顺序、断点续跑都变简单，日志可按 stage 打点 | 中等重构，需重跑 mini 验证 |
| O-02 | `--resume` / 断点续跑（按 stage 产物存在性跳过） | 大样本失败重跑不用从头 | 需要产物完整性校验 |
| O-03 | 机器可读 run manifest（JSON：版本/参数/工具版本/输入 md5/各 stage 状态） | 可复现性、排查、批量统计 | 小 |
| O-04 | analysis 统一"工具解析器" + 前置检查（与 I-03 同源） | 报错友好 | 小 |
| O-05 | 坐标换算函数补单测（正/负链 × first/second arm 四组合） | 3.1.6c 那类 star 坐标 bug 可被 CI 抓住 | 小，收益高 |
| O-06 | `--step` 组合帮助泛化（见 I-08） | 体验 | 小 |
| O-07 | 植物判据淘汰原因做成统计报告 | 调参有依据，可回答"为什么没检出" | 小 |
| O-08 | `test.sh` 增加 analysis 冒烟（Stat + DE 最小数据） | 覆盖下游 | 中（需准备最小 count/exp） |
| O-09 | `pyproject.toml` + `console_scripts` 安装入口 | 摆脱"无扩展名脚本 + chmod"的坑（I-13） | 中，会动安装方式 |
| O-10 | 启动时比对 `DATA_INDEX_RELEASE` 与已装索引版本 | 防"软件新、索引旧"的静默不一致 | 小 |

**3.1.8a 已落地**：O-01、O-02、O-03、O-08（见 §5 变更日志）。其余仍开放。

---

## 4. 运行环境备忘（不是 bug，但会绊人）

1. **三份仓库副本**：`/nfs/users/zhaojw/software/github/git/mirdeep-p3`（最新）、`../miRDeep-P3`、`../mirdp3`。三者独立，`PROJECT_ROOT` 各取自己的根 —— 改之前必须确认运行命令指向哪份。
2. **`data/index` 未下载**（I-15）。
3. **`test.sh` 走 /tmp**（I-14）。
4. **R 库并发**：`enrich_analysis.R` 会 `install.packages(..., repos=NULL, type="source")` 且**不传 `lib=`**；`Functional_analysis` 已用 `MIRDEEP_R_LIB` 私有库规避，但**手工**直接调 `enrich_analysis.R` 仍会争写共享 conda 库 —— 批量跑必须自己设 `R_LIBS_USER`。
5. **ssearch36 路径长度**：README 提到不接受过长路径（建议 <100 字符），深目录会莫名失败。

---

## 5. 变更日志（改完在这里追加）

| 日期 | 编号 | 改动 | commit | 验证方式 |
|---|---|---|---|---|
| 2026-11-18 | I-01 | 补回 `src/miRNA_consistency.py`（544 行，三种模式），`--consistency`/`--common` 恢复可用 | `ab930f9` | 参考模式与 common 模式各产出 19 行 |
| 2026-11-18 | I-02 | Target_finder 跑完默认删 `temp`，新增 `--keep-temp` | `ab930f9` | 默认无 temp / `--keep-temp` 有 temp，两次输出一致 |
| 2026-11-18 | I-03 | 新增 `src/utils/toolset.py`；六个 analysis 子命令统一解析外部工具并支持 `--rscript/--fimo/--bedtools/--ssearch36/--seqkit/--emapper/--python` | `ab930f9` | 六子命令 `-h` 均含新参数；坏路径提前退出 |
| 2026-11-18 | I-04 | `check_python_deps.py` 改为诚实检查：core 无三方依赖，numpy 标为 psrna 引擎的 optional | `ab930f9` | `test.sh` 打印 `numpy ✔ 2.4.4 (optional ...)` |
| 2026-11-18 | O-01 | identification 重构为 12 个声明式 `Stage`（新 `src/commands/stage_runner.py`），空结果改用 `StageEmpty` 表达 | `770c2bc` | mini 五产物 md5 与 3.1.7c 逐一相同 |
| 2026-11-18 | O-02 | 新增 `--resume`：按 manifest 记录 + 产物存在性跳过；签名变化整轮重跑；强制忽略 `--clean` | `770c2bc` | 全量续跑 12/12 skip；删 bed 只补 4.15；改签名触发全跑 |
| 2026-11-18 | O-03 | 每样本 `<prefix>.run.json` + 批次 `<time>.manifest.json`（版本/参数/工具版本/输入指纹/各 stage 状态与耗时，原子写） | `770c2bc` | `tests/test_output/mini/mini.run.json` 中 12 个 stage 全 status=ok |
| 2026-11-18 | O-08 | `test.sh` 第 10 段 analysis 冒烟（Stat + DE，新 fixture 180 miRNA × 6 样本），失败退出码 6 | `770c2bc` | `./test.sh -t 12 -k` → ALL TESTS PASSED |
| 2026-11-18 | I-23 | 修复 `analysis DE` / `analysis FA` 别名分发（归一化命令名） | `770c2bc` | test.sh 用 `analysis DE` 调用并通过 |
| 2026-11-18 | I-05 | `--step` 参数剥离改为按索引过滤，不再依赖固定偏移 | `ab930f9` | `-h --step 1,2` 与 `--step 1,2 -h` 输出 byte-identical |
| 2026-11-18 | I-07 | `--clean` 文案与实际一致（只删 temp）；新增 `--clean-index`，`-d` 提供的索引跳过并记日志 | `ab930f9` | 帮助正文含 `genome index is kept` |
| 2026-11-18 | I-10 | 补回 `scripts/dedup_score_pair.py`（按 (col1,col2) 保留最高分） | `ab930f9` | 脚本可执行，`docs/core-dataset.md` 流程可复现 |
| 2026-11-18 | I-13 | `mirdeep-p3`、`src/miRNA_consistency.py`、`scripts/dedup_score_pair.py` 以 100755 提交 | `ab930f9` | `git ls-files -s` |
