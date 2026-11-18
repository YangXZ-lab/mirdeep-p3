# 数据索引（data index）版本管理

> 面向维护者。目的：让「软件版本」和「数据索引版本」这两条独立的版本轴始终对得上；升级数据时知道要改哪几处。
> 代码：`src/utils/data_index.py`；真源：`data-index.env`；标记：`data/index/.index-version`。

## 1. 当前版本（2026-11-18）

| 版本轴 | 当前值 | 真源 | 查法 |
|---|---|---|---|
| 软件 | `3.1.8c` | `src/version.py` 的 `__version__` | `./mirdeep-p3 -v` |
| 数据索引（期望） | `mirdeep-p3-v3.1.4c-full` | `data-index.env:22` 的 `DATA_INDEX_RELEASE` | `grep DATA_INDEX_RELEASE data-index.env` |
| 数据索引（本机已装） | `mirdeep-p3-v3.1.4c-full` | `data/index/.index-version` | `cat data/index/.index-version` |

## 2. 为什么是两条轴

- **索引更新慢、软件更新快**：索引里是 Rfam / 已知成熟 miRNA / PmiREN 的比对库（约 592 MB），跟着数据库版本走，几个月一次；软件按功能迭代，几周一次。
- 所以**不要**用软件版本号去推索引版本号。比较只在 data 轴上做：`data/index/.index-version`（已装）vs `data-index.env`（期望）。
- 发布形式也不同：软件是 `v<version>`（源码 / conda）与 `v<version>-full`（Docker / 离线包）；索引单独一个 Release，历史名字是 `mirdeep-p3-v3.1.4c-full`，它**早于** `v<version>-full` 约定，故意与软件版本无关。

## 3. 索引里有什么

| 文件 | 用途 | 消费者 |
|---|---|---|
| `rfam_index.*.ebwt` | ncRNA（Rfam）比对，`bowtie -v 0` | identification 4.6 stage |
| `mature_index.*.ebwt` | 已知成熟 miRNA 比对，`bowtie -v 1` | identification 4.6 stage |
| `isoform-in.*`（BLAST db） | PmiREN 同源注释 | annotation 4.10（`bin/anno_miRNA.py`） |

本机 `data/index` 是从旧副本软链过来的（`data/index -> .../github/miRDeep-P3/data/index`），该副本不是 git 仓库，写标记不会脏任何仓库。

## 4. 标记文件与五种状态

`data/index/.index-version` 只有一行 tag：

```bash
echo mirdeep-p3-v3.1.4c-full > data/index/.index-version
```

- **谁写**：`test.sh:102` 解压索引后自动写；conda 包自带；**手工安装（源码方式）要自己写这一行**。
- **谁读**：`src/utils/data_index.py` 的 `index_status()` / `describe_index()`。
- **五种状态**（`index_status(root)`）：

| 状态 | 含义 | 影响运行？ |
|---|---|---|
| `missing` | `data/index` 不存在（索引没下载） | **会**，比对阶段直接失败 |
| `unrecorded` | 有索引但没写标记 | 不会，只是提示不明确 |
| `unknown` | 无标记且 `data-index.env` 也没有期望值 | 不会 |
| `match` | 已装 tag == `data-index.env` | 不会 |
| `mismatch` | 已装 tag != `data-index.env` | 不会，**只警告不阻断** |

启动时会打印一行：

```text
[index] data index: mirdeep-p3-v3.1.4c-full (matches data-index.env)
```

mismatch 时提示「结果可能与按发布索引跑出的不同」。这个状态还会写进批次 manifest 的 `params.data_index` 与样本 manifest；**索引 tag 变了会让 `--resume` 失效**（manifest 签名含索引 tag），这是刻意的——换了索引就不该续跑旧结果。

## 5. 谁在读这个 tag（改之前先看这张表）

| 文件:行 | 用法 |
|---|---|
| `data-index.env:22` | **唯一真源**，`DATA_INDEX_RELEASE=<tag>` |
| `.github/workflows/ci.yml:78` | `sed` 提取后 `gh release download` |
| `.github/workflows/docker-release.yml:66` | 同上，打进 Docker 镜像 |
| `test.sh:33` | `sed` 提取；`:35` 有硬编码兜底默认值 `mirdeep-p3-v3.1.4c-full`；`:102` 解压索引后写 `.index-version` |
| `README.md:124` | 源码安装的 `wget` 下载命令 |
| `CONTRIBUTING.md:26` | 贡献者文档的 `wget` 下载命令 |
| `docs/MIRDEEP-P3-INTERNALS.md:20,540,545` | 速查表 / 安装第 2 步 / §12.4 |
| `docs/MIRDEEP-P3-ISSUES.md:133`、`docs/RELEASE-NOTES-3.1.8b.md:36,86,116`、`CHANGELOG.md:74` | 历史记录与发布说明 |

注：上表行号会随 README / 脚本编辑而漂移，改之前先用 `grep -n` 定位；`conda.recipe/meta.yaml` 目前**不读**这个 tag（`data-index.env` 的注释里提到它，属于提醒性质）；将来若给 recipe 加索引校验，记得同步。

## 6. 更新数据索引的 checklist

换了一份新索引 / 新 PmiREN 版本时，按顺序做完这 8 步：

1. **打包并发布**：`tar czf data-index.tar.gz -C data index`，建 Release（tag 建议 `data-v<日期>`，或沿用 `...-full` 风格），上传 asset。**老 Release 与其 tag 不要删**——那是旧的 `data-index.tar.gz` 唯一托管处。
2. **改真源**：`data-index.env:22` 那一行 `DATA_INDEX_RELEASE=<新 tag>`。
3. **改两处下载命令**：`README.md:124`、`CONTRIBUTING.md:26`。
4. **写已装标记**：`echo <新 tag> > data/index/.index-version`。
5. **看启动提示**：跑一次 identification，确认打印 `[index] data index: <新 tag> (matches data-index.env)`。
6. **跑一遍测试**：`./test.sh -t 12 -k`（它按 `data-index.env` 取 tag 下载，等于 CI 预演；第 11 段还会校验 README 的版本块与本文件里的 tag 是否与真源一致，不一致直接 `exit 7`）。
7. **记一笔**：`CHANGELOG.md` 写日期 + 索引 tag 变更 + 影响（哪些注释 / 靶标结果可能变）。
8. **提醒使用者**：旧索引跑出的 annotation / target 结果与新索引不保证一致，需要复跑再比较。

只想升级软件、不动索引（**常态**）：索引相关一行都不用改，`[index]` 会继续显示 `match`。软件发版要动两处：`src/version.py` 的 `__version__`，以及 `README.md` 顶部 “Current versions” 块里的软件版本号（索引 tag 保持不变）。

## 7. 常见问题

- **忘了改 `data-index.env`？** 已装索引比 env 新 → 启动报 `mismatch`；`test.sh` / CI 仍会下旧索引。
- **软件升级了，索引还是旧的，要改吗？** 不用。这正是两条轴分开的理由，`match` 与软件版本无关。
- **删了 `.index-version` 会怎样？** 退化成 `unrecorded`，功能照跑，只是提示不明确；重写一行即可。
- **`mismatch` 会影响结果吗？** 只提示不阻断；但不同索引跑出的位点 / 注释可能不同，比较结果前先对齐索引。
- **索引换了但没换软件，算新版本吗？** 软件版本号不变，改 `data-index.env` 即可；发布记录写在 `CHANGELOG.md` 里。

## 8. 相关文件

- `src/utils/data_index.py`、`tests/test_data_index.py`
- `data-index.env`、`data/index/.index-version`
- `RELEASE.md` 第 8 节（索引 Release 的发布与不可删除）
- `docs/MIRDEEP-P3-INTERNALS.md` §12.4（实现视角）
- `docs/RELEASE-NOTES-3.1.8b.md`（O-10 的落地说明）
