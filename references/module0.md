# 模块零：运行清单

参考规范：三大部分整合文档「模块零：语言与运行时规范」。
姊妹项目 `scrna-pipeline-skill` / `spatial-pipeline-skill` 的
`references/module0.md` 是同一份文档，接口**同名同义**，
三部分的清单可以并排读。

---

## 0. 为什么要有这一层

**没有运行清单的分析结果不是结果。**

半年后拿到一份 `deg_table.csv`，如果不知道当时装的是哪个版本的 limma、
输入的表达矩阵是哪个哈希、随机种子是多少、`trend` / `robust` 开没开，
那份 CSV 就**无法被复现，也无法被质疑** —— 而不可质疑的结论没有价值。

这一层不产出任何生物学结论，它只回答一个问题：
**"这份结果是在什么条件下跑出来的？"**

产物：`results/<GSE>/run_manifest.json`

**为什么不塞进 `state.json`：** `state.json` 记的是"这一步跑没跑成"，
每步重写；manifest 记的是"本轮是在什么条件下跑出来的"，是证据，
写入后不该再变。混在一起会让后者被前者覆盖。

---

## 1. 规范条款与接口的对应

| 条款 | 要求 | 接口 | 落到 manifest 的哪个字段 |
|---|---|---|---|
| §0.3 | `sessionInfo()` / `pip freeze` **全量**输出 | `capture_versions()` | `versions` |
| §0.3 | 关键工具**逐个**记版本 | `capture_versions()` 的 `KEY_PACKAGES` | `key_versions` |
| §0.3 | 所有随机过程固定种子并记录 | `record_params()`（`params` 里含 `seed`）+ `m$seed` | `params`、`seed` |
| §0.4 | 输入数据哈希 | `record_input()` | `inputs` |
| §0.4 | 全部关键参数 | `record_params()` | `params` |
| §0.4 | Agent 决策链 | `record_decision()` | `decisions` |
| §0.4 | 人工干预记录 | `record_human_review()` | `human_review` |
| §0.2 | 跨语言转换前后维度、丢失字段 | `record_cross_language()` | `cross_language` |

---

## 2. 接口清单

| 函数 | 作用 | 备注 |
|---|---|---|
| `MANIFEST_NAME` | `"run_manifest.json"` | 三个仓库一致 |
| `KEY_PACKAGES` | 文档点名的关键工具 | 见 §4 |
| `manifest_path(cfg)` | 清单路径 | |
| `read_manifest(cfg)` | 读清单；文件不存在返回空 list/dict | |
| `init_manifest(cfg)` | **开新一轮，清掉上一轮** | 见 §3.1 |
| `capture_versions(cfg)` | 全量已装包 + `KEY_PACKAGES` 逐个 | |
| `record_input(cfg, path, ...)` | 输入文件的哈希与字节数 | |
| `record_params(cfg, params)` | 参数（`update` 合并，可多次调用） | |
| `record_decision(cfg, node, q, a, evidence)` | 决策链 | |
| `record_human_review(cfg, node, required, status, note)` | 人工复核节点 | |
| `record_cross_language(cfg, src, dst, format, before, after, lost, ...)` | 跨语言转换 | |
| `manifest_summary(cfg)` | 供验收用的摘要 | 见 §3.3 |
| `NAMED_TOOLS` / `probe_named_tools()` / `named_tools_note()` | 点名工具的缺口登记 | **仅 Python 侧**，见 §5 |

---

## 3. 四条不能省的约定

### 3.1 `init_manifest` 必须清掉上一轮

上轮的清单留在那里冒充本轮，**比没有清单更糟** —— 它看起来是证据，
实际是过期证据。这和 AGENTS 规则 14（每步开跑前先删自己的状态文件）
是同一个道理。

### 3.2 输入哈希在**步骤跑完之后**才登记

可选步骤这轮有没有产物，跑完才知道；在开头登记会把"上轮残留"记成本轮输入。

### 3.3 `manifest_summary` 同时报两个缺失数

```
inputs_missing           所有缺失的输入（可见，供人看）
inputs_missing_required  只有 required=TRUE 的缺失（供验收判 FAIL）
```

**第一版只有一个数，把可选输入缺失也算成失败**，结果每轮都红。
可选输入缺失是**正常状态**（外部队列没配、签名文件不存在），
判成 FAIL 会让红灯失去意义。但它必须**可见** —— 所以两个数都要有。

### 3.4 未装的工具要记成 `null`，不能省略键

```json
"key_versions": { "limma": "3.58.1", "GSVA": null }
```

`"GSVA": null` 和"没有 GSVA 这个键"是两件事：前者是"**查过了，没装**"，
后者是"**没查**"。省略会让读者分不清。

---

## 4. `KEY_PACKAGES`：三部分的清单

文档点名的工具**逐个列出**，装了的记版本、没装的记 `null`。

| 部分 | 数量 | 内容 |
|---|---|---|
| **Part 1（本仓库）** | 15 | `GEOquery` `limma` `WGCNA` `clusterProfiler` `GSVA` `glmnet` `survival` `survminer` `timeROC` `rms` `STRINGdb` + `TRRUST` `ChEA3` + `scTenifoldKnk` `PerturbNet` `RegVelo` |
| Part 2 | 25 | `scanpy` `anndata` `scvi-tools` `cellbender` `harmonypy` `scvelo` `celltypist` `pyscenic` `liana` `doubletdetection` `scrublet` + 拟时序 `palantir` `scfates` `cytotrace` + R 包 `monocle3` `slingshot` `cellchat` `soupx` `scdblfinder` + `scTenifoldKnk` `PerturbNet` `RegVelo` |
| Part 3 | 16 | `SpatialDE` `SpatialDE2` `spacexr` `BayesSpace` `SPARK-X` `cell2location` `STAGATE` `SpaGCN` `SpaceFlow` `stLearn` `ISORT` `Bering` `BOMS` + LIANA 等 |

**为什么 Part 1 的清单里会有 `scTenifoldKnk` / `PerturbNet` / `RegVelo`：**
规范把 §1.7/§1.8 标为保留框架、**主语言 Python**，Part 1 的职责是产出
候选靶基因 CSV（`09_export_targets.R`）。所以这三个在本仓库**本来就该是
`null`** —— 但键必须留着，否则读者不知道"是没查还是不该有"。

---

## 5. `NAMED_TOOLS`：点名工具"为什么没用上"的登记（Python 侧）

Part 2 / Part 3 各有 `NAMED_TOOLS` + `probe_named_tools()`。
**Part 1 没有这个结构** —— R 侧点名的工具都在 `KEY_PACKAGES` 里，
而 `timeROC` / `rms` 已经真的接进来了（`10_survival_diagnostics.R`）。

**判据是"理由写了没有"，不是"工具跑了没有"。**
将来某个工具能装了，验收应该依然 PASS（理由变成"已装"），
而不是因为 `available=False` 就变红 —— 那会把"如实记录"惩罚成失败。

四类 `kind`（Python 侧）：

| kind | 含义 |
|---|---|
| `r_package` | R/Bioconductor 包，CI 无 rpy2 |
| `not_on_pypi` | 真包不在 PyPI |
| `deps` / `needs_*` | PyPI 有真包，依赖链或资源跑不动 |
| `name_taken` | **PyPI 上那个名字是另一个不相干的包** |

**`name_taken` 是最危险的一类**，因为 `pip install` 会**成功**。
实测：`edgeR` 是"浏览器重定向"、`slingshot` 是"ElasticSearch 索引迁移"、
`sparkx` 是"高能物理碰撞运动学"、`ISORT` 是 Python 的 import 排序工具。
装不上会立刻报错，**装错了要到跑出结果才发现** —— 而那时结果可能
已经在图上看着挺像回事了。

反例：`SingleR` 在 PyPI 上是 **BiocPy/singler**（作者 Aaron Lun），
是 R 那个算法的官方绑定，不是顶名的。**判断依据是
summary / author / project_urls，不是"名字存不存在"。**

---

## 6. R 侧与 Python 侧必须不同的三处（不是风格问题）

### 6.1 `NA` 要写成 JSON `null`

`jsonlite` 默认把 `NA` 序列化成字符串 `"NA"`，而 Python 侧写的是 `null`。
两边不一致时，"这个工具没装"在读的人看来是"装了一个叫 NA 的工具"。
所以有独立的 `write_manifest()`，显式 `na = "null"`。

### 6.2 版本要用 `utils::installed.packages()` 全量

**不能用 `sessionInfo()$otherPkgs`。** 本仓库脚本一律 `pkg::fun()` 写全名、
不 `attach`，所以 `sessionInfo()$otherPkgs` 是**空的** —— 用它等于什么都没记。

### 6.3 哈希优先 `digest`，退回 `tools::md5sum`，并记 `hash_algo`

`digest` **不在 CI 的 R 包列表里**，所以实际走的是 md5 分支。
**算法不同的哈希不可直接比较**，不写算法等于给了个无法验证的值。
`record_input()` 会把 `hash_algo` 一起写进去。

---

## 7. 人工复核节点

| 节点 | required | 说明 |
|---|---|---|
| `geo_availability` | TRUE | 数据集是否真的可用 |
| `group_labels` | TRUE | 分组标签的推断是否正确 |
| `outlier_removal` | TRUE | 离群样本的取舍 |
| `signature_genes` | FALSE | 预后基因集的生物学合理性（只在有随访终点时才有意义） |
| `virtual_perturbation` | FALSE | 虚拟扰动靶基因的合理性（§1.7/§1.8 保留框架） |

**默认 `pending`，不算失败。** 自动化流水线不能替人签字 ——
把未确认的节点记成已确认，等于把复核节点变成摆设。
但必须**可见**：验收里作为"可见但不阻断"的项列出。

---

## 8. 清单必须和结果同时可及

`run_manifest.json` 落在 `results/<GSE>/` 下，随 artifact 一起上传 ——
**它必须和结果同时可及，否则追溯链是断的。**
只把清单写在 CI 日志里不行：日志会滚掉，文件不会。

---

## 9. 怎么读一份清单

```powershell
$m = Get-Content results/GSE42568/run_manifest.json | ConvertFrom-Json
$m.language          # "R"
$m.seed              # 随机种子
$m.inputs            # 每个输入文件的哈希与算法
$m.key_versions      # 点名工具逐个的版本（null = 查过了没装）
$m.decisions         # Agent 决策链：问题 / 结论 / 证据
$m.human_review      # 人工复核节点及状态
$m.cross_language    # 跨语言转换（本仓库只有 §1.7/§1.8 交接那一处）
```

**验收读的是 `manifest_summary(cfg)$inputs_missing_required`，不是
`inputs_missing`** —— 后者包含可选输入的缺失，那是正常状态。
