#!/usr/bin/env node
/**
 * find_dataset.mjs — GEO 数据集预检工具
 *
 * 在跑分析流水线之前，用它确认候选数据集真的满足硬约束。存在的意义：
 * GEO 的标题和二手描述经常与真实元数据不符 —— 例如 GSE197894 的摘要被描述为
 * "表达谱芯片"，其 gdstype 实际是 "Expression profiling by high throughput
 * sequencing"（RNA-seq）。用标题判断数据集类型会直接毁掉整个实验设计。
 *
 * 只依赖 Node 内置 fetch，不需要 R 或 Bioconductor。
 *
 * 用法:
 *   node scripts/find_dataset.mjs check GSE92252 [GSE...]
 *   node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10
 *   node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10 --two-groups
 *   node scripts/find_dataset.mjs samples GSE92252
 *
 * 退出码: check 模式下全部合规 -> 0，否则 -> 1
 */

const EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
const GEO_SOFT = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi'

const sleep = ms => new Promise(r => setTimeout(r, ms))

async function fetchText(url, { retries = 4 } = {}) {
  let lastErr
  for (let attempt = 0; attempt < retries; attempt++) {
    try {
      const res = await fetch(url, { headers: { 'User-Agent': 'geo-brca-pipeline/1.0' } })
      if (res.status === 429) {
        lastErr = new Error('HTTP 429 rate limited')
        await sleep(1200 * (attempt + 1))
        continue
      }
      if (!res.ok) throw new Error(`HTTP ${res.status}`)
      return await res.text()
    } catch (err) {
      lastErr = err
      await sleep(600 * (attempt + 1))
    }
  }
  throw new Error(`请求失败 ${url}: ${lastErr?.message}`)
}

async function fetchJson(url, opts) {
  const text = await fetchText(url, opts)
  if (text.trimStart().startsWith('<')) {
    throw new Error(`E-utilities 返回错误页: ${text.slice(0, 160)}`)
  }
  return JSON.parse(text)
}

/** GDS UID 与 GSE 编号的固定偏移 */
const gseToUid = gse => String(200000000 + Number(String(gse).replace(/^GSE/i, '')))

/** 抓 GEO SOFT 文本（比完整 family 文件小几个数量级） */
async function fetchSoft(acc, targ = 'self', view = 'brief') {
  const url = `${GEO_SOFT}?acc=${acc}&targ=${targ}&form=text&view=${view}`
  return await fetchText(url)
}

/** 取 series 级元数据 */
async function seriesInfo(gse) {
  const uid = gseToUid(gse)
  const json = await fetchJson(`${EUTILS}/esummary.fcgi?db=gds&retmode=json&id=${uid}`)
  return json.result?.[uid] ?? null
}

/** 把 SOFT 按 ^SAMPLE 切块并提取字段 */
function parseSamples(softText) {
  const blocks = softText.split('^SAMPLE = ').slice(1)
  return blocks.map(block => {
    const lines = block.split(/\r?\n/)
    const pick = key => lines
      .filter(l => l.startsWith(`!Sample_${key} = `))
      .map(l => l.slice(`!Sample_${key} = `.length))
    return {
      gsm: pick('geo_accession')[0] ?? '(unknown)',
      title: pick('title')[0] ?? '',
      source: pick('source_name_ch1')[0] ?? '',
      organism: pick('organism_ch1')[0] ?? '',
      characteristics: pick('characteristics_ch1'),
    }
  })
}

// ============================================================================
// 门禁常量 —— **必须与 scripts/00_validate_inputs.R 保持一致**
// ============================================================================
//
// 这里原来只有一套门禁（`n > 0 && n < maxSamples`，默认 10），也就是只实现了
// `small_sample`。后果是**预检工具和真实门禁互相矛盾**：`check GSE42568`（n=121）
// 报「样本量不合规」，而同一个数据集在 `design_mode: cohort` 下跑得好好的。
//
// 泛化到"任意疾病 + 任意芯片"之后这个矛盾更要紧 —— 任意疾病下用户经常需要
// 队列级数据集（WGCNA 要 >=15，LASSO 要有足够事件），预检却说它不合规。
//
// 两套门禁由 `design_mode` 选择，不是一个门禁换个阈值（AGENTS.md 规则 1）。
const GATES = {
  small_sample: { minSamples: 1,  maxSamples: 10, minPerGroup: 3  },
  cohort:       { minSamples: 15, maxSamples: Infinity, minPerGroup: 10 },
}

/** 对单个数据集跑全部硬约束 */
async function checkDataset(gse, opts = {}) {
  const designMode = opts.designMode ?? 'small_sample'
  const gate = GATES[designMode]
  if (!gate) {
    return { gse, ok: false, designMode,
             failures: [`design_mode 非法: ${designMode}（只能是 ${Object.keys(GATES).join(' / ')}）`],
             samples: [] }
  }
  const wantOrganism = opts.organism ?? 'Homo sapiens'
  const info = await seriesInfo(gse)
  if (!info) {
    return { gse, ok: false, designMode, failures: [`GEO 中查不到 ${gse}`], samples: [] }
  }

  const failures = []
  const warnings = []

  if (info.taxon !== wantOrganism) {
    failures.push(`物种不合规: ${info.taxon ?? '(未知)'}（要求 ${wantOrganism}）`)
  }
  const type = info.gdstype ?? ''
  if (/sequencing/i.test(type)) {
    failures.push(`类型不合规: ${type} —— 这是测序，不是基因芯片`)
  } else if (!/array/i.test(type)) {
    failures.push(`类型不合规: ${type}`)
  }
  const n = Number(info.n_samples ?? 0)
  if (designMode === 'small_sample') {
    if (!(n >= gate.minSamples && n < gate.maxSamples)) {
      failures.push(`样本量不合规: ${n}（small_sample 要求 >= ${gate.minSamples} 且 < ${gate.maxSamples}）`)
    }
  } else {
    if (!(n >= gate.minSamples)) {
      failures.push(`样本量不合规: ${n}（cohort 要求 >= ${gate.minSamples}）`)
    }
  }

  // 拉样本级信息做分组可行性判断
  let samples = []
  try {
    samples = parseSamples(await fetchSoft(gse, 'gsm'))
  } catch (err) {
    warnings.push(`无法获取样本级元数据: ${err.message}`)
  }

  const groups = groupCandidates(samples)
  if (samples.length > 0 && groups.length < 2) {
    warnings.push(`样本元数据中只识别出 ${groups.length} 个候选组，可能无法构造两组对比`)
  }
  // 每组样本数是真实门禁的一部分（R 里按 config 的 group_values 数），
  // 预检拿不到 config，只能对候选组做启发式提示 —— 所以是 warning 不是 failure。
  const thin = groups.filter(g => g.count < gate.minPerGroup)
  if (thin.length > 0) {
    warnings.push(`候选组里有 ${thin.length} 个不足 ${gate.minPerGroup} 例: ` +
                  thin.map(g => `"${g.label}"×${g.count}`).join(', '))
  }

  return {
    gse,
    designMode,
    title: info.title,
    taxon: info.taxon,
    type,
    platform: `GPL${info.gpl}`,
    nSamples: n,
    ok: failures.length === 0,
    failures,
    warnings,
    groups,
    samples,
  }
}

/** 从 characteristics 里猜分组：tissue / group / diagnosis 等字段的取值 */
function groupCandidates(samples) {
  const buckets = new Map()
  for (const s of samples) {
    const fields = s.characteristics.map(c => {
      const idx = c.indexOf(':')
      return idx >= 0 ? c.slice(0, idx).trim().toLowerCase() : c.trim().toLowerCase()
    })
    const tissue = s.characteristics
      .filter(c => /tissue|group|diagnosis|condition|source|type/i.test(c.slice(0, c.indexOf(':') + 1)))
      .map(c => c.slice(c.indexOf(':') + 1).trim())
      .join(' | ')
    const key = tissue || fields.join(' | ') || s.source || '(unknown)'
    if (!buckets.has(key)) buckets.set(key, [])
    buckets.get(key).push(s.gsm)
  }
  return [...buckets.entries()]
    .map(([label, gsms]) => ({ label, count: gsms.length, gsms }))
    .sort((a, b) => b.count - a.count)
}

async function search(opts) {
  // **`--disease` 是必填的，没有默认值。**
  //
  // 原来写的是 `opts.disease ?? 'breast cancer'`。本仓库的定位是"任意疾病 + 任意芯片"，
  // 一个静默默认在这里正是最坏的情况：用户想找前列腺癌数据集、忘了加 `--disease`，
  // 得到一份**看起来完全正常**的乳腺癌候选列表，一路做下去。
  // 这与规则 2 对 `parse_args()` 的要求是同一个理由 —— 静默默认是"跑错数据集"的来源。
  if (!opts.disease || !String(opts.disease).trim()) {
    throw new Error(
      'search 必须显式指定 --disease（没有默认值）。\n' +
      '  例: node scripts/find_dataset.mjs search --disease "prostate cancer" --max-samples 10\n' +
      '  静默默认成某个疾病，会让"找错数据集"看起来像正常结果。')
  }
  const disease = String(opts.disease).trim()
  const designMode = opts.designMode ?? 'small_sample'
  const gate = GATES[designMode]
  if (!gate) {
    throw new Error(`--design-mode 非法: ${designMode}（只能是 ${Object.keys(GATES).join(' / ')}）`)
  }
  const minSamples = Number(opts.minSamples ?? gate.minSamples)
  const maxSamples = Number(opts.maxSamples ?? (designMode === 'cohort' ? 100000 : 10))
  const limit = Number(opts.limit ?? 30)

  const term = [
    `(${disease}[Title])`,
    '"Homo sapiens"[Organism]',
    '"Expression profiling by array"[DataSet Type]',
    'gse[Entry Type]',
  ].join(' AND ')

  const search = await fetchJson(
    `${EUTILS}/esearch.fcgi?db=gds&retmode=json&retmax=500&term=${encodeURIComponent(term)}`
  )
  const ids = search.esearchresult?.idlist ?? []
  console.log(`GEO 命中 ${search.esearchresult?.count ?? 0} 个 series，检查前 ${ids.length} 个...`)

  // 批量取 series 摘要（esummary 支持逗号分隔的多 id，一次请求即可）
  const rows = []
  for (let i = 0; i < ids.length; i += 100) {
    const chunk = ids.slice(i, i + 100)
    const json = await fetchJson(`${EUTILS}/esummary.fcgi?db=gds&retmode=json&id=${chunk.join(',')}`)
    for (const uid of json.result?.uids ?? []) rows.push(json.result[uid])
    await sleep(400)
  }

  const candidates = rows
    .filter(r => r?.taxon === 'Homo sapiens')
    .filter(r => !/sequencing/i.test(r.gdstype ?? ''))
    .filter(r => /array/i.test(r.gdstype ?? ''))
    .map(r => ({ gse: r.accession, n: Number(r.n_samples ?? 0), gpl: r.gpl, title: r.title }))
    .filter(r => r.n >= minSamples && r.n < maxSamples)
    .sort((a, b) => a.n - b.n)

  console.log(`\n满足「人源 + 芯片 + ${minSamples} <= 样本数 < ${maxSamples}」的候选: ${candidates.length}\n`)

  const shown = opts.twoGroups ? candidates : candidates.slice(0, limit)
  for (const c of shown) {
    if (!opts.twoGroups) {
      console.log(`${c.gse}  n=${String(c.n).padStart(2)}  GPL${c.gpl}  ${(c.title ?? '').slice(0, 100)}`)
      continue
    }
    // --two-groups 会逐个拉样本元数据，慢但能直接筛出可做两组对比的
    try {
      const detail = await checkDataset(c.gse, { designMode })
      const viable = detail.groups.filter(g => g.count >= gate.minPerGroup)
      if (viable.length >= 2) {
        console.log(`${c.gse}  n=${c.n}  GPL${c.gpl}  ${(c.title ?? '').slice(0, 80)}`)
        for (const g of detail.groups) {
          console.log(`      ${String(g.count).padStart(2)}x  ${g.label.slice(0, 100)}`)
        }
      }
    } catch (err) {
      console.log(`${c.gse}  (样本元数据获取失败: ${err.message})`)
    }
    await sleep(350)
  }
  console.log(`\n提示: 用 "node scripts/find_dataset.mjs check <GSE>" 做完整合规校验`)
}

async function main() {
  const [command, ...rest] = process.argv.slice(2)
  const opts = {}
  const positional = []
  for (let i = 0; i < rest.length; i++) {
    if (rest[i].startsWith('--')) {
      const key = rest[i].slice(2).replace(/-([a-z])/g, (_, c) => c.toUpperCase())
      if (rest[i + 1] !== undefined && !rest[i + 1].startsWith('--')) { opts[key] = rest[i + 1]; i++ }
      else opts[key] = true
    } else positional.push(rest[i])
  }

  if (command === 'check') {
    if (positional.length === 0) { console.error('用法: check GSE92252 [GSE...] [--design-mode small_sample|cohort]'); process.exit(2) }
    const designMode = opts.designMode ?? 'small_sample'
    let allOk = true
    for (const gse of positional) {
      const r = await checkDataset(gse, { designMode })
      console.log(`\n${'='.repeat(74)}\n${r.gse}  ${r.ok ? 'PASS' : 'FAIL'}  (design_mode: ${r.designMode})\n${'='.repeat(74)}`)
      if (r.title) console.log(`title    : ${r.title}`)
      console.log(`taxon    : ${r.taxon ?? '(未知)'}`)
      console.log(`type     : ${r.type ?? '(未知)'}`)
      console.log(`platform : ${r.platform ?? '(未知)'}`)
      console.log(`n_samples: ${r.nSamples ?? '(未知)'}`)
      if (r.groups?.length) {
        console.log('groups   :')
        for (const g of r.groups) console.log(`   ${String(g.count).padStart(2)}x  ${g.label.slice(0, 110)}`)
      }
      for (const f of r.failures ?? []) console.log(`  FAIL  ${f}`)
      for (const w of r.warnings ?? []) console.log(`  WARN  ${w}`)
      if (!r.ok) allOk = false
      await sleep(350)
    }
    process.exit(allOk ? 0 : 1)
  }

  if (command === 'samples') {
    if (positional.length === 0) { console.error('用法: samples GSE92252'); process.exit(2) }
    const samples = parseSamples(await fetchSoft(positional[0], 'gsm'))
    for (const s of samples) {
      console.log(`${s.gsm}  ${s.title}\n    source: ${s.source}\n    ${s.characteristics.join('\n    ')}`)
    }
    console.log(`\n${samples.length} samples`)
    return
  }

  if (command === 'search') return await search(opts)

  console.log(`find_dataset.mjs — GEO 数据集预检

用法:
  node scripts/find_dataset.mjs check GSE92252 [GSE...] [--design-mode small_sample|cohort]
      校验物种 / 数据类型 / 样本量，并列出候选分组
      **两套门禁由 --design-mode 选择**（默认 small_sample），与 00_validate_inputs.R 一致：
        small_sample: 1 <= n < 10，每组 >= 3
        cohort:       n >= 15，每组 >= 10

  node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10
      按疾病检索人源芯片数据集，筛出样本数小于阈值的候选
      **--disease 必填**（没有默认值）；疾病名用英文，如 "prostate cancer"、"COPD"

  node scripts/find_dataset.mjs search --disease "breast cancer" --design-mode cohort --limit 20
      找**队列级**数据集（n >= 15）。cohort 模式下 --min-samples 默认 15、--max-samples 不限
      WGCNA 要 >= 15 例；LASSO-Cox 还要有随访终点（另用 check_clinical_endpoints.mjs 查）

  node scripts/find_dataset.mjs search --disease "breast cancer" --max-samples 10 --two-groups
      额外拉取样本元数据，只输出能凑出两个 >=3 样本组的数据集（慢）

  node scripts/find_dataset.mjs samples GSE92252
      打印每个样本的完整 characteristics

退出码: check 全部合规 -> 0，否则 -> 1`)
  if (command) process.exit(2)
}

main().catch(err => { console.error(`错误: ${err.message}`); process.exit(1) })
