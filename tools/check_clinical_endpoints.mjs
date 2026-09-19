#!/usr/bin/env node
/**
 * check_clinical_endpoints.mjs — 查 GEO 数据集有没有可用的临床终点
 *
 * **为什么需要这个工具。** 选数据集时最容易漏掉的一件事是：这份数据到底有没有
 * 随访。`find_dataset.mjs check` 只打印它猜出来的分组字段，看不到生存列 ——
 * 而"能不能做预后模型"完全取决于有没有随访，以及**有多少个事件**。
 *
 * 实测就是这么发现 GSE42568 和 GSE20685 都有完整的 OS / RFS / 死亡 / 转移字段，
 * 而且都是 GPL570（同平台跨队列验证，省掉跨平台校正）。
 *
 * **事件数比样本数更关键。** LASSO / Cox 模型的过拟合判据是 EPV
 * （events per variable，通行要求 >= 10）。35 个事件的队列，签名超过 3 个基因
 * 就开始过拟合 —— 文献里常见的"8 基因预后签名"配 35 个事件，是典型的过拟合。
 * 所以本工具直接把 EPV 换算成**签名基因数上限**打出来。
 *
 * 用法:
 *   node tools/check_clinical_endpoints.mjs GSE42568 [GSE...]
 *
 * 退出码: 至少一个数据集检出可用的生存终点 -> 0，否则 -> 1
 */

const GEO_SOFT = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi'
const EPV_TARGET = 10

/** 生存 / 随访相关字段的识别模式 */
const SURVIVAL_HINTS =
  /surviv|follow.?up|vital|death|dead|deceased|recurren|relapse|metastas|dfs|rfs|dss|pfs|time.?to|event/i

/** 事件列（0/1 编码）的识别模式 */
const EVENT_HINTS = /event|status|death|dead|deceased|recurren|relapse|metastas/i

/**
 * 时间列的识别模式。
 *
 * **时间列的判定必须优先于事件列。** 实测 `relapse free survival time_days`
 * 里含 "relapse"，若先按事件列的模式排除，这个时间列就被丢掉了 ——
 * 于是 RFS 只能错配到 OS 的时间列上，而两者事件数不同（48 vs 35），
 * 配错了整个模型就错了。
 */
const TIME_HINTS = /\btime\b|time_|_time|days|months|years|duration|follow.?up/i

async function soft(acc, targ = 'self') {
  const url = `${GEO_SOFT}?acc=${acc}&targ=${targ}&form=text&view=brief`
  const res = await fetch(url, { headers: { 'User-Agent': 'geo-brca-pipeline/1.0' } })
  if (!res.ok) throw new Error(`${acc} HTTP ${res.status}`)
  return await res.text()
}

/** GEO SOFT 的键名**带前缀**：!Series_title / !Series_type / !Series_platform_id */
function seriesFields(text) {
  const out = {}
  for (const line of text.split('\n')) {
    const m = /^!(\w+)\s*=\s*(.*)$/.exec(line.trim())
    if (!m) continue
    ;(out[m[1]] ??= []).push(m[2].replace(/^"|"$/g, ''))
  }
  return out
}

/** 样本级同样是 !Sample_characteristics_ch1；返回 key -> Map(value -> count) */
function characteristics(text) {
  const blocks = text.split('^SAMPLE = ').slice(1)
  const keys = new Map()
  for (const b of blocks) {
    for (const line of b.split('\n')) {
      const m = /^!Sample_characteristics_ch1\s*=\s*(.*)$/.exec(line.trim())
      if (!m) continue
      const raw = m[1].replace(/^"|"$/g, '')
      const i = raw.indexOf(':')
      const key = (i === -1 ? '(无键名)' : raw.slice(0, i)).trim().toLowerCase()
      const val = (i === -1 ? raw : raw.slice(i + 1)).trim()
      if (!keys.has(key)) keys.set(key, new Map())
      const vm = keys.get(key)
      vm.set(val, (vm.get(val) ?? 0) + 1)
    }
  }
  return { n: blocks.length, keys }
}

/** 从 0/1 编码的取值分布里数出事件数 */
function countEvents(valueMap) {
  let yes = 0
  let no = 0
  for (const [v, c] of valueMap) {
    const s = v.toLowerCase()
    if (['1', 'yes', 'true', 'dead', 'deceased', 'event', 'recurrence', 'relapse'].includes(s)) yes += c
    else if (['0', 'no', 'false', 'alive', 'censored', 'none'].includes(s)) no += c
  }
  return { yes, no, total: yes + no }
}

async function inspect(acc) {
  const self = await soft(acc, 'self')
  const sf = seriesFields(self)
  const one = k => (sf[k] ?? [''])[0].replace(/\s+/g, ' ')

  console.log('='.repeat(78))
  console.log(`  ${acc}  —  ${one('Series_title')}`)
  console.log('='.repeat(78))
  console.log(`  类型     : ${one('Series_type')}`)
  console.log(`  平台     : ${(sf.Series_platform_id ?? []).join(', ')}`)
  console.log(`  样本数   : ${(sf.Series_sample_id ?? []).length}`)
  console.log(`  PubMed   : ${(sf.Series_pubmed_id ?? []).join(', ') || '（无）'}`)

  const { n, keys } = characteristics(await soft(acc, 'gsm'))
  console.log(`  解析样本 : ${n}`)

  console.log('\n  --- characteristics 全部字段 ---')
  for (const [key, vm] of [...keys.entries()].sort()) {
    const vals = [...vm.entries()].sort((a, b) => b[1] - a[1])
    const shown = vals.slice(0, 5).map(([v, c]) => `${v}×${c}`).join('  ')
    const mark = SURVIVAL_HINTS.test(key) ? '  <<<' : ''
    console.log(`    ${key.padEnd(30)} ${String(vals.length).padStart(3)} 种  ${shown}${vals.length > 5 ? ' …' : ''}${mark}`)
  }

  // ---- 终点判定 ------------------------------------------------------------
  //
  // **生存终点必须同时有事件列和时间列。** 只匹配 "status" 会把 er_status /
  // lymph node status 这种二元临床性状也算成生存终点 —— 它们没有时间维度，
  // 做不了 Cox 模型。分开报告，因为两者的模型完全不同：
  //   事件列 + 时间列 -> Cox LASSO，EPV 按**事件数**算
  //   只有 0/1        -> logistic LASSO，EPV 按**少数类样本数**算
  //
  // **名称配对只作提示，不作断言。** 实测 GSE20685 的 `event_death` 配的是
  // `follow_up_duration (years)` —— 名字里没有任何共同词，纯靠名称匹配会漏掉。
  // 所以名称对不上时如实说"需要人工确认"，不猜。
  // **先分时间列，再分事件列** —— 见 TIME_HINTS 处的注释。
  const times = [...keys.keys()].filter(k => TIME_HINTS.test(k))
  const timeSet = new Set(times)
  const events = []
  const binaryTraits = []
  for (const [key, vm] of keys) {
    if (timeSet.has(key)) continue
    if (!EVENT_HINTS.test(key)) continue
    const { yes, no, total } = countEvents(vm)
    if (yes === 0 || no === 0 || total < 0.5 * n) continue
    const stem = key.replace(/\b(event|status)\b/gi, '').replace(/[\s_]+/g, ' ').trim()
    const words = stem.split(' ').filter(w => w.length > 3)
    // **按匹配词数取最高分，不是取第一个命中的。**
    // `relapse free survival event` 和 `overall survival event` 共享 "survival" 一词，
    // 取第一个会把 RFS 错配到 OS 的时间列上 —— 两个终点的事件数不同（48 vs 35），
    // 配错了整个模型就错了。
    let mate = null
    let bestScore = 0
    for (const t of times) {
      const tl = t.toLowerCase()
      const score = words.filter(w => tl.includes(w)).length
      if (score > bestScore) { bestScore = score; mate = t }
    }
    if (mate) events.push({ key, time: mate, events: yes, atRisk: total, score: bestScore })
    else binaryTraits.push({ key, pos: yes, neg: no, minority: Math.min(yes, no), mayPair: times.length > 0 })
  }

  console.log('\n  --- 事件列（0/1 编码）---')
  if (events.length === 0 && binaryTraits.length === 0) {
    console.log('    **无** —— 没有任何 0/1 事件列，做不了预后模型')
  }
  for (const e of events) {
    console.log(`    ${e.key.padEnd(30)} 事件 ${String(e.events).padStart(3)} / ${e.atRisk}`
      + `   时间列（${e.score} 词匹配）: ${e.time}`)
    console.log(`    ${''.padEnd(30)} -> Cox LASSO 签名上限: EPV>=10 -> ${Math.floor(e.events / EPV_TARGET)} 基因,`
      + ` EPV>=5 -> ${Math.floor(e.events / 5)} 基因`)
    if (e.events < 10) console.log(`    ${''.padEnd(30)} **事件 < 10，连单变量 Cox 都不该做**`)
    if (e.score < 2) console.log(`    ${''.padEnd(30)} **只匹配上 ${e.score} 个词，配对关系需人工确认**`)
  }
  for (const t of binaryTraits) {
    console.log(`    ${t.key.padEnd(30)} ${t.pos} vs ${t.neg}（少数类 ${t.minority}）`
      + `   时间列: 名称未配对${t.mayPair ? '，但本数据集有时间列，需人工确认' : ''}`)
    console.log(`    ${''.padEnd(30)} -> logistic LASSO 签名上限: EPV>=10 -> ${Math.floor(t.minority / EPV_TARGET)} 基因`)
  }

  console.log(`\n  --- 时间列 ---`)
  console.log(`    ${times.length ? times.join(' | ') : '**无**'}`)

  const usable = events.length > 0 || (binaryTraits.length > 0 && times.length > 0)
  console.log(`\n  结论: ${usable
    ? (events.length > 0
        ? `可做 Cox 生存模型（${events.length} 个名称可配对的事件列）`
        : '事件列与时间列的名称对不上 —— **需人工确认配对关系**后才能做生存模型')
    : '**做不了预后模型**（没有事件列，或只有事件列没有时间列）'}`)
  console.log()
  return usable
}

const accs = process.argv.slice(2).filter(a => /^GSE\d+$/i.test(a))
if (accs.length === 0) {
  console.error(`用法: node tools/check_clinical_endpoints.mjs GSE42568 [GSE...]

查 GEO 数据集有没有可用的临床终点（生存 / 复发 / 转移），
并把事件数换算成 LASSO/Cox 签名的基因数上限（EPV >= ${EPV_TARGET}）。

退出码: 至少一个数据集有可用终点 -> 0，否则 -> 1`)
  process.exit(2)
}

let any = false
for (const acc of accs) {
  try {
    if (await inspect(acc)) any = true
  } catch (err) {
    console.log(`  ${acc} 查询失败: ${err.message}\n`)
  }
  await new Promise(r => setTimeout(r, 400))   // 对 GEO 客气一点
}
process.exit(any ? 0 : 1)
