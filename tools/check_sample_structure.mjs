#!/usr/bin/env node
// ============================================================================
// tools/check_sample_structure.mjs — 样本相关结构预检（不需要 R）
// ============================================================================
// 为什么需要它：`min Pearson < 0.8 即离群` 这类 QC 规则隐含一个假设 ——
// 所有样本是同一组织的技术重复。如果分组本身对应一个全局表达位移
// （批次效应，或组织成分差异），每个样本都会与另一组的样本低相关，
// 于是 QC 报出"全部样本都是离群"，而真正的问题（分组与位移混杂）被掩盖。
//
// GSE92252 就是这样：9 个样本分成三个紧致簇，簇内 r > 0.92，
// 而 HER2− ↔ HER2+ 只有 0.31–0.38。这个结构在**原始沉积值**里就存在。
// 等 R 流水线跑完 15 分钟才发现，代价太高 —— 所以在这里先查。
//
// 用法：
//   node tools/check_sample_structure.mjs GSE92252
//   node tools/check_sample_structure.mjs GSE92252 --threshold 0.9
// ============================================================================

import { gunzipSync } from 'node:zlib'

const args = process.argv.slice(2)
const gse = (args.find(a => /^GSE\d+$/i.test(a)) || '').toUpperCase()
const thrIdx = args.indexOf('--threshold')
const THRESHOLD = thrIdx >= 0 ? Number(args[thrIdx + 1]) : 0.9

if (!gse) {
  console.error('用法: node tools/check_sample_structure.mjs <GSE编号> [--threshold 0.9]')
  process.exit(2)
}
if (!Number.isFinite(THRESHOLD) || THRESHOLD <= 0 || THRESHOLD > 1) {
  console.error('--threshold 必须是 (0, 1] 之间的数')
  process.exit(2)
}

const stub = `${gse.slice(0, -3)}nnn`
const url = `https://ftp.ncbi.nlm.nih.gov/geo/series/${stub}/${gse}/matrix/${gse}_series_matrix.txt.gz`
process.stderr.write(`下载 ${url}\n`)

const res = await fetch(url)
if (!res.ok) {
  console.error(`下载失败: HTTP ${res.status}。确认 ${gse} 存在且是 array 数据集。`)
  process.exit(1)
}
const text = gunzipSync(Buffer.from(await res.arrayBuffer())).toString('utf8')
const lines = text.split(/\r?\n/)

// ---- 解析 series matrix ----------------------------------------------------
const begin = lines.findIndex(l => /^!series_matrix_table_begin/i.test(l))
if (begin < 0) {
  console.error('series matrix 里没有表达矩阵表（!series_matrix_table_begin）。')
  process.exit(1)
}
const samples = lines[begin + 1].split('\t').slice(1).map(s => s.replace(/^"|"$/g, ''))
const ids = [], mat = []
for (let i = begin + 2; i < lines.length; i++) {
  if (/^!series_matrix_table_end/i.test(lines[i]) || !lines[i].trim()) break
  const f = lines[i].split('\t')
  ids.push(f[0].replace(/^"|"$/g, ''))
  mat.push(f.slice(1).map(Number))
}
const nS = samples.length
if (nS < 3) { console.error(`样本数只有 ${nS}，无法做相关分析。`); process.exit(1) }

// 样本标题（用于给簇起个能看懂的名字）
const titleLine = lines.find(l => /^!Sample_title/i.test(l))
const titles = titleLine
  ? titleLine.split('\t').slice(1).map(s => s.replace(/^"|"$/g, ''))
  : samples

console.log(`${gse}: ${ids.length} 探针 x ${nS} 样本\n`)

// ---- 相关性与标准化 --------------------------------------------------------
const pearson = (a, b) => {
  const n = a.length
  const ma = a.reduce((s, x) => s + x, 0) / n, mb = b.reduce((s, x) => s + x, 0) / n
  let sab = 0, sa = 0, sb = 0
  for (let i = 0; i < n; i++) { const da = a[i] - ma, db = b[i] - mb; sab += da * db; sa += da * da; sb += db * db }
  return sa === 0 || sb === 0 ? NaN : sab / Math.sqrt(sa * sb)
}
const corrOf = cols => cols.map(a => cols.map(b => pearson(a, b)))

const quantileNorm = cols => {
  const n = cols.length, len = cols[0].length
  const sorted = cols.map(c => [...c].sort((a, b) => a - b))
  const avg = Array.from({ length: len }, (_, i) => sorted.reduce((s, c) => s + c[i], 0) / n)
  return cols.map(c => {
    const order = c.map((_, i) => i).sort((a, b) => c[a] - c[b])
    const out = new Array(len)
    order.forEach((idx, rank) => { out[idx] = avg[rank] })
    return out
  })
}

const rawCols = samples.map((_, j) => mat.map(r => r[j]))
// 中位数 > 50 视为原始强度，做 log2；否则认为已经是 log 尺度
const isRaw = rawCols.every(c => {
  const s = [...c].sort((a, b) => a - b)
  return s[Math.floor(s.length / 2)] > 50
})
const logCols = isRaw ? rawCols.map(c => c.map(v => Math.log2(v + 1))) : rawCols
const normCols = quantileNorm(logCols)

console.log(`尺度判断: ${isRaw ? '原始强度（已 log2(x+1)）' : '已是 log 尺度（未再取对数）'}`)
console.log('下面的相关矩阵基于 log2 + quantile 标准化后的值。\n')

const M = corrOf(normCols)

// ---- 打印矩阵 --------------------------------------------------------------
const short = s => s.replace(/^GSM/, '').slice(-3)
console.log('样本:')
samples.forEach((s, i) => console.log(`  [${String(i + 1).padStart(2)}] ${s}  ${titles[i] ?? ''}`))
console.log('\nPearson 相关矩阵:')
console.log('      ' + samples.map((_, j) => String(j + 1).padStart(6)).join(''))
M.forEach((row, i) => {
  console.log(`  [${String(i + 1).padStart(2)}] ` + row.map(v => v.toFixed(3).padStart(6)).join(''))
})

// ---- 单链聚类（在 threshold 处切）------------------------------------------
const parent = samples.map((_, i) => i)
const find = x => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x] } return x }
const union = (a, b) => { const ra = find(a), rb = find(b); if (ra !== rb) parent[rb] = ra }
for (let i = 0; i < nS; i++) {
  for (let j = i + 1; j < nS; j++) if (M[i][j] >= THRESHOLD) union(i, j)
}
const groups = new Map()
for (let i = 0; i < nS; i++) {
  const r = find(i)
  if (!groups.has(r)) groups.set(r, [])
  groups.get(r).push(i)
}
const clusters = [...groups.values()].sort((a, b) => b.length - a.length)

console.log(`\n在 r >= ${THRESHOLD} 处切分，得到 ${clusters.length} 个相关簇:`)
clusters.forEach((c, k) => {
  console.log(`  簇 ${k + 1} (${c.length} 个): ${c.map(i => samples[i]).join(', ')}`)
  console.log(`        ${c.map(i => titles[i]).join(' | ')}`)
})

// ---- 簇内 / 簇间 -----------------------------------------------------------
const inC = [], btC = []
const clusterOf = new Array(nS)
clusters.forEach((c, k) => c.forEach(i => { clusterOf[i] = k }))
for (let i = 0; i < nS; i++) {
  for (let j = i + 1; j < nS; j++) (clusterOf[i] === clusterOf[j] ? inC : btC).push(M[i][j])
}
const mean = a => a.length ? a.reduce((s, x) => s + x, 0) / a.length : NaN
const all = []
for (let i = 0; i < nS; i++) for (let j = i + 1; j < nS; j++) all.push(M[i][j])
all.sort((a, b) => a - b)

console.log(`\n簇内平均 r = ${mean(inC).toFixed(3)}   簇间平均 r = ${mean(btC).toFixed(3)}`)
console.log(`全部配对: min=${all[0].toFixed(3)} median=${all[Math.floor(all.length / 2)].toFixed(3)} max=${all[all.length - 1].toFixed(3)}`)

// ---- 判定 ------------------------------------------------------------------
// 同组织的生物学重复通常 r > 0.95；不同组织的同类样本一般也 > 0.85。
// 簇间平均低于 0.8 说明存在一个与分簇共线的全局位移。
console.log('')
if (clusters.length === 1) {
  console.log('✓ 所有样本互相高度相关，没有明显的分簇结构。')
} else if (mean(btC) >= 0.8) {
  console.log(`△ 存在 ${clusters.length} 个相关簇，但簇间平均 r = ${mean(btC).toFixed(3)} 仍较高，`)
  console.log('  分簇可能只是真实的生物学差异。仍需确认簇的划分是否与分组一致。')
} else {
  console.log(`✗ 簇间平均 r = ${mean(btC).toFixed(3)} 偏低，存在与分簇共线的全局表达位移。`)
  console.log('  如果这些簇恰好对应你的分组，那么分组与这个位移**完全混杂**：')
  console.log('  DEG 结果分不清多少来自分组、多少来自位移，只能作为假设生成。')
  console.log('  下一步：把上面的簇与 assets/config.yml 的 group_values 对照，')
  console.log('  并在结论中显式声明这个限制（参见 EXPERIMENTAL_DESIGN.md §2.9）。')
}
