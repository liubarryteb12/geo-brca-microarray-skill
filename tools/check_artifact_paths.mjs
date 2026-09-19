// ============================================================================
// tools/check_artifact_paths.mjs — data/ 产物与 artifact 清单的一致性检查
// ============================================================================
// 为什么需要它：`check_acceptance()` 查的是 **runner 的工作目录**，
// artifact 的 `path:` 是**另一份清单**。两者不一致时，验收全绿、
// 而用户下载到的产物里缺文件 —— **只有把 artifact 下下来才会发现**。
//
// 实测踩过：`batch_assessment.json` 由 step 00 写出、验收项 PASS，
// 但 artifact 只逐个列举了 group.csv / clean_stats.json / ... 六个，
// 新文件没进去。日志和验收都看不出问题。
//
// 这个检查把两份清单交叉核对：脚本写到 data/<GSE>/ 的文件，
// 是否都在 workflow 的 artifact path 里。
//
// 用法：node tools/check_artifact_paths.mjs
// 退出码：0 = 一致，1 = 有 data/ 产物没进 artifact
// ============================================================================

import { readFileSync, readdirSync, existsSync } from 'node:fs'
import { join, dirname } from 'node:path'
import { fileURLToPath } from 'node:url'

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..')
const WORKFLOW = join(ROOT, '.github', 'workflows', 'geo_analysis.yml')
const SCRIPTS = join(ROOT, 'scripts')

// ---- 1. 从 workflow 里取 artifact path 清单 --------------------------------
const yml = readFileSync(WORKFLOW, 'utf8')
const lines = yml.split(/\r?\n/)

let inPath = false
let pathIndent = -1
const artifactPaths = []
for (const line of lines) {
  const mPath = line.match(/^(\s*)path:\s*\|\s*$/)
  if (mPath) { inPath = true; pathIndent = mPath[1].length; continue }
  if (!inPath) continue
  if (line.trim() === '') continue
  const m = line.match(/^(\s+)(\S.*)$/)
  // 条目必须比 `path:` 那一行更深；同层或更浅即块结束
  if (!m || m[1].length <= pathIndent) { inPath = false; continue }
  artifactPaths.push(m[2].trim())
}

if (artifactPaths.length === 0) {
  console.error('未能从 workflow 里解析出 artifact path 清单 —— 检查器本身失效了，不是通过')
  process.exit(1)
}

// data/ 下被显式列举的文件（去掉 data/${{ matrix.dataset }}/ 前缀）
const PREFIX = 'data/${{ matrix.dataset }}/'
const listed = new Set(
  artifactPaths
    .filter(p => p.startsWith(PREFIX))
    .map(p => p.slice(PREFIX.length))
)
const listsWholeDataDir = artifactPaths.some(p => p.trim() === 'data/${{ matrix.dataset }}/')

// ---- 2. 从 R 脚本里取写到 data_dir 的文件名 --------------------------------
// 匹配形如 file.path(cfg$output$data_dir, "x.json") / file.path(dat, "x.csv")
// 以及 sprintf 拼接的单层文件名。只取**字面量**，变量名一律忽略 ——
// 宁可漏报也不要因为猜错变量而误报。
const WRITE_RE = /file\.path\(\s*(?:cfg\$output\$data_dir|dat)\s*,\s*"([^"\/]+\.(?:json|csv|txt|rds))"/g

const written = new Map() // 文件名 -> 出现它的脚本
for (const f of readdirSync(SCRIPTS).filter(f => f.endsWith('.R'))) {
  const src = readFileSync(join(SCRIPTS, f), 'utf8')
  let m
  while ((m = WRITE_RE.exec(src)) !== null) {
    if (!written.has(m[1])) written.set(m[1], new Set())
    written.get(m[1]).add(f)
  }
}

// 也扫 lib/
const LIB = join(SCRIPTS, 'lib')
if (existsSync(LIB)) {
  for (const f of readdirSync(LIB).filter(f => f.endsWith('.R'))) {
    const src = readFileSync(join(LIB, f), 'utf8')
    let m
    while ((m = WRITE_RE.exec(src)) !== null) {
      if (!written.has(m[1])) written.set(m[1], new Set())
      written.get(m[1]).add(`lib/${f}`)
    }
  }
}

// ---- 3. 交叉核对 -----------------------------------------------------------
console.log(`artifact 里列举的 data/ 文件: ${listed.size} 个`)
if (listsWholeDataDir) console.log('  （清单包含整个 data/<GSE>/ 目录，无需逐个列举）')
console.log(`R 脚本写到 data/ 的文件: ${[...written.keys()].length} 个\n`)

if (listsWholeDataDir) {
  console.log('artifact 覆盖整个 data/ 目录，检查通过')
  process.exit(0)
}

const missing = [...written.entries()]
  .filter(([name]) => !listed.has(name))
  .sort((a, b) => a[0].localeCompare(b[0]))

if (missing.length === 0) {
  console.log('检查通过：所有 data/ 产物都在 artifact 清单里')
  process.exit(0)
}

console.error('以下文件由脚本写到 data/，但**不在 artifact 清单里**：')
console.error('（验收会 PASS，因为文件在 runner 上确实存在；但下载 artifact 拿不到）\n')
for (const [name, scripts] of missing) {
  console.error(`  ${name}`)
  console.error(`      写出位置: ${[...scripts].join(', ')}`)
}
console.error(`\n修法：在 .github/workflows/geo_analysis.yml 的 artifact path 里加：`)
for (const [name] of missing) {
  console.error(`      data/\${{ matrix.dataset }}/${name}`)
}
process.exit(1)
