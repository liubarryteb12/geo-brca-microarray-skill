#!/usr/bin/env node
/**
 * check_r_syntax.mjs — 无 R 运行时的 R 静态检查
 *
 * 本仓库的 CI 在 GitHub Actions 上跑，本地不一定装了 R。这个脚本在没有 R 的情况下
 * 抓出最常见、也最难在 CI 里定位的一类错误：括号/引号不配平、编排器引用了不存在的
 * run_XX 函数、配置字段被引用但未定义。
 *
 * 它**不是** R 解析器，不能替代 `Rscript` 的真实执行。它只是把明显错误挡在推送之前。
 *
 * 用法: node tools/check_r_syntax.mjs [目录，默认 scripts]
 */

import { readFileSync, readdirSync, statSync, existsSync } from 'node:fs'
import { join, relative, resolve, dirname } from 'node:path'
import { fileURLToPath } from 'node:url'

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..')
const targetDir = resolve(root, process.argv[2] ?? 'scripts')

const problems = []
const files = []

function walk(dir) {
  for (const entry of readdirSync(dir)) {
    const full = join(dir, entry)
    if (statSync(full).isDirectory()) walk(full)
    else if (entry.endsWith('.R')) files.push(full)
  }
}
walk(targetDir)
files.sort()

/**
 * 去掉注释与字符串字面量，只留下结构性字符。
 * 同时报告未闭合的字符串（R 里这会导致后续整段被吞掉，报错位置离真因很远）。
 */
function stripLiterals(source, file) {
  let out = ''
  let i = 0
  let line = 1
  const openers = []
  while (i < source.length) {
    const ch = source[i]
    if (ch === '\n') { line++; out += '\n'; i++; continue }

    if (ch === '#') {
      while (i < source.length && source[i] !== '\n') i++
      continue
    }
    if (ch === '"' || ch === "'") {
      const quote = ch
      const startLine = line
      i++
      let closed = false
      while (i < source.length) {
        if (source[i] === '\\') { i += 2; continue }
        if (source[i] === '\n') { line++; i++; continue }
        if (source[i] === quote) { closed = true; i++; break }
        i++
      }
      if (!closed) problems.push(`${file}:${startLine} 字符串未闭合 (${quote})`)
      out += '""'
      continue
    }
    if (ch === '`') {  // R 的反引号标识符，例如 `%||%`
      i++
      while (i < source.length && source[i] !== '`') i++
      i++
      out += 'X'
      continue
    }
    out += ch
    i++
  }
  return out
}

function checkBalance(code, file) {
  const pairs = { ')': '(', ']': '[', '}': '{' }
  const stack = []
  let line = 1
  for (let i = 0; i < code.length; i++) {
    const ch = code[i]
    if (ch === '\n') { line++; continue }
    if (ch === '(' || ch === '[' || ch === '{') stack.push({ ch, line })
    else if (ch in pairs) {
      const top = stack.pop()
      if (!top) problems.push(`${file}:${line} 多余的 '${ch}'`)
      else if (top.ch !== pairs[ch]) {
        problems.push(`${file}:${line} '${ch}' 与第 ${top.line} 行的 '${top.ch}' 不匹配`)
      }
    }
  }
  for (const item of stack) problems.push(`${file}:${item.line} '${item.ch}' 未闭合`)
}

const definedFunctions = new Set()
const calledFunctions = new Set()

for (const file of files) {
  const rel = relative(root, file).split('\\').join('/')
  const source = readFileSync(file, 'utf8')
  const code = stripLiterals(source, rel)
  checkBalance(code, rel)

  for (const m of code.matchAll(/(?:^|\n)\s*([A-Za-z_.][A-Za-z0-9_.]*)\s*<-\s*function/g)) {
    definedFunctions.add(m[1])
  }
  for (const m of code.matchAll(/\b(run_[0-9]+[a-z]?_[A-Za-z0-9_]+)\s*\(/g)) {
    calledFunctions.add(m[1])
  }
}

// 编排器引用的 run_XX 必须真的存在
for (const fn of [...calledFunctions].sort()) {
  if (!definedFunctions.has(fn)) {
    problems.push(`main_analysis.R 引用了未定义的函数 ${fn}()`)
  }
}

// 未使用的定义通常是改名后漏改的残留
const unused = [...definedFunctions]
  .filter(fn => fn.startsWith('run_') && !calledFunctions.has(fn))
if (unused.length > 0) {
  problems.push(`定义了但从未被引用的步骤函数: ${unused.join(', ')}`)
}

// 配置字段引用一致性：cfg$xxx 必须出现在**某个** assets/config.<GSE>.yml 里
//
// 一个数据集一个配置文件，所以取所有配置的字段并集。
//
// **一个配置都没找到时必须报错，不能跳过。** 原来写的是
// `if (existsSync(configPath)) { ... }` —— 配置文件改名后这个判断为假，
// 整段校验静默消失，门禁看着还在、其实已经不检查任何东西了。
// 这类"门禁自己失效"比门禁报错危险得多。
const configFiles = readdirSync(join(root, 'assets'))
  .filter(f => /^config\..*\.ya?ml$/.test(f))
  .sort()
if (configFiles.length === 0) {
  problems.push('assets/ 下没有找到任何 config.<GSE>.yml，配置字段校验无法执行')
} else {
  const configKeys = new Set()
  for (const f of configFiles) {
    const configText = readFileSync(join(root, 'assets', f), 'utf8')
    for (const m of configText.matchAll(/^([a-z_][a-z0-9_]*):/gim)) configKeys.add(m[1])
  }
  const defaults = new Set(['thresholds', 'analysis', 'enrichment', 'output', 'contrast', 'paired'])
  const referenced = new Set()
  for (const file of files) {
    const code = readFileSync(file, 'utf8')
    for (const m of code.matchAll(/\bcfg\$([a-z_][a-z0-9_]*)/g)) referenced.add(m[1])
  }
  for (const key of [...referenced].sort()) {
    if (!configKeys.has(key) && !defaults.has(key)) {
      problems.push(`代码引用了 cfg$${key}，但 ${configFiles.join(' / ')} 中都没有该字段`)
    }
  }
}

/**
 * 长副标题必须走 wrap_subtitle()。
 *
 * ggplot 的副标题**不换行**：超出画布宽度的部分被静默裁掉，不是显示成省略号。
 * 所以"字没显示全"从图上完全看不出来，只能靠静态检查挡。
 * 实测踩过：PPI 副标题 455 字符、火山图 275 字符，都被切掉了尾巴。
 *
 * 判据是**字面量总长度**（不含 sprintf 的格式化参数）超过阈值就必须包 wrap_subtitle。
 * 阈值按最窄的画布算：6.5 英寸、副标题 8.5pt，一行约容纳 96 字符；
 * 留一行余量，超过 100 字符就要求折行。
 */
const SUBTITLE_CHAR_LIMIT = 100

for (const file of files) {
  const rel = relative(root, file).split('\\').join('/')
  const source = readFileSync(file, 'utf8')
  const lines = source.split('\n')
  lines.forEach((text, idx) => {
    // 只认 labs() 的 subtitle 参数。`plot.subtitle = element_text(...)` 是主题设置，
    // 后面的 "bottom" / "horizontal" 之类字面量不是副标题文本，会算成假阳性。
    const m = /(^|[^.\w])subtitle\s*=/.exec(text)
    if (!m) return
    const at = m.index + m[0].length
    // 取这一段到下一个 labs 参数或语句结束为止，统计其中的字符串字面量
    const rest = lines.slice(idx, idx + 30).join('\n')
    const endMatch = rest.slice(at).search(/\n\s{0,20}(x|y|colour|fill|alpha|size|title|tag|shape)\s*=/)
    let seg = endMatch === -1 ? rest.slice(at) : rest.slice(at, at + endMatch)
    if (seg.includes('element_text')) return
    const literals = [...seg.matchAll(/"((?:[^"\\]|\\.)*)"/g)].map(m2 => m2[1])
    const total = literals.reduce((a, s) => a + s.length, 0)
    if (total > SUBTITLE_CHAR_LIMIT && !seg.includes('wrap_subtitle')) {
      problems.push(
        `${rel}:${idx + 1} 副标题字面量 ${total} 字符（> ${SUBTITLE_CHAR_LIMIT}）` +
        `但没有走 wrap_subtitle() —— 超出画布的部分会被静默裁掉`
      )
    }
  })
}

console.log(`检查了 ${files.length} 个 R 文件`)
console.log(`定义的函数: ${definedFunctions.size}，步骤函数调用: ${calledFunctions.size}`)

if (problems.length > 0) {
  console.log(`\n发现 ${problems.length} 个问题:`)
  for (const p of problems) console.log(`  - ${p}`)
  process.exit(1)
}
console.log('\n静态检查通过（注意：这不等于 R 能跑通，仍需真实执行验证）')
