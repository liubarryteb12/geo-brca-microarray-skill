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

// 配置字段引用一致性：cfg$xxx 必须出现在 config.yml 里
const configPath = join(root, 'assets', 'config.yml')
if (existsSync(configPath)) {
  const configText = readFileSync(configPath, 'utf8')
  const configKeys = new Set(
    [...configText.matchAll(/^([a-z_][a-z0-9_]*):/gim)].map(m => m[1])
  )
  const defaults = new Set(['thresholds', 'analysis', 'enrichment', 'output', 'contrast', 'paired'])
  const referenced = new Set()
  for (const file of files) {
    const code = readFileSync(file, 'utf8')
    for (const m of code.matchAll(/\bcfg\$([a-z_][a-z0-9_]*)/g)) referenced.add(m[1])
  }
  for (const key of [...referenced].sort()) {
    if (!configKeys.has(key) && !defaults.has(key)) {
      problems.push(`代码引用了 cfg$${key}，但 assets/config.yml 中没有该字段`)
    }
  }
}

console.log(`检查了 ${files.length} 个 R 文件`)
console.log(`定义的函数: ${definedFunctions.size}，步骤函数调用: ${calledFunctions.size}`)

if (problems.length > 0) {
  console.log(`\n发现 ${problems.length} 个问题:`)
  for (const p of problems) console.log(`  - ${p}`)
  process.exit(1)
}
console.log('\n静态检查通过（注意：这不等于 R 能跑通，仍需真实执行验证）')
