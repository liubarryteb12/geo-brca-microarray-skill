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

/**
 * 相邻字符串字面量 —— R 没有隐式字符串拼接。
 *
 *   method = ("第一段"
 *             "第二段")        # <- unexpected string constant
 *
 * Python / C 会把相邻字面量接起来，R **不会**，这是语法错误。
 * 而 `c("a", "b")` / `paste0("a", "b")` 是逗号分隔的多参数，合法。
 * 所以判据是：**两个字符串字面量之间除了空白/注释什么都没有**。
 *
 * **为什么必须在 stripLiterals 之前查：** 那个函数把每个字符串换成 `""`，
 * 于是 `("a" "b")` 变成 `("" "")`，括号配平完全正常 —— 这个错实测从
 * check_r_syntax.mjs 眼皮底下溜过去，到云端 `source()` 才炸
 * （run 35482359507，09_export_targets.R:258）。
 *
 * 跨行的字符串（R 允许字面换行）也可能让这里误报，但那种写法本身罕见，
 * 报出来人看一眼就知道是不是真的。
 */
function checkImplicitConcat(source, file) {
  let i = 0
  let line = 1
  const lineAt = () => line
  while (i < source.length) {
    const ch = source[i]
    if (ch === '\n') { line++; i++; continue }
    if (ch === '#') { while (i < source.length && source[i] !== '\n') i++; continue }
    if (ch !== '"' && ch !== "'") { i++; continue }

    // 读掉一个字符串字面量
    const quote = ch
    const startLine = lineAt()
    i++
    while (i < source.length) {
      if (source[i] === '\\') { i += 2; continue }
      if (source[i] === '\n') { line++; i++; continue }
      if (source[i] === quote) { i++; break }
      i++
    }

    // 往后跳过空白与注释，看下一个有效字符是不是又一个引号
    let j = i
    let jLine = line
    for (;;) {
      const c = source[j]
      if (c === undefined) break
      if (c === '\n') { jLine++; j++; continue }
      if (c === ' ' || c === '\t' || c === '\r') { j++; continue }
      if (c === '#') { while (j < source.length && source[j] !== '\n') j++; continue }
      break
    }
    const next = source[j]
    if (next === '"' || next === "'") {
      problems.push(
        `${file}:${startLine} 相邻字符串字面量没有逗号 —— R 没有隐式字符串拼接，` +
        `用 paste0(...) 或 c(...) 显式连接（下一段在第 ${jLine} 行）`
      )
    }
  }
}

const definedFunctions = new Set()
const calledFunctions = new Set()

// ---- 检查 A：基础绘图函数不得写成 stats::（错误台账 E-03）------------------
// 实测踩过两次：`stats::abline` / `stats::axis` —— 它们都在 **graphics** 命名空间。
// 本仓库禁止 library()，每个函数都要写全名，所以写错命名空间会直接报
// "'X' is not an exported object from 'namespace:stats'"，而那个报错
// 不会提示"应该换成 graphics"。
const GRAPHICS_FNS = new Set([
  'plot', 'image', 'axis', 'abline', 'mtext', 'par', 'layout', 'legend',
  'text', 'title', 'lines', 'points', 'box', 'grid', 'rect', 'polygon',
  'hist', 'barplot', 'pie', 'contour', 'persp', 'pairs', 'matplot',
])
const STATS_WRONGLY = /stats::([A-Za-z_.][A-Za-z0-9_.]*)/g

// ---- 检查 B：sprintf 格式串里的裸 %（错误台账 E-02）------------------------
// 实测：`sprintf("... top 2% of both ...", rho)` 报 `too few arguments` ——
// 文本里的 `% o` 被当成八进制转换符 `%o`，多吃一个参数。烧了 6 轮 CI 才定位。
//
// **判据要窄，否则误报淹没真信号**（第一版按"总 % 数 > 合法转换符数"判，
// 对跨行拼接的格式串误报 12 处 —— 那些 % 分布在多个字符串里，各自合法）。
// 真正危险的只有一种模式：**一个字符串字面量内部**，`%` 后面紧跟
// 空格 + 字母（`% o` / `% a`）或紧跟字母但不是合法转换符 —— 那才是被误读的转换符。
// 合法的 `%%`、`%d`、`%.2f`、`%s`、`%5.1f` 等一律不报。
// 真正危险的**只有一种**模式：`%` 后面紧跟**空格**、再跟字母（`% o` / `% a`）。
// R 把 `%` 后的空格当 flags、把那个字母当转换符（`o`=八进制、`a`/`e`/`f`/`g`=浮点…），
// 于是多吃一个参数。而 `%-52s`、`%5.1f%%`、`%H:%M:%S`、`%.1f%%` 这些
// **各自合法**（`%%` 是转义、`%H` 在 strftime 里不是 sprintf 格式串）。
// 第一版判据太宽 → 误报 18 处，把真信号淹没；收窄到"空格 + 字母"后只剩真问题。
const BAD_PERCENT = /(?<!%)% +[a-zA-Z]/g   // 负向后视：%% 是转义，不算裸 %
function checkSprintfPercent(source, file) {
  const lines = source.split(/\r?\n/)
  lines.forEach((line, i) => {
    if (!/sprintf\s*\(/.test(line)) return
    if (/^\s*#/.test(line)) return
    for (const m of line.matchAll(/"((?:[^"\\]|\\.)*)"/g)) {
      const lit = m[1]
      if (!lit.includes('%')) continue
      BAD_PERCENT.lastIndex = 0
      if (BAD_PERCENT.test(lit)) {
        problems.push(
          `${file}:${i + 1} sprintf 格式串里有裸 %（"${lit.slice(0, 40)}..."）—— ` +
          `文本里的 % 必须写 %%，否则 R 当成转换符多吃参数，报 too few arguments`);
      }
    }
  })
}

for (const file of files) {
  const rel = relative(root, file).split('\\').join('/')
  const source = readFileSync(file, 'utf8')
  const code = stripLiterals(source, rel)
  checkBalance(code, rel)
  // **必须在原始 source 上查**，不能用 stripLiterals 的结果（见函数注释）
  checkImplicitConcat(source, rel)
  checkSprintfPercent(source, rel)

  // stats:: 误用检查（在剥掉注释与字符串的 code 上查，避免注释里的例子误报）
  for (const m of code.matchAll(STATS_WRONGLY)) {
    if (GRAPHICS_FNS.has(m[1])) {
      problems.push(
        `${rel}: stats::${m[1]} 不存在 —— 基础绘图函数在 graphics 命名空间，` +
        `写 stats:: 会报 "not an exported object from 'namespace:stats'"`);
    }
  }

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

/**
 * 画布尺寸符号必须真的定义在 lib/common.R 里。
 *
 * 为什么需要：R 里用到未定义的名字**要到运行时才炸**，静态检查看不见。
 * 姊妹项目 Python 侧实测漏 import 一个 `W_SINGLE`，`py_compile` 照样报
 * "语法通过"，白跑了一整轮 CI。R 侧同理。
 *
 * 检查两件事：
 *   1. 用到的尺寸符号在 common.R 里有定义（挡住拼错 `W_DOUBL` 这类）
 *   2. 脚本确实 source 了 common.R（否则符号在作用域里根本不存在）
 *
 * 名单是**写死的**，不是"common.R 导出的所有名字" —— 后者会把函数参数名、
 * 局部变量当成漏定义，误报一堆（Python 侧第一版就是这么误报的）。
 */
const SIZE_SYMBOLS = ['W_SINGLE', 'W_ONE_HALF', 'W_DOUBLE', 'mm']

const commonPath = join(root, 'scripts', 'lib', 'common.R')
if (!existsSync(commonPath)) {
  problems.push('找不到 scripts/lib/common.R —— 尺寸符号检查无法执行')
} else {
  const commonSrc = readFileSync(commonPath, 'utf8')
  for (const sym of SIZE_SYMBOLS) {
    const defined =
      new RegExp(`^\\s*${sym}\\s*<-`, 'm').test(commonSrc) ||
      new RegExp(`^\\s*${sym}\\s*<-\\s*function`, 'm').test(commonSrc)
    if (!defined) problems.push(`lib/common.R 里没有定义 ${sym}`)
  }

  // 拼错检测：扫所有 `W_XXX` 形状的标识符，逐个确认 common.R 里真有定义。
  // 只看名单是不够的 —— `W_DOUBL` 不在名单里，漏掉一个字母就溜过去了，
  // 而 R 要到运行时才报 "object not found"。
  const definedW = new Set()
  for (const m of commonSrc.matchAll(/^\s*(W_[A-Z0-9_]+)\s*<-/gm)) definedW.add(m[1])

  for (const file of files) {
    const rel = relative(root, file).split('\\').join('/')
    if (rel.endsWith('lib/common.R')) continue
    const body = readFileSync(file, 'utf8').replace(/#[^\n]*/g, '')
    for (const m of body.matchAll(/(?<![\w.$])(W_[A-Z0-9_]+)\b/g)) {
      if (!definedW.has(m[1])) {
        problems.push(
          `${rel} 用了 ${m[1]}，但 lib/common.R 里没有这个符号` +
          `（已定义: ${[...definedW].sort().join(', ')}）—— 拼错了？`
        )
      }
    }
  }

  for (const file of files) {
    const rel = relative(root, file).split('\\').join('/')
    if (rel.endsWith('lib/common.R')) continue
    const source = readFileSync(file, 'utf8')
    // 去掉注释，避免"名字只出现在注释里"的误报
    const body = source.replace(/#[^\n]*/g, '')
    const used = SIZE_SYMBOLS.filter(s => new RegExp(`(?<![\\w.$])${s}\\b`).test(body))
    if (used.length === 0) continue
    // 脚本走的是候选列表模式：`cand <- c(file.path(here, "lib", "common.R"), ...)`
    // 再 `source(hit)` —— 所以 `common.R` **不在** source() 的括号里。
    // 判据是"非注释正文里出现过 common.R"，而不是在 source(...) 里找。
    const referencesCommon = /common\.R/.test(body)
    if (!referencesCommon) {
      problems.push(
        `${rel} 用了 ${used.join(', ')} 但没有加载 lib/common.R —— ` +
        `单独 Rscript 跑到这里会 "object not found"`
      )
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
