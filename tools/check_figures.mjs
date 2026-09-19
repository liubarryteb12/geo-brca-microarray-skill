#!/usr/bin/env node
/**
 * 检查 PNG 图不是空白的。
 *
 * **为什么需要这个检查：** `main_analysis.R` 的验收只判断文件**存在**，
 * 不判断图**有没有画出来**。而出图代码有一类典型的静默失败：
 * 设备开了、绘图代码没跑、设备关了 —— 文件存在、大小正常、就是一张白图。
 *
 * 本仓库真踩过一次：`save_pdf()` 想用 `force(expr)` 跑两遍绘图代码
 * （一遍 PDF、一遍 PNG），但 R 的 promise 有记忆，第二次 force 直接返回缓存值，
 * 于是 PNG 全是空白，而 PDF 正常、日志无警告、验收全过。
 *
 * 所以这里独立解码 PNG 像素，统计"非白像素占比"和"颜色数"，任一过低就报错。
 * 纯 Node 实现（zlib 内置），不引入新依赖，和 tools/ 下其他检查保持一致。
 *
 * 用法：
 *   node tools/check_figures.mjs results
 *   node tools/check_figures.mjs results --min-ink 0.01 --min-colors 8
 */
import { readFileSync, readdirSync, statSync } from 'node:fs'
import { join } from 'node:path'
import { inflateSync } from 'node:zlib'

const args = process.argv.slice(2)
const dir = args.find(a => !a.startsWith('--')) ?? 'results'
const opt = (name, dflt) => {
  const i = args.indexOf(`--${name}`)
  return i >= 0 && args[i + 1] !== undefined ? Number(args[i + 1]) : dflt
}
const MIN_INK = opt('min-ink', 0.01)
const MIN_COLORS = opt('min-colors', 8)

/** 解码 8-bit 非隔行 PNG，返回 {width, height, channels, data} */
function decodePng(buf) {
  const SIG = [137, 80, 78, 71, 13, 10, 26, 10]
  for (let i = 0; i < 8; i++) {
    if (buf[i] !== SIG[i]) throw new Error('不是 PNG（签名不匹配）')
  }

  let pos = 8
  let ihdr = null
  const idat = []
  while (pos + 8 <= buf.length) {
    const len = buf.readUInt32BE(pos)
    const type = buf.toString('ascii', pos + 4, pos + 8)
    const data = buf.subarray(pos + 8, pos + 8 + len)
    if (type === 'IHDR') {
      ihdr = {
        width: data.readUInt32BE(0),
        height: data.readUInt32BE(4),
        bitDepth: data[8],
        colorType: data[9],
        interlace: data[12],
      }
    } else if (type === 'IDAT') {
      idat.push(data)
    } else if (type === 'IEND') {
      break
    }
    pos += 12 + len // length + type + data + crc
  }
  if (!ihdr) throw new Error('缺少 IHDR')

  const CH = { 0: 1, 2: 3, 3: 1, 4: 2, 6: 4 }
  const channels = CH[ihdr.colorType]
  if (channels === undefined) throw new Error(`不支持的颜色类型 ${ihdr.colorType}`)
  if (ihdr.bitDepth !== 8) throw new Error(`只支持 8-bit，实际 ${ihdr.bitDepth}-bit`)
  if (ihdr.interlace !== 0) throw new Error('不支持隔行扫描')
  if (ihdr.colorType === 3) throw new Error('不支持调色板 PNG')

  const { width, height } = ihdr
  const bpp = channels
  const stride = width * bpp
  const raw = inflateSync(Buffer.concat(idat))
  const out = Buffer.alloc(height * stride)

  const paeth = (a, b, c) => {
    const p = a + b - c
    const pa = Math.abs(p - a), pb = Math.abs(p - b), pc = Math.abs(p - c)
    return pa <= pb && pa <= pc ? a : pb <= pc ? b : c
  }

  let rp = 0
  for (let y = 0; y < height; y++) {
    const filter = raw[rp++]
    const line = raw.subarray(rp, rp + stride)
    rp += stride
    const cur = out.subarray(y * stride, (y + 1) * stride)
    const prev = y > 0 ? out.subarray((y - 1) * stride, y * stride) : null

    for (let x = 0; x < stride; x++) {
      const a = x >= bpp ? cur[x - bpp] : 0
      const b = prev ? prev[x] : 0
      const c = prev && x >= bpp ? prev[x - bpp] : 0
      let v = line[x]
      if (filter === 1) v += a
      else if (filter === 2) v += b
      else if (filter === 3) v += (a + b) >> 1
      else if (filter === 4) v += paeth(a, b, c)
      else if (filter !== 0) throw new Error(`未知的行过滤器 ${filter}`)
      cur[x] = v & 0xff
    }
  }
  return { width, height, channels, data: out }
}

/** 抽样统计墨迹占比与颜色数 */
function measure(png) {
  const { width, height, channels, data } = png
  const step = Math.max(1, Math.floor(Math.min(width, height) / 300))
  let n = 0, ink = 0
  const colors = new Set()
  for (let y = 0; y < height; y += step) {
    for (let x = 0; x < width; x += step) {
      const i = (y * width + x) * channels
      let r, g, b
      if (channels >= 3) { r = data[i]; g = data[i + 1]; b = data[i + 2] }
      else { r = g = b = data[i] }
      n++
      if (r < 245 || g < 245 || b < 245) ink++
      if (colors.size < 4096) colors.add((r << 16) | (g << 8) | b)
    }
  }
  return { ink: ink / n, colors: colors.size }
}

let files
try {
  // 排除带版本号的副本（`PPI_network__c370bd0.png`）—— 它们是同一张图的字节拷贝，
  // 再查一遍纯属浪费，而且会让"检查了 N 个 PNG"这个数字翻倍、掩盖真实图数。
  files = readdirSync(dir)
    .filter(f => f.toLowerCase().endsWith('.png') && !/__[0-9a-f]{7}\.png$/i.test(f))
    .sort()
} catch (e) {
  console.error(`无法读取目录 ${dir}: ${e.message}`)
  process.exit(1)
}

if (files.length === 0) {
  console.error(`${dir} 里没有 PNG 图`)
  process.exit(1)
}

console.log(`检查 ${files.length} 个 PNG（墨迹 >= ${MIN_INK}，颜色数 > ${MIN_COLORS}）\n`)
const bad = []
for (const f of files) {
  const path = join(dir, f)
  const size = statSync(path).size
  try {
    const png = decodePng(readFileSync(path))
    const { ink, colors } = measure(png)
    const ok = ink >= MIN_INK && colors > MIN_COLORS
    if (!ok) bad.push(`${f}（墨迹 ${(ink * 100).toFixed(2)}%，颜色 ${colors}）`)
    console.log(
      `  ${ok ? 'OK  ' : '空白'} ${f.padEnd(30)} ${png.width}x${png.height}` +
      `  墨迹 ${(ink * 100).toFixed(2).padStart(6)}%  颜色 ${String(colors).padStart(4)}  ${(size / 1024).toFixed(0)} KB`
    )
  } catch (e) {
    bad.push(`${f}（解码失败: ${e.message}）`)
    console.log(`  错误 ${f.padEnd(30)} ${e.message}`)
  }
}

console.log()
if (bad.length > 0) {
  console.error(`以下图疑似空白或损坏：\n  - ${bad.join('\n  - ')}`)
  console.error('\n图存在但没画出来，通常意味着绘图代码没被执行到（设备开了又关）。')
  process.exit(1)
}
console.log('全部图都有实际内容')
