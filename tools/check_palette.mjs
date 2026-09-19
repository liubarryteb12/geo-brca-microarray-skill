#!/usr/bin/env node
/**
 * 调色板门禁：从 common.R 里解析出实际色值，按 better-colors 的判据验证。
 *
 * **为什么要有这道门禁：** "一个颜色一个含义"是条不变量，但没有任何东西在守着它。
 * 实测这套色值之前的状态：p 值直方图的 firebrick 与 up 红**只差 0.4°**
 * （按 15° 判据就是同一个颜色）；分类色板在 protanopia 下 up 红与 cat brown
 * 距离 0.002（几乎完全重合）；magma 序列色中段距 up 红仅 2.9°。
 * 这些都是凭眼睛看不出来的，只有算才看得见。
 *
 * 判据：
 *   1. 色相相差 < 15° 视为同一颜色 —— 承载不同含义的颜色必须拉开 15° 以上
 *   2. 三种色盲下 OKLab 距离 >= 0.05
 *   3. 序列色板亮度单调、色相跨度小
 *   4. 承载语义的颜色白底对比度 >= 2.0
 *
 * 用法：node tools/check_palette.mjs
 */
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const here = dirname(fileURLToPath(import.meta.url));
const commonPath = join(here, "..", "scripts", "lib", "common.R");

// ---- 色彩空间 -------------------------------------------------------------
const srgbToLinear = (c) => (c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4);
const linearToSrgb = (c) => {
  const v = Math.min(1, Math.max(0, c));
  return v <= 0.0031308 ? v * 12.92 : 1.055 * v ** (1 / 2.4) - 0.055;
};

function hexToRgb(h) {
  const s = h.replace("#", "");
  return [0, 2, 4].map((i) => parseInt(s.slice(i, i + 2), 16) / 255);
}

const M1 = [
  [0.4122214708, 0.5363325363, 0.0514459929],
  [0.2119034982, 0.6806995451, 0.1073969566],
  [0.0883024619, 0.2817188376, 0.6299787005],
];
const M2 = [
  [0.2104542553, 0.793617785, -0.0040720468],
  [1.9779984951, -2.428592205, 0.4505937099],
  [0.0259040371, 0.7827717662, -0.808675766],
];
const mul = (M, v) => M.map((row) => row.reduce((s, x, i) => s + x * v[i], 0));

function oklab(rgb) {
  const lin = rgb.map(srgbToLinear);
  return mul(M2, mul(M1, lin).map(Math.cbrt));
}

function oklch(rgb) {
  const [L, a, b] = oklab(rgb);
  return { L, C: Math.hypot(a, b), h: ((Math.atan2(b, a) * 180) / Math.PI + 360) % 360 };
}

const CVD = {
  protanopia: [
    [0.152286, 1.052583, -0.204868],
    [0.114503, 0.786281, 0.099216],
    [-0.003882, -0.048116, 1.051998],
  ],
  deuteranopia: [
    [0.367322, 0.860646, -0.227968],
    [0.280085, 0.672501, 0.047413],
    [-0.01182, 0.04294, 0.968881],
  ],
  tritanopia: [
    [1.255528, -0.076749, -0.178779],
    [-0.078411, 0.930809, 0.147602],
    [0.004733, 0.691367, 0.3039],
  ],
};

const simulate = (rgb, kind) => mul(CVD[kind], rgb.map(srgbToLinear)).map(linearToSrgb);
const dist = (a, b) => {
  const [x, y, z] = oklab(a);
  const [p, q, r] = oklab(b);
  return Math.hypot(x - p, y - q, z - r);
};
const hueGap = (a, b) => {
  const d = Math.abs(a - b) % 360;
  return Math.min(d, 360 - d);
};

function relLum(rgb) {
  const [r, g, b] = rgb.map(srgbToLinear);
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}
function contrast(fg, bg) {
  const l1 = relLum(fg);
  const l2 = relLum(bg);
  return (Math.max(l1, l2) + 0.05) / (Math.min(l1, l2) + 0.05);
}

// ---- 从 common.R 解析色值 --------------------------------------------------
const src = readFileSync(commonPath, "utf8");
const grab = (name) => {
  const m = src.match(new RegExp(`${name}\\s*=\\s*"(#[0-9A-Fa-f]{6})"`));
  if (!m) throw new Error(`common.R 里找不到 PAL$${name}`);
  return m[1];
};

const SEMANTIC = {
  "up / tumor": grab("up"),
  "down / normal": grab("down"),
  nominal: grab("nominal"),
  ns: grab("ns"),
  ink: grab("ink"),
};

// 分类色板与序列色板是向量/函数体，单独抓
const catMatch = src.match(/base\s*<-\s*c\(([^)]*#[0-9A-Fa-f]{6}[^)]*)\)/);
if (!catMatch) throw new Error("common.R 里找不到 pal_categorical 的 base 向量");
const MODULES = catMatch[1]
  .match(/#[0-9A-Fa-f]{6}/g)
  .reduce((acc, hx, i) => ({ ...acc, [`module ${i + 1}`]: hx }), {});

const seqMatch = src.match(/pal_sequential[\s\S]*?colorRampPalette\(\s*c\(([\s\S]*?)\)\)\(n\)/);
if (!seqMatch) throw new Error("common.R 里找不到 pal_sequential 的 ramp");
const SEQ = seqMatch[1].match(/#[0-9A-Fa-f]{6}/g);

const WHITE = hexToRgb("#FFFFFF");
const MIN_CONTRAST = 2.0;
const MIN_CVD = 0.05;
const MIN_HUE = 15;

const problems = [];
const rows = [];

// ---- 判据 1 + 2：语义色与模块色一起，两两检查 ------------------------------
const ALL = { ...SEMANTIC, ...MODULES };
const names = Object.keys(ALL);
for (let i = 0; i < names.length; i++) {
  for (let j = i + 1; j < names.length; j++) {
    const a = names[i];
    const b = names[j];
    const ca = oklch(hexToRgb(ALL[a]));
    const cb = oklch(hexToRgb(ALL[b]));
    // 近灰色不承载色相语义，跳过色相检查
    if (ca.C >= 0.04 && cb.C >= 0.04) {
      const g = hueGap(ca.h, cb.h);
      if (g < MIN_HUE) {
        problems.push(
          `色相 ${g.toFixed(1)}° < ${MIN_HUE}°：${a} ${ALL[a]} 与 ${b} ${ALL[b]} 是同一颜色`
        );
      }
    }
    for (const kind of Object.keys(CVD)) {
      const d = dist(simulate(hexToRgb(ALL[a]), kind), simulate(hexToRgb(ALL[b]), kind));
      if (d < MIN_CVD) {
        problems.push(
          `${kind} 下 Δ=${d.toFixed(3)} < ${MIN_CVD}：${a} 与 ${b} 无法分辨`
        );
      }
    }
  }
}

// ---- 判据 4：对比度 --------------------------------------------------------
for (const [n, hx] of Object.entries(ALL)) {
  const c = contrast(hexToRgb(hx), WHITE);
  rows.push(`  ${n.padEnd(16)} ${hx}  ${c.toFixed(2)}:1`);
  if (c < MIN_CONTRAST && n !== "ns") {
    problems.push(`白底对比度 ${c.toFixed(2)}:1 < ${MIN_CONTRAST}：${n} ${hx} 作为实心圆点看不清`);
  }
}

// ---- 判据 3：序列色板 ------------------------------------------------------
const seqL = SEQ.map((h) => oklch(hexToRgb(h)).L);
const seqH = SEQ.map((h) => oklch(hexToRgb(h)).h);
const seqC = SEQ.map((h) => oklch(hexToRgb(h)).C);
const monotonic = seqL.every((v, i) => i === 0 || v < seqL[i - 1]);
if (!monotonic) problems.push("序列色板亮度不是单调递减");

const hueSpan = Math.max(...seqH.map((h) => hueGap(h, seqH[0])));
if (hueSpan > 30) {
  problems.push(`序列色板色相跨度 ${hueSpan.toFixed(1)}° > 30°，不再是"保持色相"的 ramp`);
}
for (const [lbl, target] of [["up 红", oklch(hexToRgb(SEMANTIC["up / tumor"])).h],
                             ["down 蓝", oklch(hexToRgb(SEMANTIC["down / normal"])).h]]) {
  const worst = Math.min(...seqH.map((h) => hueGap(h, target)));
  if (worst < MIN_HUE) {
    problems.push(`序列色板与 ${lbl} 色相只差 ${worst.toFixed(1)}°，深色会被误读`);
  }
}

// ---- 报告 ------------------------------------------------------------------
console.log(`检查 ${commonPath}`);
console.log(`语义色 ${Object.keys(SEMANTIC).length} 个 / 模块色 ${Object.keys(MODULES).length} 个 / 序列色 ${SEQ.length} 档`);
console.log(`\n白底对比度：\n${rows.join("\n")}`);
console.log(
  `\n序列色板：亮度单调=${monotonic}  色相跨度=${hueSpan.toFixed(1)}°  ` +
    `彩度峰值=${seqC.indexOf(Math.max(...seqC))}/${SEQ.length - 1}`
);
console.log(
  `两两检查：${names.length} 个颜色，${(names.length * (names.length - 1)) / 2} 对，` +
    `判据 色相>=${MIN_HUE}° / 色盲 Δ>=${MIN_CVD}`
);

if (problems.length) {
  console.log(`\n未通过 ${problems.length} 项：`);
  for (const p of problems) console.log("  ✗ " + p);
  process.exit(1);
}
console.log("\n调色板检查通过。");
