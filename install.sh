#!/usr/bin/env bash
# 把本 skill 安装到各 agent 平台的技能目录。
#
#   ./install.sh                 安装到默认目标（~/.agents/skills 等）
#   ./install.sh --target DIR    安装到指定目录
#   ./install.sh --copy          复制而不是软链接（Windows / 跨文件系统时用）
#   ./install.sh --uninstall     从所有目标移除
set -euo pipefail

SKILL_NAME="geo-brca-microarray-skill"
SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_TARGETS=("$HOME/.agents/skills" "$HOME/.claude/skills" "$HOME/.codex/skills")

MODE="link"
TARGETS=()
UNINSTALL=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --target) TARGETS+=("$2"); shift 2 ;;
    --copy) MODE="copy"; shift ;;
    --link) MODE="link"; shift ;;
    --uninstall) UNINSTALL=1; shift ;;
    -h|--help) sed -n '2,9p' "$0"; exit 0 ;;
    *) echo "未知参数: $1" >&2; exit 2 ;;
  esac
done
[[ ${#TARGETS[@]} -eq 0 ]] && TARGETS=("${DEFAULT_TARGETS[@]}")

install_one() {
  local target_root="$1"
  local dest="$target_root/$SKILL_NAME"
  mkdir -p "$target_root"

  if [[ -e "$dest" || -L "$dest" ]]; then
    echo "  已存在，先移除: $dest"
    rm -rf "$dest"
  fi

  if [[ "$MODE" == "copy" ]]; then
    mkdir -p "$dest"
    # 排除运行产物与版本控制目录
    tar -C "$SOURCE_DIR" \
      --exclude='./.git' --exclude='./results' --exclude='./data' \
      --exclude='./.Rproj.user' -cf - . | tar -C "$dest" -xf -
    echo "  已复制 -> $dest"
  else
    ln -s "$SOURCE_DIR" "$dest"
    echo "  已链接 -> $dest"
  fi
}

uninstall_one() {
  local dest="$1/$SKILL_NAME"
  if [[ -e "$dest" || -L "$dest" ]]; then
    rm -rf "$dest"
    echo "  已移除 $dest"
  fi
}

echo "skill: $SKILL_NAME"
echo "源:    $SOURCE_DIR"
echo

for root in "${TARGETS[@]}"; do
  if [[ "$UNINSTALL" == "1" ]]; then
    echo "卸载自 $root"
    uninstall_one "$root"
  else
    echo "安装到 $root"
    install_one "$root"
  fi
done

echo
if [[ "$UNINSTALL" != "1" ]]; then
  echo "完成。验证:"
  echo "  node \"$SOURCE_DIR/scripts/find_dataset.mjs\" check GSE92252"
fi
