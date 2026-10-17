#!/bin/bash
# InoMiner 36 株回归测试：新 CLI 全流程 vs 旧 shell 流程参考输出（必须完全一致）
# 用法: bash tests/regression_36.sh
#
# 已知参考数据过时点（不影响结论）：
#  1. SKS4.fa 于 8/31 更新，而 01_FPphi_mining/outputs/dif_motif 的 SKS4 参考生成于 8/27
#     （contig 编号变化）。已单独验证：旧 awk 脚本跑当前 SKS4.fa 与新 Python 输出逐字节一致，
#     故 dif 对比跳过 SKS4。
#  2. 05_29410Vibrio/test 的 pHVR 表生成于 9/19 11:47，此后脚本把表头
#     Bakta_Gene_Annotation 改名为 Gene_Annotation。对比时统一表头，数据须完全一致。
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
GNM="$ROOT/00_36Vibrio/data/gnms"
REF="$ROOT/05_29410Vibrio/test"
REF_DIF="$ROOT/01_FPphi_mining/outputs/dif_motif"
OUT="$HERE/regression_out"
INOMINE="$HERE/../.venv/bin/inomine"

PRODIGAL=/opt/miniconda3/envs/prodigal/bin/prodigal
HMMSEARCH=/opt/miniconda3/envs/hmmer/bin/hmmsearch

rm -rf "$OUT"
"$INOMINE" run -i "$GNM" -o "$OUT" -t 8 --prodigal "$PRODIGAL" --hmmsearch "$HMMSEARCH"

fail=0
echo
echo "== tables 对比 =="
for f in FPphi_summary.tsv FPphi_gene_annotation.tsv FPphi_pHVR_annotation.tsv FPphi_pHVR.faa; do
  if [ "$f" = "FPphi_pHVR_annotation.tsv" ]; then
    d=$(diff <(sed '1s/Bakta_Gene_Annotation/Gene_Annotation/' "$REF/tables/$f") "$OUT/tables/$f")
  else
    d=$(diff "$REF/tables/$f" "$OUT/tables/$f")
  fi
  if [ -z "$d" ]; then
    echo "OK   $f"
  else
    echo "DIFF $f"; echo "$d" | head -15; fail=1
  fi
done

echo
echo "== 序列对比 =="
if diff -rq "$REF/results/FPphi" "$OUT/results/FPphi" >/dev/null 2>&1; then
  echo "OK   results/FPphi"
else
  echo "DIFF results/FPphi"; diff -rq "$REF/results/FPphi" "$OUT/results/FPphi" | head -15; fail=1
fi

echo
echo "== dif 对比（SKS4 参考过时，跳过）=="
nd=0; ndiff=0
for f in "$REF_DIF"/*-find_Vibrio_dif_motif-output.csv; do
  b=$(basename "$f")
  [ "$b" = "SKS4-find_Vibrio_dif_motif-output.csv" ] && continue
  nd=$((nd + 1))
  if ! diff -q "$f" "$OUT/dif/$b" >/dev/null 2>&1; then
    echo "DIFF dif/$b"; ndiff=$((ndiff + 1)); fail=1
  fi
done
[ "$ndiff" -eq 0 ] && echo "OK   $nd/$nd dif CSV"

nnew=$(ls "$OUT/dif"/*-find_Vibrio_dif_motif-output.csv 2>/dev/null | wc -l | tr -d ' ')
if [ "$nnew" -eq 36 ]; then
  echo "OK   dif 文件数 36"
else
  echo "DIFF dif 文件数 $nnew != 36"; fail=1
fi

echo
if [ $fail -eq 0 ]; then
  echo "== 全部一致，回归通过 =="
else
  echo "== 存在差异，回归失败 =="
fi
exit $fail
