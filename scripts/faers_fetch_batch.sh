#!/bin/bash
# 并行补齐 FAERS 季度 zip: 每次调用处理接下来 K 个未完整文件
K=${1:-8}
cd "$(dirname "$0")/../data/faers_zips" || exit 1
mapfile -t files < <(python - <<'EOF' | tr -d '\r'
import csv, os
rows=[]
with open("../faers_quarters_manifest.csv", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if r["status"]!="ok": continue
        name=os.path.basename(r["url"]); want=int(r["bytes"])
        if not (os.path.exists(name) and os.path.getsize(name)==want):
            rows.append((name,r["url"]))
for n,u in rows: print(f"{n} {u}")
EOF
)
total=${#files[@]}
echo "未完成文件数: $total"
[ "$total" -eq 0 ] && exit 0
batch=("${files[@]:0:$K}")
echo "本批并行下载 ${#batch[@]} 个"
start=$(date +%s)
pids=()
for fu in "${batch[@]}"; do
  f=${fu%% *}; u=${fu#* }
  curl -sS -C - --retry 3 --retry-delay 3 -o "$f" "$u" &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
end=$(date +%s)
echo "耗时 $((end-start)) 秒"
for fu in "${batch[@]}"; do
  f=${fu%% *}
  [ -f "$f" ] && echo "$f $(stat -c%s "$f")"
done
