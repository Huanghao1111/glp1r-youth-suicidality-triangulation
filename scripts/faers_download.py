# -*- coding: utf-8 -*-
# 下载 FAERS 季度 ASCII zip, 支持断点续传与大小校验
# 用法: python scripts/faers_download.py <start_idx> <end_idx>  (按 manifest 中 ok 行的 0 基索引)
import csv, os, subprocess, sys

os.makedirs("data/faers_zips", exist_ok=True)
rows = []
with open("data/faers_quarters_manifest.csv", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if r["status"] == "ok":
            rows.append(r)

start, end = int(sys.argv[1]), int(sys.argv[2])
todo = rows[start:end]
print(f"本批 {len(todo)} 个季度 (索引 {start}..{end-1})", flush=True)

for r in todo:
    name = os.path.basename(r["url"])
    dest = os.path.join("data/faers_zips", name)
    want = int(r["bytes"])
    if os.path.exists(dest) and os.path.getsize(dest) == want:
        print(f"跳过 {name} (已完整)", flush=True)
        continue
    cmd = ["curl", "-sS", "-C", "-", "--retry", "3", "--retry-delay", "3",
           "-o", dest, r["url"], "-w", "%{http_code} %{speed_download}"]
    p = subprocess.run(cmd, capture_output=True, text=True)
    got = os.path.getsize(dest) if os.path.exists(dest) else 0
    ok = got == want
    speed = p.stdout.split()[-1] if p.stdout else "?"
    print(f"{name}: {'OK' if ok else 'INCOMPLETE'} {got/1e6:.1f}/{want/1e6:.1f}MB 速度{float(speed)/1e6:.1f}MB/s", flush=True)
print("本批结束", flush=True)
