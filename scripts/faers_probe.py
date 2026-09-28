# -*- coding: utf-8 -*-
# 探测 FAERS 季度 ASCII 文件是否存在及大小, 生成下载清单
import subprocess, csv, re, sys

def probe(url):
    # 用 range 请求拿 Content-Range: bytes 0-0/TOTAL
    p = subprocess.run(["curl", "-s", "--max-time", "25", "-r", "0-0", "-o", "NUL",
                        "-w", "%{http_code}", url], capture_output=True, text=True)
    code = p.stdout.strip()
    if code != "206":
        return None
    p2 = subprocess.run(["curl", "-sI", "--max-time", "25", "-r", "0-0", url],
                        capture_output=True, text=True)
    m = re.search(r"content-range:\s*bytes 0-0/(\d+)", p2.stdout, re.I)
    if m:
        return int(m.group(1))
    # HEAD 不带 range 时拿不到长度, 用 range 响应头再试一次
    p3 = subprocess.run(["curl", "-s", "--max-time", "25", "-r", "0-0", "-o", "NUL",
                         "-D", "-", url], capture_output=True, text=True)
    m = re.search(r"content-range:\s*bytes 0-0/(\d+)", p3.stdout, re.I)
    return int(m.group(1)) if m else -1

rows = []
for year in range(2004, 2027):
    for q in range(1, 5):
        url = f"https://fis.fda.gov/content/Exports/faers_ascii_{year}q{q}.zip"
        size = probe(url)
        status = "ok" if size else "missing"
        rows.append((f"{year}Q{q}", url, size or 0, status))
        print(f"{year}Q{q}: {status} {size/1e6:.1f}MB" if size else f"{year}Q{q}: {status}", flush=True)

with open("data/faers_quarters_manifest.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["quarter", "url", "bytes", "status"])
    w.writerows(rows)
ok = [r for r in rows if r[3] == "ok"]
print(f"\n共 {len(ok)} 个季度可下载, 总计 {sum(r[2] for r in ok)/1e9:.1f} GB")
