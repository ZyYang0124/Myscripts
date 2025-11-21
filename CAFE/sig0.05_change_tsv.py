#sig0.05_change_tsv.py
#以 Gamma_family_results.txt 里显著家族为准，从 Gamma_change.tab 中去掉不显著家族
#输出一个新的 Gamma_change_sig0.05.tsv

#!/usr/bin/env python3
import pandas as pd

# -----------------------------
# 1. 读取显著家族列表
# -----------------------------
sig_fams = set()
with open("Gamma_family_results.txt") as f:
    next(f)  # 跳过标题
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3 and parts[2].lower() == "y":
            sig_fams.add(parts[0])

print(f"显著家族数: {len(sig_fams)}")

# -----------------------------
# 2. 读取 Gamma_change.tab
# -----------------------------
df = pd.read_csv("Gamma_change.tab", sep="\t")

# -----------------------------
# 3. 筛选显著家族
# -----------------------------
df_sig = df[df["FamilyID"].isin(sig_fams)]

# -----------------------------
# 4. 输出新的显著家族矩阵
# -----------------------------
output_file = "Gamma_change_sig0.05.tsv"
df_sig.to_csv(output_file, sep="\t", index=False)

print(f"完成，已生成 {output_file}，仅包含显著家族")