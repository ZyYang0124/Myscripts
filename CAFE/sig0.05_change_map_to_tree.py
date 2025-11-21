#sig0.05_change_map_to_tree.py
#从 Gamma_family_results.txt 读入显著家族只保留 “y” 的基因家族
#从 Gamma_change.tab 选取显著家族对应的行
#每个节点分别统计所有显著家族的扩张数和收缩数
#最后将其 map 到树上：写入 cleaned_tree_sig0.05_only.txt

#!/usr/bin/env python3
import re
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
# 2. 读取 CAFE family × node 变化矩阵
# -----------------------------
df = pd.read_csv("Gamma_change.tab", sep="\t")

node_cols = df.columns[1:]   # 第一列是 FamilyID

# 仅显著家族
df_sig = df[df["FamilyID"].isin(sig_fams)]
print(f"显著家族矩阵形状: {df_sig.shape}")

# -----------------------------
# 3. 统计显著扩张/收缩的“家族数量”
# -----------------------------
node_change = {}

for node in node_cols:
    changes = df_sig[node]

    inc = (changes > 0).sum()   # 扩张家族数量
    dec = (changes < 0).sum()   # 收缩家族数量

    node_change[node] = (int(inc), int(dec))

print("每个节点显著扩张/收缩数量统计完毕。")

# -----------------------------
# 4. 读取树
# -----------------------------
with open("cleaned_tree.txt") as f:
    tree = f.read()

# -----------------------------
# 5. 替换树中的节点名称
# -----------------------------
for node, (inc, dec) in node_change.items():

    if inc == 0 and dec == 0:
        continue    # 两者都不显著则跳过

    new_label = node
    if inc > 0:
        new_label += f"+{inc}"
    if dec > 0:
        new_label += f"-{dec}"

    tree = re.sub(re.escape(node), new_label, tree)

# -----------------------------
# 6. 输出
# -----------------------------
with open("cleaned_tree_sig0.05_only.txt", "w") as f:
    f.write(tree)

print("写入完成：cleaned_tree_sig0.05_only.txt")