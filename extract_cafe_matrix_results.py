#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_cafe_matrix_results_full.py

功能：
  从 CAFE 输出的变化矩阵（Base_change.table）和基因家族拷贝数矩阵（Orthogroups.GeneCount.tsv）
  中筛选目标物种扩张/收缩最显著的前 N 个基因家族，并整合 Orthogroups 与注释信息。
  支持注释来源：annotation table (.tsv)、protein fasta (.faa/.fa/.fasta)、或 GFF3 文件。

主要输出：
  1) <taxon>_topN_expanded.tsv   — 仅扩张家族（若 mode 包含 expand）
  2) <taxon>_topN_contracted.tsv  — 仅收缩家族（若 mode 包含 contract）
  3) <taxon>_topN_merged.tsv      — 合并的详细结果（每个基因一行，包含注释等）

作者：ChatGPT
日期：2025-11-04
"""

import argparse
import pandas as pd
from pathlib import Path
import re
import sys
import os

# -------------------------
# 参数解析（命令行）
# -------------------------
def parse_args():
    p = argparse.ArgumentParser(description="从 CAFE matrix 与 Base_change 中提取 top 扩张/收缩家族并整合注释")
    p.add_argument("--base_change", "-b", required=True, help="CAFE 输出变化矩阵（例如 Base_change.table）")
    p.add_argument("--matrix", "-m", required=True, help="家族拷贝数矩阵（例如 Orthogroups.GeneCount.tsv）")
    p.add_argument("--orthogroups", "-o", required=True, help="Orthogroups.tsv 文件（家族->基因 列表）")
    p.add_argument("--annotation", "-a", required=False, help="注释文件：family->gene->annotation 表 或 protein .faa 或 GFF3")
    p.add_argument("--gffdir", "-g", required=False, help="（可选）包含多个 gff3 的目录，用于按物种提取注释")
    p.add_argument("--taxon", "-t", required=True, help="目标物种名称（在表头中匹配，忽略大小写），如 Siler_cupreus")
    p.add_argument("--topn", "-n", type=int, default=20, help="Top N（默认 20）")
    p.add_argument("--mode", choices=["expand", "contract", "both"], default="both",
                   help="筛选模式：expand / contract / both（默认 both）")
    p.add_argument("--outdir", "-d", required=True, help="输出目录")
    return p.parse_args()

# -------------------------
# 工具函数：文件存在性检查
# -------------------------
def check_file(path, name):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"找不到 {name}：{path}")
    return p

# -------------------------
# Step1: 从 Base_change 表筛选 top family
# 注：Base_change 文件包含 FamilyID 列与每个物种的变化量列（正为扩张，负为收缩）
# -------------------------
def select_top_families(base_change_path, taxon, topn, mode):
    df = pd.read_csv(base_change_path, sep="\t", dtype=str)
    # 修正空白列名与空格
    df.columns = [c.strip() for c in df.columns]

    # 自动查找 taxon 对应列（忽略大小写、尖括号中的编号）
    match_cols = [c for c in df.columns if re.search(re.escape(taxon), c, re.IGNORECASE)]
    if not match_cols:
        # 进一步尝试用taxon中空格或下划线替换形式匹配
        alt = taxon.replace(" ", "_")
        match_cols = [c for c in df.columns if re.search(re.escape(alt), c, re.IGNORECASE)]
    if not match_cols:
        raise ValueError(f"未在 Base_change 文件中找到与目标物种匹配的列（taxon={taxon}）。可用列示例：{df.columns[:20].tolist()}")
    taxon_col = match_cols[0]

    # 将变化列转换为数值（可能包含非数字字符，errors='coerce' 将无法解析的转为 NaN）
    df["ChangeValue"] = pd.to_numeric(df[taxon_col], errors="coerce")

    # 保证 FamilyID 列存在（不同 CAFE 输出可能名为 FamilyID 或 Orthogroup）
    if "FamilyID" not in df.columns and "Orthogroup" in df.columns:
        df = df.rename(columns={"Orthogroup": "FamilyID"})
    if "FamilyID" not in df.columns:
        # 退而求其次取第一列为 ID
        df = df.rename(columns={df.columns[0]: "FamilyID"})

    # 排序并选 TopN
    df_sorted_desc = df.sort_values(by="ChangeValue", ascending=False).reset_index(drop=True)
    df_sorted_asc = df.sort_values(by="ChangeValue", ascending=True).reset_index(drop=True)

    expand_top = df_sorted_desc[df_sorted_desc["ChangeValue"] > 0].head(topn)
    contract_top = df_sorted_asc[df_sorted_asc["ChangeValue"] < 0].head(topn)

    if mode == "expand":
        selected = expand_top.copy()
    elif mode == "contract":
        selected = contract_top.copy()
    else:  # both
        selected = pd.concat([expand_top, contract_top], ignore_index=True)

    # 添加 Direction 字段（Expanded / Contracted / Stable）
    selected["Direction"] = selected["ChangeValue"].apply(lambda x: "Expanded" if pd.notna(x) and x > 0 else ("Contracted" if pd.notna(x) and x < 0 else "Stable"))
    selected = selected[["FamilyID", "ChangeValue", "Direction"]]
    return selected, taxon_col

# -------------------------
# Step2: 从 matrix（Orthogroups.GeneCount.tsv）中提取 family size（拷贝数）
# -------------------------
def read_family_sizes(matrix_path, taxon):
    df = pd.read_csv(matrix_path, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    # 匹配 taxon 列
    match_cols = [c for c in df.columns if re.search(re.escape(taxon), c, re.IGNORECASE)]
    if not match_cols:
        raise ValueError(f"未在 matrix 文件中找到与 {taxon} 匹配的列。可用列示例：{df.columns[:20].tolist()}")
    taxon_col = match_cols[0]
    # 统一列名
    if "FamilyID" not in df.columns and "Orthogroup" in df.columns:
        df = df.rename(columns={"Orthogroup": "FamilyID"})
    if "FamilyID" not in df.columns:
        df = df.rename(columns={df.columns[0]: "FamilyID"})
    sizes = df[["FamilyID", taxon_col]].rename(columns={taxon_col: "FamilySize"})
    sizes["FamilySize"] = pd.to_numeric(sizes["FamilySize"], errors="coerce")
    return sizes

# -------------------------
# Step3: 从 Orthogroups.tsv 中提取目标物种的基因ID
# -------------------------
def extract_genes_from_orthogroups(orthogroups_path, family_list, taxon):
    og = pd.read_csv(orthogroups_path, sep="\t", dtype=str)
    og.columns = [c.strip() for c in og.columns]
    if "Orthogroup" not in og.columns and "FamilyID" in og.columns:
        og = og.rename(columns={"FamilyID": "Orthogroup"})
    if "Orthogroup" not in og.columns:
        og = og.rename(columns={og.columns[0]: "Orthogroup"})

    og_sub = og[og["Orthogroup"].isin(family_list)].copy()

    # 匹配 taxon 列名（or species column）
    match_cols = [c for c in og_sub.columns if re.search(re.escape(taxon), c, re.IGNORECASE)]
    if not match_cols:
        raise ValueError(f"未在 Orthogroups.tsv 中找到 {taxon} 列。可用列示例：{og_sub.columns.tolist()[:20]}")
    taxon_col = match_cols[0]

    rows = []
    for _, row in og_sub.iterrows():
        gene_field = str(row[taxon_col]) if pd.notna(row[taxon_col]) else ""
        # 常见分隔符有 ", " 或 ";" 或 空格
        genes = re.split(r"[,;]\s*|\s+", gene_field.strip()) if gene_field.strip() else []
        for g in genes:
            if g and g.lower() != "nan":
                rows.append({"FamilyID": row["Orthogroup"], "GeneID": g})
    genes_df = pd.DataFrame(rows)
    return genes_df

# -------------------------
# Step4: 注释提取：支持 annotation table / .faa / .gff3 或 gffdir
# annotation 表格式（如果提供）：FamilyID<TAB>GeneID<TAB>Annotation  或 GeneID<TAB>Annotation
# -------------------------
def extract_annotations(annotation_path, gene_list):
    if not annotation_path:
        return pd.DataFrame(columns=["GeneID", "Annotation"])

    p = Path(annotation_path)
    if not p.exists():
        raise FileNotFoundError(f"注释文件不存在：{annotation_path}")

    suffix = p.suffix.lower()
    if suffix in [".tsv", ".txt", ".csv"]:
        # 可能是 family->gene->annotation 或 gene->annotation
        df = pd.read_csv(p, sep="\t", header=None, dtype=str)
        df.columns = [c.strip() for c in df.columns]
        # 尝试自动识别列数
        if df.shape[1] >= 3:
            # 取最后两列为 GeneID, Annotation 的可能性较大
            # 尝试匹配列名是否包含 Gene / ID / Annotation
            # 保守：假设格式 FamilyID GeneID Annotation
            df = df.rename(columns={df.columns[1]: "GeneID", df.columns[2]: "Annotation"})
            ann = df[["GeneID", "Annotation"]].dropna()
        elif df.shape[1] == 2:
            # 两列：GeneID Annotation
            df = df.rename(columns={df.columns[0]: "GeneID", df.columns[1]: "Annotation"})
            ann = df[["GeneID", "Annotation"]].dropna()
        else:
            ann = pd.DataFrame(columns=["GeneID", "Annotation"])
        ann = ann[ann["GeneID"].isin(gene_list)].drop_duplicates()
        return ann

    elif suffix in [".faa", ".fa", ".fasta"]:
        # 解析 fasta header：">geneID description..."
        ann = []
        with open(p, "r") as fh:
            cur_id = None
            cur_desc = ""
            for line in fh:
                if line.startswith(">"):
                    if cur_id:
                        ann.append((cur_id, cur_desc.strip()))
                    header = line[1:].strip()
                    parts = header.split(None, 1)
                    cur_id = parts[0]
                    cur_desc = parts[1] if len(parts) > 1 else ""
            if cur_id:
                ann.append((cur_id, cur_desc.strip()))
        ann_df = pd.DataFrame(ann, columns=["GeneID", "Annotation"])
        ann_df = ann_df[ann_df["GeneID"].isin(gene_list)].drop_duplicates()
        return ann_df

    elif suffix in [".gff3", ".gff"]:
        # 解析 GFF3 attributes，尝试提取 ID 与 product/Name/Note 等
        ann = []
        pattern = re.compile(r"ID=([^;]+).*?(?:product|Name|Note|gene)=([^;]+)", re.IGNORECASE)
        with open(p, "r") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 9:
                    continue
                attr = parts[8]
                m = pattern.search(attr)
                if m:
                    gid = m.group(1)
                    desc = m.group(2)
                    if gid in gene_list:
                        ann.append((gid, desc))
        ann_df = pd.DataFrame(ann, columns=["GeneID", "Annotation"]).drop_duplicates()
        return ann_df

    else:
        raise ValueError("不支持的注释文件类型（仅支持 .tsv/.txt/.csv/.faa/.fa/.fasta/.gff3/.gff）")

# -------------------------
# Step5: 如果提供 gffdir（目录），尝试在目录中按物种名匹配 gff 文件并从中提取注释
# -------------------------
def extract_annotations_from_dir(gffdir, taxon, gene_list):
    pdir = Path(gffdir)
    if not pdir.exists() or not pdir.is_dir():
        raise FileNotFoundError(f"GFF3 目录不存在或非目录：{gffdir}")
    # 尝试匹配带有 taxon 名称的文件
    candidates = [f for f in pdir.iterdir() if f.is_file() and re.search(re.escape(taxon), f.name, re.IGNORECASE)]
    ann_df_total = pd.DataFrame(columns=["GeneID", "Annotation"])
    for f in candidates:
        try:
            ann_df = extract_annotations(str(f), gene_list)
            ann_df_total = pd.concat([ann_df_total, ann_df], ignore_index=True)
        except Exception:
            continue
    ann_df_total = ann_df_total.drop_duplicates()
    return ann_df_total

# -------------------------
# 主流程
# -------------------------
def main():
    args = parse_args()

    # 检查输入文件
    base_change_p = check_file(args.base_change, "Base_change")
    matrix_p = check_file(args.matrix, "matrix")
    orthogroups_p = check_file(args.orthogroups, "Orthogroups")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. 从 Base_change 筛选 top 家族
    selected_df, taxon_col = select_top_families(str(base_change_p), args.taxon, args.topn, args.mode)
    selected_families = selected_df["FamilyID"].tolist()
    print(f"[INFO] 筛选到 {len(selected_families)} 个家族（含扩张与收缩，根据 mode）")

    # 2. 读取 family 拷贝数（matrix）
    sizes_df = read_family_sizes(str(matrix_p), args.taxon)

    # 3. 从 Orthogroups 提取目标物种对应的基因
    genes_df = extract_genes_from_orthogroups(str(orthogroups_p), selected_families, args.taxon)
    if genes_df.empty:
        print("[WARN] 未从 Orthogroups 中提取到任何基因（可能是命名不一致）。")
    else:
        print(f"[INFO] 从 Orthogroups 中提取到 {len(genes_df)} 条基因记录")

    # 4. 注释提取：优先使用 --annotation，其次使用 --gffdir；若都无则留空
    gene_list = genes_df["GeneID"].unique().tolist() if not genes_df.empty else []
    ann_df = pd.DataFrame(columns=["GeneID", "Annotation"])
    if args.annotation:
        try:
            ann_df = extract_annotations(args.annotation, gene_list)
            print(f"[INFO] 从注释文件提取到 {len(ann_df)} 条注释")
        except Exception as e:
            print(f"[WARN] 从注释文件提取注释失败：{e}")
    if ann_df.empty and args.gffdir:
        try:
            ann_df = extract_annotations_from_dir(args.gffdir, args.taxon, gene_list)
            print(f"[INFO] 从 gff 目录提取到 {len(ann_df)} 条注释")
        except Exception as e:
            print(f"[WARN] 从 gff 目录提取注释失败：{e}")

    # 5. 合并：Family + ChangeValue + Direction + FamilySize + GeneID + Annotation
    # 先将 selected_df 与 sizes_df 合并（按 FamilyID）
    merged_families = selected_df.merge(sizes_df, left_on="FamilyID", right_on="FamilyID", how="left")
    merged_families = merged_families.rename(columns={"ChangeValue": "DeltaCopy", "FamilySize": "CopyNumber"})
    # 输出 expanded / contracted summary 列表（只 family 层面）
    expanded_summary = merged_families[merged_families["Direction"] == "Expanded"].sort_values(by="DeltaCopy", ascending=False)
    contracted_summary = merged_families[merged_families["Direction"] == "Contracted"].sort_values(by="DeltaCopy", ascending=True)

    # 保存 family 层面输出
    if not expanded_summary.empty and args.mode in ["expand", "both"]:
        expanded_file = outdir / f"{args.taxon}_top{args.topn}_expanded_families.tsv"
        expanded_summary.to_csv(expanded_file, sep="\t", index=False)
        print(f"[OUT] 扩张家族 summary 保存为：{expanded_file}")
    if not contracted_summary.empty and args.mode in ["contract", "both"]:
        contracted_file = outdir / f"{args.taxon}_top{args.topn}_contracted_families.tsv"
        contracted_summary.to_csv(contracted_file, sep="\t", index=False)
        print(f"[OUT] 收缩家族 summary 保存为：{contracted_file}")

    # gene 层面合并（每个基因一行）
    if not genes_df.empty:
        gene_merged = genes_df.merge(merged_families[["FamilyID", "DeltaCopy", "CopyNumber", "Direction"]],
                                     left_on="FamilyID", right_on="FamilyID", how="left")
        if not ann_df.empty:
            gene_merged = gene_merged.merge(ann_df, left_on="GeneID", right_on="GeneID", how="left")
        # 输出完整的 merged 表
        merged_outfile = outdir / f"{args.taxon}_top{args.topn}_merged_genes.tsv"
        gene_merged.to_csv(merged_outfile, sep="\t", index=False)
        print(f"[OUT] 合并基因注释表保存为：{merged_outfile}")
    else:
        print("[INFO] 未生成基因层面的合并表（因为没有提取到基因）。")

    print("[DONE] 处理完成。")

if __name__ == "__main__":
    main()
