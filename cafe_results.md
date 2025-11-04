
# 基于 Cafe 5 

## 1. cafe 输出

Base_change.tab 

| FamilyID | Argiopebruennichi<0> | Argiopebruennichi<2> | Stegodyphusdumicola<3> 
| --- | --- | --- | --- |
| OG0000000 | 1 | 1 | 0 | 0 |
| OG0000001	| -7 | -7 | 1 | 5 |	
| OG0000002	| 0	| 0	| -1 | -1 	

这个表格是 CAFE 的输出结果之一。像 `OG0000001` 这样的 ID，其实是在做基因家族聚类时生成的 orthogroup ID，并不是直接对应某一个已知的基因名称。

所以要知道 `OG0000001` 代表的是什么基因家族，需要回溯到在 CAFE 之前的 基因家族构建步骤。通常流程是这样的：

- **1. 基因家族构建**
  - 一般是用 **OrthoFinder** / **OrthoMCL** / SonicParanoid 等工具，把多个物种的蛋白质序列聚类成 orthogroups（即“基因家族”）。
  - 在这一步，每个 orthogroup 会被编号，比如 `OG0000001`。
-  **2. 家族成员查询**
   -  找到 OrthoFinder 等工具的输出目录下的 **Orthogroups.tsv** 或 Orthogroups.csv 文件。
   -  在里面搜索 `OG0000001`，会看到属于这个家族的所有基因（按物种分布）。
- **3. 功能注释**
  - 拿到家族里具体的基因 ID 后，可以：
    - 对这些基因做 BLASTp，看看它们对应的已知蛋白功能；
    - 或者用 InterProScan、eggNOG-mapper 等做功能注释；
    - 也可以检查你之前的注释文件（比如 gff3 + fasta 的注释结果），看家族里的基因分别注释成了什么蛋白。
- **4. 确定家族身份**
  - 结合各物种的注释结果（最好挑模式物种的注释，如 *Drosophila*, *Homo sapiens*），判断 `OG0000001` 大概是哪个基因家族（例如 GPCR、Opsin、Cytochrome P450、Cuticular protein 等）。

换句话说：`OG0000001` 本身只是一串编号。要知道它“是什么基因家族”，必须回溯到 **聚类结果（Orthogroups.tsv）**，找到它包含的基因，再用注释或 BLAST 来确定功能。

## 2. 从 gff 注释文件提取目标注释

通过对 Base_change.tab 的筛选，找出目标物种的目标基因家族（通常是显著扩张和收缩的家族）的 ID list

回溯 Orthogroups.tsv 找到其具体对应的基因家族

再根据注释文件（.faa）提取对应的注释信息

最终输出