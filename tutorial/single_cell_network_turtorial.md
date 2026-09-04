# Single-Cell Network Example Tutorial

This tutorial explains how to prepare and upload the CSV files required for the **Single-Cell Network** workflow. It is written for users who want to quickly understand each file, check the required columns, and run the example data successfully.

## Preparation notes

1. Keep `cellgroup` names consistent across all files.
2. Make sure each file contains the required columns.
3. Use the correct numeric encoding for `weight`, CNV, SNV, expression, and methylation values.

Once all files pass the checklist, upload them in the web interface, click **Upload Files**, and then run the analysis script.
---

## 1. What this workflow does

The single-cell template is designed for **cell-group interaction network analysis**. In this workflow, each network connection describes an interaction between cell groups, often linked by ligand–receptor or gene–gene interaction information.

Typical use cases include:

- visualizing interactions between cell groups;
- showing which genes are involved in each interaction;
- adding expression, methylation, CNV, and SNV values as additional network layers;
- grouping nodes by cell type or community.

---

## 2. Required input files

Upload the following CSV files in the web interface.

| Upload field | Example file name | Required? | Purpose |
| :---: | :---: | :---: | :---: |
| Edges file | `00_singlecell_edges.csv` | Yes | Defines interactions between cell groups or nodes. |
| Nodes file | `01_singlecell_nodes.csv` | Yes | Defines nodes and the cell group each node belongs to. |
| MEG / community file | `02_singlecell_megList.csv` | Yes | Assigns each gene-cell group pair to a community. |
| Data1 | `03_singlecell_expression_matrix_data1.csv` | Yes | Provides continuous values for each cell group. |
| Data2 | `04_singlecell_methylation_matrix_data2.csv` | Yes | Provides continuous values for each cell group. |
| Data3 | `05_singlecell_cnv_data3.csv` | Yes | Provides integer values for each cell group. |
| Data4 | `06_singlecell_snv_data4.csv` | Yes | Provides integer values for each cell group. |
| Sample group | optional | Optional | Adds grouping information for visualization, if supported by your workflow. |

> **Important:** The first column of matrix files must be `cellgroup`. This tells the system to use the single-cell analysis route.

---


## 3. File-by-file guide

### 3.1 `00_singlecell_edges.csv`

This file defines the network edges. Each row represents one interaction.

#### Required columns

| Column | Required? | Description |
| :---: | :---: | :---: |
| `from` | Yes | Source cell group or source node. |
| `to` | Yes | Target cell group or target node. |
| `weight` | Yes | Numeric interaction strength. |

#### Optional columns

| Column | Required? | Description |
| :---: | :---: | :---: |
| `from_gene` | Optional | Gene involved in the source side of the interaction. |
| `to_gene` | Optional | Gene involved in the target side of the interaction. |
| `interact` | Optional | Interaction label, such as ligand-receptor pair. If provided, it can be used to color edges by interaction type. |

#### Example

```csv
from_gene,to_gene,from,to,interact,weight
FGF7,FGFR1,APOE+ FIB,APOE+ FIB,FGF7_FGFR1,0.006276
FGF7,FGFR1,FBN1+ FIB,APOE+ FIB,FGF7_FGFR1,0.004855
FGF7,FGFR1,COL11A1+ FIB,APOE+ FIB,FGF7_FGFR1,0.003576
FGF7,FGFR1,APOE+ FIB,FBN1+ FIB,FGF7_FGFR1,0.009363
```

#### Notes

- `from`, `to`, and `weight` must not be empty.
- `weight` must be numeric.
- If `interact` is not empty, the interface can use it to color directed graph edges by interaction ID.

---

### 4.2 `01_singlecell_nodes.csv`

This file defines all nodes used in the network. Each row should describe a gene and the cell group it belongs to.

#### Required columns

| Column | Required? | Description |
| :---: | :---: | :---: |
| `name` | Yes | Gene name or node name, such as `FGF7`. |
| `cellgroup` | Yes | Cell group name, such as `APOE+ FIB`. |

#### Example

```csv
name,cellgroup
FGF7,APOE+ FIB
FGF7,FBN1+ FIB
FGF7,COL11A1+ FIB
VEGFA,LC
CCL19,Inflam. FIB
CXCL12,APOE+ FIB
```

#### Notes

- The same gene can appear in multiple cell groups.
- Cell group names should be written consistently across all files.
- Nodes appearing in the edge file should also be represented in the node file.

---

### 4.3 `02_singlecell_megList.csv`

This file assigns each gene-cell group pair to a community.

#### Required columns

| Column | Required? | Description |
| :---: | :---: | :---: |
| `gene` | Yes | Gene name. |
| `cellgroup` | Yes | Cell group name. |
| `community` | Yes | Integer community ID. |

#### Example

```csv
gene,cellgroup,community
FGF7,APOE+ FIB,1
FGF7,FBN1+ FIB,1
FGF7,COL11A1+ FIB,1
VEGFA,LC,1
CCL19,Inflam. FIB,1
CXCL12,APOE+ FIB,1
```

#### Notes

- `community` should be an integer.
- If you do not have predefined communities, use a placeholder integer such as `1`.
- If all genes belong to one network without community separation, assign all rows to the same community ID.

---

### 4.4 `03_singlecell_expression_matrix_data1.csv`

This file provides gene expression values for each cell group.

#### Format

- Rows represent cell groups.
- Columns represent genes.
- The first column must be named `cellgroup`.
- Values should be numeric expression values.

#### Example

```csv
cellgroup,ANGPTL1,ANGPTL2,ANGPTL4,ANXA1,APCDD1
APOE+ FIB,0,0,0,2.078231,0
FBN1+ FIB,0.472447,0.399898,0,2.441575,0
COL11A1+ FIB,0,0,0,0.62687,0.683742
Inflam. FIB,0,0,0,1.966858,1.618042
```

#### Notes

- The column name `cellgroup` is mandatory.
- Gene names in the matrix should match gene names used in the nodes and MEG files when applicable.
- Missing values should be avoided. Use `0` only if it truly means no signal or no measured value according to your preprocessing rule.

---

### 4.5 `04_singlecell_methylation_matrix_data2.csv`

This file has the same structure as the expression matrix, but values represent methylation levels.

#### Format

| Part | Meaning |
| :---: | :---: |
| Rows | Cell groups |
| Columns | Gene names |
| First column | Must be `cellgroup` |
| Values | Numeric methylation values |

#### Example

```csv
cellgroup,ANGPTL1,ANGPTL2,ANGPTL4,ANXA1,APCDD1
APOE+ FIB,0,0,0,0.627609,0
FBN1+ FIB,0.472447,0.399898,0,0.440520,0
COL11A1+ FIB,0,0,0,0.253500,0.683742
Inflam. FIB,0,0,0,0.522068,1.618042
```

---

### 4.6 `05_singlecell_cnv_data3.csv`

This file describes copy-number variation states for each cell group.

#### Format

- Rows represent cell groups.
- Columns represent genes.
- The first column must be `cellgroup`.
- Values should be integer CNV states.

#### CNV value meaning

| Value | Meaning |
| :---: | :---: |
| `-1` | Copy-number loss |
| `0` | Copy-number neutral |
| `1` | Copy-number gain |

#### Example

```csv
cellgroup,ACKR1,ACKR2,ACKR3,ACKR4,ACVR1,ACVR1B
APOE+ FIB,1,0,0,1,1,1
FBN1+ FIB,1,1,1,0,0,0
COL11A1+ FIB,1,0,1,0,1,1
Inflam. FIB,0,1,0,0,1,1
```

---

### 4.7 `06_singlecell_snv_data4.csv`

This file describes mutation counts or mutation status for each cell group.

#### Format

- Rows represent cell groups.
- Columns represent genes.
- The first column must be `cellgroup`.
- Values should be integer mutation counts.

#### SNV value meaning

| Value | Meaning |
| :---: | :---: |
| `0` | No mutation |
| `1` | One mutation |
| `2` | Two mutations |
| `>0` | Treated as mutation present in downstream analysis |

#### Example

```csv
cellgroup,ACKR1,ACKR2,ACKR3,ACKR4,ACVR1,ACVR1B
APOE+ FIB,1,0,0,0,0,0
FBN1+ FIB,1,1,0,1,0,1
COL11A1+ FIB,1,1,1,0,1,0
Inflam. FIB,1,0,0,0,0,0
```

---

## 5. Quick checklist before uploading

Use this checklist to avoid the most common upload errors.

| Check item | Why it matters |
| :---: | :---: |
| All files are saved as `.csv` | The upload form expects CSV files. |
| Matrix files contain a `cellgroup` column | This identifies the data as single-cell input. |
| Cell group names are consistent | `APOE+ FIB` and `APOE FIB` may be treated as different groups. |
| Required columns are present | Missing columns may cause the script to fail. |
| Numeric columns contain only numeric values | Non-numeric values in `weight`, CNV, SNV, expression, or methylation columns can cause errors. |
| CNV values use `-1`, `0`, or `1` | These represent loss, neutral, and gain. |
| SNV values use `0`, `1`, `2`, or other non-negative integers | Values greater than `0` are treated as mutation present. |
| File names match the upload field | This reduces confusion during manual upload. |

---
