# Gene Network Example Tutorial

This tutorial explains how to prepare and upload CSV files for building a **Gene-Gene Interaction Network**. It is designed for first-time users who want to quickly understand which files are required, what each file means, and how to check whether the uploaded data are valid.


## Preparation notes

- Start with the provided example files first.
- Replace the example values with your own data only after you understand the format.
- Keep column names unchanged.
- Keep gene symbols and sample IDs consistent across all files.
- If your data do not include one optional layer, leave that optional file empty only if the system allows it, or skip that upload field if it is not required.
---

## 1. What this example does

The Gene Network module visualizes relationships between genes as a network.

- **Nodes** represent genes.
- **Edges** represent interactions between genes.
- **Rings around the network** show omics or clinical values, such as gene expression, methylation, CNV, SNV, or sample groups.


---

## 2. Required and optional files

Upload the files listed below. The first three files are required for the basic network structure. The data matrices are used to draw outer rings around the network.

| Upload field | Example file name | Required? | What it contains | Used for |
| :---: | :---: | :---: | :---: | :---: |
| Graph Edges | `00_graph_edges.csv` | Yes | Gene-gene interactions | Network edges |
| Graph Nodes | `01_graph_nodes.csv` | Yes | Gene node list | Network nodes |
| Node Group | `02_megList.csv` | Yes | Gene community or module assignment | Node grouping / module coloring |
| Data 1 | `03_gene_expression_data1.csv` | Recommended | Continuous gene expression values | Outermost ring |
| Data 2 | `04_methylation_data2.csv` | Recommended | Continuous DNA methylation values | Second ring |
| Data 3 | `05_cnv_data3.csv` | Recommended | CNV status values | Third ring |
| Data 4 | `06_snv_data4.csv` | Recommended | SNV / mutation values | Innermost ring |
| Sample Group | `06_stage.csv` | Optional | Patient stage or sample group information | Grouped bars / stratified display |

> Tip: Keep the file names simple and avoid spaces. CSV files should be comma-separated and saved in UTF-8 format.

---

## 3. Recommended upload order

Upload files in this order to reduce mistakes:

1. `00_graph_edges.csv`
2. `01_graph_nodes.csv`
3. `02_megList.csv`
4. `03_gene_expression_data1.csv`
5. `04_methylation_data2.csv`
6. `05_cnv_data3.csv`
7. `06_snv_data4.csv`
8. `07_stage.csv` if available


---

## 4. File format details

### 4.1 Graph Edges file

**File:** `00_graph_edges.csv`

This file defines interactions between genes. Each row is one edge.

| Column | Required? | Description | Example |
| :---: | :---: | :---: | :---: |
| `from` | Yes | Source gene symbol | `DNAJC14` |
| `to` | Yes | Target gene symbol | `AGTR1` |
| `interact` | Optional | Interaction category or label | `0` or `1` |
| `weight` | Recommended | Interaction strength. Positive values indicate positive association; negative values indicate inverse association. | `0.09132` |

Example:

```csv
from,to,interact,weight
DNAJC14,AGTR1,0,0.09132
DNAJC14,DRD1,0,0.02288
DBN1,PTEN,0,-0.1533
PTEN,MME,1,0.40162
```

Checklist:

- Every gene in `from` and `to` should also appear in `01_graph_nodes.csv`.
- Gene names should be consistent across all files.
- Do not leave empty values in `from` or `to`.

---

### 4.2 Graph Nodes file

**File:** `01_graph_nodes.csv`

This file lists all genes used as nodes in the network.

| Column | Required? | Description | Example |
| :---: | :---: | :---: | :---: |
| `name` | Yes | Gene symbol used as the node ID | `DNAJC14` |

Example:

```csv
name
DNAJC14
DBN1
NCOA3
UBE2I
DLC1
PRDX1
```

Checklist:

- The `name` column must match gene names used in the edges file.
- Avoid duplicate gene names.
- Use the same capitalization everywhere.

---

### 4.3 Node Group file

**File:** `02_megList.csv`

This file assigns each gene to a community, module, or functional group.

| Column | Required? | Description | Example |
| :---: | :---: | :---: | :---: |
| `gene` | Yes | Gene symbol | `DNAJC14` |
| `community` | Yes | Integer community ID | `1` |

Example:

```csv
gene,community
DNAJC14,1
IGFBP2,1
GGA3,1
NMT1,1
SYPL1,1
```

If you do not have community information, assign all genes to one community:

```csv
gene,community
DNAJC14,1
DBN1,1
NCOA3,1
UBE2I,1
```

Checklist:

- Every gene should appear in the node group file.
- `community` should be an integer, such as `1`, `2`, or `3`.
- If your network has no groups, use the same community value for all genes.

---

### 4.4 Data 1: Gene Expression matrix

**File:** `03_gene_expression_data1.csv`

This file contains gene expression values across patient samples.

| Part | Meaning |
| :---: | :---: |
| Rows | Patient/sample IDs, such as TCGA sample IDs |
| Columns | Gene symbols |
| Values | Log2-transformed normalized gene expression values |

Example:

```csv
,TSPAN6,DPM1,SCYL3,FGR,GCLC,NFYA
TCGA-AQ-A0Y5-01A,0.642019,1.590688,0.731065,0.327584,0.729861,1.189431
TCGA-C8-A274-01A,1.435314,1.520711,1.010347,0.295589,0.5792,1.274742
TCGA-B6-A1KC-01B,0.905893,1.46089,0.716954,0.23573,0.634024,1.395498
```

Checklist:

- The first column should contain sample IDs.
- Gene columns should use the same gene symbols as the network files.
- Values should be numeric.

---

### 4.5 Data 2: DNA Methylation matrix

**File:** `04_methylation_data2.csv`

This file has the same matrix format as gene expression.

| Part | Meaning |
| :---: | :---: |
| Rows | Patient/sample IDs |
| Columns | Gene symbols |
| Values | Normalized methylation levels |

Example:

```csv
,TSPAN6,DPM1,SCYL3,FGR,GCLC,NFYA
TCGA-AQ-A0Y5-01A,0.627609,0.048339,0.186311,0.64619,0.323994,0.105899
TCGA-C8-A274-01A,0.44052,0.048328,0.174254,0.734972,0.337164,0.113769
TCGA-B6-A1KC-01B,0.2535,0.068088,0.216154,0.69399,0.330322,0.14041
```

Checklist:

- Use numeric methylation values.
- Sample IDs should match the expression, CNV, and SNV files whenever possible.
- Gene names should be consistent with other files.

---

### 4.6 Data 3: CNV(copy number variation) matrix 

**File:** `05_cnv_data3.csv``

This file describes copy number variation states for genes across samples.

| Value | Meaning |
| :---: | :---: |
| `-1` | Copy-number loss |
| `0` | Copy-number neutral |
| `1` | Copy-number gain |

Example:

```csv
,ACAP3,ACTRT2,AGRN,ANKRD65,ATAD3A,ATAD3B
TCGA-3C-AAAU-01A,0,0,0,0,0,0
TCGA-3C-AALI-01A,-1,-1,-1,-1,-1,-1
TCGA-3C-AALJ-01A,-1,-1,-1,-1,-1,-1
TCGA-3C-AALK-01A,0,0,0,0,0,0
```

Checklist:

- Values should usually be `-1`, `0`, or `1`.
- Avoid text labels such as `gain` or `loss`; use numeric codes instead.
- Keep sample IDs consistent with other data matrices.

---

### 4.7 Data 4: SNV(single nucleotide variant) matrix

**File:** `06_snv_data4.csv`

This file records mutation counts per gene for each sample.

| Value | Meaning |
| :---: | :---: |
| `0` | No mutation |
| `1` | One mutation |
| `2` | Two mutations |
| `>0` | Treated as mutation present in downstream analysis |

Example:

```csv
,PIK3CA,DSP,PITPNM1,USP15,ACTG1,NRCAM
TCGA-EW-A2FW-01A,0,0,0,0,0,1
TCGA-AR-A0U0-01A,0,0,0,0,0,0
TCGA-AR-A24Q-01A,1,0,0,0,1,0
TCGA-D8-A27R-01A,0,0,0,0,0,0
```

Checklist:

- Values should be non-negative integers.
- Any value greater than `0` is interpreted as mutation present.
- Keep gene and sample names consistent.

---

### 4.8 Sample Group file

**File:** `07_stage.csv`

This optional file assigns each sample to a clinical stage or group. It is useful for grouped visualization and stratified analysis.

Recommended format:

```csv
sample,group
TCGA-AQ-A0Y5-01A,Stage II
TCGA-C8-A274-01A,Stage III
TCGA-B6-A1KC-01B,Stage I
```

Checklist:

- The sample ID column should match the sample IDs in the data matrices.
- Group names should be short and consistent.
- This file is optional. You can skip it if no sample group information is available.

---
## 5. How the rings are displayed

The uploaded data files are shown as rings around the network.

| Ring position | Upload field | Data type | Example file |
| :---: | :---: | :---: | :---: |
| Outermost ring | Data 1 | Continuous | `03_gene_expression_data1.csv` |
| Second ring | Data 2 | Continuous | `04_methylation_data2.csv` |
| Third ring | Data 3 | Integer / categorical numeric | `05_cnv_data3.csv` |
| Innermost ring | Data 4 | Integer / mutation count | `06_snv_data4.csv` |

---

## 6. Before uploading: quick validation checklist

Use this checklist before clicking **Upload Files**.

| Check item | Why it matters |
| :---: | :---: |
| All files are saved as `.csv` | The system expects CSV input files |
| Gene names are consistent across files | Prevents missing nodes or empty ring values |
| Sample IDs are consistent across matrix files | Ensures patient-level values align correctly |
| Edges contain valid `from` and `to` genes | Required for network construction |
| Nodes include all genes from the edge file | Prevents broken or missing nodes |
| Numeric columns contain only numbers | Prevents analysis or visualization errors |
| No empty required columns | Avoids upload failure |
| Optional sample group file matches sample IDs | Enables correct grouping in visualization |



