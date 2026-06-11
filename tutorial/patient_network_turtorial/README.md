# Patient–Patient Network Tutorial

This tutorial explains how to prepare and upload the CSV files for a **Patient–Patient interaction network**.

A Patient–Patient network is used to visualize similarity relationships between patient samples, such as TCGA samples. In this network:

- each **node** represents one patient sample;
- each **edge** represents a similarity relationship between two patients;
- multi-omics files can be displayed as circular annotation rings around the network.

---

## 1. Quick Start

Prepare the following CSV files and upload them on the **Upload CSV files** page.

| Upload field | Example file name | Required | Purpose |
| :---: | :---: | :---: | :---: |
| Graph Edges | `00_patients_edges.csv` | Yes | Defines patient-to-patient similarity edges |
| Graph Nodes | `01_patients_nodes.csv` | Yes | Defines all patient nodes in the network |
| Node Group | `02_patients_megList.csv` | Yes | Assigns each patient to a community/group |
| Data 1 | `03_patients_expression_data1.csv` | Yes | Gene expression values, shown as the outer ring |
| Data 2 | `04_patients_methylation_data2.csv` | Yes | DNA methylation values, shown as the second ring |
| Data 3 | `05_patients_cnv_data3.csv` | Yes | CNV states, shown as the third ring |
| Data 4 | `06_patients_snv_data4.csv` | Yes | SNV/mutation states, shown as the inner ring |
| Sample Group | optional sample-group CSV | Optional | Adds group labels for ring bars or stratified visualization |

> Note: some older documents may call the edge file `00_patients_edges.csv`. In the current example package, `00_patients_edges.csv` is usually the expected file name.

## Important Notes Before Editing Files

Please keep the following rules in mind before modifying the CSV files. These rules are important because the upload script reads the files according to fixed column names and matrix structures.

| Item | Do not modify | Can modify | Notes |
| :---: | :---: | :---: | :---: |
| File names | Do not change the required CSV file names unless the upload interface is also updated | You may replace the data inside each CSV | Use the exact example names to avoid upload confusion |
| Required column names | Do not rename required columns such as `from`, `to`, `weight`, `name`, `community`, or `gene` | You may add optional columns only if the workflow supports them | Column names are used by the script to identify the file structure |
| Row names / row IDs | Do not change the meaning of the first column in matrix files | You may replace the row values with your own gene symbols | In expression, methylation, CNV, and SNV files, rows represent genes |
| Column names in matrix files | Do not use inconsistent patient IDs | You may replace example patient IDs with your own patient IDs | Matrix columns must match patient names in the node file |
| Patient IDs | Do not add extra spaces, remove suffixes, or change capitalization | You may use your own patient IDs if they are consistent across all files | `TCGA-AQ-A0Y5-01A` and `TCGA-AQ-A0Y5` may be treated as different samples |
| Numeric values | Do not use text labels in numeric fields | You may replace example values with your own numeric values | `weight`, expression, methylation, CNV, and SNV values must be numeric |
| CNV / SNV coding | Do not replace numeric codes with words such as `gain`, `loss`, or `mutation` | You may update values using the accepted numeric codes | CNV should use `-1`, `0`, `1`; SNV should use non-negative integers |

---

## 2. Recommended Upload Workflow

Suggested order:

1. Prepare the patient list first.
2. Make sure every patient in the edge file also appears in the node file.
3. Prepare the node group file.
4. Prepare expression, methylation, CNV, and SNV matrices.
5. Upload all required files.
6. Click **Upload Files**.
7. Click **Run R Script**.

---

## 3. Input File Details

### 3.1 Graph Edges File

**Example file:** `00_patients_edges.csv`

This file defines similarity relationships between patients. Each row represents one edge from one patient to another patient.

The edge `weight` is a patient–patient similarity score defined by our own method. In this example, the weights are calculated from a similarity matrix generated using a Gaussian kernel function. Therefore, the `weight` column should be understood as the Gaussian-kernel-based similarity score between two patient samples.

Because the directed-graph algorithm reads each row as a directional relationship, the visualization may contain two symmetric edges: one from `source` to `target`, and another from `target` to `source`. This happens when the table contains both `source -> target` and `target -> source` records. These paired records come from the symmetric similarity matrix produced during similarity-score calculation. For clearer observation of patient–patient similarity relationships, the undirected-graph entrance is usually recommended.

| Column | Required | Description |
| :---: | :---: | :---: |
| `from` | Yes | Source patient TCGA sample ID |
| `to` | Yes | Target patient TCGA sample ID |
| `weight` | Yes | Similarity score between the two patients |

Higher `weight` values indicate stronger similarity between patients.

Example:

| from | to | weight |
| :---: | :---: | :---: |
| TCGA-E2-A1IN-01A | TCGA-AQ-A0Y5-01A | 0.001202 |
| TCGA-BH-A0H9-01A | TCGA-AQ-A0Y5-01A | 0.001773 |
| TCGA-A2-A0YT-01A | TCGA-AQ-A0Y5-01A | 0.002035 |
| TCGA-D8-A1JA-01A | TCGA-AQ-A0Y5-01A | 0.002168 |

Checklist:

- `from` and `to` must use exactly the same patient ID format as the node file.
- Do not leave `weight` empty.
- Use numeric values for `weight`.
- Do not rename `from`, `to`, or `weight`.
- The `weight` value can be replaced with your own patient–patient similarity score, but it must remain numeric.
- If the edge table contains both `A -> B` and `B -> A`, a directed graph may display two symmetric edges. Use the undirected-graph option when the goal is to observe similarity without direction.

---

### 3.2 Graph Nodes File

**Example file:** `01_patients_nodes.csv`

This file defines all patient samples that should appear as nodes in the network.

| Column | Required | Description |
| :---: | :---: | :---: |
| `name` | Yes | Patient TCGA sample ID |

Example:

- `name`
- `TCGA-3C-AAAU-01A`
- `TCGA-3C-AALI-01A`
- `TCGA-3C-AALJ-01A`
- `TCGA-3C-AALK-01A`
- `TCGA-4H-AAAK-01A`

Checklist:

- Every patient ID in `from` and `to` of the edge file should also appear in this file.
- Patient IDs must be consistent across all uploaded files.
- Avoid extra spaces before or after patient IDs.

---

### 3.3 Node Group File

**Example file:** `02_patients_megList.csv`

This file assigns each patient to a community or group.

| Column | Required | Description |
| :---: | :---: | :---: |
| `name` | Yes | Patient TCGA sample ID |
| `community` | Yes | Integer group/community ID |

Example:

| name | community |
| :---: | :---: |
| TCGA-3C-AAAU-01A | 1 |
| TCGA-3C-AALI-01A | 1 |
| TCGA-3C-AALJ-01A | 1 |
| TCGA-3C-AALK-01A | 1 |
| TCGA-4H-AAAK-01A | 1 |

If you do not have predefined patient communities, you can use a placeholder integer. For example, assign all patients to `1`.

---

## 4. Multi-Omics Data Files

The Patient–Patient network uses gene-by-patient matrices for multi-omics annotation. In these files:

- rows are **gene symbols**;
- columns are **patient TCGA sample IDs**;
- values are numeric omics measurements or mutation/CNV states.

### 4.1 Gene Expression Matrix

**Example file:** `03_patients_expression_data1.csv`

This file contains normalized gene expression values.

| Format item | Description |
| :---: | :---: |
| Rows | Gene symbols |
| Columns | Patient TCGA sample IDs |
| Values | Log2-transformed normalized expression values |

Example:

| gene | TCGA-AQ-A0Y5-01A | TCGA-C8-A274-01A | TCGA-B6-A1KC-01B |
| :---: | :---: | :---: | :---: |
| TSPAN6 | 0.642019 | 1.435314 | 0.905893 |
| DPM1 | 1.590688 | 1.520711 | 1.46089 |
| SCYL3 | 0.731065 | 1.010347 | 0.716954 |
| FGR | 0.327584 | 0.295589 | 0.23573 |

> If your file does not explicitly name the first column `gene`, make sure the first column still contains gene symbols and the remaining columns are patient IDs.

---

### 4.2 DNA Methylation Matrix

**Example file:** `04_patients_methylation_data2.csv`

This file has the same structure as the expression matrix.

| Format item | Description |
| :---: | :---: |
| Rows | Gene symbols |
| Columns | Patient TCGA sample IDs |
| Values | Normalized methylation levels |

Example:

| gene | TCGA-AQ-A0Y5-01A | TCGA-C8-A274-01A | TCGA-B6-A1KC-01B |
| :---: | :---: | :---: | :---: |
| TSPAN6 | 0.627609 | 0.44052 | 0.2535 |
| DPM1 | 0.048339 | 0.048328 | 0.068088 |
| SCYL3 | 0.186311 | 0.174254 | 0.216154 |
| FGR | 0.64619 | 0.734972 | 0.69399 |

---

### 4.3 Copy Number Variation Matrix

**Example file:** `05_patients_cnv_data3.csv`

This file describes copy number variation states for genes across patients.

| Format item | Description |
| :---: | :---: |
| Rows | Gene symbols |
| Columns | Patient TCGA sample IDs |
| Values | CNV state |

CNV value meanings:

| Value | Meaning |
| :---: | :---: |
| `0` | Copy-number neutral |
| `1` | Copy-number gain |
| `-1` | Copy-number loss |

Example:

| gene | TCGA-AQ-A0Y5-01A | TCGA-C8-A274-01A | TCGA-B6-A1KC-01B |
| :---: | :---: | :---: | :---: |
| UPF1 | -1 | 1 | 1 |
| KRIT1 | 0 | 0 | 0 |
| CREBBP | 1 | 1 | -1 |
| MPO | -1 | 0 | 1 |

---

### 4.4 Single Nucleotide Variation Matrix

**Example file:** `06_patients_snv_data4.csv`

This file records mutation information for genes across patients. Its structure is the same as the CNV matrix.

| Format item | Description |
| :---: | :---: |
| Rows | Gene symbols |
| Columns | Patient TCGA sample IDs |
| Values | Mutation count or mutation state |

SNV value meanings:

| Value | Meaning |
| :---: | :---: |
| `0` | No mutation |
| `1` | One mutation |
| `2` | Two mutations |
| `>0` | Treated as mutation present in downstream analysis |

Example:

| gene | TCGA-AQ-A0Y5-01A | TCGA-C8-A274-01A | TCGA-B6-A1KC-01B |
| :---: | :---: | :---: | :---: |
| PIK3CA | 0 | 1 | 0 |
| TP53 | 1 | 0 | 2 |
| BRCA1 | 0 | 0 | 0 |
| GATA3 | 1 | 1 | 0 |

---

## 5. Data Consistency Rules

Before uploading, check these rules carefully.

### Patient ID consistency

Patient IDs must match across files.

For example, this ID:

- `TCGA-AQ-A0Y5-01A`

should be written exactly the same way in:

- edge file;
- node file;
- node group file;
- expression matrix columns;
- methylation matrix columns;
- CNV matrix columns;
- SNV matrix columns.

Avoid problems such as:

- `TCGA-AQ-A0Y5-01A`
- `TCGA-AQ-A0Y5`
- `tcga-aq-a0y5-01a`

These may be treated as different IDs.

---

### Required columns

| File | Required columns |
| :---: | :---: |
| `00_patients_edges.csv` | `from`, `to`, `weight` |
| `01_patients_nodes.csv` | `name` |
| `02_patients_megList.csv` | `name`, `community` |
| Expression matrix | first column = gene symbols; other columns = patient IDs |
| Methylation matrix | first column = gene symbols; other columns = patient IDs |
| CNV matrix | first column = gene symbols; other columns = patient IDs |
| SNV matrix | first column = gene symbols; other columns = patient IDs |

---
## 6. Final Checklist Before Upload

- [ ] Edge file contains `from`, `to`, and `weight`, and these column names have not been renamed.
- [ ] Edge weights are numeric similarity scores; in this example they are generated from a Gaussian-kernel-based similarity matrix.
- [ ] For similarity-matrix data, symmetric pairs such as `A -> B` and `B -> A` are expected in directed-graph input; use the undirected-graph entrance for cleaner similarity visualization.
- [ ] Node file contains all patients used in the edge file.
- [ ] Node group file contains `name` and `community`.
- [ ] Expression matrix columns are patient IDs.
- [ ] Methylation matrix columns are patient IDs.
- [ ] CNV matrix uses valid CNV values: `-1`, `0`, `1`.
- [ ] SNV matrix uses valid mutation values: `0`, `1`, `2`, or other non-negative counts.
- [ ] Patient IDs are identical across all files.
- [ ] Files are saved as `.csv`.
- [ ] ZIP package includes all required CSV files.
