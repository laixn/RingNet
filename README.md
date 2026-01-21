# RingNet: an interactive platform for multi-modal data visualization in networks

RingNet is a tool for visualizing multi-modal data in networks. It supports creating networks with one or more multiple data sets. The R package export a `JSON` file that can be further customized and visualized using the JavaScript frontend at https://fip-128-214-252-149.kaj.poutavm.fi/cmt_figures/.

If you have any suggestions for improving RingNet, please feel free to contact us or report it in the issues.

## Overview
The core function is ringnet(), which:

- loads an interaction graph (edges + nodes),

- loads community membership (meglist),

- optionally loads multi-omics matrices and sample stage/index,

- intersects samples and genes across inputs,builds community subgraphs and computes layouts/degree statistics,

- selects Top-N nodes per community (optional),

- writes a compact JSON for downstream visualization.

## Repository Structure
```
ringnet/
├── R/
│   └── ringnet.R          # ringnet() implementation
├── man/
│   └── ringnet.Rd         # documentation
├── DESCRIPTION
├── NAMESPACE
└── README.md              # this file
```
## Installation
##### Local installation (GitHub)
```
install.packages("remotes")
remotes::install_github("laixn/RingNet",
                        dependencies = TRUE,
                        upgrade = "never")
library(ringnet)
```
##### Installation on a server (user folder)
```
libdir <- "PATH/TO/YOUR/USER/RLIB"   # change this
dir.create(libdir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(libdir, .libPaths()))
install.packages("remotes", lib = libdir, repos = "https://cloud.r-project.org")
remotes::install_github("laixn/RingNet",
                        lib = libdir,
                        dependencies = TRUE,
                        upgrade = "never")
library(ringnet, lib.loc = libdir)
```
## Inputs and Data Preparation
### Required inputs

RingNet requires three CSV files:

 - Graph Edges: interaction network edges

 - Graph Nodes: node list / node metadata

 - Node Group: nodes group

 - and one output out_json: output JSON path

At least one of the following omics matrices must be provided:
- Data1
- Data2
- Data3
- Data4

******1) Graph Edges CSV format (required)******

Must contain either:`from`, `to`, `weight`. 
Optional columns: interact.

******2) Graph Nodes CSV format (required)******

Must contain either: `cellgroup` which is used for `single cell data` or `name` which is used for `general data`.

RingNet will standardize node identifiers to a column named name.

******3) Node Group CSV format (required)******

Must contain: `community`,and `cellgroup` or `name` or `gene`. RingNet will build a membership vector from this table.

******4) omics matrices (optional but at least one required)******

Each omics matrix is a CSV loaded with: The samples are rows by default, genes are columns in `gene - gene network`. However, if the data is `patient - patient network`, the rows are the name of genes,column names are the sample ID of patients. Similarly, in `single cell network` type,the rows are the name of cell group, the name of columns are proteins or genes.

******5) Important alignment rules:******

RingNet will intersect samples across all provided omics matrices (and stage file if provided).RingNet will intersect genes across all provided omics matrices.The network is then restricted to genes present in the final intersection.

******6) Sample Group: CSV format (optional)******

A two-column CSV containing:a sample identifier column (name matching `sample`/`id`/`name` is auto-detected; otherwise first column),a stage/index/group column (name matching `index`/`stage`/`class`/`group` is auto-detected; otherwise second column)

If provided, samples are ordered by stage/index before export. Each sample in the first column belongs to exactly one class in the second column. Multiple samples may belong to the same class, but no sample belongs to more than one class.

## Main Function: ringnet()
#### Description

Build per-community node/edge maps from an interaction network and multi-omics matrices, then export JSON.

****Usage****
```
ringnet(graph_edges,
        graph_nodes,
        node_group,
        Data1  = NULL,
        Data2  = NULL,
        Data3  = NULL,
        Data4  = NULL,
        sample_group:  = NULL,
        out_json,
        TOP_N  = 100,
        KEEP_DIRECTED = TRUE,
        n_cores = max(1L, parallel::detectCores() - 1L))
```
****Arguments****

- `graph_edges` (required): Path to edges CSV. Must contain from/to or source/target.

- `graph_nodes` (required): Path to nodes CSV. Must contain cellgroup or name.

- `node_group` (optional): Path to community membership CSV. Must contain community and one of cellgroup/name/gene.

- `Data1` (optional): Matrix CSV. A file containing continuous values for the outermost ring.

- `Data2` (optional): Matrix CSV. A file containing continuous values for the second ring.

- `Data3` (optional): Matrix CSV. A file containing integer values for the third ring.

- `Data4` (optional): Matrix CSV. A file containing integer values for the innermost ring.

- `sample_group:` (optional): Stage/index CSV. Two columns: sample + stage/index/group.sample Group: A file containing group information for individual bars in the ring (e.g., genes in a patient similarity network or patients in a gene-gene interaction network).

- `out_json` (required): Output JSON path.

- `TOP_N` (default: 100): Keep Top-N nodes per community ranked by mean expression (if expression is provided). If fewer nodes exist, keeps all.

- `KEEP_DIRECTED` (default: TRUE): Whether to treat graph as directed.

- `n_cores` (default: detectCores()-1): Number of CPU cores for parallel community processing in Linux platform.

****Output****

A JSON file written to `out_json`

The function returns (invisibly) a list of per-community maps

`Example:` Loading the ringnet lib.
```
library(ringnet)
cwd <- getwd()
edges_path <- file.path(cwd, "data", "cellchat_edges_cell_level.csv")
nodes_path <- file.path(cwd, "data", "cellchat_nodes.csv")
nodegroup_path  <- file.path(cwd, "data", "cellchat_nodes_with_cellgroup.csv")
expression <- file.path(cwd, "data", "cell_group_expression_matrix.csv")
out_json   <- file.path(cwd, "data", "ringnet_community.json")
```
`Example:` Basic run (expression only)
```
ringnet(graph_edges   = edges_path,
        graph_nodes   = nodes_path,
        node_group = nodegroup_path,
        Data1  = expression,
        out_json = out_json,
        TOP_N  = 200,
        KEEP_DIRECTED = TRUE)
```

`Example:` Multi-omics run (expression + methylation + SNV + CNV + stage)
```
ringnet( `graph_edges`       = "data/edges.csv",
        `graph_nodes`       = "data/nodes.csv",
        `node_group`     = "data/membership.csv",
        `Data1`  = "data/expression.csv",
        `Data2` = "data/methylation.csv",
        `Data3`        = "data/cnv.csv",
        `Data4`         = "data/snv.csv",
        `sample_group`       = "data/stage.csv",
        `out_json`    = "output/ringnet_community.json",
        `TOP_N`       = 150,
        `KEEP_DIRECTED` = FALSE,
        `n_cores`     = 8)
```
### Notes on Behavior and Common Pitfalls ####

You must provide `at least one omics matrix` (`Data1`/`Data2`/`Data3`/`Data4`), otherwise the function stops.

RingNet uses the intersection of:
-  samples across all inputs (omics + stage if provided)

-  genes across all provided omics matrices

-  The network is reduced to genes retained after intersection. If your graph nodes do not match gene symbols used in matrices, you may end up with an empty or tiny graph.

-  If your nodes.csv uses `cellgroup`, RingNet renames it to name and follow the process of single cell route.

- If nodes contain cellgroup and it is non-empty, RingNet may transpose the expression matrix (the code detects a cellgroup-style mode and transposes expression accordingly). Ensure your expression file orientation matches your data conventions.

****Output JSON Schema****

The output JSON is a list of communities, each containing:

- `comm`: community id

- `max_deg`: maximum node degree in the community (after filtering)

- `nodes`: list of node objects (id, gene, optional cellgroup, layout x/y, degree, omics vectors and normalized versions, stage index)

- `edges`: list of edges (source, target, weight, w_raw, w_norm, w_z, optional interact_id),w_raw is the `original` weight value, w_norm is `min-max` weight value, and w_z is the `z-score` weight value.

-  This JSON is designed to be consumed by downstream JavaScript visualization code.

## Usage & Citation 
If you find our work useful, please consider citing it:

Liang Zhang, Xin Lai. RingNet: an interactive platform for multi-modal data visualization in networks. 
Submitted (2026).

```bash
@article{Zhang_RingNet_2026,
  title={RingNet: an interactive platform for multi-modal data visualization in networks},
  author={Liang Zhang, Fei Liu, Xin Lai},
  journal={Submitted},
  doi={},
  year={2026}
}
```

© [Lai Lab](https://sites.google.com/view/lai-lab) - This code is made available under the GPLv3 License and is available for non-commercial academic purposes.
```
