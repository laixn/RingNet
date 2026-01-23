#' RingNet: community map builder for multi-omics network
#'
#' Build per-community node/edge maps from an interaction network and multi-omics
#' matrices, then export JSON.
#'
#' Required inputs: edges, nodes, meglist, out_json.
#' At least one of expression/methylation/snv/cnv must be provided.
#' stage is optional.
#'
#' @param edges Path to edges CSV (required). Must contain columns from/to or source/target.
#' @param nodes Path to nodes CSV (required). Must contain column cellgroup or name.
#' @param meglist Path to community membership CSV (required). Must contain community and one of cellgroup/name/gene.
#' @param expression Path to expression CSV (optional). Row names are samples, columns are genes.
#' @param methylation Path to methylation CSV (optional). Row names are samples, columns are genes.
#' @param snv Path to SNV CSV (optional). Row names are samples, columns are genes.
#' @param cnv Path to CNV CSV (optional). Row names are samples, columns are genes.
#' @param stage Path to stage CSV (optional). Two columns: sample + stage/index/group.
#' @param out_json Output JSON path (required).
#' @param TOP_N Top N nodes per community by mean expression.
#' @param KEEP_DIRECTED Whether to keep the graph directed.
#' @param n_cores Number of cores for parallel computation.
#'
#' @return (invisible) A list of per-community maps; JSON is written to out_json.
#' @export
#'
#' @importFrom jsonlite write_json
#' @importFrom igraph graph_from_data_frame induced_subgraph degree vcount V E
#' @importFrom graphlayouts layout_with_stress
#' @importFrom parallel detectCores makeCluster stopCluster parLapply clusterExport clusterEvalQ
#' @importFrom dplyr rename select
ringnet <- function(graph_edges,
                    graph_nodes,
                    node_group,
                    Data1 = NULL,
                    Data2 = NULL,
                    Data3 = NULL,
                    Data4 = NULL,
                    sample_group = NULL,
                    out_json,
                    TOP_N = 100,
                    KEEP_DIRECTED = TRUE,
                    n_cores = max(1L, parallel::detectCores() - 1L)) {
  
.pkgs <- c("parallel","igraph","graphlayouts","jsonlite","tictoc","dplyr","tidygraph")
miss <- .pkgs[!vapply(.pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, dependencies = TRUE)
suppressPackageStartupMessages({
  library(parallel)
  library(igraph)
  library(graphlayouts)
  library(jsonlite)
  library(tictoc)
  library(dplyr)
  library(tidygraph)
})

## ---------- 1) Read graph (CSV) and community (megList) ----------
edges_df <- read.csv(graph_edges, check.names = FALSE)      # edges input
nodes_df_tmp <- read.csv(graph_nodes, check.names = FALSE)

if ("cellgroup" %in% names(nodes_df_tmp)) {
  nodes_df <- nodes_df_tmp["cellgroup"]
} else if ("name" %in% names(nodes_df_tmp)) {
  nodes_df <- data.frame(name = nodes_df_tmp$name, stringsAsFactors = FALSE)
} else {
  stop("❌ nodes CSV ")
}

sgcell_id <- FALSE

if ("cellgroup" %in% names(nodes_df)) {
  if (any(!is.na(nodes_df$cellgroup) & nodes_df$cellgroup != "")) {
    sgcell_id <- TRUE
    message("🔍 'cellgroup' detected in nodes → sgcell_id = TRUE")
  }
}
# Support automatic mapping of columns source/target → from/to
if (!all(c("from","to") %in% names(edges_df))) {
  if (all(c("source","target") %in% names(edges_df))) {
    names(edges_df)[match(c("source","target"), names(edges_df))] <- c("from","to")
  } else stop("edges CSV list：from, to（or source, target）")
}

if (!("cellgroup" %in% names(nodes_df))) {
  if ("name" %in% names(nodes_df)) {
    #nodes_df$cellgroup <- nodes_df$name
    message("⚠️ 'cellgroup' column not found — using 'name' instead.")
  } else {
    stop("❌ nodes CSV must contain either 'cellgroup' or 'name' column.")
  }
}

if (!("weight" %in% names(edges_df))) edges_df$weight <- NA

# --- Read community information from megList ---
a_memb <- read.csv(node_group, check.names = FALSE) 

# Standardize column names
names(a_memb) <- tolower(trimws(names(a_memb)))

# Retain only the gene and community columns (drop 'cellgroup')
if ("cellgroup" %in% names(a_memb)) {
  a_memb <- a_memb[, c("cellgroup", "community")]
} else if ("name" %in% names(a_memb)) {
  
  a_memb <- a_memb[, c("name", "community")]
  
} else if("gene" %in% names(a_memb)){
  a_memb <- a_memb[, c("gene", "community")]
  
}


# Keep unique genes (retain only one occurrence per gene)
key_col <- if ("cellgroup" %in% names(a_memb)) {
  "cellgroup"
} else if ("name" %in% names(a_memb)) {
  "name"
} else if ("gene" %in% names(a_memb)) {
  "gene"
} else {
  stop("❌ a_memb must contain a column named 'cellgroup', 'name', or 'gene'")
}
# Retain unique genes (remove duplicate gene entries)
a_memb <- a_memb[!duplicated(a_memb[[key_col]]), ]

# 构建 membership 向量（key: cellgroup/name, value: community）
mem_vec <- setNames(as.integer(a_memb$community), a_memb[[key_col]])

rm(a_memb)
#########################
# Reconstruct the entire graph (consistent with the original structure;subsequent steps still use induced_subgraph)
#graph_comp <- graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = FALSE)
#if (is.null(E(graph_comp)$weight)) E(graph_comp)$weight <- 1

nodes_df<-unique(nodes_df)

if ("cellgroup" %in% names(nodes_df)) {
  nodes_df <- nodes_df %>% rename(name = cellgroup)
  
} else {
  
}

edges_df <- edges_df %>%
  select(from, to, everything())

graph_comp <- graph_from_data_frame(
  d = edges_df,
  vertices = nodes_df,
  directed = KEEP_DIRECTED   ## <-- Previously set to FALSE; now replaced with a toggle (default = TRUE)
)
if (is.null(E(graph_comp)$weight)) E(graph_comp)$weight <- 1
if (KEEP_DIRECTED) {
  V(graph_comp)$degree_all <- degree(graph_comp, mode = "all")
  V(graph_comp)$degree_in  <- degree(graph_comp, mode = "in")
  V(graph_comp)$degree_out <- degree(graph_comp, mode = "out")
} else {
  V(graph_comp)$degree_all <- degree(graph_comp, mode = "all")
}

melanet_spg <- structure(list(
  membership = mem_vec,
  algorithm  = "csv",
  modularity = NA_real_,
  vcount     = vcount(graph_comp)
), class = "communities")

#####################################################
message("🎯 Drawing colorful multi-edge network with curvature...")

# swicth to tidygraph object
tg <- as_tbl_graph(graph_comp)

# Validate multi-edges
dup_edges <- which_multiple(graph_comp)
if (any(dup_edges)) {
  message(sprintf("✅ Found %d multi-edges in graph.", sum(dup_edges)))
} else {
  message("ℹ No duplicated (multi) edges detected.")
}
##---------- 2) Read omics matrices ----------
read_or_empty <- function(path) {
  if (is.null(path) || path %in% c("", "NA", "NULL", "null", "-") || !file.exists(path)) {
    data.frame()
  } else {
    read.csv(path, row.names = 1, check.names = FALSE)
  }
}

read_stage_vector <- function(path) {
  if (is.null(path) || path %in% c("", "NA", "NULL", "null", "-") || !file.exists(path)) {
    return(NULL)
  }
  df <- tryCatch(read.csv(path, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || ncol(df) < 2) return(NULL)
  
  # Detect columns: prioritize name matching; otherwise take the first two columns
  sample_col <- which(grepl("sample|id|name", tolower(names(df))))[1]
  index_col  <- which(grepl("index|stage|class|group", tolower(names(df))))[1]
  if (is.na(sample_col) || is.na(index_col)) {
    sample_col <- 1; index_col <- 2
  }
  s <- as.character(df[[sample_col]])
  idx <- suppressWarnings(as.integer(df[[index_col]]))
  names(idx) <- s
  idx
}


############################################
e_raw  <- read_or_empty(Data1)         # expr
m_raw  <- read_or_empty(Data2)         # meth
snv_m  <- read_or_empty(Data3)         # snv
cnv_m  <- read_or_empty(Data4)         # cnv
# If node IDs are cellgroup-style, omics matrices may be provided as genes x samples.
# We transpose non-empty matrices to ensure: rows=samples, cols=genes.
if (sgcell_id) {
  if (nrow(e_raw)  > 0 && ncol(e_raw)  > 0) {
    message("🔄 transpose e_raw")
    e_raw <- as.data.frame(t(e_raw))
  }
  if (nrow(m_raw)  > 0 && ncol(m_raw)  > 0) {
    message("🔄 transpose m_raw")
    m_raw <- as.data.frame(t(m_raw))
  }
  if (nrow(snv_m) > 0 && ncol(snv_m) > 0) {
    message("🔄 transpose snv_m")
    snv_m <- as.data.frame(t(snv_m))
  }
  if (nrow(cnv_m) > 0 && ncol(cnv_m) > 0) {
    message("🔄 transpose cnv_m")
    cnv_m <- as.data.frame(t(cnv_m))
  }
}

stage_v <- read_stage_vector(sample_group)  # stage
# At least one omics file must be non-empty
#if (ncol(e_raw)==0 && ncol(m_raw)==0 && ncol(snv_m)==0 && ncol(cnv_m)==0) {
#  stop("At least one of expr/meth/snv/cnv must be provided.")
#}
if (ncol(e_raw)==0 && ncol(m_raw)==0 && ncol(snv_m)==0 && ncol(cnv_m)==0 && is.null(stage_v)) {
  stop("At least one of expr/meth/snv/cnv/stage must be provided.")
}

## ---------- 3) Synchronize samples ----------
#samples <- Reduce(intersect, list(rownames(e_raw), rownames(m_raw), rownames(snv_m), rownames(cnv_m)))
#if (!length(samples)) stop("please check the common sample name")
rn_list <- list()
if (nrow(e_raw)  > 0) rn_list[[length(rn_list)+1]] <- rownames(e_raw)
if (nrow(m_raw)  > 0) rn_list[[length(rn_list)+1]] <- rownames(m_raw)
if (nrow(snv_m) > 0) rn_list[[length(rn_list)+1]] <- rownames(snv_m)
if (nrow(cnv_m) > 0) rn_list[[length(rn_list)+1]] <- rownames(cnv_m)
if (!is.null(stage_v)) rn_list[[length(rn_list)+1]] <- names(stage_v)
if (!length(rn_list)) stop("no sample info provided at all")
samples <- Reduce(intersect, rn_list)
if (!length(samples)) stop("please check the common sample name (no overlap among provided matrices)")
stage_idx <- NULL
if (!is.null(stage_v)) {
  stage_idx <- unname(stage_v[samples])  #  NA could exist
  ord <- order(stage_idx, na.last = TRUE)  # Sort indices in ascending order; place NA values at the end
  samples <- samples[ord]
  
  # Reorder matrices according to the new sample order
  if (nrow(e_raw)  > 0) e_raw  <- e_raw[samples,,drop=FALSE]
  if (nrow(m_raw)  > 0) m_raw  <- m_raw[samples,,drop=FALSE]
  if (nrow(snv_m) > 0) snv_m  <- snv_m[samples,,drop=FALSE]
  if (nrow(cnv_m) > 0) cnv_m  <- cnv_m[samples,,drop=FALSE]
  
  # Reorder matrices again to align with sample order
  stage_idx <- unname(stage_v[samples])
} else {
  #If no stage file is provided, fill with all-NA placeholder ）
  stage_idx <- rep(NA_integer_, length(samples))
}

if (nrow(e_raw)  > 0) e_raw <- e_raw[samples,,drop=FALSE]
if (nrow(m_raw)  > 0) m_raw <- m_raw[samples,,drop=FALSE]
if (nrow(snv_m) > 0) snv_m <- snv_m[samples,,drop=FALSE]
if (nrow(cnv_m) > 0) cnv_m <- cnv_m[samples,,drop=FALSE]

## ---------- 4) Select TOP_N genes with highest expression (intersecting genes only) ----------
#common_genes <- Reduce(intersect, list(colnames(e_raw), colnames(m_raw), colnames(snv_m), colnames(cnv_m)))
#if (!length(common_genes)) stop("no common gene")
cn_list <- list()
if (ncol(e_raw)  > 0) cn_list[[length(cn_list)+1]] <- colnames(e_raw)
if (ncol(m_raw)  > 0) cn_list[[length(cn_list)+1]] <- colnames(m_raw)
if (ncol(snv_m) > 0) cn_list[[length(cn_list)+1]] <- colnames(snv_m)
if (ncol(cnv_m) > 0) cn_list[[length(cn_list)+1]] <- colnames(cnv_m)
common_genes <- Reduce(intersect, cn_list)
if (!length(common_genes)) stop("no common gene among provided matrices")

if (ncol(e_raw)  > 0) e_raw <- e_raw[, common_genes, drop = FALSE]
if (ncol(m_raw)  > 0) m_raw <- m_raw[, common_genes, drop = FALSE]
if (ncol(snv_m) > 0) snv_m <- snv_m[, common_genes, drop = FALSE]
if (ncol(cnv_m) > 0) cnv_m <- cnv_m[, common_genes, drop = FALSE]

# — Same as the original: trim the graph by common genes first, then reorder membership accordingly —
keep_vids <- which(V(graph_comp)$name %in% common_genes)
graph_comp <- induced_subgraph(graph_comp, vids = keep_vids)

keep_names <- V(graph_comp)$name
old_mem    <- melanet_spg$membership
new_mem    <- old_mem[keep_names]
names(new_mem) <- keep_names
melanet_spg$membership <- new_mem

## ---------- 5) Parallel computation of community_map_list (same structure as original) -----
if (!dir.exists(dirname(out_json))) dir.create(dirname(out_json), recursive = TRUE, showWarnings = FALSE)

tic("build community_map")
community_ids <- sort(unique(melanet_spg$membership))  # Keep compatibility with original layout (NA values not explicitly removed)
cl <- suppressWarnings(makeCluster(max(1L, n_cores)))
clusterEvalQ(cl, {library(igraph); library(graphlayouts);library(stats)})
clusterExport(cl, varlist = c("graph_comp","melanet_spg","e_raw","m_raw","snv_m","cnv_m","TOP_N","samples","stage_idx","edges_df"), envir = environment())

community_map_list <- parLapply(cl, community_ids, function(comm) {
  vids <- which(melanet_spg$membership == comm)
  subg <- induced_subgraph(graph_comp, vids = vids)
  
  # Use stress layout (same as original implementation)
  xy  <- tryCatch(layout_with_stress(subg) * 200, error = function(e) layout_nicely(subg))
  deg <- degree(subg,mode = "all")
  max_deg <- if (length(deg)) max(deg) else 0
  build_edges <- function(sg, keep_ids = NULL) {
    ec <- ecount(sg)
    if (ec == 0L) return(list())
    
    keep_nodes <- V(sg)$name
    
    # Find all global rows corresponding to edges in this community
    rows <- which(edges_df$from %in% keep_nodes & edges_df$to %in% keep_nodes)
    if (length(rows) == 0) return(list())
    
    # 1) w_raw: directly from the original weight
    w_raw_vec <- if ("w_raw" %in% names(edges_df)) {
      as.numeric(edges_df$w_raw[rows])
    } else {
      as.numeric(edges_df$weight[rows])
    }
    
    # 2) w_norm: normalize w_raw to [-1, 1]
    if (length(w_raw_vec) &&
        is.finite(min(w_raw_vec, na.rm = TRUE)) &&
        is.finite(max(w_raw_vec, na.rm = TRUE)) &&
        min(w_raw_vec, na.rm = TRUE) < max(w_raw_vec, na.rm = TRUE)) {
      
      w_min <- min(w_raw_vec, na.rm = TRUE)
      w_max <- max(w_raw_vec, na.rm = TRUE)
      w_norm_vec <- (w_raw_vec - w_min) / (w_max - w_min)
    } else {
      w_norm_vec <- rep(0, length(w_raw_vec))
    }
    
    # 3) w_z: z-score computed from w_raw
    sd_w <- suppressWarnings(sd(w_raw_vec, na.rm = TRUE))
    if (length(w_raw_vec) && is.finite(sd_w) && sd_w > 0) {
      w_z_vec <- as.numeric(scale(w_raw_vec))
    } else {
      w_z_vec <- rep(0, length(w_raw_vec))
    }
    
    out <- vector("list", length(rows))
    idx <- 0L
    
    for (k in seq_along(rows)) {
      r   <- rows[k]
      src <- as.character(edges_df$from[r])
      tgt <- as.character(edges_df$to[r])
      
      # Top-N node filtering
      if (!is.null(keep_ids) && !(src %in% keep_ids && tgt %in% keep_ids)) {
        next
      }
      
      idx <- idx + 1L
      out[[idx]] <- list(
        source = src,
        target = tgt,
        weight = w_raw_vec[k],     # original weight
        w_raw  = w_raw_vec[k],     # original weight
        w_norm = w_norm_vec[k],    # w_raw comes from current community 
        w_z    = w_z_vec[k],       # 
        interact_id = if ("interact" %in% names(edges_df)) {
          as.character(edges_df$interact[r])
        } else if ("interact_id" %in% names(edges_df)) {
          as.character(edges_df$interact_id[r])
        } else {
          NA_character_
        }
      )
    }
    
    out <- Filter(Negate(is.null), out)
    out
  }
  
  
  
  
  ## Nodes
  nodes <- lapply(seq_len(vcount(subg)), function(i) {
    node_name <- V(subg)$name[i]
    #gene <- V(subg)$name[i]
    gene <- sub("@.*", "", node_name)
    cellgroup <- ifelse(grepl("@", node_name),
                        sub(".*@", "", node_name),
                        NA)

    exp_norm <- if (ncol(e_raw) > 0 && gene %in% colnames(e_raw)) {
      tmp <- e_raw[, gene]; rng <- range(tmp, na.rm = TRUE)
      if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] < rng[2]) -1 + 2*(tmp-rng[1])/(rng[2]-rng[1]) else rep(0, length(tmp))
    } else rep(NA_real_, length(samples))
    
    mty_norm <- if (ncol(m_raw) > 0 && gene %in% colnames(m_raw)) {
      tmp <- m_raw[, gene]; rng <- range(tmp, na.rm = TRUE)
      if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] < rng[2]) -1 + 2*(tmp-rng[1])/(rng[2]-rng[1]) else rep(0, length(tmp))
    } else rep(NA_real_, length(samples))
    
    # ----- SNV：original value + norm + zscore -----
    snv_raw <- if (ncol(snv_m) > 0 && gene %in% colnames(snv_m)) {
      as.numeric(snv_m[, gene])
    } else {
      rep(NA_real_, length(samples))
    }
    
    snv_vals <- snv_raw
    
    snv_norm <- if (all(is.na(snv_raw))) {
      rep(NA_real_, length(samples))
    } else {
      rng <- range(snv_raw, na.rm = TRUE)
      if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] < rng[2]) {
        -1 + 2 * (snv_raw - rng[1]) / (rng[2] - rng[1])
      } else {
        rep(0, length(snv_raw))
      }
    }
    
    snv_z <- if (all(is.na(snv_raw))) {
      rep(NA_real_, length(samples))
    } else {
      as.numeric(scale(snv_raw))
    }
    
    
    
    
    cnv_norm <- if (ncol(cnv_m)>0 && gene %in% colnames(cnv_m)) as.numeric(cnv_m[, gene]) else rep(NA_real_, length(samples))
    
    list(
      id        = node_name,   
      gene      = gene,        
      cellgroup = cellgroup,   
      x         = xy[i,1],
      y         = xy[i,2],
      degree    = deg[node_name],
      exp_vals  = if (ncol(e_raw)>0 && gene %in% colnames(e_raw)) as.numeric(e_raw[, gene]) else rep(NA_real_, length(samples)),
      mty_vals  = if (ncol(m_raw)>0 && gene %in% colnames(m_raw)) as.numeric(m_raw[, gene]) else rep(NA_real_, length(samples)),
      cnv_vals  = if (ncol(cnv_m)>0 && gene %in% colnames(cnv_m)) as.numeric(cnv_m[, gene]) else rep(NA_real_, length(samples)),
      exp_norm  = exp_norm,
      exp_z  = if (ncol(e_raw) > 0 && gene %in% colnames(e_raw))  as.numeric(scale(e_raw[, gene])) else rep(NA_real_, length(samples)),  
      mty_norm  = mty_norm,
      mty_z  = if (ncol(m_raw) > 0 && gene %in% colnames(m_raw))  as.numeric(scale(m_raw[, gene])) else rep(NA_real_, length(samples)), 
      cnv_norm  = cnv_norm,
      cnv_z  = if (ncol(cnv_m) > 0 && gene %in% colnames(cnv_m))  as.numeric(scale(cnv_m[, gene])) else rep(NA_real_, length(samples)),  
      snv_vals  = snv_raw,
      snv_norm  = snv_norm,
      snv_z     = snv_z,
      stage_idx = as.integer(stage_idx)
    )
  })
  
   
  ## ---------- Keep only Top-100 nodes within each community (same logic as original) ----------
  if (length(nodes) > TOP_N) {
    
    scores <- vapply(nodes, function(n) mean(n$exp_vals, na.rm = TRUE), numeric(1))
    ## Replace non-finite values with -Inf to prevent all-NA genes from entering Top100
    scores[!is.finite(scores)] <- -Inf
    
    ## === CHANGED: Defensive handling of the actual retained node count ===
    keep_n   <- min(length(scores), TOP_N)
    keep_idx <- order(scores, decreasing = TRUE)[seq_len(keep_n)]
    keep_ids <- vapply(nodes[keep_idx], `[[`, "", "id")
    nodes <- nodes[keep_idx]
    edges <- build_edges(subg, keep_ids = keep_ids)
    max_deg <- max(vapply(nodes, `[[`, 0.0, "degree"))
  } else {
    # If the community size ≤ TOP_N, retain all edges (preserve original behavior)
    
    edges <- build_edges(subg)
  }
  
  list(comm = comm, max_deg = max_deg, nodes = nodes, edges = edges)
})

stopCluster(cl)
community_map_list <- Filter(Negate(is.null), community_map_list)
toc()
## ---------- 6) Write JSON output ----------
if (!length(community_map_list)) stop("null no community is output")

out_dir <- dirname(out_json)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

write_json(community_map_list, out_json, auto_unbox = TRUE, pretty = FALSE, na = "null")
cat(sprintf("✔ Done. Saved: %s\n", out_json))
invisible(community_map_list)  # Return the result silently
}
