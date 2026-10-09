###### NPY SYSTEM GENE NETWORK ANALYSIS ######

# This script:
#
# 1. Loads the bird expression matrix and GTF annotation
# 2. Extracts the seven NPY-system genes of interest
# 3. Identifies Urban and Rural biological replicates
# 4. Calculates mean expression for each habitat
# 5. Calculates Spearman gene-gene correlations separately for Urban and Rural
# 6. Calculates network-level Psum for each habitat
# 7. Creates Cytoscape-compatible node and edge tables
# 8. Calculates network-level Psum for each habitat
# 9. Creates log2-transformed expression values for Cytoscape node coloring
# 10. Creates log2-transformed Urban vs Rural expression heatmap
#
# Cytoscape visualization:
#
# Node color  = mean expression
# Edge color  = correlation direction
# Edge width  = correlation strength
#
# Psum:
#
# Psum = sum(abs(all unique correlations)) / total possible connections
#
# For 7 genes:
#
# 7 * (7 - 1) / 2 = 21 possible connections


# 1. LOAD REQUIRED PACKAGES ########

library(dplyr)
library(stringr)
library(readr)
library(pheatmap)
library(RColorBrewer)


# 2. SET WORKING DIRECTORY ########

#setwd("~/Desktop/birds")


# 3. DEFINE NPY SYSTEM GENES ########

# Gene IDs are used to identify the correct rows in the expression
# matrix and annotation files.

npy_ids <- c(
  "EK904_000860",   # NPY
  "EK904_010633",   # NPY1R
  "EK904_002101",   # NPY2R
  "EK904_010655",   # NPY2R
  "EK904_003013",   # NPY4R
  "EK904_010631",   # NPY5R
  "EK904_002137"    # NPY6R
)


# Biological gene names corresponding to the Gene IDs

npy_names <- c(
  "NPY",
  "NPY1R",
  "NPY2R",
  "NPY2R",
  "NPY4R",
  "NPY5R",
  "NPY6R"
)


# Unique node names for Cytoscape.
#
# The two NPY2R genes have the same biological gene name, so unique
# identifiers are required to distinguish the two nodes.

node_names <- c(
  "NPY",
  "NPY1R",
  "NPY2R_1",
  "NPY2R_2",
  "NPY4R",
  "NPY5R",
  "NPY6R"
)


# 4. LOAD EXPRESSION MATRIX ########

expression_matrix <- read.delim(
  "bird_counts.noshort.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)


# Check expression matrix

head(expression_matrix)


# 5. LOAD GTF ANNOTATION ########

GTF_File <- read.delim(
  "GCA_011057915.1_Mmel_1.0_genomic.gtf.txt",
  header = FALSE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE
)


# Assign column names to the GTF

colnames(GTF_File) <- c(
  "seqname",
  "source",
  "feature",
  "start",
  "end",
  "score",
  "strand",
  "frame",
  "attribute"
)


# Check GTF

head(GTF_File)


# 6. EXTRACT GENE-LEVEL ANNOTATION ########

GTF_genes <- GTF_File %>%
  filter(feature == "gene") %>%
  mutate(
    
    gene_id = str_extract(
      attribute,
      '(?<=gene_id ")[^"]+'
    ),
    
    gene_biotype = str_extract(
      attribute,
      '(?<=gene_biotype ")[^"]+'
    ),
    
    functional_annotation = str_extract(
      attribute,
      '(?<=note ")[^"]+'
    ),
    
    InterPro = str_extract_all(
      attribute,
      "InterPro:IPR[0-9]+"
    ) %>%
      sapply(paste, collapse = "; "),
    
    PFAM = str_extract_all(
      attribute,
      "PFAM:PF[0-9]+"
    ) %>%
      sapply(paste, collapse = "; ")
  )


# 7. IDENTIFY NPY GENES IN THE GTF ########

npy_genes <- GTF_genes %>%
  filter(
    gene_id %in% npy_ids
  )


# Check whether all NPY genes were found

missing_gtf <- setdiff(
  npy_ids,
  npy_genes$gene_id
)


if (length(missing_gtf) > 0) {
  
  warning(
    paste(
      "The following NPY Gene IDs were not found in the GTF:",
      paste(
        missing_gtf,
        collapse = ", "
      )
    )
  )
}


# 8. EXTRACT NPY GENES FROM EXPRESSION MATRIX ########

# Extract only the seven NPY-system genes.

npy <- expression_matrix[
  rownames(expression_matrix) %in% npy_ids,
  ,
  drop = FALSE
]


# Add Gene_ID as an explicit column

npy$Gene_ID <- rownames(npy)


# Add biological gene names

npy$Gene <- npy_names[
  match(
    npy$Gene_ID,
    npy_ids
  )
]


# Reorder genes according to npy_ids

npy <- npy[
  match(
    npy_ids,
    npy$Gene_ID
  ),
]


# Check extracted NPY expression data

npy


# EXPORT NPY DATAFRAMES ######

# Save the extracted NPY expression dataframe.

write_csv(
  npy,
  "GHC_NPY_Expression_Data.csv"
)


# Save the NPY gene annotation dataframe.

write_csv(
  npy_genes,
  "GHC_NPY_Gene_Annotation.csv"
)


# 9. IDENTIFY URBAN AND RURAL SAMPLES ########

# Samples ending in "_U" are Urban.
# Samples ending in "_R" are Rural.

sample_cols <- grep(
  "_[UR]$",
  colnames(npy),
  value = TRUE
)


urban_cols <- grep(
  "_U$",
  sample_cols,
  value = TRUE
)


rural_cols <- grep(
  "_R$",
  sample_cols,
  value = TRUE
)


# Check sample assignments

print("Urban samples:")
print(urban_cols)

print("Rural samples:")
print(rural_cols)


# 10. CALCULATE MEAN EXPRESSION BY HABITAT ########

# Calculate mean expression across all Urban biological replicates.

mean_U <- rowMeans(
  npy[
    ,
    urban_cols,
    drop = FALSE
  ],
  na.rm = TRUE
)


# Calculate mean expression across all Rural biological replicates.

mean_R <- rowMeans(
  npy[
    ,
    rural_cols,
    drop = FALSE
  ],
  na.rm = TRUE
)


# 11. CREATE NODE EXPRESSION DATAFRAME ########

# This dataframe contains one row per NPY-system gene.
#
# This will serve as the Cytoscape node attribute table.

cyto_nodes <- data.frame(
  
  id = node_names,
  
  Gene = npy_names,
  
  Gene_ID = npy$Gene_ID,
  
  Mean_U = mean_U,
  
  Mean_R = mean_R,
  
  stringsAsFactors = FALSE
)


# Check node expression data

cyto_nodes


# 12. PREPARE URBAN EXPRESSION MATRIX ########

# Extract Urban biological replicates.

expr_U <- as.matrix(
  npy[
    ,
    urban_cols,
    drop = FALSE
  ]
)


# Assign Gene IDs to rows

rownames(expr_U) <- npy$Gene_ID


# Transpose so:
#
# rows = biological replicates
# columns = genes

expr_U <- t(expr_U)


# Ensure genes are in the predefined order

expr_U <- expr_U[
  ,
  npy_ids,
  drop = FALSE
]


# 13. PREPARE RURAL EXPRESSION MATRIX ########

# Extract Rural biological replicates.

expr_R <- as.matrix(
  npy[
    ,
    rural_cols,
    drop = FALSE
  ]
)


# Assign Gene IDs to rows

rownames(expr_R) <- npy$Gene_ID


# Transpose

expr_R <- t(expr_R)


# Ensure genes are in the predefined order

expr_R <- expr_R[
  ,
  npy_ids,
  drop = FALSE
]


# 14. CALCULATE SPEARMAN CORRELATIONS ########

# Calculate gene-gene Spearman correlations separately for each habitat.

cor_U_matrix <- cor(
  expr_U,
  method = "spearman",
  use = "pairwise.complete.obs"
)


cor_R_matrix <- cor(
  expr_R,
  method = "spearman",
  use = "pairwise.complete.obs"
)


# Check correlation matrices

cor_U_matrix

cor_R_matrix


# 15. CREATE URBAN CYTOSCAPE EDGE TABLE ########

# Because correlations are undirected, only one copy of each
# gene-gene relationship is retained using the upper triangle
# of the correlation matrix.

upper_U <- upper.tri(
  cor_U_matrix,
  diag = FALSE
)


edge_indices_U <- which(
  upper_U,
  arr.ind = TRUE
)


cyto_edges_U <- data.frame(
  
  Source = node_names[
    edge_indices_U[, "row"]
  ],
  
  Target = node_names[
    edge_indices_U[, "col"]
  ],
  
  Correlation_U = cor_U_matrix[
    upper_U
  ],
  
  stringsAsFactors = FALSE
)


# Calculate absolute Urban correlation strength

cyto_edges_U$Abs_Correlation_U <- abs(
  cyto_edges_U$Correlation_U
)


# Assign Urban correlation direction

cyto_edges_U$Direction_U <- ifelse(
  cyto_edges_U$Correlation_U >= 0,
  "Positive",
  "Negative"
)


# Keep Urban edge columns in logical order

cyto_edges_U <- cyto_edges_U[
  c(
    "Source",
    "Target",
    "Correlation_U",
    "Abs_Correlation_U",
    "Direction_U"
  )
]


# 16. CREATE RURAL CYTOSCAPE EDGE TABLE ########

# Create the Rural correlations using the same gene-pair structure.

upper_R <- upper.tri(
  cor_R_matrix,
  diag = FALSE
)


edge_indices_R <- which(
  upper_R,
  arr.ind = TRUE
)


cyto_edges_R <- data.frame(
  
  Source = node_names[
    edge_indices_R[, "row"]
  ],
  
  Target = node_names[
    edge_indices_R[, "col"]
  ],
  
  Correlation_R = cor_R_matrix[
    upper_R
  ],
  
  stringsAsFactors = FALSE
)


# Calculate absolute Rural correlation strength

cyto_edges_R$Abs_Correlation_R <- abs(
  cyto_edges_R$Correlation_R
)


# Assign Rural correlation direction

cyto_edges_R$Direction_R <- ifelse(
  cyto_edges_R$Correlation_R >= 0,
  "Positive",
  "Negative"
)


# Keep Rural edge columns in logical order

cyto_edges_R <- cyto_edges_R[
  c(
    "Source",
    "Target",
    "Correlation_R",
    "Abs_Correlation_R",
    "Direction_R"
  )
]


# 17. COMBINE URBAN AND RURAL EDGE TABLES ########

# Urban and Rural contain the same 21 possible gene-gene connections.
#
# Combine them into a single Cytoscape edge table so that the
# correlation values for Urban and Rural remain independent
# edge attributes.

cyto_edges <- cyto_edges_U %>%
  left_join(
    cyto_edges_R,
    by = c(
      "Source",
      "Target"
    )
  )


# Reorder combined edge columns

cyto_edges <- cyto_edges[
  c(
    "Source",
    "Target",
    "Correlation_U",
    "Abs_Correlation_U",
    "Direction_U",
    "Correlation_R",
    "Abs_Correlation_R",
    "Direction_R"
  )
]


# Check combined edge table

cyto_edges


# 18. CALCULATE NUMBER OF POSSIBLE CONNECTIONS ########

# For an undirected network, the number of possible connections is:
#
# n(n - 1) / 2

n_genes <- length(npy_ids)


n_connections <- (
  n_genes * (n_genes - 1)
) / 2


print(
  paste(
    "Number of genes:",
    n_genes
  )
)


print(
  paste(
    "Number of possible connections:",
    n_connections
  )
)


# 19. CALCULATE PSUM ########

# Psum represents the overall connectivity of the gene network.
#
# The absolute value is used so that both positive and negative
# correlations contribute to network connectivity.

Psum_U <- sum(
  abs(
    cyto_edges_U$Correlation_U
  ),
  na.rm = TRUE
) / n_connections


Psum_R <- sum(
  abs(
    cyto_edges_R$Correlation_R
  ),
  na.rm = TRUE
) / n_connections


# 20. CREATE PSUM DATAFRAME ########

Psum_df <- data.frame(
  
  Habitat = c(
    "Urban",
    "Rural"
  ),
  
  Psum = c(
    Psum_U,
    Psum_R
  ),
  
  stringsAsFactors = FALSE
)


# Check Psum values

Psum_df


# 21. ADD PSUM TO CYTOSCAPE NODE TABLE ########

# Psum is a network-level statistic rather than a gene-level
# statistic. The value is therefore repeated for every node.

cyto_nodes$Psum_U <- Psum_U

cyto_nodes$Psum_R <- Psum_R


# 22. LOG2-TRANSFORM EXPRESSION FOR CYTOSCAPE ########

# NPY has substantially higher expression than the receptor genes.
#
# Log2(expression + 1) compresses the large range of expression
# values while retaining information about expression magnitude.
#
# These values are intended for continuous node-color mapping
# in Cytoscape.

cyto_nodes$Log2_Mean_U <- log2(
  cyto_nodes$Mean_U + 1
)


cyto_nodes$Log2_Mean_R <- log2(
  cyto_nodes$Mean_R + 1
)


# Reorder node columns

cyto_nodes <- cyto_nodes[
  c(
    "id",
    "Gene",
    "Gene_ID",
    "Mean_U",
    "Mean_R",
    "Log2_Mean_U",
    "Log2_Mean_R",
    "Psum_U",
    "Psum_R"
  )
]


# 23. CHECK CYTOSCAPE NODE TABLE ########

cyto_nodes


# 24. CHECK COMBINED CYTOSCAPE EDGE TABLE ########

cyto_edges


# 25. VERIFY ALL NODES ARE REPRESENTED ########

# Because the edge table contains unique undirected connections,
# a node does not need to appear specifically in the Source column.
# It only needs to appear somewhere as Source or Target.

all_nodes <- unique(
  c(
    cyto_edges$Source,
    cyto_edges$Target
  )
)


print("Nodes represented in network:")
print(all_nodes)


# 26. CHECK FOR MISSING NODES ########

missing_nodes <- setdiff(
  node_names,
  all_nodes
)


print("Missing nodes:")
print(missing_nodes)


# Expected result:
#
# character(0)


# 27. EXPORT CYTOSCAPE NODE TABLE ########

write_csv(
  cyto_nodes,
  "GHC_NPY_Cytoscape_Nodes.csv"
)


# 28. EXPORT COMBINED CYTOSCAPE EDGE TABLE ########

write_csv(
  cyto_edges,
  "GHC_NPY_Cytoscape_Edges_Urban_Rural.csv"
)


# 29. EXPORT PSUM TABLE ########

write_csv(
  Psum_df,
  "GHC_NPY_Cytoscape_Psum.csv"
)


# 30. PREPARE EXPRESSION DATA FOR HEATMAP ########

# Use the same node dataframe generated above so that the heatmap
# and Cytoscape use the exact same mean-expression values.

heatmap_df <- cyto_nodes[
  ,
  c(
    "Gene_ID",
    "Gene",
    "Mean_U",
    "Mean_R"
  )
]


# Reorder genes according to npy_ids

heatmap_df <- heatmap_df[
  match(
    npy_ids,
    heatmap_df$Gene_ID
  ),
]


# 31. CREATE HEATMAP EXPRESSION MATRIX ########

heatmap_matrix <- as.matrix(
  heatmap_df[
    ,
    c(
      "Mean_U",
      "Mean_R"
    )
  ]
)


colnames(heatmap_matrix) <- c(
  "Urban",
  "Rural"
)


# 32. LOG2-TRANSFORM MEAN EXPRESSION FOR HEATMAP ########

# Log2 transformation compresses the large range of expression values
# while retaining differences in expression magnitude between genes.

heatmap_matrix <- log2(
  heatmap_matrix + 1
)


# 33. ASSIGN UNIQUE NODE NAMES TO HEATMAP ROWS ########

# Use the same unique node names as Cytoscape.
# This ensures NPY2R genes are displayed as NPY2R_1 and NPY2R_2.

rownames(heatmap_matrix) <- node_names

colnames(heatmap_matrix) <- c(
  "Urban",
  "Rural"
)



# 34. SAVE LOG2-TRANSFORMED VALUES AS DATAFRAME ########

heatmap_log2_df <- data.frame(
  
  Gene = heatmap_df$Gene,
  
  Gene_ID = heatmap_df$Gene_ID,
  
  Urban = heatmap_matrix[
    ,
    "Urban"
  ],
  
  Rural = heatmap_matrix[
    ,
    "Rural"
  ],
  
  stringsAsFactors = FALSE
)


# View heatmap dataframe

heatmap_log2_df


# Save heatmap log2 dataframe

write_csv(
  heatmap_log2_df,
  "GHC_NPY_Heatmap_Log2.csv"
)


# 35. CREATE HEATMAP COLOR PALETTE ########

# Low expression = light
# High expression = dark purple

my_palette <- colorRampPalette(
  c(
    "#F7F4F9",
    "#E0D4E8",
    "#C2A5CF",
    "#9970AB",
    "#762A83",
    "#40004B"
  )
)(100)


# 36. GENERATE NPY EXPRESSION HEATMAP ########

pheatmap(
  
  heatmap_matrix,
  
  color = my_palette,
  
  cluster_rows = FALSE,
  
  cluster_cols = FALSE,
  
  show_rownames = TRUE,
  
  show_colnames = TRUE,
  
  angle_col = "315",
  
  border_color = "grey60",
  
  legend = TRUE
)


# 37. FINAL DATAFRAMES ########

# Extracted NPY expression data:
npy


# NPY gene annotation:
npy_genes


# Cytoscape node attributes:
cyto_nodes


# Combined Cytoscape edge attributes:
cyto_edges


# Network-level Psum:
Psum_df


# Heatmap log2 dataframe:
heatmap_log2_df


### ADD SOURCE AND TARGET NODE EXPRESSION TO EDGE TABLE ########

# Add log2-transformed mean expression for the Source node.

cyto_edges$Source_Log2_Mean_U <- cyto_nodes$Log2_Mean_U[
  match(
    cyto_edges$Source,
    cyto_nodes$id
  )
]

cyto_edges$Source_Log2_Mean_R <- cyto_nodes$Log2_Mean_R[
  match(
    cyto_edges$Source,
    cyto_nodes$id
  )
]


# Add log2-transformed mean expression for the Target node.

cyto_edges$Target_Log2_Mean_U <- cyto_nodes$Log2_Mean_U[
  match(
    cyto_edges$Target,
    cyto_nodes$id
  )
]

cyto_edges$Target_Log2_Mean_R <- cyto_nodes$Log2_Mean_R[
  match(
    cyto_edges$Target,
    cyto_nodes$id
  )
]


### REORDER EDGE TABLE COLUMNS ########

cyto_edges <- cyto_edges[
  c(
    "Source",
    "Target",
    
    "Correlation_U",
    "Abs_Correlation_U",
    "Direction_U",
    
    "Correlation_R",
    "Abs_Correlation_R",
    "Direction_R",
    
    "Source_Log2_Mean_U",
    "Source_Log2_Mean_R",
    
    "Target_Log2_Mean_U",
    "Target_Log2_Mean_R"
  )
]


### CHECK EDGE TABLE ########

cyto_edges


### EXPORT UPDATED EDGE TABLE ########

write_csv(
  cyto_edges,
  "GHC_NPY_Cytoscape_Edges_Urban_Rural.csv"
)

