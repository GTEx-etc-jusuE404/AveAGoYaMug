
##### This is a script for all the heatmaps wahahahah!

### differnt working directory for ease

setwd("~/Bioinformatics Work/Meth & RNA/Meta-analysis WGCNA")

### need a cormat for the modules vs Stir meth

cormat1 <- read.table("Correlation_matrix_TCGA_meth.txt")

### need a cormat for the modules vs TCGS traits
### need a cormat for the module vs METABRIC 

## want a heat map of the methylation v modules that shows the ordering etc
### not sure about that for the other traits?



pheatmap(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE,
         cluster_cols = TRUE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         clustering_callback = identity2, cutree_rows = NA, cutree_cols = NA,
         treeheight_row = ifelse(cluster_rows, 50, 0),
         treeheight_col = ifelse(cluster_cols, 50, 0), legend = TRUE,
         legend_breaks = NA, legend_labels = NA, annotation_row = NA,
         annotation_col = NA, annotation = NA, annotation_colors = NA,
         annotation_legend = TRUE, drop_levels = TRUE, show_rownames = T,
         show_colnames = T, main = NA, fontsize = 10, fontsize_row = fontsize,
         fontsize_col = fontsize, display_numbers = F, number_format = "%.2f",
         number_color = "grey30", fontsize_number = 0.8 * fontsize,
         gaps_row = NULL, gaps_col = NULL, labels_row = NULL,
         labels_col = NULL, filename = NA, width = NA, height = NA,
         silent = FALSE, ...)
