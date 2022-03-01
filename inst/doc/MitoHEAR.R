## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, 
  fig.height=5
)



## ----setup--------------------------------------------------------------------
library(MitoHEAR)

## -----------------------------------------------------------------------------

load(system.file("extdata", "meta_data_antonio_final.Rda", package = "MitoHEAR"))

## -----------------------------------------------------------------------------
load(system.file("extdata", "output_SNP_antonio_mt.Rda", package = "MitoHEAR"))


## -----------------------------------------------------------------------------
matrix_allele_counts <- output_SNP_antonio_mt[[1]]
name_position_allele <- output_SNP_antonio_mt[[2]]
name_position <- output_SNP_antonio_mt[[3]]


## -----------------------------------------------------------------------------
row.names(meta_data_antonio_final) <- meta_data_antonio_final$antonio_array.Comment.ENA_RUN.
meta_data_antonio_final <- meta_data_antonio_final[row.names(matrix_allele_counts), ]
row.names(matrix_allele_counts) <- meta_data_antonio_final$antonio_array.Source.Name
row.names(meta_data_antonio_final) <- meta_data_antonio_final$antonio_array.Source.Name

## -----------------------------------------------------------------------------
sc_data_all <- get_heteroplasmy(matrix_allele_counts, name_position_allele, name_position, 50, 2000, filtering = 1)

## -----------------------------------------------------------------------------
sum_matrix <- sc_data_all[[1]]
sum_matrix_qc <- sc_data_all[[2]]
heteroplasmy_matrix_sc <- sc_data_all[[3]]
allele_matrix_sc <- sc_data_all[[4]]
cluster_sc <- as.character(meta_data_antonio_final[row.names(heteroplasmy_matrix_sc), ]$antonio_array.Characteristics.developmental.stage.)
index_sc <- sc_data_all[[5]]

## -----------------------------------------------------------------------------
stage_2_cells <- row.names(matrix_allele_counts)[grep("2cell_", row.names(matrix_allele_counts))]
stage_2_cells <- stage_2_cells[!grepl("32cell_", stage_2_cells)]

## -----------------------------------------------------------------------------

sc_data <- get_heteroplasmy(matrix_allele_counts[stage_2_cells, ], name_position_allele, name_position, 50, 2000, filtering = 1)


## -----------------------------------------------------------------------------
sum_matrix <- sc_data[[1]]
sum_matrix_qc <- sc_data[[2]]
heteroplasmy_matrix_sc <- sc_data[[3]]
allele_matrix_sc <- sc_data[[4]]
cluster_sc <- as.character(meta_data_antonio_final[row.names(heteroplasmy_matrix_sc), ]$antonio_array.Characteristics.developmental.stage.)
condition_sc <- rep(0, length(cluster_sc))
condition_sc[grep("2cell_1_", row.names(heteroplasmy_matrix_sc))] <- "1"
condition_sc[grep("2cell_2_", row.names(heteroplasmy_matrix_sc))] <- "2"
condition_sc[grep("2cell_3_", row.names(heteroplasmy_matrix_sc))] <- "3"
condition_sc[grep("2cell_4_", row.names(heteroplasmy_matrix_sc))] <- "4"
condition_sc[grep("2cell_5_", row.names(heteroplasmy_matrix_sc))] <- "5"
condition_sc[grep("2cell_6_", row.names(heteroplasmy_matrix_sc))] <- "6"
condition_sc[grep("2cell_7_", row.names(heteroplasmy_matrix_sc))] <- "7"
condition_sc[grep("2cell_8_", row.names(heteroplasmy_matrix_sc))] <- "8"
index_sc <- sc_data[[5]]

## -----------------------------------------------------------------------------

name_position_allele_qc <- name_position_allele[name_position%in%colnames(sum_matrix_qc)]
name_position_qc <- name_position[name_position%in%colnames(sum_matrix_qc)]


## -----------------------------------------------------------------------------

result_clustering_sc <- clustering_angular_distance(heteroplasmy_matrix_sc, allele_matrix_sc, condition_sc, length(colnames(heteroplasmy_matrix_sc)), deepSplit_param = 0, minClusterSize_param = 2, 0.2, min_value = 0.001, index = index_sc, relevant_bases = NULL)


old_new_classification <- result_clustering_sc[[1]]
dist_matrix_sc <- result_clustering_sc[[2]]
top_dist <- result_clustering_sc[[3]]
common_idx <- result_clustering_sc[[4]]

old_classification <- as.vector(old_new_classification[, 1])
new_classification <- as.vector(old_new_classification[, 2])



## -----------------------------------------------------------------------------
plot_heatmap(new_classification, old_classification, (dist_matrix_sc), cluster_columns = F, cluster_rows = F, "Euclidean distance")


## -----------------------------------------------------------------------------
q <- list()
for ( i in 1:length(top_dist[1:4])){
p <- plot_heteroplasmy(top_dist[i], heteroplasmy_matrix_sc, condition_sc, index_sc)
q <- list(q, p)
}
q




## -----------------------------------------------------------------------------
stage_8_cells <- row.names(matrix_allele_counts)[grep("8cell_", row.names(matrix_allele_counts))]


## -----------------------------------------------------------------------------

sc_data <- get_heteroplasmy(matrix_allele_counts[stage_8_cells, ], name_position_allele, name_position, 50, 2000, filtering = 1)



## -----------------------------------------------------------------------------
sum_matrix <- sc_data[[1]]
sum_matrix_qc <- sc_data[[2]]
heteroplasmy_matrix_sc <- sc_data[[3]]
allele_matrix_sc <- sc_data[[4]]
cluster_sc <- as.character(meta_data_antonio_final[row.names(heteroplasmy_matrix_sc), ]$antonio_array.Characteristics.developmental.stage.)
condition_sc <- rep(0, length(cluster_sc))
condition_sc[grep("8cell_1_", row.names(heteroplasmy_matrix_sc))] <- "1"
condition_sc[grep("8cell_2_", row.names(heteroplasmy_matrix_sc))] <- "2"
condition_sc[grep("8cell_3_", row.names(heteroplasmy_matrix_sc))] <- "3"
condition_sc[grep("8cell_4_", row.names(heteroplasmy_matrix_sc))] <- "4"

index_sc <- sc_data[[5]]

## -----------------------------------------------------------------------------

name_position_allele_qc <- name_position_allele[name_position%in%colnames(sum_matrix_qc)]
name_position_qc <- name_position[name_position%in%colnames(sum_matrix_qc)]


## -----------------------------------------------------------------------------

result_clustering_sc <- clustering_angular_distance(heteroplasmy_matrix_sc, allele_matrix_sc, condition_sc, length(colnames(heteroplasmy_matrix_sc)), deepSplit_param = 0, minClusterSize_param = 8, 0.2, min_value = 0.001, index = index_sc, relevant_bases = NULL)


old_new_classification <- result_clustering_sc[[1]]
dist_matrix_sc <- result_clustering_sc[[2]]
top_dist <- result_clustering_sc[[3]]
common_idx <- result_clustering_sc[[4]]

old_classification <- as.vector(old_new_classification[, 1])
new_classification <- as.vector(old_new_classification[, 2])


## -----------------------------------------------------------------------------
plot_heatmap(new_classification, old_classification, (dist_matrix_sc), cluster_columns = F, cluster_rows = F, "Euclidean distance")


## -----------------------------------------------------------------------------
utils::sessionInfo()

