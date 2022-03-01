
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MitoHEAR

<!-- badges: start -->

<!-- badges: end -->

MitoHEAR (*Mito*chondrial *HE*teroplasmy Analyze*R*) is an R package
that allows the estimation as well as downstream statistical analysis of
the mtDNA heteroplasmy calculated from single-cell datasets.

## Installation

You can install the released version of MitoHEAR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MitoHEAR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/ScialdoneLab/MitoHEAR/tree/master")
library(MitoHEAR)
```

## Getting started

The package has two main functions: get\_raw\_counts\_allele and
get\_heteroplasmy.

``` r
library(MitoHEAR)

load(system.file("extdata", "meta_data_antonio_final.Rda", package = "MitoHEAR"))
#cell_names <- meta_data_antonio_final$antonio_array.Comment.ENA_RUN.
#path_to_bam <- "full_path_to_bam_files"
#bam_input <- paste(path_to_bam,cell_names, ".unique.bam", sep = "")
#path_fasta <- "full_path_to_fasta_file"
#output_SNP_antonio_mt <- get_raw_counts_allele(bam_input, path_fasta, cell_names, cores_number = 1 )
```

``` r
load(system.file("extdata", "output_SNP_antonio_mt.Rda", package = "MitoHEAR"))
matrix_allele_counts <- output_SNP_antonio_mt[[1]]
# In this example we have 124 cells and 65196 columns (4 possible alleles for the 16299 bases in the mouse MT genome)
name_position_allele <- output_SNP_antonio_mt[[2]]
name_position <- output_SNP_antonio_mt[[3]]
```

``` r
row.names(meta_data_antonio_final) <- meta_data_antonio_final$antonio_array.Comment.ENA_RUN.
meta_data_antonio_final <- meta_data_antonio_final[row.names(matrix_allele_counts), ]
row.names(matrix_allele_counts) <- meta_data_antonio_final$antonio_array.Source.Name
row.names(meta_data_antonio_final) <- meta_data_antonio_final$antonio_array.Source.Name
```

We select only the cells for the 2-cells stage for down-stream
analysis.

``` r
stage_2_cells <- row.names(matrix_allele_counts)[grep("2cell_", row.names(matrix_allele_counts))]
stage_2_cells <- stage_2_cells[!grepl("32cell_", stage_2_cells)]
```

The next step is to obtain a matrix with allele frequencies and a matrix
with heteroplasmy values for each pair of cell-base. This is obtained
with the function *get\_heteroplasmy*. This function performs a two step
filtering procedure, the first on the cells and the second on the bases.
The aim is to keep only the cells that have more than *number\_reads*
counts in more than *number\_positions* bases and to keep only the bases
that are covered by more than *number\_reads* counts in all the
remaining cells (*filtering*=1) or in at least 50% of cells in each
cluster (*filtering*=2).

``` r

sc_data <- get_heteroplasmy(matrix_allele_counts[stage_2_cells, ], name_position_allele, name_position, 50, 2000, filtering = 1)
```

Among the output of *get\_heteroplasmy* there are the matrix with
heteroplasmy values and the matrix with allele frequencies, for all the
cells and bases that pass the two steps filtering procedure. The
heteroplasmy is computed as *1-max(f)*, where *f* are the frequencies of
the four alleles for every cell-base pair. For more info about the
output see *?get\_heteroplasmy*.

``` r
sum_matrix <- sc_data[[1]]
sum_matrix_qc <- sc_data[[2]]
heteroplasmy_matrix_sc <- sc_data[[3]]
allele_matrix_sc <- sc_data[[4]]
cluster_sc <- as.character(meta_data_antonio_final[row.names(heteroplasmy_matrix_sc), ]$antonio_array.Characteristics.developmental.stage.)
index_sc <- sc_data[[5]]
```

## Down-stream analysis

*MitoHEAR* offers several ways to extrapolate relevant information from
heteroplasmy measurement. For the identification of most different bases
according to heteroplasmy between two group of cells (i.e. two
clusters), an unpaired two-samples Wilcoxon test is performed with the
function *get\_wilcox\_test*. The heteroplasmy and the corresponding
allele frequencies for a specific base can be plotted with
*plot\_heteroplasmy* and *plot\_allele\_frequency*. If for each sample a
diffusion pseudo time information is available, then it is possible to
detect the bases whose heteroplasmy changes in a significant way along
pseudo-time with dpt\_test and to plot the trend with *plot\_dpt*. It is
also possible to perform a cluster analysis on the samples based on
distance matrix obtained from allele frequencies with
*clustering\_angular\_distance* and to visualize an heatmap of the
distance matrix with samples sorted according to the cluster result with
*plot\_heatmap*. This approach could be useful for lineage tracing
analysis. For more exhaustive information about the functions offered by
MitoHEAR see *Vignettes* section below and the help page of the single
functions. (*?function\_name*).

## Vignettes

The following vignette is provided within the package MitoHEAR and is
accessible within R:

``` r
#utils::vignette("MitoHEAR")
```
