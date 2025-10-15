
# GenomeScaleEmbeddings

We explore the SNP-level LLM embeddings from the paper “Incorporating
LLM Embeddings for Variation Across the Human Genome”. We use `duckdb`
to download the embeddings from huggingFace, save them in
[`houba`](github.com/HervePerdry/houba) memory-mapped matrices, and use
[`bigPCACpp`](https://github.com/fbertran/bigPCAcpp/) to perform PCA and
some explorations.

## Load the packages

``` r
library(GenomeScaleEmbeddings)
library(knitr)
library(httr2)
library(jsonlite)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
```

## 1. Peek into remote parquet files

``` r
# Use OpenRemoteParquetView to inspect the first few rows
OpenRemoteParquetView()
#> # Source:   table<embeddings> [?? x 6]
#> # Database: DuckDB 1.4.0 [root@Linux 6.8.0-78-generic:R 4.5.1//tmp/Rtmps6ucDx/file8a5ff195dec19.duckdb]
#>    chrom pos       ref_UKB alt_UKB rsid       embedding    
#>    <chr> <chr>     <chr>   <chr>   <chr>      <list>       
#>  1 5     148899362 T       G       rs4705280  <dbl [3,072]>
#>  2 5     148899764 C       T       rs6872985  <dbl [3,072]>
#>  3 5     148900624 C       A       rs1181141  <dbl [3,072]>
#>  4 5     148901339 C       T       rs1181139  <dbl [3,072]>
#>  5 5     148902146 G       T       rs4705282  <dbl [3,072]>
#>  6 5     148902250 G       A       rs10036926 <dbl [3,072]>
#>  7 5     148903759 C       T       rs17108911 <dbl [3,072]>
#>  8 5     148904179 C       A       rs1181137  <dbl [3,072]>
#>  9 5     148904381 T       G       rs4705283  <dbl [3,072]>
#> 10 5     148906034 C       T       rs1181135  <dbl [3,072]>
#> # ℹ more rows
```

## 2. Copy remote parquet files into local DuckDB

``` r
system.time(

CopyParquetToDuckDB(db_path = "local_embeddings.duckdb", overwrite = FALSE)
)
#> Database file 'local_embeddings.duckdb' already exists. Skipping copy.
#>    user  system elapsed 
#>   0.024   0.007   0.028
file.info("local_embeddings.duckdb")$size
#> [1] 12106084352
```

## 3. Write embeddings to houba mmatrix

``` r
system.time(

houba <- writeEmbeddingsHoubaFromDuckDB(dbPath = "local_embeddings.duckdb", 
overwrite = FALSE)
)
#> Warning in mk.descriptor.file(object@file, object@dim[1], object@dim[2], :
#> local_embeddings.houba.desc already exists.
#> Using existing houba file and info data.
#>    user  system elapsed 
#>   0.581   0.033   0.570
houba
#> Houba mmatrix file: local_embeddings.houba 
#> Embeddings (houba::mmatrix):
#> A read-only mmatrix with 616386 rows and 3072 cols
#> data type:  double 
#> File: local_embeddings.houba 
#> --- excerpt
#>              [,1]       [,2]        [,3]         [,4]         [,5]
#> [1,] -0.008660397 0.02568399 -0.01839902 -0.001554275 -0.014016987
#> [2,] -0.003067377 0.01635006 -0.02215753  0.003889058 -0.013320981
#> [3,] -0.002486568 0.04113422 -0.01771097  0.024944620 -0.022576459
#> [4,] -0.002763281 0.02626457 -0.01791935  0.010160015 -0.011453237
#> [5,] -0.013722803 0.02826897 -0.01959616 -0.005454814 -0.009612823
#> 
#> Info data.frame (first 10 rows):
#>    chrom       pos ref_UKB alt_UKB       rsid
#> 1      5 148899362       T       G  rs4705280
#> 2      5 148899764       C       T  rs6872985
#> 3      5 148900624       C       A  rs1181141
#> 4      5 148901339       C       T  rs1181139
#> 5      5 148902146       G       T  rs4705282
#> 6      5 148902250       G       A rs10036926
#> 7      5 148903759       C       T rs17108911
#> 8      5 148904179       C       A  rs1181137
#> 9      5 148904381       T       G  rs4705283
#> 10     5 148906034       C       T  rs1181135
file.info("local_embeddings.houba")$size
#> [1] 15148302336
```

## 4. Run PCA on houba mmatrix

``` r
system.time(
  pca_res <- houbaPCA("local_embeddings.houba")
)
#> Attached big.matrix from descriptor file: local_embeddings.houba.desc
#> Dimensions: 616386 x 3072
#> Running PCA with center=TRUE, scale=TRUE, ncomp=15
#>    user  system elapsed 
#> 358.806 112.452  54.725
```

## 5. Get PCA scores

``` r
pc_scores <- getPcaScores(pca_res)
object.size(pc_scores)
#> 73966536 bytes
```

## 6. Plot PCA dimensions (PC1 vs PC2, colored by chromosome)

``` r
plotPcaDims(pc_scores, houba$info, annotation_col = "chrom", dim1 = 1, dim2 = 2)
```

<img src="docs/README-unnamed-chunk-8-1.png" width="100%" />

## 7. Annotate variants with coffee consumption GWAS and plot PC1 vs PC2 by gene

``` r
# Helper function to fetch and extract annotation
gen_phenotype_anno <- function(phenotype_term) {
  query <- utils::URLencode(paste0("https://rest.ensembl.org/phenotype/term/homo_sapiens/", phenotype_term))
  resp <- request(query) |>
    req_headers("Content-Type" = "application/json") |>
    req_perform()
  json <- resp_body_json(resp)
  lapply(json, function(x) {
    if (!is.null(x$Variation)) {
      data.frame(
        rsid = x$Variation,
        gen = if (!is.null(x$attributes$associated_gene)) x$attributes$associated_gene else NA,
        stringsAsFactors = FALSE
      )
    }
  }) |> bind_rows()
}

# Fetch annotation data frames using helper function
coffee_anno <- gen_phenotype_anno("coffee consumption")
preeclampsia_anno <- gen_phenotype_anno("preeclampsia")
ckd_anno <- gen_phenotype_anno("chronic kidney disease") |> rename(gen_ckd = gen)

# Merge PCA scores and houba info
plotDf <- cbind(pc_scores,houba$info)
names(plotDf)[1:15] <- paste0("PC",1:15)

# Annotate plotDf with coffee, preeclampsia, and CKD
plotDf <- plotDf |>
  left_join(coffee_anno, by = c("rsid" = "rsid")) |>
  left_join(preeclampsia_anno, by = c("rsid" = "rsid"), suffix = c("_coffee", "_preeclampsia")) |>
  left_join(ckd_anno, by = c("rsid" = "rsid"), suffix = c("", "_ckd"))

# Subset for coffee, preeclampsia, and CKD SNPs
coffee_snps <- subset(plotDf, !is.na(gen_coffee))|>
        unique()
preeclampsia_snps <- subset(plotDf, !is.na(gen_preeclampsia))|>
        unique()
ckd_snps <- subset(plotDf, !is.na(gen_ckd)) |>
        unique()
# Plot: coffee SNPs as circles, preeclampsia SNPs as triangles, CKD SNPs as squares
ggplot() +
  geom_point(data = coffee_snps, aes(x = PC1, y = PC2, color = gen_coffee), shape = 16, size = 2, alpha = 0.7) +
  geom_point(data = preeclampsia_snps, aes(x = PC1, y = PC2, color = gen_preeclampsia), shape = 17, size = 3, alpha = 0.7) +
  geom_point(data = ckd_snps, aes(x = PC1, y = PC2, color = gen_ckd), shape = 15, size = 3, alpha = 0.7) +
  labs(title = "PC1 vs PC2: coffee (circles), preeclampsia (triangles), CKD (squares)") +
  theme_minimal() +
  theme(legend.position = "none")
```

<img src="docs/README-unnamed-chunk-9-1.png" width="100%" />

## Notes

- All file sizes are displayed after each step for transparency.
- Remote parquet files are previewed before local processing.
- For large datasets, consider subsampling for faster plotting.
- Use `overwrite = FALSE` for DuckDB and houba steps to avoid slow
  download and reuse existing files.

## References

- The embeddings explored here were published by the group in:
  - **Incorporating LLM Embeddings for Variation Across the Human
    Genome** ([arXiv:2509.20702v1](https://arxiv.org/html/2509.20702v1))
