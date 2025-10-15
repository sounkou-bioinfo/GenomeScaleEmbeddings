
# GenomeScaleEmbeddings

## Load the package

``` r
library(GenomeScaleEmbeddings)
library(patchwork)
```

## 1. Peek into remote parquet files

``` r
# Use OpenRemoteParquetView to inspect the first few rows
OpenRemoteParquetView()
#> # Source:   table<embeddings> [?? x 6]
#> # Database: DuckDB 1.4.0 [root@Linux 6.8.0-78-generic:R 4.5.1//tmp/RtmpvN9HWX/file82b0073fce3a7.duckdb]
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

## 2. Copy remote parquet files into local DuckDB (overwrite = FALSE recommended for speed)

``` r
CopyParquetToDuckDB(db_path = "local_embeddings.duckdb", overwrite = FALSE)
#> Database file 'local_embeddings.duckdb' already exists. Skipping copy.
file.info("local_embeddings.duckdb")$size
#> [1] 12128628736
```

## 3. Write embeddings to houba mmatrix (overwrite = FALSE recommended for speed)

``` r
houba <- writeEmbeddingsHoubaFromDuckDB(dbPath = "local_embeddings.duckdb", overwrite = FALSE)
#> Warning in writeEmbeddingsHoubaFromDuckDB(dbPath = "local_embeddings.duckdb", :
#> Embedding file already exists: local_embeddings.houba. Use overwrite = TRUE to
#> overwrite.
#> Warning in mk.descriptor.file(object@file, object@dim[1], object@dim[2], :
#> local_embeddings.houba.desc already exists.
#> Writing 616386 rows in batches of 100000...
#> Warning in writeEmbeddingsHoubaFromDuckDB(dbPath = "local_embeddings.duckdb", :
#> Embedding file already exists: local_embeddings.houba. Use overwrite = TRUE to
#> overwrite.
#> Done writing embeddings and info to houba mmatrix.
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
#> 355.975 116.238  54.191
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

<img src="README_files/figure-gfm/unnamed-chunk-8-1.png" width="100%" />

## 7. Plot spatial correlation (PC1 vs genomic position, faceted by chromosome)

``` r
plots <- lapply(unique(houba$info$chrom), function(chr) {
  idx <- houba$info$chrom == chr
  ggplot2::ggplot(data.frame(
    PC1 = pc_scores[idx, 1],
    Position = houba$info$pos[idx]
  ), ggplot2::aes(x = Position, y = PC1)) +
    ggplot2::geom_point(size = 1.5, alpha = 0.8) +
    ggplot2::labs(title = paste("PC1 vs Position (chr", chr, ")"), x = paste("Genomic Position (chr", chr, ")"), y = "PC1") +
    ggplot2::theme_minimal()
})
wrap_plots(plots, ncol = 2)
```

<img src="README_files/figure-gfm/unnamed-chunk-9-1.png" width="100%" />

------------------------------------------------------------------------

- All file sizes are displayed after each step for transparency.
- Remote parquet files are previewed before local processing.
- For large datasets, consider subsampling for faster plotting.
- Use `overwrite = FALSE` for DuckDB and houba steps to avoid slow
  download and reuse existing files.
