
# GenomeScaleEmbeddings

## Load the package

``` r
library(GenomeScaleEmbeddings)
library(knitr)
```

## 1. Peek into remote parquet files

``` r
# Use OpenRemoteParquetView to inspect the first few rows
OpenRemoteParquetView()
#> # Source:   table<embeddings> [?? x 6]
#> # Database: DuckDB 1.4.0 [root@Linux 6.8.0-78-generic:R 4.5.1//tmp/RtmpizOfOU/file8797818176bec.duckdb]
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
#>   0.031   0.003   0.029
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
#>   0.555   0.050   0.562
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
#> 351.989 117.295  54.023
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

## 7. Print spatial correlations (PC1 vs genomic position, by chromosome)

``` r
cor_table <- data.frame(
  Chromosome = unique(houba$info$chrom),
  Correlation = NA,
  P_value = NA
)
for (i in seq_along(cor_table$Chromosome)) {
  chr <- cor_table$Chromosome[i]
  idx <- houba$info$chrom == chr
  pos_numeric <- as.numeric(as.character(houba$info$pos[idx]))
  pc1 <- pc_scores[idx, 1]
  # Use differences
  d_pc1 <- diff(pc1)
  d_pos <- diff(pos_numeric)
  ct <- cor.test(d_pc1, d_pos, use = "complete.obs")
  cor_table$Correlation[i] <- ct$estimate
  cor_table$P_value[i] <- ct$p.value
}
kable(cor_table, digits = 4, caption = "Correlation (and p-value) between PC1 and Genomic Position differences by Chromosome")
```

| Chromosome | Correlation | P_value |
|:-----------|------------:|--------:|
| 5          |     -0.0001 |  0.9786 |
| 6          |      0.0056 |  0.0724 |
| 7          |      0.0054 |  0.1119 |
| 8          |      0.0023 |  0.5001 |
| 9          |     -0.0038 |  0.3210 |
| 1          |     -0.0053 |  0.4246 |
| 18         |      0.0001 |  0.9828 |
| 19         |     -0.0088 |  0.1108 |
| 20         |     -0.0015 |  0.7599 |
| 2          |     -0.0022 |  0.7113 |
| 21         |      0.0046 |  0.5550 |

Correlation (and p-value) between PC1 and Genomic Position differences
by Chromosome

------------------------------------------------------------------------

- All file sizes are displayed after each step for transparency.
- Remote parquet files are previewed before local processing.
- For large datasets, consider subsampling for faster plotting.
- Use `overwrite = FALSE` for DuckDB and houba steps to avoid slow
  download and reuse existing files.
