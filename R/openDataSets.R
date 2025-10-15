#' List of the huggingface datasets for the paper
#' "Incorporating LLM Embeddings for Variation Across the Human Genome"
#' @return A character vector of dataset names
#' @export
DatasetParquetUrlList <- function() {
    # curl -X GET "https://huggingface.co/api/datasets/hong-niu/Variant-Foundation-Embeddings/parquet/default/train"
    urlList <- readLines(
        paste0(
            "https://huggingface.co/api/datasets/hong-niu/",
            "Variant-Foundation-Embeddings/parquet/default/train"
        ),
        warn = FALSE
    ) |>
        jsonlite::fromJSON(simplifyVector = FALSE) |>
        unlist()
    stopifnot(length(urlList) == 2)
    urlList
}

#' Copy remote parquet files into a local DuckDB database file
#' @param db_path Path to the local DuckDB database file
#' @param urlList Character vector of parquet URLs
#' @param table_name Name of the table to create in the DuckDB database
#' @param overwrite Whether to overwrite the existing database
#' @export
CopyParquetToDuckDB <- function(
    db_path = "local_embeddings.duckdb",
    urlList = DatasetParquetUrlList(),
    table_name = "embeddings",
    overwrite = FALSE) {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path, overwrite = FALSE)
    table_names <- DBI::dbListTables(con)
    if (!table_name %in% table_names || overwrite) {
        sql <- sprintf(
            "INSTALL httpfs; load httpfs;
             CREATE TABLE %s AS SELECT * FROM read_parquet(%s)",
            table_name,
            paste0("['", paste(urlList, collapse = "','"), "']")
        )

        DBI::dbExecute(con, sql)
        message(
            "Copied parquet files to DuckDB table '",
            table_name,
            "' in database '",
            db_path,
            "'."
        )
        # peak at the data
        DBI::dbDisconnect(con)
    } else {
        message("Database file '", db_path, "' already exists. Skipping copy.")
    }
    invisible(db_path)
}

#' Iterate over embeddings as matrix batches from a local DuckDB file
#' @param chunk_size Number of rows per batch
#' @param db_path Path to the local DuckDB database file
#' @param table_name Name of the table in the DuckDB database
#' @return A list of embedding batches
#' @export
IterateEmbeddingsMatrixBatches <- function(
    chunk_size = 100000,
    db_path = "local_embeddings.duckdb",
    table_name = "embeddings") {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
    meta <- DBI::dbGetQuery(
        con,
        sprintf(
            "SELECT COUNT(*) AS n FROM %s",
            table_name
        )
    )
    n <- meta$n
    message("Total number of rows: ", n)
    batches <- list()
    for (offset in seq(0, n - 1, by = chunk_size)) {
        query <- sprintf(
            "SELECT embedding FROM %s LIMIT %d OFFSET %d",
            table_name,
            chunk_size,
            offset
        )
        message(
            "Processing rows ",
            offset + 1,
            " to ",
            min(offset + chunk_size, n)
        )
        batch <- DBI::dbGetQuery(con, query)
        batches[[length(batches) + 1]] <- batch$embedding
    }
    DBI::dbDisconnect(con)
    batches
}

#' Open remote parquet files as a DuckDB VIEW and return as tibble (minimal, http(s) only)
#' @param urlList Character vector of parquet URLs
#' @param view_name Name of the DuckDB view to create
#' @param db_path Path to DuckDB database file (default: temporary file, null for in-memory)
#' @param unify_schemas Whether to unify schemas across files
#' @return dplyr tibble referencing the DuckDB view
#' @export
OpenRemoteParquetView <- function(
    urlList = DatasetParquetUrlList(),
    view_name = "embeddings",
    db_path = tempfile(fileext = ".duckdb"),
    unify_schemas = FALSE) {
    if (is.null(db_path)) {
        conn <- DBI::dbConnect(duckdb::duckdb())
    } else {
        conn <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
    }
    sql <- sprintf(
        "CREATE OR REPLACE TEMPORARY VIEW %s AS SELECT * FROM parquet_scan(%s, UNION_BY_NAME=%s)",
        view_name,
        paste0("['", paste(urlList, collapse = "','"), "']"),
        tolower(unify_schemas)
    )
    DBI::dbExecute(conn, sql)
    dplyr::tbl(conn, view_name)
}

#' Write embeddings to houba mmatrix and return info as data.frame from a local DuckDB file
#' @param dbPath Path to the local DuckDB database file
#' @param tableName Name of the table in the DuckDB database
#' @param embeddingCol Name of the embeddings column
#' @param batchSize Number of rows per batch
#' @param embeddingDim Dimension of each embedding vector
#' @param embeddingFile Path for houba mmatrix file
#' @param overwrite Whether to overwrite existing houba file
#' @return An object of class 'HoubaEmbeddings' with mmatrix, info data.frame, and houba file path
#' @export
writeEmbeddingsHoubaFromDuckDB <- function(
    dbPath = "local_embeddings.duckdb",
    tableName = "embeddings",
    embeddingCol = "embedding",
    batchSize = 100000,
    embeddingDim = 3072,
    embeddingFile = gsub("\\.duckdb$", ".houba", dbPath),
    overwrite = FALSE) {
    if (!requireNamespace("houba", quietly = TRUE)) stop("houba package required.")
    if (file.exists(embeddingFile) && !overwrite) {
        stop("Embedding file already exists: ", embeddingFile, ". Use overwrite = TRUE to overwrite.")
    }
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = dbPath)
    n <- DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM %s", tableName))$n
    # Get info column names from table
    infoColNames <- setdiff(DBI::dbListFields(con, tableName), embeddingCol)
    embMat <- houba::mmatrix(datatype = "double", nrow = n, ncol = embeddingDim, filename = embeddingFile, readonly = FALSE)
    # Create descriptor file for bigmemory compatibility
    dscFile <- houba::descriptor.file(embMat)
    infoDf <- data.frame(matrix(nrow = n, ncol = length(infoColNames)))
    colnames(infoDf) <- infoColNames
    idx <- 1
    message(sprintf("Writing %d rows in batches of %d...", n, batchSize))
    while (idx <= n) {
        rows <- min(batchSize, n - idx + 1)
        query <- sprintf(
            "SELECT %s, %s FROM %s LIMIT %d OFFSET %d",
            embeddingCol,
            paste(infoColNames, collapse = ", "),
            tableName,
            rows,
            idx - 1
        )
        batch <- DBI::dbGetQuery(con, query)
        # Assign each embedding vector directly to the mmatrix row
        for (row in seq_len(rows)) {
            embMat[idx + row - 1, ] <- batch[[embeddingCol]][[row]]
        }
        # Directly assign info columns to infoDf
        for (j in seq_along(infoColNames)) {
            infoDf[idx:(idx + rows - 1), j] <- batch[[infoColNames[j]]]
        }
        message(sprintf("Processed rows %d to %d", idx, idx + rows - 1))
        idx <- idx + rows
    }
    DBI::dbDisconnect(con)
    message("Done writing embeddings and info to houba mmatrix.")
    structure(
        list(embeddings = embMat, info = infoDf, houba_file = embeddingFile, descriptor_file = dscFile),
        class = "HoubaEmbeddings"
    )
}

#' @export
print.HoubaEmbeddings <- function(x, ...) {
    cat("Houba mmatrix file:", x$houba_file, "\n")
    cat("Embeddings (houba::mmatrix):\n")
    print(x$embeddings)
    cat("\nInfo data.frame (first 10 rows):\n")
    print(utils::head(x$info, 10))
    invisible(x)
}

#' Quick summary for houba mmatrix
#' @param embMat houba mmatrix
#' @return List with dim, colMeans, rowMeans
#' @export
embeddingSummary <- function(embMat) {
    list(
        dim = dim(embMat),
        colMeans = houba::colMeans(embMat),
        rowMeans = houba::rowMeans(embMat)
    )
}

#' Quick summary for houba info mmatrix
#' @param infoMat houba mmatrix
#' @return List with dim, unique chroms, and example rsids
#' @export
infoSummary <- function(infoMat) {
    list(
        dim = base::dim(infoMat),
        chroms = base::unique(infoMat[, 1]),
        rsids = base::head(infoMat[, 5])
    )
}


#' PCA using bigPCAcpp on houba mmatrix or bigmemory::big.matrix
#' @param embMat houba mmatrix path or bigmemory::big.matrix
#' @param center logical, whether to center columns
#' @param scale logical, whether to scale columns
#' @param ncomp number of principal components to compute
#' @return List with PCA result object and houbaM big.matrix
#' @export
houbaPCA <- function(embMat = "local_embeddings.houba", center = TRUE, scale = TRUE, ncomp = 15) {
    if (is.character(embMat)) {
        descFile <- paste0(embMat, ".desc")
        if (!file.exists(descFile)) {
            stop("Descriptor file does not exist: ", descFile)
        }
        houbaM <- bigmemory::attach.big.matrix(descFile)
        message("Attached big.matrix from descriptor file: ", descFile)
        message("Dimensions: ", paste(dim(houbaM), collapse = " x "))
    } else {
        if (!inherits(embMat, "big.matrix")) {
            stop("embMat must be a bigmemory::big.matrix or path to descriptor file.")
        }
        houbaM <- embMat
    }
    message("Running PCA with center=", center, ", scale=", scale, ", ncomp=", ncomp)
    pca <- bigPCAcpp::pca_bigmatrix(houbaM, center = center, scale = scale, ncomp = ncomp)
    list(pca = pca, houbaM = houbaM)
}

#' Get PCA scores from houbaPCA result
#' @param houbaPCA_res List returned by houbaPCA (with 'pca' and 'houbaM')
#' @return Matrix of principal component scores (variants x PCs)
#' @export
getPcaScores <- function(houbaPCA_res) {
    bigPCAcpp::pca_scores_bigmatrix(
        xpMat = houbaPCA_res$houbaM,
        rotation = houbaPCA_res$pca$rotation,
        center = houbaPCA_res$pca$center,
        scale = houbaPCA_res$pca$scale
    )
}

#' Plot spatial correlation between PC scores and genomic position, faceted by chromosome
#' @param pc_scores Matrix of principal component scores (variants x PCs)
#' @param info_df Data frame with variant info (must contain 'chrom' and 'pos')
#' @param pc Which principal component to plot (default: 1)
#' @export
plotPCSpatialCorrelation <- function(pc_scores, info_df, pc = 1) {
    chroms <- as.factor(info_df$chrom)
    positions <- info_df$pos
    oldpar <- par(mfrow = c(ceiling(sqrt(length(levels(chroms)))), ceiling(sqrt(length(levels(chroms))))))
    on.exit(par(oldpar))
    for (chr in levels(chroms)) {
        idx <- which(chroms == chr)
        plot(positions[idx], pc_scores[idx, pc],
            col = chr, pch = 16,
            xlab = paste0("Genomic Position (chr ", chr, ")"),
            ylab = paste0("PC", pc),
            main = paste0("PC", pc, " vs Position (chr ", chr, ")")
        )
    }
}

#' Compute correlation between PC scores and genomic position, per chromosome
#' @param pc_scores Matrix of principal component scores (variants x PCs)
#' @param info_df Data frame with variant info (must contain 'chrom' and 'pos')
#' @param pc Which principal component to correlate (default: 1)
#' @param method Correlation method (default: 'spearman')
#' @return Named vector of correlation values per chromosome
#' @export
correlatePCWithPosition <- function(pc_scores, info_df, pc = 1, method = "spearman") {
    chroms <- as.factor(info_df$chrom)
    positions <- info_df$pos
    corrs <- sapply(levels(chroms), function(chr) {
        idx <- which(chroms == chr)
        if (length(idx) > 1) {
            stats::cor(pc_scores[idx, pc], positions[idx], method = method)
        } else {
            NA
        }
    })
    names(corrs) <- levels(chroms)
    corrs
}

#' Plot PCA dimensions using ggplot2, colored by annotation
#' @param pc_scores Matrix of principal component scores (variants x PCs)
#' @param info_df Data frame with variant info (must contain annotation column)
#' @param annotation_col Name of column in info_df to color by (e.g. 'gwas')
#' @param dim1 First PC dimension to plot (default: 1)
#' @param dim2 Second PC dimension to plot (default: 2)
#' @export
plotPcaDims <- function(pc_scores, info_df, annotation_col = "chrom", dim1 = 1, dim2 = 2) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package required.")
    df <- data.frame(
        PC1 = pc_scores[, dim1],
        PC2 = pc_scores[, dim2],
        annotation = as.factor(info_df[[annotation_col]])
    )
    ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, color = annotation)) +
        ggplot2::geom_point(alpha = 0.7) +
        ggplot2::labs(x = paste0("PC", dim1), y = paste0("PC", dim2), color = annotation_col) +
        ggplot2::theme_minimal()
}

#' Attach a houba file and return the bigmemory::big.matrix
#' @param houba_file Path to houba file (without .desc extension)
#' @return bigmemory::big.matrix object
#' @export
attachHoubaBigMatrix <- function(houba_file) {
    descFile <- paste0(houba_file, ".desc")
    if (!file.exists(descFile)) {
        stop("Descriptor file does not exist: ", descFile)
    }
    bigmemory::attach.big.matrix(descFile)
}
