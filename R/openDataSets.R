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
    table_name = "embeddings") {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path, overwrite = FALSE)
    if (!file.exists(db_path) || overwrite) {
        sql <- sprintf(
            "CREATE TABLE %s AS SELECT * FROM read_parquet(%s)",
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
#' @return An object of class 'HoubaEmbeddings' with mmatrix, info data.frame, and houba file path
#' @export
writeEmbeddingsHoubaFromDuckDB <- function(
    dbPath = "local_embeddings.duckdb",
    tableName = "embeddings",
    embeddingCol = "embedding",
    batchSize = 100000,
    embeddingDim = 3072,
    embeddingFile = gsub("\\.duckdb$", ".houba", dbPath)) {
    if (!requireNamespace("houba", quietly = TRUE)) stop("houba package required.")
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
