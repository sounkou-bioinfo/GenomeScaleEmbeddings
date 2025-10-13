#' List of the huggingface datasets for the paper
#' "Incorporating LLM Embeddings for Variation Across the Human Genome"
#' @return A character vector of dataset names
#' @export
datasetParquetUrlList <- function() {
    # curl -X GET "https://huggingface.co/api/datasets/hong-niu/Variant-Foundation-Embeddings/parquet/default/train"
    urlList <- readLines(paste0(
        "https://huggingface.co/api/datasets/hong-niu/",
        "Variant-Foundation-Embeddings/parquet/default/train"
    ), warn = FALSE) |>
        jsonlite::fromJSON(simplifyVector = FALSE) |>
        unlist()
    stopifnot(length(urlList) == 2)
    urlList
}

#' Open the parquet files from the huggingface dataset
#'  using duckdbfs
#' @param urlList A character vector of dataset names
#' @param unify_schemas Whether to unify the schemas of the datasets
#' @param ... Additional arguments passed to duckdbfs::open_dataset
#' @return A duckdbfs connection to the parquet files
#' @export
openParquetFiles <- function(
    urlList = datasetParquetUrlList(),
    unify_schemas = FALSE,
    ...) {
    stopifnot(length(urlList) == 2)
    duckdbfs::open_dataset(
        urlList,
        unify_schemas = unify_schemas,
        ...
    )
}


#' Copy remote parquet files into a local DuckDB database file
#' @param db_path Path to the local DuckDB database file
#' @param urlList Character vector of parquet URLs
#' @param table_name Name of the table to create in the DuckDB database
#' @param overwrite Whether to overwrite the existing database
#' @export
copyParquetToDuckDB <- function(db_path = "local_embeddings.duckdb", urlList = datasetParquetUrlList(), table_name = "embeddings") {
    con <- DBI::dbConnect(duckdb::duckdb(),
        dbdir = db_path,
        overwrite = FALSE
    )
    if (!file.exists(db_path) || overwrite) {
        sql <- sprintf(
            "CREATE TABLE %s AS SELECT * FROM read_parquet(%s)",
            table_name,
            paste0("['", paste(urlList, collapse = "','"), "']")
        )

        DBI::dbExecute(con, sql)
        message("Copied parquet files to DuckDB table '", table_name, "' in database '", db_path, "'.")
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
iterateEmbeddingsMatrixBatches <- function(
    chunk_size = 100000, db_path = "local_embeddings.duckdb",
    table_name = "embeddings") {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = db_path)
    meta <- DBI::dbGetQuery(con, sprintf(
        "SELECT COUNT(*) AS n FROM %s",
        table_name
    ))
    n <- meta$n
    message("Total number of rows: ", n)
    batches <- list()
    for (offset in seq(0, n - 1, by = chunk_size)) {
        query <- sprintf(
            "SELECT embedding FROM %s LIMIT %d OFFSET %d",
            table_name, chunk_size, offset
        )
        message(
            "Processing rows ", offset + 1, " to ",
            min(offset + chunk_size, n)
        )
        batch <- DBI::dbGetQuery(con, query)
        batches[[length(batches) + 1]] <- batch$embedding
    }
    DBI::dbDisconnect(con)
    batches
}
