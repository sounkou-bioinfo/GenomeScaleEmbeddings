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
