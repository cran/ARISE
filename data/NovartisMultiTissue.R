NovartisMultiTissue <- local({
  candidates <- c(
    file.path("..", "inst", "extdata", "novartis_multi_tissue.csv.gz"),
    file.path("inst", "extdata", "novartis_multi_tissue.csv.gz"),
    system.file("extdata", "novartis_multi_tissue.csv.gz", package = "ARISE")
  )
  file <- candidates[file.exists(candidates)]
  if (!length(file)) {
    stop("The bundled Novartis dataset file could not be found.", call. = FALSE)
  }
  connection <- gzfile(file[1L], open = "rt")
  values <- tryCatch(
    scan(connection, what = double(), sep = ",", skip = 1L, quiet = TRUE),
    finally = close(connection)
  )
  n_columns <- 5566L
  if (length(values) != 102L * n_columns) {
    stop("Unexpected size of the bundled Novartis dataset.", call. = FALSE)
  }
  dat <- matrix(values, ncol = n_columns, byrow = TRUE)
  X <- dat[, -ncol(dat), drop = FALSE]
  storage.mode(X) <- "double"
  colnames(X) <- paste0("Feature", seq_len(ncol(X)))
  rownames(X) <- paste0("Sample", seq_len(nrow(X)))
  Y <- factor(dat[, ncol(dat)])
  list(X = X, Y = Y)
})
