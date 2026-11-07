library(ARISE)

run_regression_tests <- function() {
  helper <- function(name) {
    get(name, envir = environment(ARISEFS), inherits = FALSE)
  }
  expect_error <- function(expr, pattern = NULL) {
    result <- tryCatch({
      force(expr)
      NULL
    }, error = identity)
    stopifnot(inherits(result, "error"))
    if (!is.null(pattern)) {
      stopifnot(grepl(pattern, conditionMessage(result), fixed = TRUE))
    }
    invisible(NULL)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    original_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  ## A large location shift must not turn a variable predictor into a constant.
  variance <- helper(".arise_fast_variance")
  shifted <- cbind(
    shifted = 1e12 + seq_len(12L),
    original = seq_len(12L),
    constant = rep(1e12, 12L)
  )
  stopifnot(isTRUE(all.equal(
    unname(variance(shifted)), c(13, 13, 0), tolerance = 1e-12
  )))

  ## Discrete predictors retain their information, including binary features.
  mutual_information <- helper(".arise_mutual_information")
  labels <- factor(rep(c("a", "b"), each = 4L))
  scores <- mutual_information(cbind(
    perfect = rep(0:1, each = 4L),
    independent = rep(0:1, 4L),
    constant = rep(0, 8L)
  ), labels)
  stopifnot(
    abs(scores["perfect"] - log(2)) < 1e-12,
    abs(scores["independent"]) < 1e-12,
    abs(scores["constant"]) < 1e-12
  )

  ## Invalid inputs fail before expensive tuning and without coercion warnings.
  small_x <- cbind(
    signal = c(-3, -2, -1, 1, 2, 3),
    noise = c(1, 3, 2, 3, 1, 2)
  )
  small_y <- factor(rep(c("a", "b"), each = 3L))
  for (bad_top_k in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), TRUE,
                         1 + 1i, matrix(1))) {
    expect_error(ARISEFS(small_x, small_y, bad_top_k), "'topK'")
  }
  expect_error(ARISEFS(small_x + 1i, small_y), "real numeric")
  expect_error(ARISEFS(data.frame(a = seq_len(6L), b = letters[1:6]), small_y))
  matrix_column <- data.frame(a = seq_len(6L))
  matrix_column$b <- I(matrix(seq_len(12L), nrow = 6L))
  expect_error(ARISEFS(matrix_column, small_y), "Invalid columns")
  list_column <- data.frame(a = seq_len(6L))
  list_column$b <- I(as.list(seq_len(6L)))
  expect_error(ARISEFS(list_column, small_y), "Invalid columns")
  expect_error(ARISEFS(array(seq_len(12L), c(3L, 2L, 2L)), small_y))
  expect_error(ARISEFS(small_x, as.complex(seq_len(6L))), "real class labels")
  expect_error(ARISEFS(small_x, matrix(small_y, ncol = 1L)))
  expect_error(ARISEFS(small_x, small_y[-1L]), "length(y)")
  expect_error(ARISEFS(small_x, c("a", "a", "a", "b", "b", NA)),
               "missing or non-finite")
  expect_error(ARISEFS(small_x, rep("a", 6L)), "At least two")
  expect_error(ARISEFS(small_x, c("a", "a", "a", "a", "b", "b")),
               "at least 3")
  duplicate_x <- small_x
  colnames(duplicate_x) <- c("a", "a")
  expect_error(ARISEFS(duplicate_x, small_y), "unique")

  ## Imputation and feature removal use only training data.
  control <- helper(".arise_control")(tune = FALSE, stability_subsamples = 2L)
  preprocess <- helper(".arise_fit_preprocessor")
  apply_preprocessor <- helper(".arise_apply_preprocessor")
  training <- cbind(
    signal = c(1, 2, 3, 4, 5, NA_real_),
    constant = rep(1, 6L),
    unavailable = rep(NA_real_, 6L)
  )
  prepared <- preprocess(training, control)
  stopifnot(
    identical(colnames(prepared$x), "signal"),
    prepared$x[6L, "signal"] == 3
  )
  validation <- rbind(c(100, 1, NA_real_), c(NA_real_, 1, NA_real_))
  colnames(validation) <- colnames(training)
  transformed <- apply_preprocessor(prepared$preprocessor, validation)
  stopifnot(transformed[1L, 1L] == 100, transformed[2L, 1L] == 3)
  expect_error(preprocess(matrix(1, nrow = 6L, ncol = 2L), control),
               "No variable feature")
  expect_error(preprocess(matrix(NA_real_, nrow = 6L, ncol = 2L), control),
               "entirely non-finite")
  stopifnot(ncol(preprocess(shifted, control)$x) == 2L)

  ## Repeated feature selection is deterministic and preserves the caller's RNG.
  fit <- helper(".arise_fit")
  set.seed(4821)
  rng_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  first <- fit(small_x, small_y, 2L, control, seed = 13L)
  second <- fit(small_x, small_y, 2L, control, seed = 13L)
  stopifnot(
    identical(first, second),
    length(first) == 2L,
    !anyDuplicated(first),
    identical(rng_before, get(".Random.seed", envir = .GlobalEnv))
  )
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(helper(".arise_with_seed")(19L, function() stats::runif(3L)))
  stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  ## Native SVM pair order must work for labels containing ambiguous slashes.
  svm_x <- cbind(
    first = rep(c(3, 1, -1, -3), each = 4L) + rep(c(-.1, 0, .1, .2), 4L),
    second = rep(c(-3, 1, 3, -1), each = 4L) + rep(c(.2, -.1, 0, .1), 4L)
  )
  simple_labels <- paste0("class", seq_len(4L))
  slash_labels <- c("a", "b/c", "a/b", "c")
  class_index <- rep(4:1, each = 4L)
  simple_y <- factor(simple_labels[class_index], levels = simple_labels)
  slash_y <- factor(slash_labels[class_index], levels = slash_labels)
  predict_fixed <- helper(".arise_fit_predict")
  simple <- predict_fixed("svm", svm_x, simple_y, svm_x, control, 7L)
  slash <- predict_fixed("svm", svm_x, slash_y, svm_x, control, 7L)
  stopifnot(
    isTRUE(all.equal(unname(simple$scores), unname(slash$scores),
                     tolerance = 1e-10)),
    identical(as.integer(simple$class), as.integer(slash$class)),
    all(is.finite(slash$scores)),
    max(abs(rowSums(slash$scores) - 1)) < 1e-12
  )

  ## One complete public call exercises both tuning stages and all learners.
  ## A request exceeding the available predictors returns each available index.
  set.seed(8273)
  rng_before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  selected <- ARISEFS(small_x, small_y, topK = 10L)
  stopifnot(
    is.integer(selected),
    is.null(names(selected)),
    length(selected) == ncol(small_x),
    identical(sort(selected), seq_len(ncol(small_x))),
    identical(rng_before, get(".Random.seed", envir = .GlobalEnv))
  )

  ## The bundled dataset must load from the installed package.
  data_environment <- new.env(parent = baseenv())
  utils::data("NovartisMultiTissue", package = "ARISE", envir = data_environment)
  benchmark <- get("NovartisMultiTissue", envir = data_environment)
  stopifnot(
    is.list(benchmark),
    is.matrix(benchmark$X),
    is.numeric(benchmark$X),
    identical(dim(benchmark$X), c(102L, 5565L)),
    nrow(benchmark$X) == length(benchmark$Y),
    all(is.finite(benchmark$X)),
    identical(levels(benchmark$Y), as.character(seq_len(4L))),
    identical(as.integer(table(benchmark$Y)), c(25L, 26L, 28L, 23L)),
    identical(colnames(benchmark$X), paste0("Feature", seq_len(5565L))),
    identical(rownames(benchmark$X), paste0("Sample", seq_len(102L)))
  )
  invisible(TRUE)
}

run_regression_tests()
