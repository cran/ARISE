ARISEFS <- local({
  .arise_stop <- function(...) {
    stop(..., call. = FALSE)
  }


  .arise_scalar_logical <- function(x, name) {
    if (!is.logical(x) || !is.null(dim(x)) || length(x) != 1L || is.na(x)) {
      .arise_stop("'", name, "' must be TRUE or FALSE.")
    }
    x
  }

  .arise_scalar_integer <- function(x, name, minimum = -Inf, maximum = Inf) {
    valid <- is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
      length(x) == 1L && !is.na(x) &&
      is.finite(x) && x == floor(x) && x >= minimum && x <= maximum &&
      x >= -.Machine$integer.max && x <= .Machine$integer.max
    if (!valid) {
      bounds <- if (is.finite(minimum) && is.finite(maximum)) {
        paste0(" between ", minimum, " and ", maximum)
      } else if (is.finite(minimum)) {
        paste0(" >= ", minimum)
      } else if (is.finite(maximum)) {
        paste0(" <= ", maximum)
      } else {
        ""
      }
      .arise_stop("'", name, "' must be a whole-number integer", bounds, ".")
    }
    as.integer(x)
  }

  .arise_integer_vector <- function(x, name, minimum = 1L) {
    valid <- is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
      length(x) > 0L && !anyNA(x) &&
      all(is.finite(x)) && all(x == floor(x)) && all(x >= minimum) &&
      all(x <= .Machine$integer.max)
    if (!valid) {
      .arise_stop(
        "'", name, "' must contain whole-number integers >= ", minimum, "."
      )
    }
    sort(unique(as.integer(x)))
  }

  .arise_scalar_numeric <- function(x, name, lower = -Inf, upper = Inf,
                                    lower_open = FALSE, upper_open = FALSE) {
    valid <- is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
      length(x) == 1L && !is.na(x) && is.finite(x)
    if (valid) {
      valid <- if (lower_open) x > lower else x >= lower
    }
    if (valid) {
      valid <- if (upper_open) x < upper else x <= upper
    }
    if (!valid) {
      .arise_stop("Invalid value supplied for '", name, "'.")
    }
    as.numeric(x)
  }

  .arise_feature_cap <- function(x) {
    valid_type <- is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
      length(x) == 1L
    valid_inf <- valid_type && !is.na(x) && is.infinite(x) && x > 0
    valid_integer <- valid_type && !is.na(x) &&
      is.finite(x) && x == floor(x) && x >= 1 && x <= .Machine$integer.max
    if (!(valid_inf || valid_integer)) {
      .arise_stop("'feature_cap' must be a positive whole-number integer or Inf.")
    }
    if (valid_inf) Inf else as.integer(x)
  }

  .arise_seed <- function(seed, ...) {
    values <- c(as.double(seed), as.double(unlist(list(...), use.names = FALSE)))
    values[!is.finite(values)] <- 0
    index <- seq_along(values)
    value <- sum((index * 104729 + 1009) * (abs(values) + 1))
    as.integer(value %% 2147483629 + 1)
  }

  .arise_with_seed <- function(seed, fun) {
    if (!is.function(fun)) {
      .arise_stop("'fun' must be a function.")
    }
    seed <- as.integer(seed)[1L]
    if (!is.finite(seed)) {
      seed <- 1L
    }
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    fun()
  }

  .arise_fast_variance <- function(x) {
    x <- as.matrix(x)
    n <- nrow(x)
    if (n < 2L) {
      out <- rep(0, ncol(x))
    } else {
      out <- matrixStats::colVars(x)
      out[!is.finite(out) | out < 0] <- 0
    }
    stats::setNames(out, colnames(x))
  }

  .arise_order_desc <- function(x, n = length(x), tie_index = seq_along(x)) {
    x <- as.numeric(x)
    tie_index <- as.numeric(tie_index)
    if (length(x) != length(tie_index)) {
      .arise_stop("Tie-breaking index has the wrong length.")
    }
    n <- min(length(x), max(0L, as.integer(n)[1L]))
    if (n == 0L) {
      return(integer())
    }
    x[!is.finite(x)] <- -Inf
    utils::head(order(-x, tie_index, method = "radix"), n)
  }

  .arise_percentile_rank <- function(x) {
    x <- as.numeric(x)
    x[!is.finite(x)] <- 0
    n <- length(x)
    if (n <= 1L || diff(range(x)) <= sqrt(.Machine$double.eps)) {
      return(rep(0, n))
    }
    (rank(x, ties.method = "average") - 1) / (n - 1)
  }

  .arise_percentile_matrix <- function(x) {
    x <- as.matrix(x)
    if (ncol(x) == 0L) {
      return(x)
    }
    out <- vapply(
      seq_len(ncol(x)),
      function(j) .arise_percentile_rank(x[, j]),
      numeric(nrow(x))
    )
    out <- matrix(
      out,
      nrow = nrow(x),
      ncol = ncol(x),
      dimnames = dimnames(x)
    )
    out
  }

  .arise_numeric_matrix <- function(x, name = "x") {
    if (is.data.frame(x)) {
      ok <- vapply(x, function(z) {
        (is.numeric(z) || is.logical(z)) && !is.complex(z) && is.null(dim(z))
      }, logical(1))
      if (!all(ok)) {
        .arise_stop(
          "All columns of '", name,
          "' must be real numeric or logical vectors. Invalid columns: ",
          paste(names(x)[!ok], collapse = ", "), "."
        )
      }
      x <- as.matrix(x)
    } else if (is.matrix(x)) {
      if (!(is.numeric(x) || is.logical(x)) || is.complex(x)) {
        .arise_stop("'", name, "' must be a real numeric matrix or data frame.")
      }
    } else {
      .arise_stop("'", name, "' must be a real numeric matrix or data frame.")
    }
    storage.mode(x) <- "double"
    if (nrow(x) == 0L || ncol(x) == 0L) {
      .arise_stop("'", name, "' must contain at least one row and one column.")
    }
    if (is.null(colnames(x))) {
      colnames(x) <- paste0("Feature_", seq_len(ncol(x)))
    }
    colnames(x) <- as.character(colnames(x))
    if (anyNA(colnames(x)) || any(!nzchar(colnames(x)))) {
      .arise_stop("All feature names must be non-missing and non-empty.")
    }
    if (anyDuplicated(colnames(x))) {
      .arise_stop("Feature names must be unique.")
    }
    x
  }

  .arise_validate_input <- function(x, y, tuning = TRUE, inner_folds = 3L) {
    x <- .arise_numeric_matrix(x, "x")
    if (!is.atomic(y) || !is.null(dim(y)) || is.complex(y)) {
      .arise_stop("'y' must be an atomic outcome vector of real class labels.")
    }
    if (length(y) != nrow(x)) {
      .arise_stop("length(y) must equal nrow(x).")
    }
    if (anyNA(y) || (is.numeric(y) && any(!is.finite(y)))) {
      .arise_stop("'y' must not contain missing or non-finite values.")
    }
    y <- droplevels(factor(y))
    if (nlevels(y) < 2L) {
      .arise_stop("At least two outcome classes are required.")
    }
    counts <- table(y)
    minimum <- if (isTRUE(tuning)) {
      folds <- .arise_scalar_integer(inner_folds, "inner_folds", minimum = 2L)
      max(3L, folds)
    } else {
      2L
    }
    if (any(counts < minimum)) {
      suffix <- if (isTRUE(tuning)) {
        paste0(
          " when tune = TRUE and inner_folds = ", inner_folds,
          ", so every inner training split retains at least two observations per class."
        )
      } else {
        " when tune = FALSE."
      }
      .arise_stop(
        "Every class must contain at least ", minimum,
        " observations", suffix
      )
    }
    list(x = x, y = y)
  }

  .arise_validate_processed <- function(x, y) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    y <- droplevels(factor(y))
    if (nrow(x) != length(y)) {
      .arise_stop("Processed feature matrix and outcome have incompatible sizes.")
    }
    if (ncol(x) == 0L) {
      .arise_stop("No feature remains after preprocessing.")
    }
    if (any(!is.finite(x))) {
      .arise_stop("The processed feature matrix contains non-finite values.")
    }
    if (any(table(y) < 2L)) {
      .arise_stop("Every class needs at least two observations for feature scoring.")
    }
    list(x = x, y = y)
  }

  .arise_effective_budgets <- function(budgets, maximum) {
    budgets <- sort(unique(as.integer(budgets)))
    budgets <- budgets[is.finite(budgets) & budgets > 0L]
    if (length(budgets) == 0L) {
      .arise_stop("At least one positive feature budget is required.")
    }
    pmin(budgets, as.integer(maximum))
  }

  .arise_inverse_frequency_weights <- function(y) {
    counts <- table(y)
    weights <- length(y) / (length(counts) * as.numeric(counts))
    stats::setNames(weights, names(counts))
  }

  .arise_normalize_classifier <- function(x) {
    key <- gsub("[^a-z]", "", tolower(as.character(x)[1L]))
    if (key %in% c("knn", "knearestneighbors", "knearestneighbours")) {
      return("knn")
    }
    if (key %in% c("svm", "supportvectormachine")) {
      return("svm")
    }
    if (key %in% c("rf", "randomforest")) {
      return("rf")
    }
    .arise_stop("Unknown fixed classifier: ", x, ".")
  }

  .arise_log <- function(control, ...) {
    if (isTRUE(control$verbose)) {
      message(...)
    }
    invisible(NULL)
  }


  .arise_weight_columns <- function() {
    c(
      POS = "W_POS",
      Fisher = "W_Fisher",
      SNR = "W_SNR",
      Welch = "W_Welch",
      Bhattacharyya = "W_Bhattacharyya",
      RankAUC = "W_RankAUC",
      MI = "W_MI"
    )
  }

  .arise_profiles <- function(screen_max = 700L) {
    screen_max <- .arise_scalar_integer(
      screen_max, "screen_max", minimum = 1L
    )
    mixed <- rbind(
      balanced = c(.16, .14, .12, .12, .10, .12, .12),
      location_shift = c(.08, .20, .18, .22, .10, .08, .06),
      distributional = c(.25, .10, .06, .08, .20, .18, .07),
      nonlinear = c(.06, .08, .05, .07, .08, .08, .28),
      stability_first = rep(.125, 7L),
      diversity_first = c(.12, .12, .10, .12, .10, .12, .14),
      low_redundancy = c(.18, .15, .12, .12, .10, .13, .10),
      maximin = c(.15, .12, .10, .12, .12, .14, .13)
    )
    mixed <- mixed / rowSums(mixed)
    pure <- diag(7L)
    rownames(pure) <- c(
      "pure_POS", "pure_Fisher", "pure_SNR", "pure_Welch",
      "pure_Bhattacharyya", "pure_rank_AUC", "pure_MI"
    )
    weights <- rbind(mixed, pure)
    colnames(weights) <- unname(.arise_weight_columns())

    profile_names <- rownames(weights)
    mixed_count <- nrow(mixed)
    screen_factor <- c(.90, .75, 1.00, 1.00, .90, 1.00, .75, 1.00, rep(1, 7L))
    out <- data.frame(
      Config_ID = sprintf("ARISE_%02d", seq_len(nrow(weights))),
      Profile = profile_names,
      Priority = seq_len(nrow(weights)),
      weights,
      Stability_Weight = c(
        .25, .28, .25, .20, .42, .32, .22, .35, rep(0, 7L)
      ),
      Redundancy_Lambda = c(
        .18, .12, .22, .08, .18, .28, .05, .15, rep(0, 7L)
      ),
      Redundancy_Decay = c(
        .55, .55, .50, .65, .55, .45, .70, .55, rep(0, 7L)
      ),
      Pair_Coverage_Weight = c(
        .08, .06, .14, .08, .12, .12, .04, .18, rep(0, 7L)
      ),
      Screen_Max = pmax(1L, as.integer(round(screen_factor * screen_max))),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    rownames(out) <- NULL
    if (nrow(out) != mixed_count + 7L) {
      .arise_stop("Internal profile construction failed.")
    }
    out
  }

  .arise_weight_matrix <- function(profiles) {
    columns <- .arise_weight_columns()
    missing <- setdiff(unname(columns), names(profiles))
    if (length(missing)) {
      .arise_stop("Profile table is missing: ", paste(missing, collapse = ", "), ".")
    }
    weights <- as.matrix(profiles[, unname(columns), drop = FALSE])
    storage.mode(weights) <- "double"
    weights[!is.finite(weights) | weights < 0] <- 0
    totals <- rowSums(weights)
    if (any(totals <= 0)) {
      .arise_stop("Every ARISE profile must have positive component weight.")
    }
    weights <- weights / totals
    colnames(weights) <- names(columns)
    rownames(weights) <- profiles$Config_ID
    weights
  }


  .arise_control <- function(
      tune = TRUE,
      inner_folds = 3L,
      budgets = c(1L, 5L, 10L, 15L, 20L, 25L, 30L, 35L),
      stage1_budgets = c(5L, 15L),
      stability_subsamples = 8L,
      stability_fraction = 0.70,
      pos_whisker = 1.5,
      mi_bins = 3L,
      screen_max = 700L,
      screen_multiplier = 15L,
      stage2_profiles = 3L,
      softmax_temperature = 0.02,
      max_regret_weight = 0.70,
      mean_regret_weight = 0.30,
      one_se_multiplier = 1,
      knn_k = 5L,
      svm_cost = 1,
      rf_trees = 120L,
      class_weights = TRUE,
      variance_tolerance = 1e-12,
      feature_cap = Inf,
      untuned_profile = "balanced",
      store_components = TRUE,
      verbose = FALSE) {
    tune <- .arise_scalar_logical(tune, "tune")
    inner_folds <- .arise_scalar_integer(
      inner_folds, "inner_folds", minimum = 2L
    )
    budgets <- .arise_integer_vector(budgets, "budgets", minimum = 1L)
    stage1_budgets <- .arise_integer_vector(
      stage1_budgets, "stage1_budgets", minimum = 1L
    )
    stability_subsamples <- .arise_scalar_integer(
      stability_subsamples, "stability_subsamples", minimum = 1L
    )
    stability_fraction <- .arise_scalar_numeric(
      stability_fraction, "stability_fraction", 0, 1,
      lower_open = TRUE, upper_open = FALSE
    )
    pos_whisker <- .arise_scalar_numeric(pos_whisker, "pos_whisker", 0, Inf)
    mi_bins <- .arise_scalar_integer(mi_bins, "mi_bins", minimum = 2L)
    screen_max <- .arise_scalar_integer(screen_max, "screen_max", minimum = 1L)
    screen_multiplier <- .arise_scalar_integer(
      screen_multiplier, "screen_multiplier", minimum = 1L
    )
    stage2_profiles <- .arise_scalar_integer(
      stage2_profiles, "stage2_profiles", minimum = 1L, maximum = 15L
    )
    softmax_temperature <- .arise_scalar_numeric(
      softmax_temperature, "softmax_temperature", 0, Inf, lower_open = TRUE
    )
    max_regret_weight <- .arise_scalar_numeric(
      max_regret_weight, "max_regret_weight", 0, Inf
    )
    mean_regret_weight <- .arise_scalar_numeric(
      mean_regret_weight, "mean_regret_weight", 0, Inf
    )
    if (max_regret_weight + mean_regret_weight <= 0) {
      .arise_stop("The two regret weights cannot both be zero.")
    }
    total <- max_regret_weight + mean_regret_weight
    max_regret_weight <- max_regret_weight / total
    mean_regret_weight <- mean_regret_weight / total
    one_se_multiplier <- .arise_scalar_numeric(
      one_se_multiplier, "one_se_multiplier", 0, Inf
    )
    knn_k <- .arise_scalar_integer(knn_k, "knn_k", minimum = 1L)
    svm_cost <- .arise_scalar_numeric(
      svm_cost, "svm_cost", 0, Inf, lower_open = TRUE
    )
    rf_trees <- .arise_scalar_integer(rf_trees, "rf_trees", minimum = 1L)
    class_weights <- .arise_scalar_logical(class_weights, "class_weights")
    variance_tolerance <- .arise_scalar_numeric(
      variance_tolerance, "variance_tolerance", 0, Inf
    )
    feature_cap <- .arise_feature_cap(feature_cap)
    if (!is.character(untuned_profile) || length(untuned_profile) != 1L ||
        is.na(untuned_profile) || !nzchar(untuned_profile)) {
      .arise_stop("'untuned_profile' must be one non-empty profile name.")
    }
    valid_profiles <- .arise_profiles(screen_max)$Profile
    if (!untuned_profile %in% valid_profiles) {
      .arise_stop(
        "Unknown 'untuned_profile'. Available profiles: ",
        paste(valid_profiles, collapse = ", "), "."
      )
    }
    store_components <- .arise_scalar_logical(
      store_components, "store_components"
    )
    verbose <- .arise_scalar_logical(verbose, "verbose")

    structure(list(
      tune = tune,
      inner_folds = inner_folds,
      budgets = budgets,
      stage1_budgets = stage1_budgets,
      stage1_classifiers = c("knn", "svm"),
      stage2_classifiers = c("knn", "svm", "rf"),
      stability_subsamples = stability_subsamples,
      stability_fraction = stability_fraction,
      pos_whisker = pos_whisker,
      mi_bins = mi_bins,
      screen_max = screen_max,
      screen_multiplier = screen_multiplier,
      stage2_profiles = stage2_profiles,
      softmax_temperature = softmax_temperature,
      max_regret_weight = max_regret_weight,
      mean_regret_weight = mean_regret_weight,
      one_se_multiplier = one_se_multiplier,
      knn_k = knn_k,
      svm_cost = svm_cost,
      rf_trees = rf_trees,
      class_weights = class_weights,
      variance_tolerance = variance_tolerance,
      feature_cap = feature_cap,
      untuned_profile = untuned_profile,
      store_components = store_components,
      verbose = verbose
    ), class = ".arise_control")
  }

  .arise_revalidate_control <- function(control) {
    if (!inherits(control, ".arise_control") || !is.list(control)) {
      .arise_stop("'control' must be created by .arise_control().")
    }
    argument_names <- names(formals(.arise_control))
    missing <- setdiff(argument_names, names(control))
    if (length(missing)) {
      .arise_stop(
        "The ARISE control object is missing: ", paste(missing, collapse = ", "), "."
      )
    }
    arguments <- unclass(control)[argument_names]
    do.call(.arise_control, arguments)
  }


  .arise_fit_preprocessor <- function(x_train, control) {
    x_train <- .arise_numeric_matrix(x_train, "x_train")
    input_features <- colnames(x_train)
    finite_count <- colSums(is.finite(x_train))
    finite_features <- input_features[finite_count > 0L]
    if (length(finite_features) == 0L) {
      .arise_stop("Every input feature is entirely non-finite in the training data.")
    }
    x <- x_train[, finite_features, drop = FALSE]
    medians <- vapply(seq_len(ncol(x)), function(j) {
      value <- x[, j]
      stats::median(value[is.finite(value)], na.rm = TRUE)
    }, numeric(1))
    medians[!is.finite(medians)] <- 0
    names(medians) <- finite_features
    for (j in seq_len(ncol(x))) {
      bad <- !is.finite(x[, j])
      if (any(bad)) {
        x[bad, j] <- medians[j]
      }
    }
    variance <- .arise_fast_variance(x)
    variable <- is.finite(variance) & variance > control$variance_tolerance
    variable_features <- names(variance)[variable]
    if (length(variable_features) == 0L) {
      .arise_stop("No variable feature remains after training-only preprocessing.")
    }
    x <- x[, variable_features, drop = FALSE]
    variance <- variance[variable_features]
    if (is.finite(control$feature_cap) && ncol(x) > control$feature_cap) {
      index <- .arise_order_desc(
        variance,
        as.integer(control$feature_cap),
        match(names(variance), input_features)
      )
      selected_features <- names(variance)[index]
      x <- x[, selected_features, drop = FALSE]
      variance <- variance[selected_features]
    } else {
      selected_features <- colnames(x)
    }
    preprocessor <- structure(list(
      input_features = input_features,
      finite_features = finite_features,
      selected_features = selected_features,
      medians = medians[selected_features],
      variances = variance,
      dropped_all_nonfinite = setdiff(input_features, finite_features),
      dropped_nonvariable = setdiff(finite_features, variable_features),
      variance_tolerance = control$variance_tolerance,
      feature_cap = control$feature_cap
    ), class = "arise_preprocessor")
    list(x = x, preprocessor = preprocessor)
  }

  .arise_apply_preprocessor <- function(preprocessor, newdata) {
    if (!inherits(preprocessor, "arise_preprocessor")) {
      .arise_stop("Invalid ARISE preprocessor.")
    }
    if ((is.matrix(newdata) || is.data.frame(newdata)) &&
        is.null(colnames(newdata)) &&
        ncol(newdata) == length(preprocessor$input_features)) {
      colnames(newdata) <- preprocessor$input_features
    }
    x <- .arise_numeric_matrix(newdata, "newdata")
    needed <- preprocessor$selected_features
    missing <- setdiff(needed, colnames(x))
    if (length(missing)) {
      .arise_stop(
        "'newdata' is missing required feature", if (length(missing) > 1L) "s" else "",
        ": ", paste(missing, collapse = ", "), "."
      )
    }
    x <- x[, needed, drop = FALSE]
    medians <- preprocessor$medians[needed]
    for (j in seq_len(ncol(x))) {
      bad <- !is.finite(x[, j])
      if (any(bad)) {
        x[bad, j] <- medians[j]
      }
    }
    storage.mode(x) <- "double"
    x
  }

  .arise_prepare_split <- function(x_train, x_valid, control) {
    fitted <- .arise_fit_preprocessor(x_train, control)
    list(
      x_train = fitted$x,
      x_valid = .arise_apply_preprocessor(fitted$preprocessor, x_valid),
      preprocessor = fitted$preprocessor
    )
  }

  .arise_standardize_split <- function(x_train, x_test) {
    means <- colMeans(x_train)
    scales <- apply(x_train, 2L, stats::sd)
    scales[!is.finite(scales) | scales < 1e-12] <- 1
    list(
      x_train = sweep(sweep(x_train, 2L, means, "-"), 2L, scales, "/"),
      x_test = sweep(sweep(x_test, 2L, means, "-"), 2L, scales, "/"),
      center = means,
      scale = scales
    )
  }


  .arise_col_quantiles <- function(x, probs = c(.25, .75)) {
    x <- as.matrix(x)
    q <- matrixStats::colQuantiles(x, probs = probs, drop = FALSE)
    q <- t(q)
    rownames(q) <- paste0("Q", probs)
    colnames(q) <- colnames(x)
    q
  }

  .arise_class_summaries <- function(x, y) {
    classes <- levels(y)
    p <- ncol(x)
    k <- length(classes)
    template <- matrix(
      0,
      nrow = k,
      ncol = p,
      dimnames = list(classes, colnames(x))
    )
    means <- variances <- minima <- maxima <- q25 <- q75 <- template
    counts <- integer(k)
    names(counts) <- classes
    for (i in seq_len(k)) {
      xi <- x[y == classes[i], , drop = FALSE]
      counts[i] <- nrow(xi)
      means[i, ] <- colMeans(xi)
      variances[i, ] <- .arise_fast_variance(xi)
      minima[i, ] <- matrixStats::colMins(xi)
      maxima[i, ] <- matrixStats::colMaxs(xi)
      quantiles <- .arise_col_quantiles(xi, c(.25, .75))
      q25[i, ] <- quantiles[1L, ]
      q75[i, ] <- quantiles[2L, ]
    }
    list(
      classes = classes,
      counts = counts,
      means = means,
      variances = variances,
      minima = minima,
      maxima = maxima,
      q25 = q25,
      q75 = q75
    )
  }

  .arise_pair_index <- function(classes) {
    pair <- utils::combn(seq_along(classes), 2L)
    data.frame(
      first = pair[1L, ],
      second = pair[2L, ],
      Pair = paste(classes[pair[1L, ]], classes[pair[2L, ]], sep = "__vs__"),
      stringsAsFactors = FALSE
    )
  }

  .arise_rank_auc_pair <- function(x_pair, positive) {
    n_positive <- sum(positive)
    n_negative <- length(positive) - n_positive
    if (n_positive == 0L || n_negative == 0L) {
      return(rep(0, ncol(x_pair)))
    }
    ranks <- matrixStats::colRanks(
      x_pair,
      ties.method = "average",
      preserveShape = TRUE
    )
    u <- colSums(ranks[positive, , drop = FALSE]) -
      n_positive * (n_positive + 1) / 2
    abs(2 * u / (n_positive * n_negative) - 1)
  }

  .arise_pairwise_components <- function(x, y, summaries, pos_whisker = 1.5) {
    pairs <- .arise_pair_index(summaries$classes)
    p <- ncol(x)
    pair_names <- pairs$Pair
    component_names <- c(
      "POS", "Fisher", "SNR", "Welch", "Bhattacharyya", "RankAUC"
    )
    raw <- stats::setNames(lapply(component_names, function(unused) {
      matrix(
        0,
        nrow = p,
        ncol = nrow(pairs),
        dimnames = list(colnames(x), pair_names)
      )
    }), component_names)

    for (h in seq_len(nrow(pairs))) {
      a <- pairs$first[h]
      b <- pairs$second[h]
      variance_a <- pmax(summaries$variances[a, ], 1e-12)
      variance_b <- pmax(summaries$variances[b, ], 1e-12)
      delta <- summaries$means[b, ] - summaries$means[a, ]

      raw$Fisher[, h] <- delta^2 / (variance_a + variance_b + 1e-12)
      raw$SNR[, h] <- abs(delta) /
        (sqrt(variance_a) + sqrt(variance_b) + 1e-12)
      raw$Welch[, h] <- abs(delta) / sqrt(
        variance_a / summaries$counts[a] +
          variance_b / summaries$counts[b] + 1e-12
      )
      raw$Bhattacharyya[, h] <-
        .25 * log(.25 * (variance_a / variance_b + variance_b / variance_a + 2)) +
        .25 * delta^2 / (variance_a + variance_b)

      iqr_a <- summaries$q75[a, ] - summaries$q25[a, ]
      iqr_b <- summaries$q75[b, ] - summaries$q25[b, ]
      lower_a <- pmax(
        summaries$minima[a, ],
        summaries$q25[a, ] - pos_whisker * iqr_a
      )
      upper_a <- pmin(
        summaries$maxima[a, ],
        summaries$q75[a, ] + pos_whisker * iqr_a
      )
      lower_b <- pmax(
        summaries$minima[b, ],
        summaries$q25[b, ] - pos_whisker * iqr_b
      )
      upper_b <- pmin(
        summaries$maxima[b, ],
        summaries$q75[b, ] + pos_whisker * iqr_b
      )
      overlap <- pmax(0, pmin(upper_a, upper_b) - pmax(lower_a, lower_b))
      hull <- pmax(upper_a, upper_b) - pmin(lower_a, lower_b)
      raw$POS[, h] <- ifelse(hull > 1e-12, 1 - overlap / hull, 0)

      keep <- y %in% summaries$classes[c(a, b)]
      raw$RankAUC[, h] <- .arise_rank_auc_pair(
        x[keep, , drop = FALSE],
        y[keep] == summaries$classes[b]
      )
    }

    raw <- lapply(raw, function(value) {
      value[!is.finite(value)] <- 0
      value
    })
    percentiles <- lapply(raw, .arise_percentile_matrix)
    aggregate_percentiles <- vapply(percentiles, rowMeans, numeric(p))
    aggregate_percentiles <- matrix(
      aggregate_percentiles,
      nrow = p,
      ncol = length(component_names),
      dimnames = list(colnames(x), component_names)
    )
    list(
      pairs = pairs,
      raw = raw,
      percentiles = percentiles,
      aggregate_percentiles = aggregate_percentiles
    )
  }

  .arise_mutual_information <- function(x, y, bins = 3L) {
    bins <- max(2L, as.integer(bins)[1L])
    y <- factor(y)
    score_one <- function(value) {
      support <- sort(unique(value))
      if (length(support) <= bins) {
        discretized <- match(value, support)
      } else {
        breaks <- unique(stats::quantile(
          value,
          probs = seq(0, 1, length.out = bins + 1L),
          names = FALSE,
          type = 7
        ))
        if (length(breaks) < 2L) {
          return(0)
        }
        discretized <- cut(
          value,
          breaks = breaks,
          include.lowest = TRUE,
          labels = FALSE
        )
      }
      table_xy <- table(discretized, y)
      probability_xy <- table_xy / sum(table_xy)
      expected <- outer(rowSums(probability_xy), colSums(probability_xy))
      valid <- probability_xy > 0 & expected > 0
      sum(probability_xy[valid] * log(probability_xy[valid] / expected[valid]))
    }
    out <- vapply(seq_len(ncol(x)), function(j) score_one(x[, j]), numeric(1))
    out[!is.finite(out)] <- 0
    stats::setNames(out, colnames(x))
  }

  .arise_component_cache <- function(x, y, control) {
    checked <- .arise_validate_processed(x, y)
    x <- checked$x
    y <- checked$y
    summaries <- .arise_class_summaries(x, y)
    pairwise <- .arise_pairwise_components(
      x,
      y,
      summaries,
      pos_whisker = control$pos_whisker
    )
    mutual_information <- .arise_mutual_information(x, y, control$mi_bins)
    component_percentiles <- cbind(
      pairwise$aggregate_percentiles,
      MI = .arise_percentile_rank(mutual_information)
    )
    component_percentiles <- component_percentiles[
      , names(.arise_weight_columns()), drop = FALSE
    ]
    rownames(component_percentiles) <- colnames(x)
    list(
      classes = levels(y),
      summaries = summaries,
      pairwise = pairwise,
      mutual_information = mutual_information,
      component_percentiles = component_percentiles
    )
  }


  .arise_profile_relevance <- function(component_percentiles, profiles) {
    component_percentiles <- as.matrix(component_percentiles)
    weights <- .arise_weight_matrix(profiles)
    missing <- setdiff(colnames(weights), colnames(component_percentiles))
    if (length(missing)) {
      .arise_stop("Component cache is missing: ", paste(missing, collapse = ", "), ".")
    }
    out <- component_percentiles[, colnames(weights), drop = FALSE] %*% t(weights)
    colnames(out) <- profiles$Config_ID
    rownames(out) <- rownames(component_percentiles)
    out
  }

  .arise_screen_indices <- function(cache, profiles, target, control) {
    relevance <- .arise_profile_relevance(cache$component_percentiles, profiles)
    required_per_profile <- min(nrow(relevance), max(1L, as.integer(target)))
    required <- unique(unlist(lapply(seq_len(ncol(relevance)), function(j) {
      .arise_order_desc(relevance[, j], required_per_profile)
    }), use.names = FALSE))
    mean_relevance <- rowMeans(relevance)
    screen_cap <- min(
      nrow(relevance),
      max(
        length(required),
        as.integer(target),
        as.integer(control$screen_max),
        as.integer(control$screen_multiplier) * as.integer(target)
      )
    )
    fill <- .arise_order_desc(mean_relevance, nrow(relevance))
    fill <- fill[!fill %in% required]
    c(required, utils::head(fill, max(0L, screen_cap - length(required))))
  }

  .arise_balanced_subsample_indices <- function(y, fraction, seed) {
    y <- factor(y)
    counts <- table(y)
    take <- min(
      min(counts),
      max(2L, as.integer(floor(min(counts) * as.numeric(fraction))))
    )
    .arise_with_seed(seed, function() {
      unlist(lapply(levels(y), function(class) {
        sample(which(y == class), take, replace = FALSE)
      }), use.names = FALSE)
    })
  }

  .arise_stability_matrix <- function(x, y, profiles, control, seed) {
    weights <- .arise_weight_matrix(profiles)
    accumulator <- matrix(
      0,
      nrow = ncol(x),
      ncol = nrow(profiles),
      dimnames = list(colnames(x), profiles$Config_ID)
    )
    for (b in seq_len(control$stability_subsamples)) {
      index <- .arise_balanced_subsample_indices(
        y,
        control$stability_fraction,
        .arise_seed(seed, 3001, b)
      )
      local_cache <- .arise_component_cache(x[index, , drop = FALSE], y[index], control)
      combined <- local_cache$component_percentiles[, colnames(weights), drop = FALSE] %*%
        t(weights)
      combined <- .arise_percentile_matrix(combined)
      colnames(combined) <- profiles$Config_ID
      accumulator <- accumulator + combined
    }
    accumulator / control$stability_subsamples
  }

  .arise_within_class_residual <- function(x, y) {
    residual <- x
    for (class in levels(y)) {
      index <- which(y == class)
      residual[index, ] <- sweep(
        x[index, , drop = FALSE],
        2L,
        colMeans(x[index, , drop = FALSE]),
        "-"
      )
    }
    degrees_freedom <- max(1L, nrow(x) - nlevels(y))
    scale <- sqrt(colSums(residual * residual) / degrees_freedom)
    scale[!is.finite(scale) | scale < 1e-12] <- 1
    residual <- sweep(residual, 2L, scale, "/")
    residual[!is.finite(residual)] <- 0
    list(x = residual, degrees_freedom = degrees_freedom, scale = scale)
  }

  .arise_build_workspace <- function(
      x,
      y,
      profiles,
      target,
      control,
      seed,
      full_cache = NULL) {
    checked <- .arise_validate_processed(x, y)
    x <- checked$x
    y <- checked$y
    target <- min(ncol(x), max(1L, as.integer(target)))
    if (is.null(full_cache)) {
      full_cache <- .arise_component_cache(x, y, control)
    }
    if (!identical(rownames(full_cache$component_percentiles), colnames(x))) {
      .arise_stop("The supplied component cache does not match the feature universe.")
    }
    screen_index <- .arise_screen_indices(full_cache, profiles, target, control)
    screen_features <- colnames(x)[screen_index]
    x_screen <- x[, screen_index, drop = FALSE]
    pairwise_percentiles <- lapply(
      full_cache$pairwise$percentiles,
      function(value) value[screen_index, , drop = FALSE]
    )
    stability <- .arise_stability_matrix(
      x_screen,
      y,
      profiles,
      control,
      .arise_seed(seed, 5003)
    )
    residual <- .arise_within_class_residual(x_screen, y)
    relevance <- .arise_profile_relevance(
      full_cache$component_percentiles[screen_index, , drop = FALSE],
      profiles
    )
    list(
      x = x,
      y = y,
      profiles = profiles,
      target = target,
      control = control,
      seed = seed,
      full_cache = full_cache,
      screen_features = screen_features,
      component_percentiles = full_cache$component_percentiles[screen_index, , drop = FALSE],
      pairwise_percentiles = pairwise_percentiles,
      pair_index = full_cache$pairwise$pairs,
      relevance = relevance,
      stability = stability,
      residual_z = residual$x,
      residual_df = residual$degrees_freedom
    )
  }

  .arise_config_row <- function(config, profiles = NULL) {
    if (is.character(config) && length(config) == 1L) {
      if (is.null(profiles)) {
        .arise_stop("A profile table is required when a configuration ID is supplied.")
      }
      config <- profiles[
        profiles$Config_ID == config | profiles$Profile == config,
        ,
        drop = FALSE
      ]
    }
    config <- as.data.frame(config, stringsAsFactors = FALSE)
    if (nrow(config) != 1L) {
      .arise_stop("Exactly one ARISE configuration is required.")
    }
    config
  }

  .arise_pair_profile <- function(workspace, config, features) {
    if (nrow(workspace$pair_index) <= 1L) {
      return(NULL)
    }
    pair_components <- c(
      "POS", "Fisher", "SNR", "Welch", "Bhattacharyya", "RankAUC"
    )
    weights <- .arise_weight_matrix(config)[1L, pair_components]
    if (sum(weights) <= 0) {
      return(NULL)
    }
    weights <- weights / sum(weights)
    index <- match(features, workspace$screen_features)
    out <- matrix(
      0,
      nrow = length(features),
      ncol = nrow(workspace$pair_index),
      dimnames = list(features, workspace$pair_index$Pair)
    )
    for (component in pair_components) {
      out <- out + weights[component] *
        workspace$pairwise_percentiles[[component]][index, , drop = FALSE]
    }
    out
  }

  .arise_rank_path <- function(workspace, config, target = workspace$target) {
    config <- .arise_config_row(config, workspace$profiles)
    config_id <- config$Config_ID[1L]
    target <- min(
      max(1L, as.integer(target)),
      length(workspace$screen_features),
      ncol(workspace$x)
    )
    relevance <- workspace$relevance[, config_id]
    stability <- workspace$stability[, config_id]
    stability_weight <- min(1, max(0, as.numeric(config$Stability_Weight[1L])))
    base_score <- (1 - stability_weight) * relevance +
      stability_weight * stability

    local_screen_n <- min(
      length(base_score),
      max(target, as.integer(config$Screen_Max[1L]))
    )
    universe_index <- match(workspace$screen_features, colnames(workspace$x))
    local_index <- .arise_order_desc(
      base_score,
      local_screen_n,
      universe_index
    )
    features <- workspace$screen_features[local_index]
    base_score <- base_score[local_index]
    relevance <- relevance[local_index]
    stability <- stability[local_index]
    residual_z <- workspace$residual_z[, local_index, drop = FALSE]
    pair_profile <- .arise_pair_profile(workspace, config, features)

    selected <- integer(target)
    available <- rep(TRUE, length(features))
    max_abs_correlation <- numeric(length(features))
    selected_criterion <- numeric(target)
    selected_redundancy <- numeric(target)
    selected_coverage_gain <- numeric(target)
    selected_lambda <- numeric(target)
    pair_coverage <- if (is.null(pair_profile)) {
      numeric()
    } else {
      rep(0, ncol(pair_profile))
    }
    lambda0 <- max(0, as.numeric(config$Redundancy_Lambda[1L]))
    decay <- min(1, max(0, as.numeric(config$Redundancy_Decay[1L])))
    coverage_weight <- if (is.null(pair_profile)) {
      0
    } else {
      max(0, as.numeric(config$Pair_Coverage_Weight[1L]))
    }

    for (step in seq_len(target)) {
      fraction <- if (target <= 1L) 0 else (step - 1) / (target - 1)
      lambda_step <- lambda0 * (1 - decay * fraction)
      coverage_gain <- if (coverage_weight > 0) {
        rowMeans(pmax(sweep(pair_profile, 2L, pair_coverage, "-"), 0))
      } else {
        rep(0, length(features))
      }
      criterion <- base_score + coverage_weight * coverage_gain -
        lambda_step * max_abs_correlation
      criterion[!available] <- -Inf
      chosen <- which.max(criterion)
      selected[step] <- chosen
      selected_criterion[step] <- criterion[chosen]
      selected_redundancy[step] <- max_abs_correlation[chosen]
      selected_coverage_gain[step] <- coverage_gain[chosen]
      selected_lambda[step] <- lambda_step
      available[chosen] <- FALSE

      correlation <- abs(drop(crossprod(
        residual_z,
        residual_z[, chosen, drop = FALSE]
      ))) / max(1L, workspace$residual_df)
      correlation[!is.finite(correlation)] <- 0
      max_abs_correlation <- pmax(max_abs_correlation, pmin(1, correlation))
      if (coverage_weight > 0) {
        pair_coverage <- pmax(pair_coverage, pair_profile[chosen, ])
      }
    }

    chosen_features <- features[selected]
    names(selected_criterion) <- chosen_features
    list(
      features = chosen_features,
      scores = selected_criterion,
      config = config,
      details = data.frame(
        Position = seq_len(target),
        Feature = chosen_features,
        Relevance = as.numeric(relevance[selected]),
        Stability = as.numeric(stability[selected]),
        Pair_Coverage_Gain = selected_coverage_gain,
        Redundancy_Lambda = selected_lambda,
        Max_Residual_Correlation_When_Selected = selected_redundancy,
        Greedy_Criterion = selected_criterion,
        stringsAsFactors = FALSE
      ),
      screen_n = local_screen_n
    )
  }

  .arise_rank_all_profiles <- function(workspace, target = workspace$target) {
    stats::setNames(lapply(workspace$profiles$Config_ID, function(id) {
      .arise_rank_path(workspace, id, target)
    }), workspace$profiles$Config_ID)
  }

  .arise_soft_consensus_config <- function(profiles, profile_weights) {
    if (is.null(names(profile_weights))) {
      if (length(profile_weights) != nrow(profiles)) {
        .arise_stop("Unnamed profile weights must match the number of profiles.")
      }
      names(profile_weights) <- profiles$Config_ID
    }
    weights <- as.numeric(profile_weights[profiles$Config_ID])
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) {
      .arise_stop("At least one consensus weight must be positive.")
    }
    weights <- weights / sum(weights)
    numeric_columns <- c(
      unname(.arise_weight_columns()),
      "Stability_Weight",
      "Redundancy_Lambda",
      "Redundancy_Decay",
      "Pair_Coverage_Weight"
    )
    out <- data.frame(
      Config_ID = "ARISE_SOFT_CONSENSUS",
      Profile = "one_se_soft_consensus",
      Priority = 1L,
      stringsAsFactors = FALSE
    )
    for (name in numeric_columns) {
      out[[name]] <- sum(weights * as.numeric(profiles[[name]]))
    }
    out$Screen_Max <- max(
      100L,
      as.integer(round(sum(weights * profiles$Screen_Max)))
    )
    component_columns <- unname(.arise_weight_columns())
    total <- sum(as.numeric(out[1L, component_columns, drop = TRUE]))
    out[, component_columns] <- out[, component_columns, drop = FALSE] / total
    attr(out, "profile_weights") <- stats::setNames(weights, profiles$Config_ID)
    out
  }


  .arise_prediction <- function(class, scores, classifier) {
    scores <- as.matrix(scores)
    storage.mode(scores) <- "double"
    if (nrow(scores) != length(class) || ncol(scores) < 2L ||
        is.null(colnames(scores))) {
      .arise_stop("A fixed learner returned an invalid score matrix.")
    }
    if (any(!is.finite(scores))) {
      .arise_stop("A fixed learner returned non-finite class scores.")
    }
    class <- factor(as.character(class), levels = colnames(scores))
    if (anyNA(class)) {
      .arise_stop("A fixed learner returned a class outside its score columns.")
    }
    list(class = class, scores = scores, classifier = classifier)
  }

  .arise_predict_knn <- function(x_train, y_train, x_test, k, class_weights) {
    classes <- levels(y_train)
    number_classes <- length(classes)
    k <- min(nrow(x_train), max(1L, as.integer(k)))
    train_norm <- rowSums(x_train * x_train)
    test_norm <- rowSums(x_test * x_test)
    distance_squared <- outer(test_norm, train_norm, "+") -
      2 * tcrossprod(x_test, x_train)
    distance_squared[distance_squared < 0 & distance_squared > -1e-9] <- 0
    distance_squared <- pmax(distance_squared, 0)
    y_code <- as.integer(y_train)
    scores <- matrix(
      0,
      nrow = nrow(x_test),
      ncol = number_classes,
      dimnames = list(rownames(x_test), classes)
    )
    if (is.null(class_weights)) {
      class_weights <- rep(1, number_classes)
    } else {
      class_weights <- as.numeric(class_weights[classes])
    }
    for (i in seq_len(nrow(x_test))) {
      neighbours <- order(distance_squared[i, ], method = "radix")[seq_len(k)]
      distance <- sqrt(distance_squared[i, neighbours])
      zero_distance <- distance <= 1e-12
      if (any(zero_distance)) {
        neighbours <- neighbours[zero_distance]
        vote_weight <- rep(1, length(neighbours))
      } else {
        vote_weight <- 1 / pmax(distance, 1e-12)
      }
      vote_weight <- vote_weight * class_weights[y_code[neighbours]]
      scores[i, ] <- vapply(seq_len(number_classes), function(class_index) {
        sum(vote_weight[y_code[neighbours] == class_index])
      }, numeric(1))
    }
    scores <- scores / pmax(rowSums(scores), .Machine$double.eps)
    predicted <- classes[max.col(scores, ties.method = "first")]
    .arise_prediction(predicted, scores, "kNN")
  }

  .arise_svm_margin_scores <- function(fit, x_train, y_train, x_test,
                                       test_prediction = NULL) {
    classes <- levels(y_train)
    number_classes <- length(classes)
    train_prediction <- stats::predict(fit, x_train, decision.values = TRUE)
    if (is.null(test_prediction)) {
      test_prediction <- stats::predict(fit, x_test, decision.values = TRUE)
    }
    train_margin_raw <- attr(train_prediction, "decision.values")
    test_margin_raw <- attr(test_prediction, "decision.values")
    train_margin <- matrix(train_margin_raw, nrow = nrow(x_train))
    test_margin <- matrix(test_margin_raw, nrow = nrow(x_test))
    # e1071 emits pairwise margins in its native label order. Class names may
    # themselves contain '/', so the printed margin names are not unambiguous.
    class_order <- match(fit$levels[fit$labels], classes)
    if (length(class_order) != number_classes || anyNA(class_order) ||
        anyDuplicated(class_order)) {
      .arise_stop("The SVM did not return the expected class ordering.")
    }
    pairs <- utils::combn(class_order, 2L)
    if (ncol(train_margin) != ncol(pairs) ||
        ncol(test_margin) != ncol(pairs)) {
      .arise_stop("The SVM did not return the expected pairwise margins.")
    }
    scores <- matrix(
      0,
      nrow = nrow(x_test),
      ncol = number_classes,
      dimnames = list(rownames(x_test), classes)
    )
    for (j in seq_len(ncol(test_margin))) {
      a <- pairs[1L, j]
      b <- pairs[2L, j]
      train_j <- train_margin[, j]
      test_j <- test_margin[, j]
      mean_a <- mean(train_j[y_train == classes[a]])
      mean_b <- mean(train_j[y_train == classes[b]])
      if (is.finite(mean_a) && is.finite(mean_b) && mean_a < mean_b) {
        train_j <- -train_j
        test_j <- -test_j
      }
      scale_j <- stats::sd(train_j[y_train %in% classes[c(a, b)]])
      if (!is.finite(scale_j) || scale_j < 1e-8) {
        scale_j <- 1
      }
      probability_a <- stats::plogis(test_j / scale_j)
      scores[, a] <- scores[, a] + probability_a
      scores[, b] <- scores[, b] + 1 - probability_a
    }
    scores / ncol(pairs)
  }

  .arise_predict_svm <- function(x_train, y_train, x_test, control, seed,
                                 class_weights) {
    gamma <- 1 / max(1L, ncol(x_train))
    fit <- .arise_with_seed(seed, function() {
      e1071::svm(
        x = x_train,
        y = y_train,
        type = "C-classification",
        kernel = "radial",
        cost = control$svm_cost,
        gamma = gamma,
        scale = FALSE,
        probability = FALSE,
        shrinking = TRUE,
        class.weights = class_weights,
        fitted = FALSE,
        cachesize = 200
      )
    })
    test_prediction <- stats::predict(
      fit,
      x_test,
      probability = FALSE,
      decision.values = TRUE
    )
    scores <- .arise_svm_margin_scores(
      fit,
      x_train,
      y_train,
      x_test,
      test_prediction = test_prediction
    )
    classes <- levels(y_train)
    predicted <- classes[max.col(scores, ties.method = "first")]
    .arise_prediction(predicted, scores, "SVM")
  }

  .arise_predict_rf <- function(x_train, y_train, x_test, control, seed,
                                class_weights) {
    q <- ncol(x_train)
    fit <- .arise_with_seed(seed, function() {
      randomForest::randomForest(
        x = x_train,
        y = y_train,
        xtest = x_test,
        ntree = control$rf_trees,
        mtry = max(1L, floor(sqrt(q))),
        nodesize = 1L,
        replace = TRUE,
        classwt = class_weights,
        importance = FALSE,
        proximity = FALSE,
        keep.forest = FALSE,
        norm.votes = TRUE
      )
    })
    scores <- fit$test$votes
    if (is.null(scores)) {
      .arise_stop("The random-forest learner returned no validation votes.")
    }
    scores <- as.matrix(scores)
    classes <- levels(y_train)
    if (!all(classes %in% colnames(scores))) {
      .arise_stop("The random-forest score matrix is missing an outcome class.")
    }
    scores <- scores[, classes, drop = FALSE]
    scores[!is.finite(scores)] <- 0
    scores <- scores / pmax(rowSums(scores), .Machine$double.eps)
    predicted <- classes[max.col(scores, ties.method = "first")]
    .arise_prediction(predicted, scores, "Random Forest")
  }

  .arise_fit_predict <- function(classifier, x_train, y_train, x_test,
                                 control, seed) {
    classifier <- .arise_normalize_classifier(classifier)
    standardized <- .arise_standardize_split(x_train, x_test)
    x_train <- standardized$x_train
    x_test <- standardized$x_test
    class_weights <- if (isTRUE(control$class_weights)) {
      .arise_inverse_frequency_weights(y_train)
    } else {
      NULL
    }
    if (classifier == "knn") {
      return(.arise_predict_knn(
        x_train,
        y_train,
        x_test,
        control$knn_k,
        class_weights
      ))
    }
    if (classifier == "svm") {
      return(.arise_predict_svm(
        x_train,
        y_train,
        x_test,
        control,
        seed,
        class_weights
      ))
    }
    .arise_predict_rf(
      x_train,
      y_train,
      x_test,
      control,
      seed,
      class_weights
    )
  }

  .arise_metrics <- function(y_true, predicted) {
    y_true <- factor(y_true)
    predicted <- factor(predicted, levels = levels(y_true))
    classes <- levels(y_true)
    recall <- f1 <- numeric(length(classes))
    for (k in seq_along(classes)) {
      class <- classes[k]
      true_positive <- sum(y_true == class & predicted == class)
      false_negative <- sum(y_true == class & predicted != class)
      false_positive <- sum(y_true != class & predicted == class)
      recall[k] <- true_positive / max(1, true_positive + false_negative)
      f1[k] <- 2 * true_positive /
        max(1, 2 * true_positive + false_positive + false_negative)
    }
    confusion <- table(y_true, predicted)
    total <- sum(confusion)
    accuracy <- sum(diag(confusion)) / max(1, total)
    expected <- if (total > 0) {
      sum(rowSums(confusion) * colSums(confusion)) / (total * total)
    } else {
      NA_real_
    }
    kappa <- if (is.finite(expected) &&
        abs(1 - expected) > sqrt(.Machine$double.eps)) {
      (accuracy - expected) / (1 - expected)
    } else {
      -1
    }
    c(
      Balanced_Accuracy = mean(recall),
      Macro_F1 = mean(f1),
      Kappa = kappa
    )
  }


  .arise_stratified_folds <- function(y, v, seed) {
    y <- droplevels(factor(y))
    v <- min(as.integer(v), min(table(y)))
    if (v < 2L) {
      .arise_stop("At least two inner folds are required.")
    }
    validation <- vector("list", v)
    validation <- lapply(validation, function(unused) integer())
    for (class_index in seq_along(levels(y))) {
      class <- levels(y)[class_index]
      index <- which(y == class)
      shuffled <- .arise_with_seed(
        .arise_seed(seed, class_index),
        function() sample(index, length(index), replace = FALSE)
      )
      assignment <- rep(seq_len(v), length.out = length(shuffled))
      for (fold in seq_len(v)) {
        validation[[fold]] <- c(validation[[fold]], shuffled[assignment == fold])
      }
    }
    folds <- lapply(seq_len(v), function(fold) {
      valid <- sort(unique(validation[[fold]]))
      train <- setdiff(seq_along(y), valid)
      if (any(table(droplevels(y[train])) < 2L)) {
        .arise_stop(
          "An inner training fold contains fewer than two observations in a class."
        )
      }
      list(
        Fold_ID = paste0("F", fold),
        train = train,
        valid = valid
      )
    })
    names(folds) <- vapply(folds, `[[`, character(1), "Fold_ID")
    folds
  }

  .arise_evaluate_path <- function(path, fold_data, config_id, stage,
                                   budgets, classifiers, control, seed) {
    rows <- vector("list", length(budgets) * length(classifiers))
    cursor <- 0L
    for (b in seq_along(budgets)) {
      effective_budget <- min(as.integer(budgets[b]), length(path$features))
      features <- path$features[seq_len(effective_budget)]
      for (h in seq_along(classifiers)) {
        cursor <- cursor + 1L
        prediction <- .arise_fit_predict(
          classifiers[h],
          fold_data$x_train[, features, drop = FALSE],
          fold_data$y_train,
          fold_data$x_valid[, features, drop = FALSE],
          control,
          .arise_seed(seed, b, h, effective_budget)
        )
        metric <- .arise_metrics(fold_data$y_valid, prediction$class)
        rows[[cursor]] <- data.frame(
          Fold_ID = fold_data$Fold_ID,
          Stage = stage,
          Config_ID = config_id,
          Feature_Budget = as.integer(budgets[b]),
          Effective_Feature_Budget = effective_budget,
          Classifier = classifiers[h],
          Balanced_Accuracy = unname(metric["Balanced_Accuracy"]),
          Macro_F1 = unname(metric["Macro_F1"]),
          Kappa = unname(metric["Kappa"]),
          N_Train = nrow(fold_data$x_train),
          N_Validation = nrow(fold_data$x_valid),
          N_Features_After_Preprocessing = ncol(fold_data$x_train),
          stringsAsFactors = FALSE
        )
      }
    }
    do.call(rbind, rows[seq_len(cursor)])
  }

  .arise_fold_metric_means <- function(evaluations, candidates) {
    metric_names <- c("Balanced_Accuracy", "Macro_F1", "Kappa")
    evaluations <- evaluations[
      evaluations$Config_ID %in% candidates,
      ,
      drop = FALSE
    ]
    if (nrow(evaluations) == 0L) {
      .arise_stop("No tuning evaluations are available.")
    }
    stats::aggregate(
      evaluations[, metric_names, drop = FALSE],
      evaluations[, c("Config_ID", "Fold_ID"), drop = FALSE],
      function(value) mean(value[is.finite(value)], na.rm = TRUE)
    )
  }

  .arise_regret_summary <- function(evaluations, candidates, profiles, control) {
    metric_names <- c("Balanced_Accuracy", "Macro_F1", "Kappa")
    fold_means <- .arise_fold_metric_means(evaluations, candidates)
    for (metric in metric_names) {
      value <- fold_means[[metric]]
      value[!is.finite(value)] <- if (metric == "Kappa") -1 else 0
      if (metric == "Kappa") {
        value <- (value + 1) / 2
      }
      fold_means[[metric]] <- pmin(1, pmax(0, value))
    }
    fold_ids <- unique(fold_means$Fold_ID)
    regret_rows <- vector("list", length(fold_ids))
    for (i in seq_along(fold_ids)) {
      local <- fold_means[fold_means$Fold_ID == fold_ids[i], , drop = FALSE]
      local <- local[match(candidates, local$Config_ID), , drop = FALSE]
      if (anyNA(local$Config_ID)) {
        .arise_stop("A profile is missing from one or more inner folds.")
      }
      regrets <- matrix(
        0,
        nrow = nrow(local),
        ncol = length(metric_names),
        dimnames = list(NULL, metric_names)
      )
      for (metric in metric_names) {
        regrets[, metric] <- max(local[[metric]]) - local[[metric]]
      }
      fold_loss <- control$max_regret_weight * apply(regrets, 1L, max) +
        control$mean_regret_weight * rowMeans(regrets)
      regret_rows[[i]] <- data.frame(
        Config_ID = local$Config_ID,
        Fold_ID = local$Fold_ID,
        Regret_Balanced_Accuracy = regrets[, "Balanced_Accuracy"],
        Regret_Macro_F1 = regrets[, "Macro_F1"],
        Regret_Kappa = regrets[, "Kappa"],
        Fold_Loss = fold_loss,
        stringsAsFactors = FALSE
      )
    }
    fold_regrets <- do.call(rbind, regret_rows)
    summary_rows <- lapply(candidates, function(id) {
      local <- fold_regrets[fold_regrets$Config_ID == id, , drop = FALSE]
      mean_regret <- c(
        Balanced_Accuracy = mean(local$Regret_Balanced_Accuracy),
        Macro_F1 = mean(local$Regret_Macro_F1),
        Kappa = mean(local$Regret_Kappa)
      )
      mean_loss <- mean(local$Fold_Loss)
      se_loss <- if (nrow(local) > 1L) {
        stats::sd(local$Fold_Loss) / sqrt(nrow(local))
      } else {
        0
      }
      data.frame(
        Config_ID = id,
        Mean_Loss = mean_loss,
        SE_Loss = se_loss,
        Worst_Mean_Regret = max(mean_regret),
        Mean_Regret = mean(mean_regret),
        Regret_Balanced_Accuracy = mean_regret["Balanced_Accuracy"],
        Regret_Macro_F1 = mean_regret["Macro_F1"],
        Regret_Kappa = mean_regret["Kappa"],
        N_Folds = nrow(local),
        stringsAsFactors = FALSE
      )
    })
    summary <- do.call(rbind, summary_rows)
    priority <- profiles$Priority[match(summary$Config_ID, profiles$Config_ID)]
    summary <- summary[
      order(summary$Mean_Loss, summary$SE_Loss, priority, method = "radix"),
      ,
      drop = FALSE
    ]
    rownames(summary) <- NULL
    list(
      summary = summary,
      fold_regrets = fold_regrets,
      fold_means = fold_means
    )
  }

  .arise_stage1_survivors <- function(summary, profiles, control) {
    priority <- profiles$Priority[match(summary$Config_ID, profiles$Config_ID)]
    ordered <- summary[
      order(summary$Mean_Loss, summary$SE_Loss, priority, method = "radix"),
      ,
      drop = FALSE
    ]
    threshold <- ordered$Mean_Loss[1L] +
      control$one_se_multiplier * ordered$SE_Loss[1L]
    preferred <- ordered$Config_ID[
      ordered$Mean_Loss <= threshold + 1e-15
    ]
    padded <- c(preferred, setdiff(ordered$Config_ID, preferred))
    survivors <- utils::head(unique(padded), control$stage2_profiles)
    list(ids = survivors, threshold = threshold, preferred = preferred)
  }

  .arise_final_profile_weights <- function(summary, profiles, control) {
    priority <- profiles$Priority[match(summary$Config_ID, profiles$Config_ID)]
    ordered <- summary[
      order(summary$Mean_Loss, summary$SE_Loss, priority, method = "radix"),
      ,
      drop = FALSE
    ]
    threshold <- ordered$Mean_Loss[1L] +
      control$one_se_multiplier * ordered$SE_Loss[1L]
    eligible <- ordered[
      ordered$Mean_Loss <= threshold + 1e-15,
      ,
      drop = FALSE
    ]
    if (nrow(eligible) == 0L) {
      eligible <- ordered[1L, , drop = FALSE]
    }
    temperature <- control$softmax_temperature
    raw <- exp(-(eligible$Mean_Loss - min(eligible$Mean_Loss)) / temperature) /
      (1 + eligible$SE_Loss / temperature)
    weights <- raw / sum(raw)
    list(
      weights = stats::setNames(weights, eligible$Config_ID),
      eligible = eligible,
      threshold = threshold
    )
  }

  .arise_tune <- function(x_raw, y, profiles, target, control, seed) {
    folds <- .arise_stratified_folds(y, control$inner_folds, .arise_seed(seed, 101))
    fold_cache <- vector("list", length(folds))
    names(fold_cache) <- names(folds)
    maximum_budget <- max(c(
      target,
      control$budgets,
      control$stage1_budgets
    ))

    for (i in seq_along(folds)) {
      fold <- folds[[i]]
      prepared <- .arise_prepare_split(
        x_raw[fold$train, , drop = FALSE],
        x_raw[fold$valid, , drop = FALSE],
        control
      )
      y_train <- droplevels(factor(y[fold$train], levels = levels(y)))
      y_valid <- factor(y[fold$valid], levels = levels(y))
      effective_target <- min(maximum_budget, ncol(prepared$x_train))
      .arise_log(
        control,
        "ARISE inner fold ", i, "/", length(folds),
        ": constructing all profile paths."
      )
      workspace <- .arise_build_workspace(
        prepared$x_train,
        y_train,
        profiles,
        effective_target,
        control,
        .arise_seed(seed, 201, i)
      )
      fold_cache[[i]] <- list(
        Fold_ID = fold$Fold_ID,
        x_train = prepared$x_train,
        x_valid = prepared$x_valid,
        y_train = y_train,
        y_valid = y_valid,
        paths = .arise_rank_all_profiles(workspace, effective_target),
        train_indices = fold$train,
        validation_indices = fold$valid
      )
    }

    stage1_rows <- list()
    cursor <- 0L
    for (i in seq_along(fold_cache)) {
      fold_data <- fold_cache[[i]]
      for (id in profiles$Config_ID) {
        cursor <- cursor + 1L
        stage1_rows[[cursor]] <- .arise_evaluate_path(
          fold_data$paths[[id]],
          fold_data,
          id,
          "stage1",
          control$stage1_budgets,
          control$stage1_classifiers,
          control,
          .arise_seed(seed, 301, i)
        )
      }
    }
    stage1_evaluations <- do.call(rbind, stage1_rows)
    stage1_regret <- .arise_regret_summary(
      stage1_evaluations,
      profiles$Config_ID,
      profiles,
      control
    )
    stage1_selection <- .arise_stage1_survivors(
      stage1_regret$summary,
      profiles,
      control
    )

    stage2_rows <- list()
    cursor <- 0L
    for (i in seq_along(fold_cache)) {
      fold_data <- fold_cache[[i]]
      for (id in stage1_selection$ids) {
        cursor <- cursor + 1L
        stage2_rows[[cursor]] <- .arise_evaluate_path(
          fold_data$paths[[id]],
          fold_data,
          id,
          "stage2",
          control$budgets,
          control$stage2_classifiers,
          control,
          .arise_seed(seed, 401, i)
        )
      }
    }
    stage2_evaluations <- do.call(rbind, stage2_rows)
    final_regret <- .arise_regret_summary(
      stage2_evaluations,
      stage1_selection$ids,
      profiles,
      control
    )
    final_selection <- .arise_final_profile_weights(
      final_regret$summary,
      profiles,
      control
    )
    eligible_profiles <- profiles[
      profiles$Config_ID %in% names(final_selection$weights),
      ,
      drop = FALSE
    ]
    eligible_profiles <- eligible_profiles[
      match(names(final_selection$weights), eligible_profiles$Config_ID),
      ,
      drop = FALSE
    ]
    consensus <- .arise_soft_consensus_config(
      eligible_profiles,
      final_selection$weights
    )

    list(
      consensus_config = consensus,
      profile_weights = final_selection$weights,
      eligible_profiles = final_selection$eligible,
      survivors = stage1_selection$ids,
      stage1_preferred = stage1_selection$preferred,
      stage1_threshold = stage1_selection$threshold,
      final_threshold = final_selection$threshold,
      stage1_summary = stage1_regret$summary,
      final_summary = final_regret$summary,
      stage1_evaluations = stage1_evaluations,
      stage2_evaluations = stage2_evaluations,
      fold_regrets = final_regret$fold_regrets,
      folds = lapply(fold_cache, function(value) {
        list(
          Fold_ID = value$Fold_ID,
          train_indices = value$train_indices,
          validation_indices = value$validation_indices
        )
      }),
      objective = list(
        metrics = c("Balanced_Accuracy", "Macro_F1", "Kappa"),
        worst_regret_weight = control$max_regret_weight,
        mean_regret_weight = control$mean_regret_weight,
        one_se_multiplier = control$one_se_multiplier,
        softmax_temperature = control$softmax_temperature
      )
    )
  }


  .arise_fit <- function(x, y, n_features, control, seed = 1L) {
    control <- .arise_revalidate_control(control)
    n_features <- .arise_scalar_integer(
      n_features, "n_features", minimum = 1L
    )
    seed <- .arise_scalar_integer(seed, "seed", minimum = 0L)
    checked <- .arise_validate_input(
      x,
      y,
      tuning = control$tune,
      inner_folds = control$inner_folds
    )
    x_raw <- checked$x
    y <- checked$y
    requested_maximum <- if (isTRUE(control$tune)) {
      max(c(n_features, control$budgets, control$stage1_budgets))
    } else {
      n_features
    }
    if (requested_maximum > control$screen_max) {
      .arise_stop(
        "The maximum requested feature budget exceeds the internal screen size."
      )
    }

    profiles <- .arise_profiles(control$screen_max)

    if (isTRUE(control$tune)) {
      tuning <- .arise_tune(
        x_raw,
        y,
        profiles,
        requested_maximum,
        control,
        .arise_seed(seed, 701)
      )
      final_config <- tuning$consensus_config
    } else {
      final_config <- profiles[
        profiles$Profile == control$untuned_profile,
        ,
        drop = FALSE
      ]
    }

    preprocessing <- .arise_fit_preprocessor(x_raw, control)
    x_processed <- preprocessing$x
    final_target <- min(requested_maximum, ncol(x_processed))
    workspace <- .arise_build_workspace(
      x_processed,
      y,
      final_config,
      final_target,
      control,
      .arise_seed(seed, 801)
    )
    path <- .arise_rank_path(workspace, final_config, final_target)
    path$features[seq_len(min(n_features, length(path$features)))]
  }

  function(X, Y, topK = 10) {
    X <- .arise_numeric_matrix(X, "X")
    topK <- .arise_scalar_integer(topK, "topK", minimum = 1L)
    control <- .arise_control(
      screen_max = max(700L, topK),
      store_components = FALSE,
      verbose = FALSE
    )
    selected <- .arise_fit(
      x = X,
      y = Y,
      n_features = topK,
      control = control,
      seed = 1L
    )
    unname(as.integer(match(selected, colnames(X))))
  }
})
