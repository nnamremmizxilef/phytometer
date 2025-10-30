# Avoid NSE notes
utils::globalVariables(c(".data", "PC", "PC1", "PC2", "Variance",
                         "color_var", "correlation", "env", "pheno"))

#' Merge multi-omics datasets
#'
#' Intelligently merge environmental, phenotypic, and genotypic data with
#' flexible join keys and options.
#'
#' @param enviro Environmental data frame (optional)
#' @param pheno Phenotypic data frame (optional)
#' @param geno Genotypic data frame (optional)
#' @param by.env Character vector of columns to join environmental data (default: "site_id")
#' @param by.pheno Character vector of columns to join phenotypic data (default: c("site_id", "tree_id"))
#' @param by.geno Character vector of columns to join genotypic data (default: c("site_id", "tree_id"))
#' @param join.type Type of join: "inner", "left", "full" (default: "inner")
#' @return A merged data frame
#' @export
#' @examples
#' env <- envirodata(n_sites = 5)
#' pheno <- phenotydata(n_trees = 25)
#' geno <- genotydata(n_trees = 25)
#' merged <- mergeomics(enviro = env, pheno = pheno, geno = geno)
mergeomics <- function(enviro = NULL,
                       pheno = NULL,
                       geno = NULL,
                       by.env = "site_id",
                       by.pheno = c("site_id", "tree_id"),
                       by.geno = c("site_id", "tree_id"),
                       join.type = c("inner", "left", "full")) {

  join.type <- match.arg(join.type)

  # Determine join function
  join_fun <- switch(join.type,
                     "inner" = dplyr::inner_join,
                     "left" = dplyr::left_join,
                     "full" = dplyr::full_join)

  result <- NULL

  # Start with first available dataset
  if (!is.null(pheno)) {
    result <- pheno
  } else if (!is.null(geno)) {
    result <- geno
  } else if (!is.null(enviro)) {
    result <- enviro
  } else {
    stop("At least one dataset must be provided")
  }

  # Join phenotypic data
  if (!is.null(pheno) && is.null(result)) {
    result <- pheno
  } else if (!is.null(pheno) && !identical(result, pheno)) {
    result <- join_fun(result, pheno, by = intersect(by.pheno, names(result)))
  }

  # Join genotypic data
  if (!is.null(geno)) {
    join_cols <- intersect(by.geno, names(result))
    if (length(join_cols) > 0) {
      result <- join_fun(result, geno, by = join_cols)
    }
  }

  # Join environmental data
  if (!is.null(enviro)) {
    join_cols <- intersect(by.env, names(result))
    if (length(join_cols) > 0) {
      result <- join_fun(result, enviro, by = join_cols)
    }
  }

  return(result)
}

#' Calculate composite stress indices
#'
#' Compute stress indices from environmental or physiological variables.
#'
#' @param data Data frame containing stress-related variables
#' @param stress.type Type of stress: "drought", "heat", "cold", "herbivory", or "custom"
#' @param variables Character vector of variable names to use (required for "custom")
#' @param weights Numeric vector of weights for each variable (default: equal weights)
#' @param scale.vars Logical, scale variables before combining (default: TRUE)
#' @param invert Logical vector indicating which variables should be inverted (default: NULL)
#' @return Data frame with added stress_index column
#' @export
#' @importFrom stats sd
#' @examples
#' env <- envirodata()
#' env_stress <- calculatestress(env, stress.type = "drought")
#' head(env_stress)
calculatestress <- function(data,
                           stress.type = c("drought", "heat", "cold", "herbivory", "custom"),
                           variables = NULL,
                           weights = NULL,
                           scale.vars = TRUE,
                           invert = NULL) {

  stress.type <- match.arg(stress.type)

  # Define default variables for each stress type
  if (stress.type == "drought" && is.null(variables)) {
    variables <- c("drought_index", "soil_moisture", "vpd")
    if (is.null(invert)) invert <- c(FALSE, TRUE, FALSE)
  } else if (stress.type == "heat" && is.null(variables)) {
    variables <- c("temp_mean", "temp_max", "vpd")
    if (is.null(invert)) invert <- c(FALSE, FALSE, FALSE)
  } else if (stress.type == "cold" && is.null(variables)) {
    variables <- c("temp_mean", "temp_min")
    if (is.null(invert)) invert <- c(TRUE, TRUE)
  } else if (stress.type == "herbivory" && is.null(variables)) {
    variables <- c("herbivory_aboveground", "herbivory_belowground")
    if (is.null(invert)) invert <- c(FALSE, FALSE)
  } else if (stress.type == "custom" && is.null(variables)) {
    stop("For custom stress type, 'variables' must be specified")
  }

  # Check variables exist
  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Extract variables
  stress_data <- data[, variables, drop = FALSE]

  # Invert specified variables
  if (!is.null(invert) && length(invert) == length(variables)) {
    for (i in seq_along(variables)) {
      if (invert[i]) {
        stress_data[[i]] <- -stress_data[[i]]
      }
    }
  }

  # Scale variables
  if (scale.vars) {
    stress_data <- as.data.frame(scale(stress_data))
  }

  # Set weights
  if (is.null(weights)) {
    weights <- rep(1 / length(variables), length(variables))
  } else if (length(weights) != length(variables)) {
    stop("Length of weights must match length of variables")
  } else {
    weights <- weights / sum(weights)
  }

  # Calculate stress index
  data$stress_index <- as.matrix(stress_data) %*% weights

  return(data)
}

#' Compare conditions across groups
#'
#' Perform statistical comparisons of variables across groups with multiple
#' test options.
#'
#' @param data Data frame
#' @param group.by Character, column name to group by
#' @param variables Character vector of variables to compare
#' @param test Statistical test: "anova", "kruskal", "ttest" (default: "anova")
#' @param p.adjust Method for p-value adjustment: "none", "bonferroni", "holm", "fdr" (default: "fdr")
#' @param plot Logical, create comparison plots (default: TRUE)
#' @return List containing test results and optionally plots
#' @export
#' @importFrom stats aov kruskal.test t.test p.adjust
#' @examples
#' env <- envirodata(n_sites = 5)
#' results <- compareconditions(env, group.by = "site_id",
#'                              variables = c("temp_mean", "precipitation"))
compareconditions <- function(data,
                             group.by,
                             variables,
                             test = c("anova", "kruskal", "ttest"),
                             p.adjust = c("fdr", "bonferroni", "holm", "none"),
                             plot = TRUE) {

  test <- match.arg(test)
  p.adjust <- match.arg(p.adjust)

  # Check inputs
  if (!group.by %in% names(data)) {
    stop("group.by column not found in data")
  }

  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Initialize results
  results <- data.frame(
    variable = variables,
    test = test,
    statistic = NA,
    p.value = NA,
    p.adjusted = NA
  )

  # Perform tests
  for (i in seq_along(variables)) {
    var <- variables[i]
    formula_obj <- stats::as.formula(paste(var, "~", group.by))

    if (test == "anova") {
      fit <- stats::aov(formula_obj, data = data)
      summary_fit <- summary(fit)
      results$statistic[i] <- summary_fit[[1]]$`F value`[1]
      results$p.value[i] <- summary_fit[[1]]$`Pr(>F)`[1]
    } else if (test == "kruskal") {
      fit <- stats::kruskal.test(formula_obj, data = data)
      results$statistic[i] <- fit$statistic
      results$p.value[i] <- fit$p.value
    } else if (test == "ttest") {
      groups <- unique(data[[group.by]])
      if (length(groups) != 2) {
        warning("t-test requires exactly 2 groups. Skipping variable: ", var)
        next
      }
      fit <- stats::t.test(formula_obj, data = data)
      results$statistic[i] <- fit$statistic
      results$p.value[i] <- fit$p.value
    }
  }

  # Adjust p-values
  if (p.adjust != "none") {
    results$p.adjusted <- stats::p.adjust(results$p.value, method = p.adjust)
  } else {
    results$p.adjusted <- results$p.value
  }

  # Create plots if requested
  plots <- NULL
  if (plot) {
    plots <- list()
    for (var in variables) {
      p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[group.by]], y = .data[[var]])) +
        ggplot2::geom_boxplot(fill = "steelblue", alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Comparison of", var),
                     x = group.by,
                     y = var)
      plots[[var]] <- p
    }
  }

  return(list(results = results, plots = plots))
}

#' Aggregate data over time periods
#'
#' Aggregate temporal data by specified time periods with multiple aggregation
#' functions.
#'
#' @param data Data frame with temporal data
#' @param time.col Character, name of time column (must be numeric for months/weeks)
#' @param variables Character vector of variables to aggregate
#' @param group.by Character vector of grouping columns (default: NULL)
#' @param agg.fun Aggregation function: "mean", "sum", "max", "min", "sd" (default: "mean")
#' @param period Time period: "week", "month", "season", "year" (default: "month")
#' @return Aggregated data frame
#' @export
#' @importFrom stats sd
#' @examples
#' env <- envirodata(n_months = 24)
#' env_monthly <- aggregatetemporal(env, time.col = "month",
#'                                  variables = c("temp_mean", "precipitation"),
#'                                  group.by = "site_id")
aggregatetemporal <- function(data,
                             time.col,
                             variables,
                             group.by = NULL,
                             agg.fun = c("mean", "sum", "max", "min", "sd"),
                             period = c("month", "week", "season", "year")) {

  agg.fun <- match.arg(agg.fun)
  period <- match.arg(period)

  # Check inputs
  if (!time.col %in% names(data)) {
    stop("time.col not found in data")
  }

  missing_vars <- setdiff(variables, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # Create period column
  if (period == "week") {
    data$period <- ceiling(data[[time.col]] / 7)
  } else if (period == "month") {
    data$period <- data[[time.col]]
  } else if (period == "season") {
    data$period <- ceiling(data[[time.col]] / 3)
  } else if (period == "year") {
    data$period <- ceiling(data[[time.col]] / 12)
  }

  # Select aggregation function
  agg_func <- switch(agg.fun,
                     "mean" = mean,
                     "sum" = sum,
                     "max" = max,
                     "min" = min,
                     "sd" = stats::sd)

  # Prepare grouping variables
  group_vars <- c("period", group.by)

  # Aggregate
  result <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(variables),
                                   ~agg_func(.x, na.rm = TRUE)),
                     .groups = "drop")

  return(as.data.frame(result))
}

#' Detect outliers in data
#'
#' Identify outliers using multiple detection methods.
#'
#' @param data Data frame
#' @param variables Character vector of variables to check (default: all numeric)
#' @param method Detection method: "iqr", "zscore", "mad" (default: "iqr")
#' @param threshold Threshold value (default: 3 for zscore/mad, 1.5 for iqr)
#' @param mark.only Logical, only mark outliers without removing (default: TRUE)
#' @return Data frame with outlier information
#' @export
#' @importFrom stats mad
#' @examples
#' env <- envirodata()
#' env_outliers <- detectoutliers(env, variables = c("temp_mean", "precipitation"))
#' table(env_outliers$is_outlier)
detectoutliers <- function(data,
                          variables = NULL,
                          method = c("iqr", "zscore", "mad"),
                          threshold = NULL,
                          mark.only = TRUE) {

  method <- match.arg(method)

  # Select numeric variables if not specified
  if (is.null(variables)) {
    variables <- names(data)[sapply(data, is.numeric)]
  }

  # Set default threshold
  if (is.null(threshold)) {
    threshold <- if (method == "iqr") 1.5 else 3
  }

  # Initialize outlier marker
  data$is_outlier <- FALSE
  data$outlier_variable <- NA_character_

  for (var in variables) {
    x <- data[[var]]
    outliers <- rep(FALSE, length(x))

    if (method == "iqr") {
      q1 <- stats::quantile(x, 0.25, na.rm = TRUE)
      q3 <- stats::quantile(x, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      outliers <- x < (q1 - threshold * iqr) | x > (q3 + threshold * iqr)
    } else if (method == "zscore") {
      z <- abs((x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE))
      outliers <- z > threshold
    } else if (method == "mad") {
      med <- stats::median(x, na.rm = TRUE)
      mad_val <- stats::mad(x, na.rm = TRUE)
      outliers <- abs(x - med) / mad_val > threshold
    }

    outliers[is.na(outliers)] <- FALSE
    data$is_outlier <- data$is_outlier | outliers
    data$outlier_variable[outliers] <- paste(data$outlier_variable[outliers], var, sep = ";")
  }

  # Clean up outlier_variable column
  data$outlier_variable <- gsub("^NA;", "", data$outlier_variable)
  data$outlier_variable <- gsub("^;", "", data$outlier_variable)

  if (!mark.only) {
    data <- data[!data$is_outlier, ]
  }

  return(data)
}

#' Correlate phenotype with environment
#'
#' Calculate correlations between phenotypic traits and environmental variables.
#'
#' @param pheno.data Phenotypic data frame
#' @param env.data Environmental data frame
#' @param pheno.vars Character vector of phenotypic variables (default: all numeric)
#' @param env.vars Character vector of environmental variables (default: all numeric)
#' @param by Character vector of columns to join datasets (default: "site_id")
#' @param method Correlation method: "pearson", "spearman" (default: "pearson")
#' @param plot Logical, create correlation heatmap (default: TRUE)
#' @return List with correlation matrix and optionally plot
#' @export
#' @importFrom stats cor
#' @examples
#' pheno <- phenotydata(n_trees = 30)
#' env <- envirodata(n_sites = 6)
#' cor_results <- correlatephenoenv(pheno, env)
correlatephenoenv <- function(pheno.data,
                             env.data,
                             pheno.vars = NULL,
                             env.vars = NULL,
                             by = "site_id",
                             method = c("pearson", "spearman"),
                             plot = TRUE) {

  method <- match.arg(method)

  # Merge datasets
  merged <- dplyr::inner_join(pheno.data, env.data, by = by)

  # Select variables
  if (is.null(pheno.vars)) {
    pheno.vars <- names(pheno.data)[sapply(pheno.data, is.numeric)]
  }
  if (is.null(env.vars)) {
    env.vars <- names(env.data)[sapply(env.data, is.numeric)]
  }

  # Calculate correlations
  cor_matrix <- stats::cor(merged[, pheno.vars],
                          merged[, env.vars],
                          method = method,
                          use = "pairwise.complete.obs")

  # Create plot
  cor_plot <- NULL
  if (plot) {
    # Reshape for plotting
    cor_long <- as.data.frame(cor_matrix)
    cor_long$pheno <- rownames(cor_long)
    cor_long <- tidyr::pivot_longer(cor_long,
                                    cols = -pheno,
                                    names_to = "env",
                                    values_to = "correlation")

    cor_plot <- ggplot2::ggplot(cor_long,
                                ggplot2::aes(x = env, y = pheno, fill = correlation)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                   midpoint = 0, limits = c(-1, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = paste0("Phenotype-Environment Correlations (", method, ")"),
                   x = "Environmental Variables",
                   y = "Phenotypic Traits",
                   fill = "Correlation")
  }

  return(list(correlation_matrix = cor_matrix, plot = cor_plot))
}

#' Run PCA on omics data
#'
#' Perform Principal Component Analysis with visualization options.
#'
#' @param data Data frame with numeric variables
#' @param variables Character vector of variables to include (default: all numeric)
#' @param scale Logical, scale variables before PCA (default: TRUE)
#' @param center Logical, center variables before PCA (default: TRUE)
#' @param n.components Number of components to return (default: 5)
#' @param color.by Character, column name for coloring points in plot (default: NULL)
#' @param plot Logical, create PCA plots (default: TRUE)
#' @return List with PCA results, variance explained, and optionally plots
#' @export
#' @importFrom stats prcomp
#' @examples
#' geno <- genotydata(n_trees = 50)
#' pca_results <- runpca(geno, n.components = 3)
runpca <- function(data,
                  variables = NULL,
                  scale = TRUE,
                  center = TRUE,
                  n.components = 5,
                  color.by = NULL,
                  plot = TRUE) {

  # Select numeric variables
  if (is.null(variables)) {
    variables <- names(data)[sapply(data, is.numeric)]
  }

  # Run PCA
  pca_result <- stats::prcomp(data[, variables],
                             scale. = scale,
                             center = center)

  # Calculate variance explained
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

  # Extract scores
  scores <- as.data.frame(pca_result$x[, 1:min(n.components, ncol(pca_result$x))])

  # Add color variable if specified
  if (!is.null(color.by) && color.by %in% names(data)) {
    scores$color_var <- data[[color.by]]
  }

  # Create plots
  plots <- NULL
  if (plot) {
    plots <- list()

    # Scree plot
    scree_data <- data.frame(
      PC = 1:length(var_explained),
      Variance = var_explained
    )
    plots$scree <- ggplot2::ggplot(scree_data, ggplot2::aes(x = PC, y = Variance)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::geom_line(color = "red") +
      ggplot2::geom_point(color = "red") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Scree Plot",
                   x = "Principal Component",
                   y = "Variance Explained (%)")

    # Biplot PC1 vs PC2
    if (!is.null(color.by)) {
      plots$biplot <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2, color = color_var)) +
        ggplot2::geom_point(size = 2, alpha = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "PCA Biplot",
                     x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
                     y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
                     color = color.by)
    } else {
      plots$biplot <- ggplot2::ggplot(scores, ggplot2::aes(x = PC1, y = PC2)) +
        ggplot2::geom_point(size = 2, alpha = 0.7, color = "steelblue") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "PCA Biplot",
                     x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
                     y = paste0("PC2 (", round(var_explained[2], 1), "%)"))
    }
  }

  return(list(
    pca = pca_result,
    scores = scores,
    variance_explained = var_explained,
    plots = plots
  ))
}

#' Summarize data by groups
#'
#' Calculate summary statistics for variables grouped by one or more factors.
#'
#' @param data Data frame
#' @param group.vars Character vector of grouping variables
#' @param summary.vars Character vector of variables to summarize
#' @param funs Character vector of summary functions: "mean", "sd", "median", "min", "max", "n" (default: c("mean", "sd", "n"))
#' @param na.rm Logical, remove NA values (default: TRUE)
#' @return Data frame with summary statistics
#' @export
#' @importFrom stats median sd
#' @examples
#' env <- envirodata()
#' summary_stats <- summarizebygroup(env,
#'                                   group.vars = "site_id",
#'                                   summary.vars = c("temp_mean", "precipitation"))
summarizebygroup <- function(data,
                            group.vars,
                            summary.vars,
                            funs = c("mean", "sd", "n"),
                            na.rm = TRUE) {

  # Check inputs
  missing_groups <- setdiff(group.vars, names(data))
  if (length(missing_groups) > 0) {
    stop("Grouping variables not found: ", paste(missing_groups, collapse = ", "))
  }

  missing_vars <- setdiff(summary.vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Summary variables not found: ", paste(missing_vars, collapse = ", "))
  }

  # Create summary functions
  summary_list <- list()
  for (fun in funs) {
    if (fun == "n") {
      summary_list[[fun]] <- function(x) sum(!is.na(x))
    } else if (fun == "mean") {
      summary_list[[fun]] <- function(x) mean(x, na.rm = na.rm)
    } else if (fun == "sd") {
      summary_list[[fun]] <- function(x) stats::sd(x, na.rm = na.rm)
    } else if (fun == "median") {
      summary_list[[fun]] <- function(x) stats::median(x, na.rm = na.rm)
    } else if (fun == "min") {
      summary_list[[fun]] <- function(x) min(x, na.rm = na.rm)
    } else if (fun == "max") {
      summary_list[[fun]] <- function(x) max(x, na.rm = na.rm)
    }
  }

  # Perform summary
  result <- data |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group.vars))) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(summary.vars),
                   summary_list,
                   .names = "{.col}_{.fn}"),
      .groups = "drop"
    )

  return(as.data.frame(result))
}
