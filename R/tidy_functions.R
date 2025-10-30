#' Explore data distribution and correlations
#'
#' Converts data to long format, creates distribution plots for all numeric
#' variables, and generates a correlation heatmap.
#'
#' @param dataframe A data frame to analyze
#' @param cor_method Correlation method: "pearson" or "spearman" (default: "pearson")
#' @return A list containing the long-format data, distribution plot, and correlation plot
#' @export
#' @importFrom stats cor sd complete.cases
#' @examples
#' \dontrun{
#' df <- data.frame(temp = rnorm(100, 20, 5), 
#'                  precip = rnorm(100, 50, 10),
#'                  gdd = rnorm(100, 1500, 200))
#' result <- showtidydata(df)
#' }
showtidydata <- function(dataframe, cor_method = "pearson") {
  
  # Avoid NSE notes
  id <- value <- NULL
  
  # Select only numeric columns
  numeric_df <- dataframe |>
    dplyr::select(dplyr::where(is.numeric))
  
  if (ncol(numeric_df) == 0) {
    stop("No numeric columns found in dataframe")
  }
  
  # Remove columns with zero variance (constant values)
  var_cols <- sapply(numeric_df, function(x) stats::sd(x, na.rm = TRUE) > 0)
  numeric_df <- numeric_df[, var_cols, drop = FALSE]
  
  if (ncol(numeric_df) == 0) {
    stop("All numeric columns have zero variance")
  }
  
  if (ncol(numeric_df) < 2) {
    warning("Only one variable with variance found. Correlation plot skipped.")
    cor_matrix <- NULL
  }
  
  # Convert to long format
  long_data <- numeric_df |>
    dplyr::mutate(id = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -id, 
                        names_to = "variable", 
                        values_to = "value")
  
  # Create distribution plots
  dist_plot <- ggplot2::ggplot(long_data, 
                                ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Distribution of Variables",
                  x = "Value", 
                  y = "Count")
  
  print(dist_plot)
  
  # Calculate and plot correlation matrix only if multiple variables
  if (ncol(numeric_df) >= 2) {
    cor_matrix <- stats::cor(numeric_df, 
                             method = cor_method, 
                             use = "pairwise.complete.obs")
    
    # Check for valid correlation matrix
    if (any(is.na(cor_matrix)) || any(is.infinite(cor_matrix))) {
      warning("Correlation matrix contains NA or infinite values. Using complete cases only.")
      numeric_df <- numeric_df[stats::complete.cases(numeric_df), ]
      cor_matrix <- stats::cor(numeric_df, method = cor_method)
    }
    
    # Create correlation heatmap
    corrplot::corrplot(cor_matrix, 
                       method = "color",
                       type = "upper",
                       order = "original",
                       addCoef.col = "black",
                       number.cex = 0.7,
                       tl.col = "black",
                       tl.srt = 45,
                       title = paste0("Correlation Heatmap (", cor_method, ")"),
                       mar = c(0, 0, 2, 0))
  } else {
    cor_matrix <- NULL
  }
  
  # Return results
  return(invisible(list(
    long_data = long_data,
    distribution_plot = dist_plot,
    correlation_matrix = cor_matrix
  )))
}

#' Create tidy long format data
#'
#' Converts a dataframe to tidy long format without saving. Returns the 
#' transformed dataframe for further analysis or manipulation.
#'
#' @param dataframe A data frame to convert
#' @param id_col Name for the ID column in long format (default: "id")
#' @param keep_cols Character vector of column names to keep as identifiers 
#'   (e.g., site_id, tree_id). These won't be pivoted. (default: NULL)
#' @param var_types Types of variables to include: "numeric", "character", 
#'   "factor", or "all" (default: "numeric")
#' @param var_name Name for the variable column (default: "variable")
#' @param value_name Name for the value column (default: "value")
#' @param add_timestamp Logical, add a timestamp column (default: FALSE)
#' @param na_rm Logical, remove rows with NA values (default: FALSE)
#' @param arrange_by Character vector of column names to arrange by (default: NULL)
#' @return A tidy dataframe in long format
#' @export
#' @importFrom rlang :=
#' @examples
#' \dontrun{
#' data(envirodata)
#' tidy_env <- createtidydata(envirodata, keep_cols = "site_id")
#' 
#' data(phenotydata)
#' tidy_pheno <- createtidydata(phenotydata, 
#'                              keep_cols = c("tree_id", "site_id"),
#'                              var_name = "trait",
#'                              value_name = "measurement",
#'                              arrange_by = c("tree_id", "trait"))
#' }
createtidydata <- function(dataframe, 
                           id_col = "id",
                           keep_cols = NULL,
                           var_types = c("numeric", "character", "factor", "all"),
                           var_name = "variable",
                           value_name = "value",
                           add_timestamp = FALSE,
                           na_rm = FALSE,
                           arrange_by = NULL) {
  
  # Match arguments
  var_types <- match.arg(var_types)
  
  # Select columns based on var_types
  if (var_types == "numeric") {
    data_cols <- dataframe |> dplyr::select(dplyr::where(is.numeric))
  } else if (var_types == "character") {
    data_cols <- dataframe |> dplyr::select(dplyr::where(is.character))
  } else if (var_types == "factor") {
    data_cols <- dataframe |> dplyr::select(dplyr::where(is.factor))
  } else {
    data_cols <- dataframe
  }
  
  # Separate keep_cols from data to pivot
  if (!is.null(keep_cols)) {
    if (!all(keep_cols %in% names(dataframe))) {
      stop("Some keep_cols not found in dataframe: ", 
           paste(setdiff(keep_cols, names(dataframe)), collapse = ", "))
    }
    id_data <- dataframe[, keep_cols, drop = FALSE]
    data_cols <- data_cols |> 
      dplyr::select(-dplyr::any_of(keep_cols))
  } else {
    id_data <- NULL
  }
  
  if (ncol(data_cols) == 0) {
    stop("No columns to pivot after applying filters")
  }
  
  # Add row ID
  data_cols <- data_cols |> 
    dplyr::mutate(!!id_col := dplyr::row_number())
  
  # Convert to long format
  tidy_data <- data_cols |>
    tidyr::pivot_longer(cols = -dplyr::all_of(id_col),
                        names_to = var_name,
                        values_to = value_name)
  
  # Add back identifier columns
  if (!is.null(id_data)) {
    id_data <- id_data |> 
      dplyr::mutate(!!id_col := dplyr::row_number())
    tidy_data <- tidy_data |>
      dplyr::left_join(id_data, by = id_col)
  }
  
  # Add timestamp if requested
  if (add_timestamp) {
    tidy_data <- tidy_data |>
      dplyr::mutate(timestamp = Sys.time())
  }
  
  # Remove NAs if requested
  if (na_rm) {
    tidy_data <- tidy_data |> 
      tidyr::drop_na(dplyr::all_of(value_name))
  }
  
  # Arrange if requested
  if (!is.null(arrange_by)) {
    if (!all(arrange_by %in% names(tidy_data))) {
      warning("Some arrange_by columns not found in tidy data. Skipping arrangement.")
    } else {
      tidy_data <- tidy_data |>
        dplyr::arrange(dplyr::across(dplyr::all_of(arrange_by)))
    }
  }
  
  return(tidy_data)
}

#' Save data in tidy long format
#'
#' Converts a dataframe to tidy long format and saves it to a file. Optionally
#' includes metadata columns (ID, timestamp, grouping variables) and allows
#' filtering of variable types.
#'
#' @param dataframe A data frame to convert and save
#' @param filename Output filename (default: "tidy_data.csv")
#' @param format Output format: "csv", "tsv", "rds", "xlsx" (default: "csv")
#' @param id_col Name for the ID column in long format (default: "id")
#' @param keep_cols Character vector of column names to keep as identifiers 
#'   (e.g., site_id, tree_id). These won't be pivoted. (default: NULL)
#' @param var_types Types of variables to include: "numeric", "character", 
#'   "factor", or "all" (default: "numeric")
#' @param add_timestamp Logical, add a timestamp column (default: FALSE)
#' @param na_rm Logical, remove rows with NA values (default: FALSE)
#' @return Invisibly returns the tidy dataframe and prints save confirmation
#' @export
#' @importFrom rlang :=
#' @examples
#' \dontrun{
#' data(envirodata)
#' savetidydata(envirodata, "environment_long.csv", keep_cols = "site_id")
#' savetidydata(phenotydata, "phenotype_long.xlsx", format = "xlsx", 
#'              keep_cols = c("tree_id", "site_id"))
#' }
savetidydata <- function(dataframe, 
                         filename = "tidy_data.csv",
                         format = c("csv", "tsv", "rds", "xlsx"),
                         id_col = "id",
                         keep_cols = NULL,
                         var_types = c("numeric", "character", "factor", "all"),
                         add_timestamp = FALSE,
                         na_rm = FALSE) {
  
  # Match arguments
  format <- match.arg(format)
  var_types <- match.arg(var_types)
  
  # Avoid NSE notes
  variable <- value <- NULL
  
  # Select columns based on var_types
  if (var_types == "numeric") {
    data_cols <- dataframe |> dplyr::select(dplyr::where(is.numeric))
  } else if (var_types == "character") {
    data_cols <- dataframe |> dplyr::select(dplyr::where(is.character))
  } else if (var_types == "factor") {
    data_cols <- dataframe |> dplyr::select(dplyr::where(is.factor))
  } else {
    data_cols <- dataframe
  }
  
  # Separate keep_cols from data to pivot
  if (!is.null(keep_cols)) {
    if (!all(keep_cols %in% names(dataframe))) {
      stop("Some keep_cols not found in dataframe: ", 
           paste(setdiff(keep_cols, names(dataframe)), collapse = ", "))
    }
    id_data <- dataframe[, keep_cols, drop = FALSE]
    data_cols <- data_cols |> 
      dplyr::select(-dplyr::any_of(keep_cols))
  } else {
    id_data <- NULL
  }
  
  if (ncol(data_cols) == 0) {
    stop("No columns to pivot after applying filters")
  }
  
  # Add row ID
  data_cols <- data_cols |> 
    dplyr::mutate(!!id_col := dplyr::row_number())
  
  # Convert to long format
  tidy_data <- data_cols |>
    tidyr::pivot_longer(cols = -dplyr::all_of(id_col),
                        names_to = "variable",
                        values_to = "value")
  
  # Add back identifier columns
  if (!is.null(id_data)) {
    id_data <- id_data |> 
      dplyr::mutate(!!id_col := dplyr::row_number())
    tidy_data <- tidy_data |>
      dplyr::left_join(id_data, by = id_col)
  }
  
  # Add timestamp if requested
  if (add_timestamp) {
    tidy_data <- tidy_data |>
      dplyr::mutate(timestamp = Sys.time())
  }
  
  # Remove NAs if requested
  if (na_rm) {
    n_before <- nrow(tidy_data)
    tidy_data <- tidy_data |> tidyr::drop_na(value)
    n_removed <- n_before - nrow(tidy_data)
    message(sprintf("Removed %d rows with NA values", n_removed))
  }
  
  # Save based on format
  if (format == "csv") {
    readr::write_csv(tidy_data, filename)
  } else if (format == "tsv") {
    readr::write_tsv(tidy_data, filename)
  } else if (format == "rds") {
    saveRDS(tidy_data, filename)
  } else if (format == "xlsx") {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Package 'writexl' needed for xlsx format. Install with: install.packages('writexl')")
    }
    writexl::write_xlsx(tidy_data, filename)
  }
  
  message(sprintf("Tidy data saved to: %s", filename))
  message(sprintf("Dimensions: %d rows x %d columns", nrow(tidy_data), ncol(tidy_data)))
  
  return(invisible(tidy_data))
}