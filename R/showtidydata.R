#' Explore data distribution and correlations
#'
#' Converts data to long format, creates distribution plots for all numeric
#' variables, and generates a correlation heatmap.
#'
#' @param dataframe A data frame to analyze
#' @param cor_method Correlation method: "pearson" or "spearman" (default: "pearson")
#' @return A list containing the long-format data, distribution plot, and correlation plot
#' @export
#' @importFrom stats cor
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
  
  # Calculate correlation matrix
  cor_matrix <- stats::cor(numeric_df, 
                           method = cor_method, 
                           use = "pairwise.complete.obs")
  
  # Create correlation heatmap
  corrplot::corrplot(cor_matrix, 
                     method = "color",
                     type = "upper",
                     order = "hclust",
                     addCoef.col = "black",
                     tl.col = "black",
                     tl.srt = 45,
                     title = paste0("Correlation Heatmap (", cor_method, ")"),
                     mar = c(0, 0, 2, 0))
  
  # Return results
  return(list(
    long_data = long_data,
    distribution_plot = dist_plot,
    correlation_matrix = cor_matrix
  ))
}