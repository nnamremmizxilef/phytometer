#' Calculate growing degree days
#'
#' @param temp_min Vector of minimum daily temperatures
#' @param temp_max Vector of maximum daily temperatures
#' @param base_temp Base temperature threshold (default 5°C)
#' @return Vector of growing degree days
#' @export
#' @examples
#' calculate_gdd(c(10, 12), c(20, 22), base_temp = 5)
calculate_gdd <- function(temp_min, temp_max, base_temp = 5) {
  temp_mean <- (temp_min + temp_max) / 2
  gdd <- pmax(temp_mean - base_temp, 0)
  return(gdd)
}