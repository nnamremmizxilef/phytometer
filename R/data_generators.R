#' Generate example environmental data
#'
#' Creates a customizable environmental monitoring dataset suitable for testing
#' holobiont response analyses.
#'
#' @param n_sites Number of monitoring sites (default: 10)
#' @param n_months Number of months per site (default: 12)
#' @param seed Random seed for reproducibility (default: 123)
#' @param temp_range Temperature range in °C as c(min, max) (default: c(5, 25))
#' @param precip_mean Mean precipitation in mm (default: 60)
#' @return A data frame with environmental measurements
#' @export
#' @importFrom stats rnorm runif rpois
#' @examples
#' # Default dataset
#' env <- envirodata()
#' 
#' # Larger dataset with more sites
#' env_large <- envirodata(n_sites = 20, n_months = 24)
#' 
#' # Custom temperature range
#' env_cold <- envirodata(temp_range = c(-5, 15))
envirodata <- function(n_sites = 10, 
                       n_months = 12, 
                       seed = 123,
                       temp_range = c(5, 25),
                       precip_mean = 60) {
  
  set.seed(seed)
  n_obs <- n_sites * n_months
  temp_mean <- mean(temp_range)
  temp_amplitude <- diff(temp_range) / 2
  
  data.frame(
    site_id = rep(paste0("Site_", 1:n_sites), each = n_months),
    month = rep(1:n_months, times = n_sites),
    year = 2024,
    temp_mean = rnorm(n_obs, temp_mean, 3) + 
      rep(sin(seq(0, 2*pi, length.out = n_months)) * temp_amplitude, n_sites),
    temp_max = rnorm(n_obs, temp_mean + 7, 4) + 
      rep(sin(seq(0, 2*pi, length.out = n_months)) * (temp_amplitude + 2), n_sites),
    temp_min = rnorm(n_obs, temp_mean - 7, 3) + 
      rep(sin(seq(0, 2*pi, length.out = n_months)) * (temp_amplitude - 2), n_sites),
    precipitation = abs(rnorm(n_obs, precip_mean, 30)),
    soil_moisture = runif(n_obs, 10, 40),
    vpd = abs(rnorm(n_obs, 1.2, 0.6)),
    par = abs(rnorm(n_obs, 400, 150)),
    drought_index = runif(n_obs, 0, 1),
    herbivory_aboveground = rpois(n_obs, 3),
    herbivory_belowground = rpois(n_obs, 2)
  )
}

#' Generate example phenotypic data
#'
#' Creates customizable phenotypic measurements from oak saplings.
#'
#' @param n_trees Number of trees (default: 50)
#' @param n_sites Number of sites (default: 10)
#' @param clone Clone designation (default: "DF159")
#' @param seed Random seed for reproducibility (default: 123)
#' @param height_range Height range in cm as c(min, max) (default: c(150, 210))
#' @return A data frame with phenotypic measurements
#' @export
#' @examples
#' # Default dataset
#' pheno <- phenotydata()
#' 
#' # More trees
#' pheno_large <- phenotydata(n_trees = 100, n_sites = 20)
#' 
#' # Smaller saplings
#' pheno_small <- phenotydata(height_range = c(80, 120))
phenotydata <- function(n_trees = 50,
                        n_sites = 10,
                        clone = "DF159",
                        seed = 123,
                        height_range = c(150, 210)) {
  
  set.seed(seed)
  height_mean <- mean(height_range)
  height_sd <- diff(height_range) / 6
  
  data.frame(
    tree_id = paste0("Tree_", 1:n_trees),
    site_id = rep(paste0("Site_", 1:n_sites), length.out = n_trees),
    clone = clone,
    height_cm = rnorm(n_trees, height_mean, height_sd),
    dbh_mm = rnorm(n_trees, 45, 10),
    leaf_area_cm2 = rnorm(n_trees, 120, 25),
    root_biomass_g = abs(rnorm(n_trees, 85, 20)),
    shoot_biomass_g = abs(rnorm(n_trees, 150, 35)),
    stomatal_conductance = abs(rnorm(n_trees, 0.3, 0.1)),
    photosynthesis_rate = abs(rnorm(n_trees, 12, 3)),
    water_use_efficiency = abs(rnorm(n_trees, 4.5, 1.2)),
    stress_tolerance_index = runif(n_trees, 0, 1)
  )
}

#' Generate example genotypic data
#'
#' Creates customizable multi-omics data including SNPs, epigenetic marks,
#' gene expression, and microbiome diversity.
#'
#' @param n_trees Number of trees (default: 50)
#' @param n_sites Number of sites (default: 10)
#' @param n_snps Number of SNP markers (default: 5)
#' @param seed Random seed for reproducibility (default: 123)
#' @return A data frame with genotypic and molecular data
#' @export
#' @examples
#' # Default dataset
#' geno <- genotydata()
#' 
#' # More SNPs
#' geno_dense <- genotydata(n_snps = 20)
#' 
#' # Larger dataset
#' geno_large <- genotydata(n_trees = 100, n_sites = 20)
genotydata <- function(n_trees = 50,
                       n_sites = 10,
                       n_snps = 5,
                       seed = 123) {
  
  set.seed(seed)
  
  # Create SNP columns
  snp_data <- as.data.frame(
    matrix(sample(c(0, 1, 2), n_trees * n_snps, replace = TRUE),
           nrow = n_trees, ncol = n_snps)
  )
  names(snp_data) <- paste0("snp_", sprintf("%03d", 1:n_snps))
  
  # Combine with other data
  data.frame(
    tree_id = paste0("Tree_", 1:n_trees),
    site_id = rep(paste0("Site_", 1:n_sites), length.out = n_trees),
    snp_data,
    epigenetic_mark_1 = rnorm(n_trees, 0, 1),
    epigenetic_mark_2 = rnorm(n_trees, 0, 1),
    gene_expression_drought = abs(rnorm(n_trees, 100, 30)),
    gene_expression_herbivory = abs(rnorm(n_trees, 80, 25)),
    microbiome_diversity = abs(rnorm(n_trees, 3.5, 0.8))
  )
}