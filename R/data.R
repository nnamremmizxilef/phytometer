#' MOCK environmental monitoring data from oak holobiont study
#'
#' Monthly environmental measurements from 10 sites monitoring Quercus robur
#' clone DF159 across environmental gradients.
#'
#' @format A data frame with 120 rows and 13 variables:
#' \describe{
#'   \item{site_id}{Site identifier}
#'   \item{month}{Month (1-12)}
#'   \item{year}{Year of measurement}
#'   \item{temp_mean}{Mean temperature (°C)}
#'   \item{temp_max}{Maximum temperature (°C)}
#'   \item{temp_min}{Minimum temperature (°C)}
#'   \item{precipitation}{Precipitation (mm)}
#'   \item{soil_moisture}{Soil moisture content (%)}
#'   \item{vpd}{Vapor pressure deficit (kPa)}
#'   \item{par}{Photosynthetically active radiation (μmol/m²/s)}
#'   \item{drought_index}{Drought stress index (0-1)}
#'   \item{herbivory_aboveground}{Aboveground herbivory count}
#'   \item{herbivory_belowground}{Belowground herbivory count}
#' }
#' @source PhytOakmeter research unit
"envirodata"

#' MOCK Phenotypic measurements from oak saplings
#'
#' Physiological and morphological measurements from Quercus robur clone DF159
#' saplings across 10 sites.
#'
#' @format A data frame with 50 rows and 12 variables:
#' \describe{
#'   \item{tree_id}{Tree identifier}
#'   \item{site_id}{Site identifier}
#'   \item{clone}{Clone designation (DF159)}
#'   \item{height_cm}{Tree height (cm)}
#'   \item{dbh_mm}{Diameter at breast height (mm)}
#'   \item{leaf_area_cm2}{Total leaf area (cm²)}
#'   \item{root_biomass_g}{Root biomass (g)}
#'   \item{shoot_biomass_g}{Shoot biomass (g)}
#'   \item{stomatal_conductance}{Stomatal conductance (mol/m²/s)}
#'   \item{photosynthesis_rate}{Photosynthesis rate (μmol CO₂/m²/s)}
#'   \item{water_use_efficiency}{Water use efficiency (μmol CO₂/mmol H₂O)}
#'   \item{stress_tolerance_index}{Stress tolerance index (0-1)}
#' }
#' @source PhytOakmeter research unit
"phenotydata"

#' MOCK Genotypic and molecular data from oak holobiont
#'
#' Multi-omics data including SNPs, epigenetic marks, gene expression, and
#' microbiome diversity from Quercus robur clone DF159.
#'
#' @format A data frame with 50 rows and 12 variables:
#' \describe{
#'   \item{tree_id}{Tree identifier}
#'   \item{site_id}{Site identifier}
#'   \item{snp_001 to snp_005}{SNP genotypes (0=homozygous ref, 1=heterozygous, 2=homozygous alt)}
#'   \item{epigenetic_mark_1, epigenetic_mark_2}{Epigenetic modification levels}
#'   \item{gene_expression_drought}{Drought response gene expression level}
#'   \item{gene_expression_herbivory}{Herbivory response gene expression level}
#'   \item{microbiome_diversity}{Shannon diversity index of associated microbiome}
#' }
#' @source PhytOakmeter research unit
"genotydata"