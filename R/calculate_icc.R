library(lme4)
#' Calculate Intraclass Correlation Coefficient (ICC) from a Null Mixed-Effects Model
#'
#' This function fits a null mixed-effects model (a random intercept only 
#' model without any predictors) using lmer and calculates the ICC, 
#' which is the ratio of between-cluster variance to total variance.
#'
#' @param data The input data frame.
#' @param dependent_var The name of the dependent (outcome) variable.
#' @param cluster_id_col The name of the column containing the cluster IDs (e.g., school IDs).
#' @param min_icc The minimum ICC value to determine whether GLMM is recommended or not.
#' 
#' @return A list containing the ICC value and the lmer model object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'data' exists with 'score' and 'school_id' columns
#' result <- calculate_icc(
#'   data = score, 
#'   dependent_var = "score", 
#'   cluster_id_col = "IDSCHOOL",
#'   min_icc = 0.1
#' )
#' print(result$icc)
#' summary(result$model)
#' }
calculate_icc <- function(data, dependent_var, cluster_id_col, min_icc = 0.1) {
  
  # Ensure lme4 is available
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("The 'lme4' package is required but not installed. Please install it.")
  }
  
  # 1. Construct the Model Formula
  # Formula structure: outcome ~ 1 + (1 | cluster_id_col)
  # This is the null model: fixed intercept + random intercept by cluster
  formula_str <- paste0(
    dependent_var, 
    " ~ 1 + (1 | ", 
    cluster_id_col, 
    ")"
  )
  null_formula <- as.formula(formula_str)
  
  # 2. Fit the Null Model
  null_model <- lme4::lmer(null_formula, data = data)
  
  # 3. Extract Variance Components
  # VarCorr() returns a list of variance/covariance components
  var_components <- as.data.frame(lme4::VarCorr(null_model))
  tau_sq <- var_components[1, "vcov"] # Between-Cluster Variance
  sigma_sq <- var_components[2, "vcov"] # Residual Variance
  
  # 4. Calculate ICC
  # ICC = Between-Cluster Variance / (Between-Cluster Variance + Within-Cluster Variance)
  ICC <- tau_sq / (tau_sq + sigma_sq)
  
  # 5. Provide Interpretation
  if (ICC > min_icc) {
    message(paste("ICC is", round(ICC, 3), ". Since ICC > ", min_icc, ", GLMM is recommended."))
  } else {
    message(paste("ICC is", round(ICC, 3), ". Since ICC <= ", min_icc, ", GLMM may not be necessary."))
  }
  
  # 6. Return the ICC and the model object
  return(list(
    icc = ICC,
    model = null_model
  ))
}