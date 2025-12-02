library(LongituRF)
library(randomForest)
library(lme4)
library(rlang)

#' Performs MERF-based Optimal Model Selection for Linear Mixed Models (LMM)
#'
#' This function integrates Mixed-Effects Random Forest (MERF) for predictor ranking
#' with a forward selection process to identify optimal linear mixed models (LMMs).
#' This sequentially performs predictor and model selection steps.
#'
#' @param data A data frame containing all variables.
#' @param dep_col_name A string for the dependent variable column name (e.g., "dep").
#' @param id_col_name A string for the cluster ID column name (e.g., "schID").
#' @param L1_cols A character vector of Level 1 predictor column names.
#' @param L2_cols A character vector of Level 2 predictor column names.
#' @param dep_is_continuous Logical. TRUE (default) if the dependent variable is continuous (Regression);
#'   FALSE if it is categorical (Classification). 
#' @param add_slopes Logical. If TRUE (default), significant random slopes are added.
#'   If FALSE, random slope selection is skipped.
#' @param alpha The significance level for the LMM forward selection (default is 0.05).
#' @param tol The tolerance for checking model singularity (default is 1e-4).
#' @param merf_iter The number of iterations for the MERF fitting process (default is 100).
#' @param ntree The number of trees in the MERF forest (default is 1000).
#' @return A list containing the results:
#'   \itemize{
#'     \item \code{optimal_model}: The final LMM identified.
#'     \item \code{final_merf}: The final MERF model fit with the optimal random effects structure.
#'     \item \code{included_predictors}: Character vector of significant fixed effects.
#'     \item \code{significant_interactions}: Character vector of significant interaction terms.
#'     \item \code{included_slopes}: Character vector of significant random slopes (empty if add_slopes = FALSE).
#'   }
#'
#' @examples
#' optimal_model <- optimalLMM(
#'   data = data,
#'   dep_col_name = "dep",
#'   id_col_name = "schID",
#'   L1_cols = paste0("Stu",1:50),
#'   L2_cols = paste0("Sch",1:50)
#' )
#'
#' optimal_model$optimal_model
#' optimal_model$final_merf
#' optimal_model$included_predictors
#' optimal_model$significant_interactions
#' optimal_model$included_slopes
optimalLMM <- function(
    data, dep_col_name, id_col_name, L1_cols, L2_cols,
    dep_is_continuous = TRUE, add_slopes = TRUE,
    alpha = 0.05, tol = 1e-4, merf_iter = 100, ntree = 1000
) {
  
  # --- 1. Setup and MERF Ranking ---
  cat("## 1. MERF Ranking for Predictor Importance\n")
  
  # Determine importance type
  importance_type <- ifelse(dep_is_continuous, 1, 2)
  cat(paste("Using importance type =", importance_type, "for", ifelse(dep_is_continuous, "Regression (IncMSE)", "Classification (MeanDecreaseAccuracy)"), ".\n"))
  
  # Prepare MERF input
  X_merf_cols <- c(L1_cols, L2_cols)
  X_merf <- data[, X_merf_cols, drop = FALSE]
  Y_merf <- data[[dep_col_name]]
  Z <- model.matrix(~ data[[id_col_name]] - 1)
  time <- rep(1, dim(data)[1])
  mtry_merf <- max(1, floor(ncol(X_merf) / 3))
  
  # Fit the initial MERF model
  merf.fit <- LongituRF::MERF(
    X = X_merf, Y = Y_merf, id = data[[id_col_name]], Z = Z, iter = merf_iter,
    mtry = mtry_merf, ntree = ntree, time = time, sto = "none", delta = 0.001
  )
  
  # Get predictor importance
  imp <- randomForest::importance(merf.fit$forest, type = importance_type)
  
  # Sort importance
  importance_name <- if (importance_type == 1) "IncMSE" else "MeanDecreaseAccuracy"
  full_sorted <- data.frame(
    Predictor = row.names(imp),
    Importance = if (importance_type == 1) as.vector(imp) else imp[, importance_name],
    stringsAsFactors = FALSE
  )
  full_sorted <- full_sorted[order(full_sorted$Importance, decreasing = TRUE), ]
  
  # Get sorted predictor names
  imp.col <- full_sorted$Predictor
  imp.col_L1 <- imp.col[imp.col %in% L1_cols]
  imp.col_L2 <- imp.col[imp.col %in% L2_cols]
  
  # Define interaction candidates automatically from top-3 predictors
  
  if (length(imp.col) >= 3) {
    top_3_predictors <- imp.col[1:3]
    interaction_candidates <- c(
      paste(top_3_predictors[1], top_3_predictors[2], sep = ":"),
      paste(top_3_predictors[1], top_3_predictors[3], sep = ":"),
      paste(top_3_predictors[2], top_3_predictors[3], sep = ":")
    )
    cat("\nInteractions among the top-3 MERF predictors:\n")
    print(top_3_predictors)
  } else {
    interaction_candidates <- character(0)
    cat("\nFewer than 3 total predictors. Skipping two-way interaction search.\n")
  }
  
  # 2. Forward Predictor Selection
  cat("\n## 2. Forward Predictor Selection (Fixed Effects)\n")
  
  # Fit the null (a random intercept) model
  null_formula <- as.formula(paste(dep_col_name, "~ 1 + (1 |", id_col_name, ")"))
  null_model <- lme4::lmer(null_formula, data = data, REML = FALSE)
  
  cat("\n-- Initial Model: Null LMM (Random Intercept Only) --\n")
  print(summary(null_model))
  
  simpler_model <- null_model
  included_predictors <- character(0)
  sig_inter <- character(0)
  included_slope <- character(0)
  
  check_model_conv <- function(model) {
    if (lme4::isSingular(model, tol = tol)) {
      message("Model is Singular (tol = ", tol, "). Skipping...")
      return(TRUE)
    }
    return(FALSE)
  }
  
  # Search for significant Level-1 predictors
  cat("\n-- Searching for Significant Level-1 Predictors --\n")
  for (q in 1:length(imp.col_L1)) {
    p_name <- imp.col_L1[q]
    trial_predictors <- c(included_predictors, p_name)
    trial_formula <- as.formula(paste(dep_col_name, "~", paste(trial_predictors, collapse = "+"), "+ (1 |", id_col_name, ")"))
    
    trial_model <- tryCatch({
      lme4::lmer(trial_formula, data = data, REML = FALSE)
    }, error = function(e) {
      message("Error fitting model with ", p_name, ". Skipping...")
      return(NULL)
    })
    
    if (is.null(trial_model) || check_model_conv(trial_model)) next
    
    pvalue <- tryCatch({
      as.numeric(anova(simpler_model, trial_model)[2, "Pr(>Chisq)"])
    }, error = function(e) {
      message("Error running ANOVA. Skipping...")
      return(NaN)
    })
    
    if (is.finite(pvalue) && pvalue < alpha) {
      cat(sprintf("L1 Predictor '%s' is significant (p = %.3g). Adding to model.\n", p_name, pvalue))
      simpler_model <- trial_model
      included_predictors <- trial_predictors
    } else {
      cat(sprintf("L1 Predictor '%s' is not significant (p = %.3g). Stopping L1 search.\n", p_name, pvalue))
      break
    }
  }
  
  # Search for significant Level-2 predictors
  cat("\n-- Searching for Significant Level-2 Predictors --\n")
  for (r in 1:length(imp.col_L2)) {
    p_name <- imp.col_L2[r]
    trial_predictors <- c(included_predictors, p_name)
    trial_formula <- as.formula(paste(dep_col_name, "~", paste(trial_predictors, collapse = "+"), "+ (1 |", id_col_name, ")"))
    
    trial_model <- tryCatch({
      lme4::lmer(trial_formula, data = data, REML = FALSE)
    }, error = function(e) {
      message("Error fitting model with ", p_name, ". Skipping...")
      return(NULL)
    })
    
    if (is.null(trial_model) || check_model_conv(trial_model)) next
    
    pvalue <- tryCatch({
      as.numeric(anova(simpler_model, trial_model)[2, "Pr(>Chisq)"])
    }, error = function(e) {
      message("Error running ANOVA. Skipping...")
      return(NaN)
    })
    
    if (is.finite(pvalue) && pvalue < alpha) {
      cat(sprintf("L2 Predictor '%s' is significant (p = %.3g). Adding to model.\n", p_name, pvalue))
      simpler_model <- trial_model
      included_predictors <- trial_predictors
    } else {
      cat(sprintf("L2 Predictor '%s' is not significant (p = %.3g). Stopping L2 search.\n", p_name, pvalue))
      break
    }
  }
  
  # Search for significant Interactions
  cat("\n-- Searching for Significant Interaction Terms --\n")
  
  if (length(interaction_candidates) > 0) {
    for (k in seq_along(interaction_candidates)) {
      inter_term <- interaction_candidates[k]
      
      trial_inter <- c(sig_inter, inter_term)
      
      trial_formula <- as.formula(paste(dep_col_name, "~", paste(c(included_predictors, trial_inter), collapse = "+"), "+ (1 |", id_col_name, ")"))
      
      trial_model <- tryCatch({
        lme4::lmer(trial_formula, data = data, REML = FALSE)
      }, error = function(e) {
        message("Error fitting model with interaction ", inter_term, ". Skipping...")
        return(NULL)
      })
      
      if (is.null(trial_model) || check_model_conv(trial_model)) next
      
      pvalue <- tryCatch({
        as.numeric(anova(simpler_model, trial_model)[2, "Pr(>Chisq)"])
      }, error = function(e) {
        message("Error running ANOVA. Skipping...")
        return(NaN)
      })
      
      if (is.finite(pvalue) && pvalue < alpha) {
        cat(sprintf("Interaction term '%s' is significant (p = %.3g). Adding to model.\n", inter_term, pvalue))
        simpler_model <- trial_model
        sig_inter <- trial_inter
      } else {
        cat(sprintf("Interaction term '%s' is not significant (p = %.3g). Stopping interaction search.\n", inter_term, pvalue))
        break
      }
    }
  } else {
    cat("No interaction candidates to test.\n")
  }
  
  
  # Search for significant Random Slopes
  if (add_slopes == TRUE) {
    cat("\n## 2d. Searching for Significant Random Slopes (Step Enabled)\n")
    
    # Candidates are all currently included main effect predictors
    slope_candidates <- included_predictors
    
    for (current_pred in slope_candidates) {
      trial_slopes <- c(included_slope, current_pred)
      
      random_part <- paste("(1 +", paste(trial_slopes, collapse = "+"), "|", id_col_name, ")")
      fixed_part <- paste(c(included_predictors, sig_inter), collapse = "+")
      
      trial_formula <- as.formula(paste(dep_col_name, "~", fixed_part, "+", random_part))
      
      message(sprintf("Testing random slope candidate: '%s'", current_pred))
      
      trial_model <- tryCatch({
        lme4::lmer(trial_formula, data = data, REML = FALSE)
      }, error = function(e) {
        message("Error fitting model with random slope ", current_pred, ". Skipping...")
        return(NULL)
      })
      
      if (is.null(trial_model) || check_model_conv(trial_model)) next
      
      # LRT for random effects using a mixture Chi-square test
      stat <- as.numeric(anova(simpler_model, trial_model)[2, "Chisq"])
      df_k <- length(trial_slopes) - length(included_slope)
      
      if (!is.finite(stat) || stat < 0) {
        message("Invalid Chi-squared statistic. Skipping.")
        next
      }
      
      pvalue <- 0.5 * (pchisq(stat, df = df_k, lower.tail = FALSE) +
                         pchisq(stat, df = df_k + 1, lower.tail = FALSE))
      
      if (is.finite(pvalue) && pvalue < alpha) {
        cat(sprintf("Random slope for '%s' is significant (mixture p = %.3g). Adding to model.\n", current_pred, pvalue))
        simpler_model <- trial_model
        included_slope <- trial_slopes
      } else {
        cat(sprintf("Random slope for '%s' is not significant (mixture p = %.3g). Stopping random slope search.\n", current_pred, pvalue))
        break
      }
    }
  } else {
    cat("\n--Random Slope Selection SKIPPED--\n")
  }
  
  optimal_model <- simpler_model
  cat("\n## 3. Final Linear Mixed Model (LMM) Summary\n")
  print(summary(optimal_model))
  
  # 4. Final MERF Fit with Optimal Structure
  cat("\n## 4. Final MERF Fit with Optimal Random Effects\n")
  
  # Reconstruct X for MERF (fixed effects only)
  fixed_effects <- c(included_predictors, sig_inter)
  
  if (rlang::is_empty(fixed_effects)) {
    merf_form <- as.formula("~ 0")
  } else {
    merf_form <- as.formula(paste("~ 0 + ", paste(fixed_effects, collapse = "+")))
  }
  
  X_final <- model.matrix(merf_form, data = data)
  
  # Reconstruct Z for MERF (random intercept + included random slopes)
  if (rlang::is_empty(included_slope)) {
    # Random intercept only
    Z_final <- as.matrix(rep(1, dim(data)[1]))
  } else {
    # Random intercept and significant random slopes
    Z_final <- as.matrix(cbind(1, data[, included_slope, drop = FALSE]))
  }
  
  mtry_final <- max(1, floor(ncol(X_final) / 3))
  
  final_merf <- LongituRF::MERF(
    X = X_final, Y = data[[dep_col_name]], id = data[[id_col_name]], Z = Z_final,
    iter = merf_iter, mtry = mtry_final, ntree = ntree, time = time, sto = "none", delta = 0.001
  )
  
  cat("\n-- Final MERF Model OOB Error --\n")
  print(final_merf$OOB[length(final_merf$OOB)])
  
  cat("\n-- Final MERF Variable Importance --\n")
  print(randomForest::importance(final_merf$forest))

  return(list(
    optimal_model = optimal_model,
    final_merf = final_merf,
    included_predictors = included_predictors,
    significant_interactions = sig_inter,
    included_slopes = included_slope
  ))
}