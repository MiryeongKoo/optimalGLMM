library(LongituRF)
library(randomForest)
library(lme4)

#' Performs MERF-based Forward Predictor Selection for Linear Mixed Models (LMM)
#'
#' This function integrates Mixed-Effects Random Forest (MERF) for predictor ranking
#' with a forward selection process to identify optimal linear mixed models (LMMs).
#' This alternates predictor and model selection steps.
#'
#' @param data A data frame containing all variables.
#' @param dep_col_name A string for the dependent variable column name (e.g., "dep").
#' @param id_col_name A string for the cluster ID column name (e.g., "schID").
#' @param L1_cols A character vector of Level 1 predictor column names.
#' @param L2_cols A character vector of Level 2 predictor column names.
#' @param dep_is_continuous Logical. TRUE (default) if the dependent variable is continuous.
#' @param add_slopes Logical. TRUE (default) to include the random slope search after each significant L1 predictor.
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
#' optimal_model <- optimalLMM2(
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

optimalLMM2 <- function(
    data, dep_col_name, id_col_name, L1_cols, L2_cols,
    dep_is_continuous = TRUE, add_slopes = TRUE,
    alpha = 0.05, tol = 1e-4, merf_iter = 100, ntree = 1000
) {
  
  # 1. Setup and MERF Ranking
  cat("## 1. MERF Ranking for Predictor Importance\n")
  
  # Renaming for convenience inside the function
  dep_var <- data[[dep_col_name]]
  id_var_name <- id_col_name
  
  # Determine importance type
  importance_type <- ifelse(dep_is_continuous, 1, 2)
  
  # Prepare MERF input
  X_merf_cols <- c(L1_cols, L2_cols)
  X_merf <- data[, X_merf_cols, drop = FALSE]
  Y_merf <- data[[dep_col_name]]
  Z <- model.matrix(~ data[[id_col_name]] - 1)
  time <- rep(1, dim(data)[1])
  mtry_merf <- max(1, floor(ncol(X_merf) / 3))
  
  # Fit the initial MERF model (using default random intercept structure)
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
  
  # Get sorted predictor names (replacing imp.col, imp.col_L1, imp.col_L2)
  imp_col_full <- full_sorted$Predictor
  imp_col_L1 <- imp_col_full[imp_col_full %in% L1_cols]
  imp_col_L2 <- imp_col_full[imp_col_full %in% L2_cols]
  
  cat("MERF Ranking Complete. Starting LMEM Forward Selection.\n")
  cat("-- Top 3 Predictors for Interaction Search:", paste(imp_col_full[1:3], collapse = ", "), "\n")
  
  make_formula <- function(fixed_terms, random_slopes) {
    fixed_part <- paste(fixed_terms, collapse = "+")
    
    if (is.null(random_slopes) || length(random_slopes) == 0) {
      random_part <- paste("(1 |", id_var_name, ")")
    } else {
      random_part <- paste("(1 +", paste(random_slopes, collapse = "+"), "|", id_var_name, ")")
    }
    
    if (nchar(fixed_part) == 0) fixed_part <- "1"
    
    return(as.formula(paste(dep_col_name, "~", fixed_part, "+", random_part)))
  }
  
  check_model_fit <- function(formula, data) {
    conv_warn <- FALSE
    trial_model <- tryCatch({
      withCallingHandlers(
        lmer(formula, data = data, REML = FALSE),
        warning = function(w) {
          cc <- conditionCall(w); if (is.null(cc)) return()
          fun <- if (is.call(cc[[1]]) && as.character(cc[[1]][[1]]) %in% c("::",":::")) as.character(cc[[1]][[3]]) else as.character(cc[[1]])
          is_lmerish <- (fun %in% c("lmer","glmer","checkConv")) || grepl("\\b(checkConv|lme4)\\b", paste(deparse(cc), collapse = " "), ignore.case = TRUE)
          if (is_lmerish && grepl("failed to converge|max\\|grad\\||singular|nearly unidentifiable|Hessian|positive-definite", conditionMessage(w), ignore.case = TRUE)) {
            conv_warn <<- TRUE
            invokeRestart("muffleWarning")
          }
        },
        message = function(m) {
          cc <- conditionCall(m); if (is.null(cc)) return()
          fun <- if (is.call(cc[[1]]) && as.character(cc[[1]][[1]]) %in% c("::",":::")) as.character(cc[[1]][[3]]) else as.character(cc[[1]])
          is_lmerish <- (fun %in% c("lmer","glmer","checkConv")) || grepl("\\b(checkConv|lme4)\\b", paste(deparse(cc), collapse = " "), ignore.case = TRUE)
          if (is_lmerish && grepl("convergence code|failed to converge|max\\|grad\\||Hessian|positive-definite", conditionMessage(m), ignore.case = TRUE)) {
            conv_warn <<- TRUE
            invokeRestart("muffleMessage")
          }
        }
      )
    }, error = function(e) {
      message("Error fitting model: ", conditionMessage(e))
      return(NULL)
    })
    
    if (is.null(trial_model)) return(list(model = NULL, skip = TRUE))
    if (isTRUE(conv_warn) || lme4::isSingular(trial_model, tol = tol)) {
      message("Model skipped due to convergence warning or singularity.")
      return(list(model = NULL, skip = TRUE))
    }
    return(list(model = trial_model, skip = FALSE))
  }
  
  # 2. Forward Predictor Selection
  
  # Fit the null (a random intercept) model
  null_model_formula <- make_formula(fixed_terms = character(0), random_slopes = character(0))
  null_model <- lmer(null_model_formula, data = data, REML = FALSE)
  
  cat("## 2. Forward Selection\n")
  print(summary(null_model))
  
  simpler_model <- null_model
  included_predictors <- character(0) # Fixed effects (main + interactions)
  included_slope <- character(0)      # Random slopes
  sig_inter <- character(0)           # Only interaction terms
  

  # Search for significant Level-1 predictors
  cat("\n## Searching for Significant Level-1 Predictors\n")
  
  for (q in 1:length(imp_col_L1)) {
    p_name <- imp_col_L1[q]
    trial_predictors_fixed <- c(included_predictors, p_name)
    
    # 1. Test adding new L1 predictor as FIXED EFFECT
    trial_formula_fixed <- make_formula(fixed_terms = trial_predictors_fixed, random_slopes = included_slope)
    
    fit_result <- check_model_fit(trial_formula_fixed, data)
    if (fit_result$skip) next
    trial_model_fixed <- fit_result$model
    
    # Calculate LRT p-value
    pvalue <- tryCatch(
      as.numeric(anova(simpler_model, trial_model_fixed)[2,"Pr(>Chisq)"]),
      error = function(e) { message("Error running ANOVA. Skipping..."); return(NaN) }
    )
    
    if (is.finite(pvalue) && pvalue < alpha) {
      print(paste(p_name, " is significant. Testing for random slope..."))
      
      # Random Slope Test
      if (add_slopes == TRUE) {
        # Test adding the significant predictor p_name as a random slope
        trial_slopes <- c(included_slope, p_name) 
        
        trial_formula_rslope <- make_formula(fixed_terms = trial_predictors_fixed, random_slopes = trial_slopes)
        
        fit_result_rslope <- check_model_fit(trial_formula_rslope, data)
        
        if (!fit_result_rslope$skip) {
          trial_model_rslope <- fit_result_rslope$model
          
          # Mixture chi-square p-value comparing trial_model_fixed vs trial_model_rslope
          stat <- as.numeric(anova(trial_model_fixed, trial_model_rslope)[2, "Chisq"])
          df_k <- length(trial_slopes) - length(included_slope) 
          
          if (is.finite(stat) && stat >= 0) {
            pvalue_rs <- 0.5 * (pchisq(stat, df = df_k, lower.tail = FALSE) + 
                                  pchisq(stat, df = df_k + 1, lower.tail = FALSE))
          } else {
            pvalue_rs <- NaN
          }
          
          if (is.finite(pvalue_rs) && pvalue_rs < alpha) {
            message(sprintf("The random slope for '%s' is significant (mixture p = %.3g).", p_name, pvalue_rs))
            simpler_model <- trial_model_rslope # Update to model with new slope
            included_slope <- trial_slopes
          } else {
            message(sprintf("Random slope for '%s' not significant (p = %.3g). Keeping fixed effect only.", p_name, pvalue_rs))
            simpler_model <- trial_model_fixed # Update to model with fixed effect only
          }
        } else {
          message(sprintf("Random slope model for '%s' failed to fit or was singular. Keeping fixed effect only.", p_name))
          simpler_model <- trial_model_fixed # Update to model with fixed effect only
        }
      } else {
        # Skip random slope test if add_slopes == FALSE
        message(sprintf("Random slope test for '%s' skipped (add_slopes = FALSE).", p_name))
        simpler_model <- trial_model_fixed # Update to model with fixed effect only
      }
      # Final update of fixed predictors after all checks
      included_predictors <- trial_predictors_fixed
      
    } else {
      print(paste(p_name, " is not significant. Stop L1 search."))
      break
    }
  }
  
  # Search for significant Level-2 predictors
  cat("\n## Searching for Significant Level-2 Predictors\n")
  
  for (r in 1:length(imp_col_L2)) {
    p_name <- imp_col_L2[r]
    trial_predictors_fixed <- c(included_predictors, p_name)
    
    trial_formula <- make_formula(fixed_terms = trial_predictors_fixed, random_slopes = included_slope)
    
    fit_result <- check_model_fit(trial_formula, data)
    if (fit_result$skip) next
    trial_model <- fit_result$model
    
    # Calculate LRT p-value
    pvalue <- tryCatch(
      as.numeric(anova(simpler_model, trial_model)[2,"Pr(>Chisq)"]),
      error = function(e) { message("Error running ANOVA. Skipping..."); return(NaN) }
    )
    
    if (is.finite(pvalue) && pvalue < alpha) {
      print(paste(p_name, " is significant. Continue searching for significant L2 predictors."))
      simpler_model <- trial_model
      included_predictors <- trial_predictors_fixed
    } else {
      message(sprintf("L2 predictor '%s' not significant (p = %.3g). Stop L2 search.", p_name, pvalue))
      break
    }
  }
  
  # Search for significant Interaction terms
  cat("\n## Searching for Significant Interaction Terms\n")
  
  if (length(imp_col_full) >= 3) {
    top_3 <- imp_col_full[1:3]
    interaction_candidates <- c(
      paste(top_3[1], top_3[2], sep = ":"),
      paste(top_3[1], top_3[3], sep = ":"),
      paste(top_3[2], top_3[3], sep = ":")
    )
    cat("Testing interactions:", paste(interaction_candidates, collapse = ", "), "\n")
  } else {
    interaction_candidates <- character(0)
    cat("Fewer than 3 total predictors. Skipping interaction search.\n")
  }
  
  if (length(interaction_candidates) > 0) {
    for (k in 1:length(interaction_candidates)) {
      inter_term <- interaction_candidates[k]
      trial_interaction <- c(sig_inter, inter_term)
      
      # Fixed predictors include main effects and current interactions
      trial_predictors_fixed <- c(included_predictors, trial_interaction)
      
      trial_formula <- make_formula(fixed_terms = trial_predictors_fixed, random_slopes = included_slope)
      
      fit_result <- check_model_fit(trial_formula, data)
      if (fit_result$skip) next
      trial_model <- fit_result$model
      
      # Calculate LRT p-value
      pvalue <- tryCatch(
        as.numeric(anova(simpler_model, trial_model)[2,"Pr(>Chisq)"]),
        error = function(e) { message("Error running ANOVA. Skipping..."); return(NaN) }
      )
      
      if (is.finite(pvalue) && pvalue < alpha) {
        print(paste(inter_term, " is significant. Continue searching for interaction terms."))
        simpler_model <- trial_model
        included_predictors <- trial_predictors_fixed
        sig_inter <- trial_interaction
      } else {
        message(sprintf("Interaction term '%s' not significant (p = %.3g). Stop interaction search.", inter_term, pvalue))
        break
      }
    }
  }
  
  # 3. Final LMM Summary
  optimal_model <- simpler_model
  cat("\n## 3. Final Optimal Linear Mixed Model (LMM) Summary\n")
  print(summary(optimal_model))
  
  # 4. Final MERF Fit with Optimal Structure
  cat("\n## 4. Final MERF Fit with Optimal Random Effects\n")
  
  # Reconstruct X for MERF (fixed effects only)
  fixed_effects_final <- c(included_predictors)
  fixed_formula_terms <- c(imp_col_full[imp_col_full %in% fixed_effects_final], sig_inter)
  
  if(rlang::is_empty(fixed_formula_terms)) {
    merf_form <- as.formula("~ 0")
  } else {
    merf_form <- as.formula(paste("~ 0 + ", paste(fixed_formula_terms, collapse = "+")))
  }
  
  X_final <- model.matrix(merf_form, data = data)
  
  # Reconstruct Z for MERF (random intercept + included random slopes)
  G <- model.matrix(~ 0 + factor(data[[id_var_name]]), data = data)
  
  if (length(included_slope) == 0) {
    # Random intercept only
    Z_final <- G
  } else {
    # Random intercept and significant random slopes
    Z_final <- cbind(G, do.call(cbind, lapply(included_slope, function(s) G * data[[s]])))
  }
  
  time <- rep(1, dim(data)[1])
  mtry_final <- max(1, floor(ncol(X_final) / 3))
  
  final_merf <- LongituRF::MERF(
    X = X_final, Y = data[[dep_var]], id = data[[id_var_name]], Z = Z_final,
    iter = merf_iter, mtry = mtry_final, ntree = ntree, time = time, sto = "none", delta = 0.001
  )
  
  cat("\n-- Final MERF Model OOB Error --\n")
  print(final_merf$OOB[length(final_merf$OOB)])
  
  cat("\n-- Final MERF Variable Importance --\n")
  print(randomForest::importance(final_merf$forest))
  
  return(list(
    optimal_model = optimal_model,
    final_merf = final_merf,
    included_predictors = fixed_effects_final[fixed_effects_final %in% c(L1_cols, L2_cols)], 
    significant_interactions = sig_inter,
    included_slopes = included_slope
  ))
}