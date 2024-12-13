#' Compute Propensity Score
#'
#' This function computes the propensity score for each individual in a dataset 
#' based on a logistic regression model. It calculates the probability of receiving 
#' the treatment based on the specified covariates.
#'
#' @param data A data frame containing the variables.
#' @param treatment A character string specifying the treatment variable (binary: 0 or 1).
#' @param covariates A character vector of covariate names used in the logistic regression.
#'
#' @return A numeric vector of propensity scores.
#' @export
#' @name compute_propensity_score
#' @import stats
#' @import microbenchmark

compute_propensity_score <- function(data, treatment, covariates) {
  # 构建公式
  formula <- as.formula(paste(treatment, "~", paste(covariates, collapse = "+")))
  # 逻辑回归拟合
  model <- glm(formula, data = data, family = binomial)
  # 返回倾向性得分
  return(predict(model, type = "response"))
}

#' Nearest Neighbor Matching
#'
#' This function performs nearest neighbor matching based on propensity scores, 
#' matching treated units to control units. Optionally, a caliper can be specified to 
#' limit the maximum allowable distance for matching.
#'
#' @param data A data frame containing the treatment variable and covariates.
#' @param treatment A character string specifying the treatment variable (binary: 0 or 1).
#' @param covariates A character vector of covariate names used for propensity score calculation.
#' @param caliper A numeric value specifying the maximum allowable distance for matching. Default is 0.05.
#' @param replace A logical value indicating whether to allow replacement of matched controls. Default is FALSE.
#'
#' @return A data frame with treated and matched control units.
#' @export
#' @name nearest_neighbor_matching
#' @import stats
nearest_neighbor_matching <- function(data, treatment, covariates, caliper = 0.05, replace = FALSE) {
  prop_scores <- compute_propensity_score(data, treatment, covariates)
  
  treated <- data[data[[treatment]] == 1, ]
  control <- data[data[[treatment]] == 0, ]
  
  matched_control_indices <- numeric(nrow(treated))
  
  for (i in 1:nrow(treated)) {
    dist <- abs(prop_scores[data[[treatment]] == 0] - prop_scores[data[[treatment]] == 1][i])
    if (!is.null(caliper)) {
      dist <- dist[dist <= caliper]
    }
    if (length(dist) > 0) {
      matched_control_indices[i] <- which.min(dist)
    } else {
      matched_control_indices[i] <- NA
    }
    if (!replace) {
      prop_scores[data[[treatment]] == 0][matched_control_indices[i]] <- NA
    }
  }
  
  matched_control <- control[matched_control_indices, ]
  matched_data <- rbind(treated, matched_control)
  
  return(matched_data)
}

#' Full Matching using MatchIt
#'
#' This function uses `matchit` package to perform full matching based on propensity scores.
#'
#' @param data A data frame containing the treatment variable and covariates.
#' @param treatment A character string specifying the treatment variable (binary: 0 or 1).
#' @param covariates A character vector of covariate names used for propensity score calculation.
#'
#' @return A data frame with matched data from `matchit` using full matching.
#' @export
#' @import MatchIt
#' @import stats
#' @name full_matchit_matching
NULL
full_matchit_matching <- function(data, treatment, covariates) {
  formula <- as.formula(paste(treatment, "~", paste(covariates, collapse = "+")))
  m.out <- matchit(formula, data = data, method = "full", ratio = 1, caliper = 0.05)
  matched_data <- match.data(m.out)
  return(matched_data)
}

#' Optimal Matching using MatchIt
#'
#' This function uses `matchit` package to perform optimal matching based on propensity scores.
#'
#' @param data A data frame containing the treatment variable and covariates.
#' @param treatment A character string specifying the treatment variable (binary: 0 or 1).
#' @param covariates A character vector of covariate names used for propensity score calculation.
#'
#' @return A data frame with matched data from `matchit` using optimal matching.
#' @export
#' @name optimal_matchit_matching
#' @import stats
optimal_matchit_matching <- function(data, treatment, covariates) {
  formula <- as.formula(paste(treatment, "~", paste(covariates, collapse = "+")))
  m.out <- matchit(formula, data = data, method = "optimal", ratio = 1, caliper = 0.05)
  matched_data <- match.data(m.out)
  return(matched_data)
}

#' Standardized Mean Difference
#'
#' This function computes the standardized mean difference (SMD) between the treatment 
#' and control groups for each covariate before and after matching.
#'
#' @param data A data frame containing the treatment variable and covariates.
#' @param treatment A character string specifying the treatment variable (binary: 0 or 1).
#' @param covariates A character vector of covariate names used to calculate the SMD.
#'
#' @return A named numeric vector of standardized mean differences for each covariate.
#' @export
#' @name standardized_mean_diff
standardized_mean_diff <- function(data, treatment, covariates) {
  treated <- data[data[[treatment]] == 1, ]
  control <- data[data[[treatment]] == 0, ]
  
  treated_means <- colMeans(treated[, covariates], na.rm = TRUE)
  control_means <- colMeans(control[, covariates], na.rm = TRUE)
  
  treated_sds <- apply(treated[, covariates], 2, sd, na.rm = TRUE)
  control_sds <- apply(control[, covariates], 2, sd, na.rm = TRUE)
  
  standardized_diff <- (treated_means - control_means) / sqrt((treated_sds^2 + control_sds^2) / 2)
  
  return(standardized_diff)
}

#' Compare Matching Methods
#'
#' This function compares the performance of custom nearest neighbor matching, optimal matching, 
#' full matching, and matching using the `matchit` package. It compares the matching time and 
#' standardized mean differences before and after matching.
#'
#' @param data A data frame containing the treatment variable and covariates.
#' @param treatment A character string specifying the treatment variable (binary: 0 or 1).
#' @param covariates A character vector of covariate names used for propensity score calculation.
#'
#' @return A summary of the matching performance including time and standardized mean differences.
#' @export
#' @name compare_matching
compare_matching <- function(data, treatment, covariates) {
  # 1. 测试自定义匹配函数
  start_time <- Sys.time()
  custom_matched_data <- nearest_neighbor_matching(data, treatment, covariates)
  custom_time <- Sys.time() - start_time
  
  # 2. 测试 Matchit 的 Full Matching
  start_time <- Sys.time()
  matchit_full_matched_data <- full_matchit_matching(data, treatment, covariates)
  matchit_full_time <- Sys.time() - start_time
  
  # 3. 测试 Matchit 的 Optimal Matching
  start_time <- Sys.time()
  matchit_optimal_matched_data <- optimal_matchit_matching(data, treatment, covariates)
  matchit_optimal_time <- Sys.time() - start_time
  
  # 4. 计算匹配前后的标准化均值差异
  pre_matching_diff <- standardized_mean_diff(data, treatment, covariates)
  post_custom_matching_diff <- standardized_mean_diff(custom_matched_data, treatment, covariates)
  post_matchit_full_matching_diff <- standardized_mean_diff(matchit_full_matched_data, treatment, covariates)
  post_matchit_optimal_matching_diff <- standardized_mean_diff(matchit_optimal_matched_data, treatment, covariates)
  
  # 5. 输出比较结果
  cat("Custom Matching Time: ", custom_time, "\n")
  cat("Matchit Full Matching Time: ", matchit_full_time, "\n")
  cat("Matchit Optimal Matching Time: ", matchit_optimal_time, "\n")
  
  cat("\nStandardized Mean Differences (Before Matching):\n")
  print(pre_matching_diff)
  
  cat("\nStandardized Mean Differences (After Custom Matching):\n")
  print(post_custom_matching_diff)
  
  cat("\nStandardized Mean Differences (After Matchit Full Matching):\n")
  print(post_matchit_full_matching_diff)
  
  cat("\nStandardized Mean Differences (After Matchit Optimal Matching):\n")
  print(post_matchit_optimal_matching_diff)
}


