
######################## Get mixing weights with validation set

rr = function(x,digit=10) return(round(x,digit))

#' Evaluate prediction accuracy of a single PRS 
#'
#' This function get prediction accuracy for a single PRS
#'
#' @param data_df Data to assess prediction accuracy
#' @param prs_name PGS list of the trait
#' @param covar_list Array of covariates
#' @param isbinary TRUE if binary and FALSE otherwise (default = FALSE)
#' @param liabilityR2 TRUE if liability R2 should be reported (default = FALSE)
#' @param alpha Significance level to estimate power (default = 0.05)
#' @return A dataframe for prediction accuracy of a single PRS and their power with R2, R2 for output form, standard error, lower 95% CI, upper 95% CI, P-value and power
#' @export
eval_single_PRS = function(data_df, pheno = "trait", prs_name, covar_list, isbinary=F, liabilityR2=F, alpha=0.05, regression_output=F) {

	# if (var(data_df[,prs_name],na.rm=T)==0) {
	# 	writeLines(paste0("WARNING: Variance of ", prs_name, "=0. If this is a combined PRS, then covariates already explained the phenotype"))
	# 	res = data.frame(pgs=prs_name, R2=0, se=0, lowerCI=0, upperCI=0, pval=0, power=0)
	# 	return(res)
	# }

	colnames(data_df)[which(colnames(data_df)==pheno)] = "trait"

	if (isbinary & !liabilityR2) {
		data_df$trait = as.numeric(data_df$trait)
		formula = as.formula(paste0(pheno, " ~ scale(", prs_name, ") + ", paste0(covar_list, collapse="+")))
		model_full = glm(formula, data=data_df, family="binomial")
		r_full = suppressWarnings(logLik(model_full, REML=FALSE))[1]
		
		formula = as.formula(paste0(pheno, " ~ ", paste0(covar_list, collapse="+")))
		model_null = glm(formula, data=data_df, family="binomial")
		r_null = suppressWarnings(logLik(model_null, REML=FALSE))[1]
		
		m = r_full
		n = r_null
		N = nobs(model_full)
		cs = 1 - exp(-2/N * (m - n))
		nk = cs/(1 - exp(2/N * n))
		R2 = nk
		
	} else {
		data_df$trait = as.numeric(data_df$trait)
		
		formula = as.formula(paste0("trait ~ scale(", prs_name, ") + ", paste0(covar_list, collapse="+")))
		model_full = lm(formula, data=data_df)
		r_full = summary(model_full)$r.squared
		
		formula = as.formula(paste0("trait ~ ", paste0(covar_list, collapse="+")))
		model_null = lm(formula, data=data_df)
		r_null = summary(model_null)$r.squared
		
		N = nobs(model_full)
		R2 = r_full - r_null		
	}
	if (isbinary & liabilityR2) {
		N = nrow(data_df)
		K = mean(data_df$trait, na.rm=T)
		R2 = R2 * K * (1-K) / (dnorm(qnorm(p=1-K, lower.tail=T))^2)
	}
	
	
	NCP = N * R2 / (1-R2)
	power = 1-pnorm(qnorm(1-alpha/2)-NCP^0.5) + pnorm(qnorm(alpha/2)-NCP^0.5)
	
	vv = (4*R2*(1-R2)^2 *(N-2)^2) / ((N^2-1)*(N+3))
	
	se = sqrt(vv)
	lower_r2 = R2 - 1.97*se
	upper_r2 = R2 + 1.97*se
	pval = pchisq((R2/se)^2, df=1, lower.tail=F)
	
	r2_out = paste0(rr(R2,3), " (", rr(lower_r2,3), "-", rr(upper_r2,3), ")")

	if (regression_output) {
		return(data.frame(
			pgs=prs_name, 
			R2=R2, se=se, lowerCI=lower_r2, upperCI=upper_r2, pval=pval, power=power, 
			coef_regression=coef(summary(model_full))[2,1], 
			se_regression=coef(summary(model_full))[2,2], 
			pval_regression=coef(summary(model_full))[2,4]))
	} else {
		return(data.frame(pgs=prs_name, R2=R2, R2_out=r2_out, se=se, lowerCI=lower_r2, upperCI=upper_r2, pval=pval, power=power))
	}
}


#' Evaluate multiple PRSs from the given vector of PRS
#'
#' This function get prediction accuracy of the list of PRS
#'
#' @param data_df A dataframe where rows for samples, columns for PRS
#' @param pgs_list PGS list to evaluate, must be exist as columns in the dataframe
#' @param covar_list Array of covariates, must be exist as columns in the dataframe
#' @param liabilityR2 TRUE if liability R2 should be reported (default = FALSE)
#' @param alpha Significance level to estimate power (default = 0.05)
#' @param isbinary TRUE if binary and FALSE otherwise (default = FALSE)
#' @return A dataframe for prediction accuracy of multiple PRSss and their power with R2, R2 for output form, standard error, lower 95% CI, upper 95% CI, P-value and power
#'
#' @export
eval_multiple_PRS = function(data_df, pgs_list, covar_list, liabilityR2=F, alpha=0.05, isbinary=F, ncores=1, regression_output=F, pheno="trait") {

	colnames(data_df)[which(colnames(data_df)==pheno)] = "trait"	

	writeLines("Eval all PGS")
	if (isbinary) {
		writeLines("Case - control numbers:")
		print(table(data_df$trait))
	}

	idx = which(!pgs_list %in% colnames(data_df))
	if (length(idx)>0) {
		writeLines(paste0(length(idx), " scores not found or sum equal to 0 in data"))
		pgs_list = pgs_list[-idx]
	}

	# pred_acc_test = NULL
	tmp = mclapply(1:length(pgs_list), function(prs_i) {
	# for (prs_i in 1:length(pgs_list)) {
		
		if (prs_i %% 100 == 0) {
			writeLines(paste0("Evaluated ",  prs_i, " scores"))
		}
		      
		prs_name = pgs_list[prs_i]

		# print(prs_i)
		# print(prs_name)
		# print(var(data_df[,prs_name]))

		pred_acc_test_tmp = eval_single_PRS(data_df=data_df, pheno="trait", prs_name=prs_name, covar_list=covar_list, isbinary=isbinary, liabilityR2=liabilityR2, alpha=alpha, regression_output=regression_output)
		# pred_acc_test = rbind(pred_acc_test, pred_acc_test_tmp)

		return(pred_acc_test_tmp)

	}, mc.cores=ncores)
	
	pred_acc_test = do.call(rbind, tmp)
	pred_acc_test = pred_acc_test[order(pred_acc_test$R2, decreasing=T),]
	
	return(pred_acc_test)
}


######################## Get mixing weights with validation set

#' Evaluate prediction accuracy of a single PRS (no covariates; R2 = cor^2)
#'
#' @param data_df Data to assess prediction accuracy
#' @param prs_name Name of the PRS column
#' @param covar_list Ignored (kept for compatibility)
#' @param isbinary TRUE if binary phenotype
#' @param liabilityR2 TRUE to convert observed-scale R2 to liability scale (binary only)
#' @param alpha Significance level for power (default 0.05)
#' @param regression_output If TRUE, also returns regression coef/se/pval from a PRS-only model
#' @return Data.frame with pgs, R2 (cor^2), CI/se, pval (from cor test), and power
#' @export
eval_single_PRS_v3 = function(data_df, pheno = "trait", prs_name, covar_list,
                           isbinary = FALSE, liabilityR2 = FALSE, alpha = 0.05,
                           regression_output = FALSE) {

  # Standardize column name like upstream code expects
  colnames(data_df)[which(colnames(data_df) == pheno)] = "trait"

  # Pull vectors and ensure numeric
  x <- suppressWarnings(as.numeric(data_df[[prs_name]]))
  y <- suppressWarnings(as.numeric(data_df[["trait"]]))
  use <- stats::complete.cases(x, y)

  if (!any(use)) {
    out <- data.frame(pgs = prs_name, R2 = NA_real_, R2_out = NA_character_,
                      se = NA_real_, lowerCI = NA_real_, upperCI = NA_real_,
                      pval = NA_real_, power = NA_real_)
    if (regression_output) {
      out$coef_regression <- NA_real_
      out$se_regression   <- NA_real_
      out$pval_regression <- NA_real_
    }
    return(out)
  }

  x <- x[use]; y <- y[use]
  N <- length(x)

  # ---- R2 by correlation ----
  r  <- suppressWarnings(stats::cor(x, y))
  r  <- if (is.finite(r)) r else 0
  R2 <- r^2

  # Optional liability transform for binary phenotypes
  if (isbinary && liabilityR2) {
    K <- mean(y, na.rm = TRUE)             # prevalence in sample (expects 0/1 coding)
    # Convert observed-scale R2 to liability scale (Lee et al.)
    thr <- stats::qnorm(1 - K)
    R2 <- R2 * K * (1 - K) / (stats::dnorm(thr)^2)
  }

  # ---- p-value for correlation (two-sided) ----
  if (abs(r) < 1 && N > 2) {
    tval <- r * sqrt((N - 2) / (1 - r^2))
    pval <- 2 * stats::pt(-abs(tval), df = N - 2)
  } else {
    pval <- NA_real_
  }

  # ---- "power" definition kept as in PRSMix (from R2 via NCP) ----
  R2_clamped <- pmin(pmax(R2, 0), 1 - 1e-12)
  NCP <- N * R2_clamped / (1 - R2_clamped)
  z   <- stats::qnorm(1 - alpha/2)
  power <- 1 - (stats::pnorm(z - sqrt(NCP)) - stats::pnorm(-z - sqrt(NCP)))

  # ---- SE and CI for R2 (same formula as PRSMix) ----
  vv <- (4 * R2_clamped * (1 - R2_clamped)^2 * (N - 2)^2) /
        ((N^2 - 1) * (N + 3))
  se <- sqrt(vv)
  lower_r2 <- R2 - 1.97 * se
  upper_r2 <- R2 + 1.97 * se
  r2_out   <- paste0(rr(R2, 3), " (", rr(lower_r2, 3), "-", rr(upper_r2, 3), ")")

  out <- data.frame(
    pgs     = prs_name,
    R2      = R2,
    R2_out  = r2_out,
    se      = se,
    lowerCI = lower_r2,
    upperCI = upper_r2,
    pval    = pval,
    power   = power
  )

  if (regression_output) {
    # For completeness: PRS-only model on scaled PRS (matches prior API shape)
    z <- as.numeric(scale(x))
    if (isbinary) {
      fit <- stats::glm(y ~ z, family = "binomial")
    } else {
      fit <- stats::lm(y ~ z)
    }
    s <- summary(fit)
    out$coef_regression <- s$coefficients[2, 1]
    out$se_regression   <- s$coefficients[2, 2]
    out$pval_regression <- s$coefficients[2, 4]
  }

  return(out)
}


#' Evaluate multiple PRSs (no covariates; uses cor^2)
#'
#' @export
eval_multiple_PRS_v3 = function(data_df, pgs_list, covar_list, liabilityR2 = FALSE,
                             alpha = 0.05, isbinary = FALSE, ncores = 1,
                             regression_output = FALSE, pheno = "trait") {

  colnames(data_df)[which(colnames(data_df) == pheno)] = "trait"

  if (isbinary) {
    writeLines("Case - control numbers:")
    print(table(data_df$trait))
  }

  # Keep only PRS columns that exist and vary
  missing_idx <- which(!pgs_list %in% colnames(data_df))
  if (length(missing_idx) > 0) {
    writeLines(paste0(length(missing_idx), " scores not found in data; skipping"))
    pgs_list <- pgs_list[-missing_idx]
  }

  # Parallel eval (covariates ignored)
  tmp <- parallel::mclapply(seq_along(pgs_list), function(prs_i) {
    if (prs_i %% 100 == 0) writeLines(paste0("Evaluated ", prs_i, " scores"))
    prs_name <- pgs_list[prs_i]
    eval_single_PRS(
      data_df = data_df,
      pheno   = "trait",
      prs_name = prs_name,
      covar_list = NULL,                # ignored
      isbinary = isbinary,
      liabilityR2 = liabilityR2,
      alpha = alpha,
      regression_output = regression_output
    )
  }, mc.cores = ncores)

  pred_acc_test <- do.call(rbind, tmp)
  pred_acc_test <- pred_acc_test[order(pred_acc_test$R2, decreasing = TRUE), ]
  return(pred_acc_test)
}
