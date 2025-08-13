#' Perform linear combination of the polygenic risk scores (no covariates; R2 by correlation)
#'
#' This function performs a linear combination of scores. It no longer requires covariates.
#'
#' @param pheno_file Directory to the phenotype file
#' @param covariate_file (OPTIONAL) Path to covariate file; if NULL, covariates are ignored (DEFAULT = NULL)
#' @param score_files_list A vector of paths to the PGS files to read
#' @param trait_specific_score_file A filename containing PGS IDs of trait-specific scores to combine (one per line)
#' @param pheno_name Column name of the phenotype in pheno_file
#' @param isbinary TRUE if binary trait
#' @param out Prefix of output
#' @param metascore Meta-information from PGS Catalog (unused here)
#' @param liabilityR2 Ignored (kept for compatibility)
#' @param IID_pheno Column name of IID in phenotype file
#' @param covar_list Ignored (kept for compatibility; DEFAULT = NULL)
#' @param cat_covar_list Ignored (kept for compatibility; DEFAULT = NULL)
#' @param ncores Number of CPU cores for parallel processing (DEFAULT = 1)
#' @param is_extract_adjSNPeff TRUE to extract adjSNPeff (unchanged)
#' @param original_beta_files_list Vector of paths to SNP effect sizes (unchanged)
#' @param train_ids_file A file containing the IIDs for training split
#' @param power_thres_list Vector of power thresholds to select scores (DEFAULT = 0.95); here power == R2
#' @param pval_thres_list Vector of P-value thresholds to select scores (DEFAULT = 0.05)
#' @param read_pred_training Read precomputed training results (unchanged hooks)
#' @param read_pred_testing Read precomputed testing results (unchanged hooks)
#' @return Writes the same outputs as before except any "NULL model" AUC files; returns 0 on success
#'
#' @importFrom stats cor qnorm pt as.formula predict sd
#' @importFrom data.table fread fwrite
#' @importFrom dplyr bind_rows select all_of mutate group_by summarise rowwise filter
#' @importFrom magrittr %>%
#' @importFrom parallel makePSOCKcluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom caret train trainControl
#' @importFrom utils head read.table
#' @export
combine_PRS_v2 = function(
  pheno_file,
  covariate_file = NULL,
  score_files_list,
  trait_specific_score_file,
  pheno_name,
  isbinary,
  out,
  allPGS_list = NULL,
  metascore = NULL,
  liabilityR2 = FALSE,
  IID_pheno = "IID",
  covar_list = NULL,
  cat_covar_list = NULL,
  ncores = 1,
  is_extract_adjSNPeff = FALSE,
  original_beta_files_list = NULL,
  train_ids_file = NULL,
  training_result_file = NULL,
  power_thres_list = c(0.95),
  pval_thres_list = c(0.05),
  nfold_cv = 3,
  read_pred_training = FALSE,
  read_pred_testing = FALSE,
  debug = FALSE
) {

  options(datatable.fread.datatable = FALSE)

  # ------------------------- helpers (new) -------------------------
  irnt <- function(x) return(qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
  rr   <- function(x, d = 3) round(x, d)

  # Replace your eval_single_PRS_nocov with this version
    # --- replace this helper ---
    eval_single_PRS_nocov <- function(df,
                                    pheno = "trait",
                                    prs_name,
                                    isbinary = FALSE,
                                    liabilityR2 = FALSE,
                                    alpha = 0.05,
                                    regression_output = FALSE) {
    # match original: rename column to "trait"
    colnames(df)[which(colnames(df) == pheno)] <- "trait"

    # vectors
    x <- suppressWarnings(as.numeric(df[[prs_name]]))
    y <- suppressWarnings(as.numeric(df[["trait"]]))
    keep <- stats::complete.cases(x, y)
    x <- x[keep]; y <- y[keep]
    N <- length(x)

    if (!length(x)) {
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

    # R2 from correlation
    r  <- suppressWarnings(stats::cor(x, y))
    r  <- if (is.finite(r)) r else 0
    R2 <- r^2

    # liability-scale (same place as original)
    if (isbinary && liabilityR2) {
        K <- mean(y, na.rm = TRUE)  # assumes 0/1 coding
        thr <- stats::qnorm(1 - K)
        R2 <- R2 * K * (1 - K) / (stats::dnorm(thr)^2)
    }

    # power (same formula as PRSMix)
    R2c <- pmin(pmax(R2, 0), 1 - 1e-12)
    NCP <- N * R2c / (1 - R2c)
    z   <- stats::qnorm(1 - alpha/2)
    power <- 1 - (stats::pnorm(z - sqrt(NCP)) - stats::pnorm(-z - sqrt(NCP)))

    # SE/CI/pval (same as PRSMix)
    vv <- (4 * R2c * (1 - R2c)^2 * (N - 2)^2) / ((N^2 - 1) * (N + 3))
    se <- sqrt(vv)
    lower_r2 <- R2 - 1.97 * se
    upper_r2 <- R2 + 1.97 * se
    pval <- stats::pchisq((R2 / se)^2, df = 1, lower.tail = FALSE)
    r2_out <- paste0(rr(R2, 3), " (", rr(lower_r2, 3), "-", rr(upper_r2, 3), ")")

    out <- data.frame(pgs = prs_name, R2 = R2, R2_out = r2_out,
                        se = se, lowerCI = lower_r2, upperCI = upper_r2,
                        pval = pval, power = power)

    if (regression_output) {
        # optional: PRS-only model (no covariates), matching original return shape
        zprs <- as.numeric(scale(x))
        fit <- if (isbinary) stats::glm(y ~ zprs, family = "binomial") else stats::lm(y ~ zprs)
        s <- summary(fit)
        out$coef_regression <- s$coefficients[2, 1]
        out$se_regression   <- s$coefficients[2, 2]
        out$pval_regression <- s$coefficients[2, 4]
    }

    out
    }

    # --- replace this helper ---
    eval_multiple_PRS_nocov <- function(df,
                                        pgs_vec,
                                        isbinary = FALSE,
                                        ncores = 1L,
                                        liabilityR2 = FALSE,
                                        alpha = 0.05,
                                        regression_output = FALSE,
                                        pheno = "trait") {
    # match original: rename column to "trait"
    colnames(df)[which(colnames(df) == pheno)] <- "trait"

    if (isbinary) {
        writeLines("Case - control numbers:")
        print(table(df$trait))
    }

    # drop missing PRS columns
    missing_idx <- which(!pgs_vec %in% colnames(df))
    if (length(missing_idx) > 0) {
        writeLines(paste0(length(missing_idx), " scores not found in data; skipping"))
        pgs_vec <- pgs_vec[-missing_idx]
    }

    # optional: drop zero-variance PRS to avoid NA cor
    if (length(pgs_vec)) {
        v <- sapply(pgs_vec, function(p) stats::var(df[[p]], na.rm = TRUE))
        pgs_vec <- pgs_vec[v > 0]
    }

    # parallel map (no covariates)
    res_list <- parallel::mclapply(seq_along(pgs_vec), function(i) {
        if (i %% 100 == 0) writeLines(paste0("Evaluated ", i, " scores"))
        eval_single_PRS_nocov(
        df, pheno = "trait",
        prs_name = pgs_vec[i],
        isbinary = isbinary,
        liabilityR2 = liabilityR2,
        alpha = alpha,
        regression_output = regression_output
        )
    }, mc.cores = ncores)

    out <- do.call(rbind, res_list)
    out <- out[order(out$R2, decreasing = TRUE), ]
    out
    }


  # ---------------------- read scores & pheno ----------------------
  writeLines("--- Reading all polygenic risk scores ---")
  writeLines("Ensure score column names end with _SUM")
  all_scores <- NULL
  for (score_file_i in seq_along(score_files_list)) {
    score_file <- score_files_list[score_file_i]
    dd <- data.table::fread(score_file)
    idx <- which(endsWith(colnames(dd), "_SUM") & colnames(dd) != "NAMED_ALLELE_DOSAGE_SUM")
    idx2 <- which(colnames(dd) == "IID")
    dd_sub <- dd[, c(idx2, idx)]
    if (is.null(all_scores)) {
      all_scores <- dd_sub
    } else {
      all_scores <- merge(all_scores, dd_sub, by = "IID")
    }
  }
  colnames(all_scores)[2:ncol(all_scores)] <-
    substring(colnames(all_scores)[2:ncol(all_scores)], 1, nchar(colnames(all_scores)[2:ncol(all_scores)]) - 4)

  if (!is.null(allPGS_list)) {
    pgs_extract <- intersect(colnames(all_scores)[2:ncol(all_scores)], allPGS_list)
    all_scores <- all_scores[, c(1, match(pgs_extract, colnames(all_scores)))]
  }
  score_names <- colnames(all_scores)[2:ncol(all_scores)]

  # Which PRS are trait-specific?
  pgs_list <- NULL
  for (ff_i in seq_along(trait_specific_score_file)) {
    writeLines(paste0("Reading: ", trait_specific_score_file[ff_i]))
    pgs_list_tmp <- data.table::fread(trait_specific_score_file[ff_i], header = FALSE)[, 1]
    pgs_list <- c(pgs_list, pgs_list_tmp)
  }
  pgs_list <- intersect(pgs_list, colnames(all_scores))

  # Phenotype
  pheno <- data.table::fread(pheno_file)
  idx <- which(colnames(pheno) %in% c(IID_pheno, pheno_name))
  pheno <- pheno[, idx]
  colnames(pheno) <- c("IID", "trait")

  writeLines("--- Merging Phenotype and PRS files ---")
  pheno_prs <- merge(pheno, all_scores, by = "IID")
  pheno_prs <- pheno_prs[!is.na(pheno_prs$trait), ]

  # (Covariates intentionally ignored)
  if (!is.null(covariate_file)) {
    writeLines("NOTE: covariate_file provided but will be ignored in modeling/metrics per your request.")
  }

  # --------------------------- split ---------------------------
  if (!is.null(train_ids_file)) {
    train_iids <- data.table::fread(train_ids_file)
    train_idx  <- unique(stats::na.omit(match(train_iids$IID, pheno_prs$IID)))
    if (!length(train_idx)) stop("No training IIDs matched merged IID")
  } else {
    writeLines("Not using custom data split!")
    set.seed(1)
    train_idx <- sample(seq_len(nrow(pheno_prs)), size = floor(0.7 * nrow(pheno_prs)))
  }
  if (isbinary) data.table::fwrite(as.data.frame(table(pheno_prs$trait)),
                                   paste0(out, "_case_counts.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

  train_df <- pheno_prs[train_idx, ]
  test_df  <- pheno_prs[-train_idx, ]

  if (!isbinary) {
    train_df$trait <- irnt(train_df$trait)
    test_df$trait  <- irnt(test_df$trait)
  }

  data.table::fwrite(train_df[, c("IID", "trait")], paste0(out, "_train_df.txt"),
                     row.names = FALSE, quote = FALSE, sep = "\t")
  data.table::fwrite(test_df[, c("IID", "trait")], paste0(out, "_test_df.txt"),
                     row.names = FALSE, quote = FALSE, sep = "\t")

  # -------------------- training evaluation --------------------
  writeLines("--- Evaluating PRS in training set (R2 = cor^2) ---")
  if (!read_pred_training || is.null(training_result_file)) {
    training_file <- paste0(out, "_train_allPRS.txt")
    read_pred_training_1 <- (read_pred_training && file.exists(training_file))
  } else {
    training_file <- training_result_file
    read_pred_training_1 <- read_pred_training && all(file.exists(training_file))
    if (!read_pred_training_1) {
      writeLines("Declared:")
      writeLines(paste0("read_pred_training = ", read_pred_training))
      writeLines(paste0("training_result_file = ", paste(training_result_file, collapse = ";")))
      stop("Missing at least one training result file; set read_pred_training=FALSE or fix paths.")
    }
    writeLines("Reading all files in training_result_file")
  }

  if (!read_pred_training_1) {
    sumscore <- apply(as.data.frame(train_df[, score_names, drop = FALSE]), 2, stats::var)
    idx0 <- which(sumscore == 0)
    if (length(idx0) > 0) train_df <- train_df[, -match(names(idx0), colnames(train_df))]
    pgs_list_all <- intersect(colnames(train_df), colnames(all_scores)[2:ncol(all_scores)])
    pred_acc_train_allPGS_summary <- eval_multiple_PRS_nocov(train_df, pgs_list_all, isbinary, ncores)
    data.table::fwrite(pred_acc_train_allPGS_summary, training_file, row.names = FALSE, sep = "\t", quote = FALSE)
  } else {
    pred_acc_train_allPGS_summary <- NULL
    for (file_i in seq_along(training_file)) {
      writeLines(paste0("Reading training file: ", training_file[file_i]))
      tmp <- data.table::fread(training_file[file_i])
      pred_acc_train_allPGS_summary <- rbind(pred_acc_train_allPGS_summary, tmp)
    }
    pred_acc_train_allPGS_summary <- pred_acc_train_allPGS_summary[!duplicated(pred_acc_train_allPGS_summary$pgs), ]
  }
  pred_acc_train_allPGS_summary <- as.data.frame(pred_acc_train_allPGS_summary)

  pred_acc_train_trait_summary <- pred_acc_train_allPGS_summary
  pred_acc_train_trait_summary <- pred_acc_train_trait_summary[order(as.numeric(pred_acc_train_trait_summary$pval), decreasing = FALSE), ]

  pred_acc_train_allPGS_summary1 <- dplyr::filter(pred_acc_train_allPGS_summary, pgs %in% pgs_list)
  pred_acc_train_allPGS_summary1 <- pred_acc_train_allPGS_summary1[order(as.numeric(pred_acc_train_allPGS_summary1$R2), decreasing = TRUE), ]
  bestPRS <- pred_acc_train_allPGS_summary1$pgs[1]
  writeLines(paste0("The best single trait-specific score in the training set is ", bestPRS))

  bestPRS_acc <- eval_single_PRS_nocov(test_df, pheno = "trait", prs_name = bestPRS, isbinary = isbinary)
  data.table::fwrite(bestPRS_acc, paste0(out, "_best_acc.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

  # -------------------- testing evaluation ---------------------
  writeLines("--- Evaluating PRS in testing set (R2 = cor^2) ---")
  testing_file <- paste0(out, "_test_allPRS.txt")
  read_pred_testing_1 <- (read_pred_testing && file.exists(testing_file))

  if (!read_pred_testing_1) {
    sumscore <- apply(as.data.frame(test_df[, score_names, drop = FALSE]), 2, stats::var)
    idx0 <- which(sumscore == 0)
    if (length(idx0) > 0) test_df <- test_df[, -match(names(idx0), colnames(test_df))]

    pred_acc_test_trait <- eval_multiple_PRS_nocov(test_df, pgs_list, isbinary, ncores)
    pred_acc_test_trait_summary <- pred_acc_test_trait
    pred_acc_test_trait_summary <- pred_acc_test_trait_summary[order(as.numeric(pred_acc_test_trait_summary$pval), decreasing = FALSE), ]
    data.table::fwrite(pred_acc_test_trait_summary, testing_file, row.names = FALSE, sep = "\t", quote = FALSE)
  } else {
    writeLines(paste0("Reading testing file: ", testing_file))
    pred_acc_test_trait_summary <- data.table::fread(testing_file)
  }

  # If binary, (optionally) write AUC for bestPRS (no covariates)
  if (isbinary) {
    suppressMessages({
      ctrl <- caret::trainControl(method = "repeatedcv", allowParallel = TRUE, number = nfold_cv,
                                  returnData = FALSE, trim = TRUE, verboseIter = TRUE,
                                  classProbs = TRUE, summaryFunction = twoClassSummary)
    })
    # Recode binary trait to factor with levels 0/1 -> "X0"/"X1" to satisfy caret
    recode_bin <- function(v) { factor(ifelse(v == 1, "X1", "X0"), levels = c("X0", "X1")) }
    train_tmp <- data.frame(trait = recode_bin(train_df$trait),
                            z = scale(train_df[[bestPRS]]))
    set.seed(123)
    model_best <- caret::train(trait ~ z, data = train_tmp, method = "glmnet",
                               trControl = ctrl, family = "binomial", tuneLength = 50, metric = "ROC")
    test_tmp  <- data.frame(z = scale(test_df[[bestPRS]]), trait = recode_bin(test_df$trait))
    probs <- stats::predict(model_best, newdata = test_tmp, type = "prob")[, "X1"]
    # Use pROC for CI on AUC if available
    if (requireNamespace("pROC", quietly = TRUE)) {
      auc_ci <- pROC::ci.auc(test_tmp$trait, probs)
      auc_out <- data.frame(method = "bestPGS", auc = as.numeric(auc_ci[2]),
                            lowerCI = as.numeric(auc_ci[1]), upperCI = as.numeric(auc_ci[3]))
    } else {
      roc_obj <- pROC::roc(test_tmp$trait, probs, quiet = TRUE)
      auc_out <- data.frame(method = "bestPGS", auc = as.numeric(pROC::auc(roc_obj)),
                            lowerCI = NA_real_, upperCI = NA_real_)
    }
    data.table::fwrite(auc_out, paste0(out, "_auc_BestPGS.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

    # OR for bestPRS (PRS only)
    model1 <- stats::glm(trait ~ scale(z), data = transform(test_df, z = test_df[[bestPRS]]), family = "binomial")
    model1s <- summary(model1)
    mm <- exp(model1s$coefficients[2, 1])
    ll <- exp(model1s$coefficients[2, 1] - 1.97 * model1s$coefficients[2, 2])
    uu <- exp(model1s$coefficients[2, 1] + 1.97 * model1s$coefficients[2, 2])
    pval <- format.pval(model1s$coefficients[2, 4])
    writeLines(paste0("OR(bestPRS) = ", rr(mm), " (", rr(ll), "-", rr(uu), "); P-value=", pval))
    data.table::fwrite(data.frame(mm, ll, uu, pval), paste0(out, "_OR_BestPGS.txt"),
                       row.names = FALSE, sep = "\t", quote = FALSE)
  }

  pred_acc_train_trait_summary <- dplyr::filter(pred_acc_train_allPGS_summary, pgs %in% pgs_list)
  pred_acc_test_trait_summary_out <- pred_acc_test_trait_summary

  writeLines("--- Iterating power and p-value parameters ---")
  for (power_thres in power_thres_list) {
    for (pval_thres in pval_thres_list) {

      writeLines(paste0("P = ", pval_thres))
      writeLines(paste0("Power (R2) >= ", power_thres))
      writeLines("PRSmix:")

      topprs <- pred_acc_train_trait_summary %>%
        dplyr::filter(pval <= pval_thres & R2 >= power_thres) %>%
        dplyr::pull(pgs)
      topprs <- intersect(topprs, colnames(train_df))
      start_time <- Sys.time()

      if (length(topprs) == 0) {
        print("No high power trait-specific PRS for PRSmix")
        ww_raw <- 1
      } else {

        if (!isbinary) {
          x_train <- dplyr::select(train_df, dplyr::all_of(topprs))
          if (length(topprs) > 1) sd_train <- apply(as.data.frame(x_train[, topprs, drop = FALSE]), 2, stats::sd, na.rm = TRUE)
          x_train[, topprs] <- scale(x_train[, topprs])
          y_train <- as.vector(train_df$trait)
          train_data <- data.frame(x_train, trait = y_train)

          x_test <- dplyr::select(test_df, dplyr::all_of(topprs))
          y_test <- as.vector(test_df$trait)
          test_data <- data.frame(x_test, trait = y_test)

          formula <- stats::as.formula(paste0("trait ~ ", paste0(topprs, collapse = "+")))
          train_tmp <- train_data[, c("trait", topprs)]

          if (length(topprs) == 1) {
            ww <- ww_raw <- c(1); names(ww) <- names(ww_raw) <- topprs
          } else {
            ctrl <- caret::trainControl(method = "repeatedcv", allowParallel = TRUE,
                                        number = nfold_cv, verboseIter = TRUE)
            cl <- parallel::makePSOCKcluster(ncores)
            doParallel::registerDoParallel(cl)
            set.seed(123)
            model_prsmix <- caret::train(formula, data = train_tmp, method = "glmnet",
                                         trControl = ctrl, tuneLength = 50, verbose = TRUE)
            parallel::stopCluster(cl)
            ww <- coef(model_prsmix$finalModel, model_prsmix$bestTune$lambda)[, 1][-1]
            ww_raw <- ww
            if (all(ww == 0)) {
              writeLines("No weight for PRS; falling back to bestPRS")
              ww <- c(1); names(ww) <- bestPRS_acc$pgs; topprs <- bestPRS_acc$pgs
            } else {
              ww <- ww / sd_train[match(names(ww), names(sd_train))]
            }
          }

          test_df1 <- cbind(test_data, IID = test_df$IID)
          test_df1$newprs <- as.matrix(test_df1[, topprs, drop = FALSE]) %*% as.vector(ww)

          res_lm1 <- eval_single_PRS_nocov(test_df1, pheno = "trait", prs_name = "newprs", isbinary = isbinary, alpha = pval_thres)
          res_lm1$pgs <- "PRSmix"

        } else {
          x_train <- dplyr::select(train_df, dplyr::all_of(topprs))
          if (length(topprs) > 1) sd_train <- apply(as.data.frame(x_train[, topprs, drop = FALSE]), 2, stats::sd, na.rm = TRUE)
          x_train[, topprs] <- scale(x_train[, topprs])
          y_train <- as.vector(train_df$trait)
          train_data <- data.frame(x_train, trait = y_train)

          x_test <- dplyr::select(test_df, dplyr::all_of(topprs))
          y_test <- as.vector(test_df$trait)
          test_data <- data.frame(x_test, trait = y_test)

          formula <- stats::as.formula(paste0("trait ~ ", paste0(topprs, collapse = "+")))
          train_tmp <- train_data[, c("trait", topprs)]
          train_tmp$trait <- as.factor(train_tmp$trait)

          if (length(topprs) == 1) {
            ww <- ww_raw <- c(1); names(ww) <- names(ww_raw) <- topprs
          } else {
            ctrl <- caret::trainControl(method = "repeatedcv", allowParallel = TRUE,
                                        number = nfold_cv, verboseIter = TRUE)
            cl <- parallel::makePSOCKcluster(ncores)
            doParallel::registerDoParallel(cl)
            set.seed(123)
            model_prsmix <- caret::train(formula, data = train_tmp, method = "glmnet",
                                         trControl = ctrl, family = "binomial", tuneLength = 50, verbose = TRUE)
            parallel::stopCluster(cl)
            ww <- coef(model_prsmix$finalModel, model_prsmix$bestTune$lambda)[, 1][-1]
            ww_raw <- ww
            if (all(ww == 0)) {
              writeLines("No weight for PRS; falling back to bestPRS")
              ww <- c(1); names(ww) <- bestPRS_acc$pgs; topprs <- bestPRS_acc$pgs
            } else {
              ww <- ww / sd_train[match(names(ww), names(sd_train))]
            }

            # Optional PRSmix AUC (no covariates)
            test_data1 <- test_data
            test_data1[, topprs] <- scale(test_data1[, topprs])
            if (requireNamespace("pROC", quietly = TRUE)) {
              ctrl2 <- caret::trainControl(method = "repeatedcv", allowParallel = TRUE, number = nfold_cv,
                                           classProbs = TRUE, summaryFunction = twoClassSummary)
              recode_bin <- function(v) { factor(ifelse(v == 1, "X1", "X0"), levels = c("X0", "X1")) }
              train_tmp2 <- train_tmp; train_tmp2$trait <- recode_bin(as.integer(as.character(train_tmp2$trait)))
              set.seed(123)
              model_prsmix_auc <- caret::train(formula, data = train_tmp2, method = "glmnet",
                                               trControl = ctrl2, family = "binomial", tuneLength = 50, metric = "ROC")
              test_tmp2 <- data.frame(test_data1, trait = recode_bin(test_df$trait))
              probs <- stats::predict(model_prsmix_auc, newdata = test_tmp2, type = "prob")[, "X1"]
              auc_ci <- pROC::ci.auc(test_tmp2$trait, probs)
              auc_out <- data.frame(method = "prsmix", auc = as.numeric(auc_ci[2]),
                                    lowerCI = as.numeric(auc_ci[1]), upperCI = as.numeric(auc_ci[3]))
              data.table::fwrite(auc_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_auc_PRSmix.txt"),
                                 row.names = FALSE, sep = "\t", quote = FALSE)
            }
          }

          test_df1 <- cbind(test_data, IID = test_df$IID)
          test_df1$newprs <- as.matrix(test_df1[, topprs, drop = FALSE]) %*% as.vector(ww)
          res_lm1 <- eval_single_PRS_nocov(test_df1, pheno = "trait", prs_name = "newprs", isbinary = isbinary, alpha = pval_thres)
          res_lm1$pgs <- "PRSmix"

          # OR (PRS only)
          model1 <- stats::glm(trait ~ scale(newprs), data = test_df1, family = "binomial")
          model1s <- summary(model1)
          mm <- exp(model1s$coefficients[2, 1])
          ll <- exp(model1s$coefficients[2, 1] - 1.97 * model1s$coefficients[2, 2])
          uu <- exp(model1s$coefficients[2, 1] + 1.97 * model1s$coefficients[2, 2])
          pval <- format.pval(model1s$coefficients[2, 4])
          writeLines(paste0("OR(PRSmix) = ", rr(mm), " (", rr(ll), "-", rr(uu), "); P-value=", pval))
          data.table::fwrite(data.frame(mm, ll, uu, pval),
                             paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_OR_PRSmix.txt"),
                             row.names = FALSE, sep = "\t", quote = FALSE)
        }

        data.table::fwrite(data.frame(c(topprs), ww),
                           paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_PRSmix.txt"),
                           sep = "\t", quote = FALSE, row.names = FALSE)

        if (is_extract_adjSNPeff) {
          mixing_weight_file <- paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_PRSmix.txt")
          outfile <- paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_adjSNPeff_PRSmix.txt")
          extract_adjSNPeff(mixing_weight_file, original_beta_files_list, outfile)
        }

        res_lm1_summary <- res_lm1
        res_lm1_summary$pgs <- "PRSmix"
        pred_acc_test_trait_summary_out <- dplyr::bind_rows(res_lm1, pred_acc_test_trait_summary)
        data.table::fwrite(pred_acc_test_trait_summary_out,
                           paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_test_summary_traitPRS_withPRSmix.txt"),
                           row.names = FALSE, sep = "\t", quote = FALSE)

        prs_out <- test_df1 %>% dplyr::select(IID, newprs)
        colnames(prs_out) <- c("IID", "prsmix")
        data.table::fwrite(prs_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_prsmix.txt"),
                           row.names = FALSE, sep = "\t", quote = FALSE)
      }

      end_time <- Sys.time()
      timerunning <- difftime(end_time, start_time, units = "secs")[[1]]
      timedf <- data.frame(pgs = "PRSmix", npgs = length(ww_raw), time = timerunning)
      data.table::fwrite(timedf, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_time_PRSmix.txt"),
                         row.names = FALSE, sep = "\t", quote = FALSE)

      # --------------------------- PRSmix+ ---------------------------
      writeLines("PRSmix+:")

      topprs <- pred_acc_train_allPGS_summary %>%
        dplyr::filter(pval <= pval_thres & R2 >= power_thres) %>%
        dplyr::pull(pgs)
      topprs <- intersect(topprs, colnames(train_df))
      start_time <- Sys.time()

      if (length(topprs) == 0) {
        print("No high power PRS for PRSmix+")
        ww_raw <- 1
      } else {

        if (!isbinary) {
          x_train <- dplyr::select(train_df, dplyr::all_of(topprs))
          if (length(topprs) > 1) sd_train <- apply(as.data.frame(x_train[, topprs, drop = FALSE]), 2, stats::sd, na.rm = TRUE)
          x_train[, topprs] <- scale(x_train[, topprs])
          y_train <- as.vector(train_df$trait)
          train_data <- data.frame(x_train, trait = y_train)

          x_test <- dplyr::select(test_df, dplyr::all_of(topprs))
          y_test <- as.vector(test_df$trait)
          test_data <- data.frame(x_test, trait = y_test)

          formula <- stats::as.formula(paste0("trait ~ ", paste0(topprs, collapse = "+")))
          train_tmp <- train_data[, c("trait", topprs)]

          if (length(topprs) == 1) {
            ww <- c(1); names(ww) <- topprs
          } else {
            ctrl <- caret::trainControl(method = "repeatedcv", allowParallel = TRUE,
                                        number = nfold_cv, verboseIter = TRUE)
            cl <- parallel::makePSOCKcluster(ncores)
            doParallel::registerDoParallel(cl)
            set.seed(123)
            model_prsmix <- caret::train(formula, data = train_tmp, method = "glmnet",
                                         trControl = ctrl, tuneLength = 50, verbose = TRUE)
            parallel::stopCluster(cl)
            ww <- coef(model_prsmix$finalModel, model_prsmix$bestTune$lambda)[, 1][-1]
            ww_raw <- ww
            if (all(ww == 0)) {
              writeLines("No weight for PRS; falling back to bestPRS")
              ww <- c(1); names(ww) <- bestPRS_acc$pgs; topprs <- bestPRS_acc$pgs
            } else {
              ww <- ww / sd_train[match(names(ww), names(sd_train))]
            }
          }

          test_df1 <- cbind(test_data, IID = test_df$IID)
          test_df1$newprs <- as.matrix(test_df1[, topprs, drop = FALSE]) %*% as.vector(ww)
          res_lm <- eval_single_PRS_nocov(test_df1, pheno = "trait", prs_name = "newprs", isbinary = isbinary, alpha = pval_thres)
          res_lm$pgs <- "PRSmix+"
          nonzero_w <- names(ww[ww != 0])

        } else {
          x_train <- dplyr::select(train_df, dplyr::all_of(topprs))
          if (length(topprs) > 1) sd_train <- apply(as.data.frame(x_train[, topprs, drop = FALSE]), 2, stats::sd, na.rm = TRUE)
          x_train[, topprs] <- scale(x_train[, topprs])
          y_train <- as.vector(train_df$trait)
          train_data <- data.frame(x_train, trait = y_train)

          x_test <- dplyr::select(test_df, dplyr::all_of(topprs))
          y_test <- as.vector(test_df$trait)
          test_data <- data.frame(x_test, trait = y_test)

          formula <- stats::as.formula(paste0("trait ~ ", paste0(topprs, collapse = "+")))
          train_tmp <- train_data[, c("trait", topprs)]
          train_tmp$trait <- as.factor(train_tmp$trait)

          if (length(topprs) == 1) {
            ww <- c(1); names(ww) <- topprs
          } else {
            ctrl <- caret::trainControl(method = "repeatedcv", allowParallel = TRUE,
                                        number = nfold_cv, verboseIter = TRUE)
            cl <- parallel::makePSOCKcluster(ncores)
            doParallel::registerDoParallel(cl)
            set.seed(123)
            model_prsmix <- caret::train(formula, data = train_tmp, method = "glmnet",
                                         trControl = ctrl, family = "binomial", tuneLength = 50, verbose = TRUE)
            parallel::stopCluster(cl)
            ww <- coef(model_prsmix$finalModel, model_prsmix$bestTune$lambda)[, 1][-1]
            ww <- ww[ww != 0]
            ww_raw <- ww
            if (all(ww == 0)) {
              writeLines("No weight for PRS; falling back to bestPRS")
              ww <- c(1); names(ww) <- bestPRS_acc$pgs; topprs <- bestPRS_acc$pgs
            } else {
              ww <- ww / sd_train[match(names(ww), names(sd_train))]
            }

            # Optional PRSmix+ AUC (no covariates)
            test_data1 <- test_data
            test_data1[, topprs] <- scale(test_data1[, topprs])
            if (requireNamespace("pROC", quietly = TRUE)) {
              ctrl2 <- caret::trainControl(method = "repeatedcv", allowParallel = TRUE, number = nfold_cv,
                                           classProbs = TRUE, summaryFunction = twoClassSummary)
              recode_bin <- function(v) { factor(ifelse(v == 1, "X1", "X0"), levels = c("X0", "X1")) }
              train_tmp2 <- train_tmp; train_tmp2$trait <- recode_bin(as.integer(as.character(train_tmp2$trait)))
              set.seed(123)
              model_prsmix_auc <- caret::train(formula, data = train_tmp2, method = "glmnet",
                                               trControl = ctrl2, family = "binomial", tuneLength = 50, metric = "ROC")
              test_tmp2 <- data.frame(test_data1, trait = recode_bin(test_df$trait))
              probs <- stats::predict(model_prsmix_auc, newdata = test_tmp2, type = "prob")[, "X1"]
              auc_ci <- pROC::ci.auc(test_tmp2$trait, probs)
              auc_out <- data.frame(method = "prsmixP", auc = as.numeric(auc_ci[2]),
                                    lowerCI = as.numeric(auc_ci[1]), upperCI = as.numeric(auc_ci[3]))
              data.table::fwrite(auc_out, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_auc_PRSmixPlus.txt"),
                                 row.names = FALSE, sep = "\t", quote = FALSE)
            }
          }

          test_df1 <- cbind(test_data, IID = test_df$IID)
          test_df1$newprs <- as.matrix(test_df1[, topprs, drop = FALSE]) %*% as.vector(ww)
          res_lm <- eval_single_PRS_nocov(test_df1, pheno = "trait", prs_name = "newprs", isbinary = isbinary, alpha = pval_thres)
          res_lm$pgs <- "PRSmix+"
          nonzero_w <- names(ww[ww != 0])

          # OR (PRS only)
          model <- stats::glm(trait ~ scale(newprs), data = test_df1, family = "binomial")
          model1s <- summary(model)
          mm <- exp(model1s$coefficients[2, 1])
          ll <- exp(model1s$coefficients[2, 1] - 1.97 * model1s$coefficients[2, 2])
          uu <- exp(model1s$coefficients[2, 1] + 1.97 * model1s$coefficients[2, 2])
          pval <- format.pval(model1s$coefficients[2, 4])
          writeLines(paste0("OR(PRSmix+) = ", rr(mm), " (", rr(ll), "-", rr(uu), "); P-value=", pval))
          data.table::fwrite(data.frame(mm, ll, uu, pval),
                             paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_OR_PRSmixPlus.txt"),
                             row.names = FALSE, sep = "\t", quote = FALSE)
        }

        pred_acc_test_trait_summary_out <- dplyr::bind_rows(res_lm, pred_acc_test_trait_summary_out)
        data.table::fwrite(pred_acc_test_trait_summary_out,
                           paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_test_summary_traitPRS_withPRSmixPlus.txt"),
                           row.names = FALSE, sep = "\t", quote = FALSE)

        prsmixplus <- test_df1 %>% dplyr::select(IID, newprs)
        data.table::fwrite(prsmixplus, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_prsmixPlus.txt"),
                           row.names = FALSE, sep = "\t", quote = FALSE)

        data.table::fwrite(data.frame(topprs, ww),
                           paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_PRSmixPlus.txt"),
                           row.names = FALSE, sep = "\t", quote = FALSE)

        if (is_extract_adjSNPeff) {
          mixing_weight_file <- paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_PRSmixPlus.txt")
          outfile <- paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_adjSNPeff_PRSmixPlus.txt")
          extract_adjSNPeff(mixing_weight_file, original_beta_files_list, outfile)
        }
        data.table::fwrite(data.frame(topprs, ww_raw),
                           paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_weight_raw_PRSmixPlus.txt"),
                           row.names = FALSE, sep = "\t", quote = FALSE)
      }

      end_time <- Sys.time()
      timerunning <- difftime(end_time, start_time, units = "secs")[[1]]
      timedf <- data.frame(pgs = "PRSmix+", npgs = length(ww_raw), time = timerunning)
      data.table::fwrite(timedf, paste0(out, "_power.", power_thres, "_pthres.", pval_thres, "_time_PRSmixPlus.txt"),
                         row.names = FALSE, sep = "\t", quote = FALSE)
    }
  }

  writeLines("Finished")
  return(0)
}
