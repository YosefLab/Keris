#' @title Fitting the keris model
#'
#' @description Fits the TCA model for an input matrix of features by observations that are coming from a mixture of \code{k} sources, under the assumption that each observation is a mixture of unique (unobserved) source-specific values (in each feature in the data). This function further allows to statistically test the effect of covariates on source-specific values. For example, in the context of tissue-level bulk DNA methylation data coming from a mixture of cell types (i.e. the input is methylation sites by individuals), \code{tca} allows to model the methylation of each individual as a mixture of cell-type-specific methylation levels that are unique to the individual. In addition, it allows to statistically test the effects of covariates and phenotypes on methylation at the cell-type level.
#'
#' @param X An \code{n} by \code{m} matrix of bulk gene expression measurements of \code{m} genes for \code{n} observations (samples). Each row in \code{X} is assumed to be a mixture of \code{k} sources. Note that \code{X} must include row names and column names and that NA values are not supported.
#' @param W An \code{n} by \code{k} matrix of cell-type proportions. Each row, which corresponds to a specific sample in \code{X}, includes the cell-type proportions of \code{k} cell types in one sample (nonnegative and sum up to 1). Note that \code{W} must include row names and column names and that NA values are not supported.
#' @param y an \code{n}-length vector of binary grouping of the \code{n} samples in \code{X} into two groups (values should be coded as 0 or 1). NA values are not supported.
#' @param C1 An \code{n} by \code{p1} design matrix of covariates that may have cell-type-specific effects. Note that \code{C1} must include row names and column names and should not include an intercept term. NA values are not supported. Note that each covariate in \code{C1} results in \code{k} additional parameters in the model of each feature, therefore, in order to alleviate the possibility of model overfitting, it is advised to be mindful of the balance between the size of \code{C1} and the number of samples \code{n}in \code{X}.
#' @param C2 An \code{n} by \code{p2} design matrix of covariates that may affect the mixture (i.e., rather than have cell-type-specific effects). Note that \code{C2} must include row names and column names and should not include an intercept term. NA values are not supported.
#' @param covar_pairs A matrix with 2 columns, where each row specifies the names of two genes for which keris will learn their cell-type level gene-gene covariances. Each particular gene should be indicated by its name as it appears in the column names of \code{X}. If \code{covar_pairs == NULL} then no covariances will be estimated.
#' @param num_perms The number of permutations to perform in evaluate permutation-based statistical significance of differential moments. If \code{num_perms == 0} then no permutations will be performed.
#' @param fdr A desired false discovery rate for the permutation-based evaluation of differential moments (applied only if \code{num_perms != NULL}).
#' @param max_iters A numeric value indicating the maximal number of iterations to use in the GMM optimization of the Keris model (\code{max_iters} iterations will be used as long as the optimization has not converged earlier).
#' @param max_sds A numeric value for determining which observations should be considered as outliers and therefore should not be taken into account for estimating moments. if \code{max_sds == NULL} then all observations will be used for estimating moments, and otherwise only observations that are within \code{max_sds} standard deviations from the mean will be considered.
#' @param penalty An L2 penalty to be used in the optimization for numeric stability. The actual penalty value used in a given GMM optimization is the provided value scaled by the standard deviation of the data.
#' @param learn_tau A logical value indicating whether to learn an i.i.d. component of variance for each gene.
#' @param constrain_mus A logical value indicating whether to constrain the estimates of the mean parameters.
#' @param scale_data A logical value indicating whether to scale the data before estimating moments. If \code{scale_data == FALSE} then estimates will be based on non-scaled data, however, if \code{num_perms != NULL} data will eventually be scaled and moments will be further estimated again in order to allow more stable statistical testing.
#' @param rand_seed A numeric value for setting random seed.
#' @param parallel A logical value indicating whether to use parallel computing (possible when using a multi-core machine). Note that parallel computing is only used for for running permutation (i.e., if \code{num_perms != 0}).
#' @param num_cores A numeric value indicating the number of cores to use (activated only if \code{parallel == TRUE}). If \code{num_cores == NULL} then all available cores except for one will be used.
#' @param log_file A path to an output log file. Note that if the file \code{log_file} already exists then logs will be appended to the end of the file. Set \code{log_file} to \code{NULL} to prevent output from being saved into a file; note that if \code{verbose == FALSE} then no output file will be generated regardless of the value of \code{log_file}.
#' @param verbose A logical value indicating whether to print logs.
#' @param debug_mode A logical value indicating whether to set the logger to a more detailed debug level; set \code{debug} to \code{TRUE} before reporting issues.
#'
#' @details ...
#'
#' @return A list with two items. The first item is a list with the estimated parameters of the model under \code{y=0} and the second item is a list with the estimated parameters of the model under \code{y=1}. Each of these lists include estimates of moments. In addition, the second list includes results of the statistical testing for evaluating whether each moment is different between the two conditions (i.e. \code{y=0} and \code{y=1}).
#'
#' @export keris
keris <- function(X, W, y, C1 = NULL, C2 = NULL, covar_pairs = NULL,
                  num_perms = 0, fdr = 0.05, max_iters = 3, max_sds = 4, penalty = 1e-4,
                  learn_tau = FALSE, constrain_mus = TRUE, scale_data = TRUE, rand_seed = NULL,
                  parallel = FALSE, num_cores = NULL, log_file = "keris.log", verbose = TRUE,
                  debug_mode = FALSE){

  start_logger(log_file, debug_mode, verbose)

  flog.info("Starting keris...")

  # global parameters
  globals <- list(
    "MIN_SQRT_INV_WEIGHT" = 1, # the minimal inverse of the square root of the weight allowed for a single observation in the optimization; required for numeric stability so that a small number of data points don't become too dominant in the estimation.
    "MIN_VAR" = 1e-10, # minimum variance; for numeric stability
    "EPSILON" = 1e-5, # criterion for convergence ; this number will be multiplied by the number of parameters in the data
    "MIN_MEAN" = -10e6, # minimum mean; for numeric stability
    "MAX_MEAN" = 10e6 # maximum mean; for numeric stability
  )

  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(W)

  # input validations and warnings
  flog.info("Validating input...")
  valid_input <- keris.validate_input(X, W, y, C1, C2, covar_pairs, num_perms, fdr, max_iters, max_sds, penalty, learn_tau, constrain_mus, scale_data, rand_seed, parallel, num_cores, debug)
  if (!valid_input){
    flog.error("There was at least one critical warning that must be adressed; terminating execution.")
    return (NULL);
  }

  if (parallel && num_perms > 0){
    flog.debug("Initiating cluster...")
    cl <- init_cluster(num_cores)
    clusterEvalQ(cl, library("keris"))
  }else{
    cl <- NULL
  }

  flog.info("Detect extreme data points to exclude from the model optimization...")
  if (is.null(max_sds)) max_sds <- Inf
  valid_data <- abs(scale(X, center = TRUE, scale = FALSE)) <= repmat(colSds(X)*max_sds,n,1)

  flog.info("Looking for constant genes...")
  sds <- unlist(lapply(1:m, function(j) sd(X[valid_data[,j],j])))
  if (sum(sds==0)){
    flog.info(sprintf("Found constant genes; that is, genes with no variance after excluding outlier data points: %s.\nProceeding without these genes...",paste(colnames(X)[sds==0], collapse=", ")))
    X <- X[,sds>0]
    valid_data <- valid_data[,sds>0]
    m <- ncol(X)
    if (m==0){
      flog.error("No genes left after exclusion; terminating execution.")
      return (NULL);
    }
  }
  sds <- unlist(lapply(1:ncol(X), function(j) sd(X[valid_data[,j],j])))

  flog.info("Initializing parameters...")
  if (!is.null(rand_seed)) set.seed(rand_seed)
  if (!is.null(covar_pairs)){
    covar_pairs.indices <- get_covar_pairs_indices(covar_pairs, colnames(X))
    rownames(covar_pairs) <- rownames(covar_pairs.indices)
  }

  if (is.null(C1)){
    C1 <- matrix(0, nrow=nrow(X), ncol=0)
    rownames(C1) <- rownames(X)
  }
  if (is.null(C2)){
    C2 <- matrix(0, nrow=nrow(X), ncol=0)
    rownames(C2) <- rownames(X)
  }
  y <- as.matrix(y)
  colnames(y) <- "y"
  p1 <- ncol(C1)
  p2 <- ncol(C2)
  C1_ <- calc_C1_W_interactions(W,C1)
  W_squared <- W**2
  W_cross_prods <- calc_W_cross_products(W)
  W.inter.y <- if(is.null(y)) NULL else calc_C1_W_interactions(W,y)
  W.inter.0 <- W*repmat(1-y,1,k)
  gammas.y_indices <- seq(p1+1,k*(p1+1),p1+1)

  profiles <- list()
  for (i in 1:(1+!is.null(y))){
    profiles[[i]] <- list("means" = list(), "vars" = list(), "covars" = list())
    profiles[[i]]$means$features <- colnames(X)
  }

  if (scale_data){
    scale_factors <- unlist(lapply(1:m, function(j) attr(scale(X[valid_data[,j],j],center=FALSE),"scaled:scale")))
    X <- X/repmat(scale_factors,n,1)
  }
  mus_range <- if (constrain_mus) c(min(X),max(X)) else NULL

  flog.info("Learn means (mu, delta, and gamma parameters)...")
  means.lst <- keris.fit_means(X, W, W_squared, C2, C1_, y, W.inter.0, W.inter.y, valid_data, mus_range, max_iters, penalty, cl = NULL, globals)
  flog.debug("Updating profiles...")
  for (i in 1:length(profiles)){
    profiles[[i]]$means$features <- colnames(X)
  }
  profiles[[1]]$means$mus_hat <- means.lst$mus_hat
  profiles[[1]]$means$betas_hat <- means.lst$betas_hat
  profiles[[1]]$means$gammas_hat <- means.lst$gammas_hat[,(1*(p1>0)):(p1*k),drop=F]
  if (!is.null(y)){
    profiles[[2]]$means$mus_hat <- means.lst$gammas_hat[,(p1*k+1):ncol(means.lst$gammas_hat)]
    profiles[[2]]$means$betas_hat <- means.lst$betas_hat
    profiles[[2]]$means$gammas_hat <- means.lst$gammas_hat[,(1*(p1>0)):(p1*k),drop=F]
  }
  X_means_hat <- means.lst$X_means_hat


  # Learn variances
  squared_residuals <- calc_squared_residuals(X, X_means_hat)
  W_.inter.y <- if (is.null(y)) NULL else calc_C1_W_interactions(W_squared,y)
  W_.inter.0 <- if (is.null(y)) NULL else calc_C1_W_interactions(W_squared,1-y)
  learn_tau_msg <- if (learn_tau) "and tau parameters" else "parameters"
  flog.info(sprintf("Learn variances (sigma %s)...",learn_tau_msg))
  for (i in 1:(1+!is.null(y))){
    profiles[[i]]$vars <- list()
  }
  vars.lst <- keris.fit_vars(X, X_means_hat, W, W_squared, y, W_.inter.0, W_.inter.y, W_cross_prods, squared_residuals, valid_data, max_iters, penalty, learn_tau, globals, cl = NULL)
  flog.debug("Updating profiles...")
  for (i in 1:length(profiles)){
    profiles[[i]]$vars$features <- colnames(X)
  }
  tau_squared_hat <- vars.lst$tau_squared_hat
  sigmas_squared_hat <- vars.lst$sigmas_squared_hat
  sigmas_squared_hat[sigmas_squared_hat < globals[["MIN_VAR"]]] <- 0
  profiles[[1]]$vars$sigmas_squared_hat <- sigmas_squared_hat
  profiles[[1]]$vars$tau_squared_hat <- vars.lst$tau_squared_hat
  if (!is.null(y)){
    sigmas_squared_hat.delta <- vars.lst$sigmas_squared_hat.delta
    rownames(sigmas_squared_hat.delta) <- colnames(X)
    profiles[[2]]$vars$sigmas_squared_hat <- sigmas_squared_hat + sigmas_squared_hat.delta
    profiles[[2]]$vars$sigmas_squared_hat[profiles[[2]]$vars$sigmas_squared_hat < globals[["MIN_VAR"]]] <- 0
    profiles[[2]]$vars$tau_squared_hat <- vars.lst$tau_squared_hat
  }

  if (!is.null(covar_pairs)){
    flog.info("Learn covars...")
    squared_residuals.covars <- calc_squared_residuals(X, X_means_hat, covar_pairs = covar_pairs.indices)
    covars.lst <- keris.fit_covars(X, W, W_squared, covar_pairs.indices, y, W_.inter.0, W_.inter.y, W_cross_prods, squared_residuals.covars, valid_data, sigmas_squared_hat, sigmas_squared_hat.delta, tau_squared_hat, max_iters, penalty, globals, cl = NULL)
    profiles[[1]]$covars <- list("features" = covar_pairs,
                                 "sigmas_squared_hat" = covars.lst$sigmas_squared_hat.covars)
    if (!is.null(y)){
      profiles[[2]]$covars <- list("features" = covar_pairs,
                                   "sigmas_squared_hat" = covars.lst$sigmas_squared_hat.covars + covars.lst$sigmas_squared_hat.covars.delta)
    }
  }

  if (!is.null(y)){
    if (scale_data){
      flog.info("Calculate t ratio statistics for the means...")
      profiles[[2]]$means$mus_hat.stats <- calc_t_statistics(X, W, S = cbind(W,W.inter.y,C1_,C2),
                                                             U_hat_sqrt_inv = calc_U_sqrt_diag(X, W, W_squared, X_means_hat, iter = max_iters, globals = globals),
                                                             theta_hat = profiles[[2]]$means$mus_hat-profiles[[1]]$means$mus_hat,
                                                             valid_data, globals, cl = NULL)
      flog.info("Calculate t ratio statistics for the vars...")
      S.vars <- cbind(W_squared,W_.inter.y)
      if (learn_tau) S.vars <- cbind(S.vars,numeric(n)+1)
      profiles[[2]]$vars$sigmas_squared_hat.stats <- calc_t_statistics(X, W, S = S.vars,
                                                                       U_hat_sqrt_inv = calc_V_sqrt_diag(X, W_.inter.0, W_.inter.y, sigmas_squared_hat, sigmas_squared_hat.delta, tau_squared_hat, squared_residuals, iter = max_iters, globals = globals),
                                                                       theta_hat = profiles[[2]]$vars$sigmas_squared_hat-profiles[[1]]$vars$sigmas_squared_hat,
                                                                       valid_data, globals, cl = NULL)
      flog.info("Calculate t ratio statistics for the covars...")
      if (!is.null(covar_pairs)){
        profiles[[2]]$covars$sigmas_squared_hat.stats <- calc_t_statistics(X, W, S = cbind(W_squared,W_.inter.y),
                                                                           U_hat_sqrt_inv = calc_V_sqrt_diag(X, W_.inter.0, W_.inter.y, profiles[[1]]$covars$sigmas_squared_hat, profiles[[2]]$covars$sigmas_squared_hat-profiles[[1]]$covars$sigmas_squared_hat, tau_squared_hat = 0, squared_residuals = squared_residuals.covars, iter = max_iters, globals = globals),
                                                                           theta_hat = profiles[[2]]$covars$sigmas_squared_hat-profiles[[1]]$covars$sigmas_squared_hat,
                                                                           valid_data, globals, cl = NULL, covar_pairs.indices)
      }
    }else{

      flog.info("Rerunning Keris on a scaled version of the data in order to allow calculation of t ratio statistics...")
      keris.scaled.mdl <- keris(X, W, y, C1, C2, covar_pairs, num_perms = 0, fdr = fdr,
                                max_iters = max_iters, max_sds = max_sds, penalty = penalty, learn_tau = learn_tau,
                                constrain_mus = constrain_mus, scale_data = TRUE,
                                rand_seed = NULL, parallel = parallel, num_cores = num_cores, log_file = NULL, verbose = verbose, debug_mode = debug_mode)
      profiles[[2]]$means$mus_hat.stats <- keris.scaled.mdl[[2]]$means$mus_hat.stats
      profiles[[2]]$vars$sigmas_squared_hat.stats <- keris.scaled.mdl[[2]]$vars$sigmas_squared_hat.stats
      profiles[[2]]$covars$sigmas_squared_hat.stats <- keris.scaled.mdl[[2]]$covars$sigmas_squared_hat.stats
    }
  }

  if (!is.null(y) & (num_perms > 0)){
    moments <- c("means","vars")
    num_features <- c(m,m)
    if (!is.null(covar_pairs)){
      moments <- c(moments,"covars")
      num_features <- c(num_features,nrow(covar_pairs))
    }
    null_statistics <- list()
    for (i in 1:length(moments)){
      moment = moments[i]
      null_statistics[[moment]] <- vector("list", length = k)
      for (h in 1:k){
        null_statistics[[moment]][[h]] <- matrix(0,num_perms,num_features[i])
      }
    }

    # permutations
    flog.info(sprintf("Running  %d permutations...", num_perms))
    res <- pblapply(1:num_perms,function(p){
      pboptions(type = "none")
      y.perm <- y[sample(length(y))]
      keris.mdl.perm <- keris(X, W, y.perm, C1, C2, covar_pairs, num_perms = 0, fdr = fdr,
                              max_iters = max_iters, max_sds = max_sds, penalty = penalty, learn_tau = learn_tau,
                              constrain_mus = constrain_mus, scale_data = TRUE,
                              rand_seed = rand_seed+p, parallel = FALSE, num_cores = num_cores, log_file = log_file, verbose = FALSE, debug_mode = debug_mode)
      #flog.info(sprintf("Permutation %d completed", p))
      print(p)
      return(keris.mdl.perm)
    }, cl = cl)
    pboptions(type = "txt")

    if (parallel & !is.null(cl)){
      flog.debug("Stopping cluster...")
      stop_cluster(cl)
    }

    for (p in 1:num_perms){
      for (h in 1:k){
        null_statistics$means[[h]][p,] <- res[[p]][[2]]$means$mus_hat.stats$tstat[,h]
        null_statistics$vars[[h]][p,] <- res[[p]][[2]]$vars$sigmas_squared_hat.stats$tstat[,h]
        if (!is.null(covar_pairs)){
          null_statistics$covars[[h]][p,] <- res[[p]][[2]]$covars$sigmas_squared_hat.stats$tstat[,h]
        }
      }
    }

    profiles[[2]]$means$mus_hat.stats$is_significant <- matrix(0,m,k)
    rownames(profiles[[2]]$means$mus_hat.stats$tstat) <- colnames(X)
    colnames(profiles[[2]]$means$mus_hat.stats$tstat) <- colnames(W)
    profiles[[2]]$vars$sigmas_squared_hat.stats$is_significant <- matrix(0,m,k)
    rownames(profiles[[2]]$vars$sigmas_squared_hat.stats$is_significant) <- colnames(X)
    colnames(profiles[[2]]$vars$sigmas_squared_hat.stats$is_significant) <- colnames(W)
    if (!is.null(covar_pairs)){
      profiles[[2]]$covars$sigmas_squared_hat.stats$is_significant <- matrix(0,nrow(covar_pairs),k)
      rownames(profiles[[2]]$covars$sigmas_squared_hat.stats$is_significant) <- rownames(covar_pairs)
      colnames(profiles[[2]]$covars$sigmas_squared_hat.stats$is_significant) <- colnames(W)
    }
    win_size <- round(10*sqrt(num_perms))
    for (h in 1:k){
      # means
      citical_vals <- get_critical_values(S = profiles[[2]]$means$mus_hat.stats$tstat[,h], S.perm = null_statistics$means[[h]], fdr = fdr, window_size = win_size)
      profiles[[2]]$means$mus_hat.stats$is_significant[,h] <- (profiles[[2]]$means$mus_hat.stats$tstat[,h] <= citical_vals[1]) |
        (profiles[[2]]$means$mus_hat.stats$tstat[,h] >= citical_vals[2])
      # vars
      citical_vals <- get_critical_values(S = profiles[[2]]$vars$sigmas_squared_hat.stats$tstat[,h], S.perm = null_statistics$vars[[h]], fdr = fdr, window_size = win_size)
      profiles[[2]]$vars$sigmas_squared_hat.stats$is_significant[,h] <- (profiles[[2]]$vars$sigmas_squared_hat.stats$tstat[,h] <= citical_vals[1]) |
        (profiles[[2]]$vars$sigmas_squared_hat.stats$tstat[,h] >= citical_vals[2])
      if (!is.null(covar_pairs)){
        # covars
        citical_vals <- get_critical_values(S = profiles[[2]]$covars$sigmas_squared_hat.stats$tstat[,h], S.perm = null_statistics$covars[[h]], fdr = fdr, window_size = win_size)
        profiles[[2]]$covars$sigmas_squared_hat.stats$is_significant[,h] <- (profiles[[2]]$covars$sigmas_squared_hat.stats$tstat[,h] <= citical_vals[1]) |
          (profiles[[2]]$covars$sigmas_squared_hat.stats$tstat[,h] >= citical_vals[2])
      }
    }
  }

  return(profiles)
}



