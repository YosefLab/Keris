
assert <- function (expr, error) {
  if (! expr) stop(error, call. = FALSE)
}

# covar_pairs is a matrix with two columns, each row has the names of a pair of features.
get_covar_pairs_indices <- function(covar_pairs, feature_names){
  # convert gene names to indices in X
  covar_pairs.indices <- covar_pairs
  for (i in 1:nrow(covar_pairs)){
    covar_pairs.indices[i,1] <- which(feature_names == covar_pairs[i,1])
    covar_pairs.indices[i,2] <- which(feature_names == covar_pairs[i,2])
  }
  covar_pairs.indices <- Reshape(as.numeric(covar_pairs.indices),nrow(covar_pairs.indices),2)
  rownames(covar_pairs.indices) <- paste0(covar_pairs[,1],"@",covar_pairs[,2])
  return(covar_pairs.indices)
}

calc_C1_W_interactions <- function(W,C1){
  n <- nrow(W)
  k <- ncol(W)
  p1 <- ncol(C1)
  if (p1){
    return( hadamard.prod(Reshape(Reshape(apply(W, 2, function(v) repmat(v,1,p1)), n*p1*k,1), n,p1*k), repmat(C1, 1, k)) )
  }else{
    return(matrix(0,n,0))
  }
}

calc_W_cross_products <- function(W){
  l <- vector("list", length = nrow(W))
  for (i in 1:nrow(W)){
    l[[i]] <- crossprod(t(W[i,]),t(W[i,]))
  }
  return(l)
}


# covar_pairs - if not null then values will be calculated for the pairs of features in covar_pairs
calc_squared_residuals <- function(X, X_means_hat, covar_pairs = NULL){
  n <- nrow(X)
  if (is.null(covar_pairs)) return ((X-X_means_hat)**2)
  squared_residuals <- matrix(0,n,nrow(covar_pairs))
  for (j in 1:nrow(covar_pairs)){
    squared_residuals[,j] <- (X[,covar_pairs[j,1]]-X_means_hat[,covar_pairs[j,1]])*(X[,covar_pairs[j,2]]-X_means_hat[,covar_pairs[j,2]])
  }
  return(squared_residuals)
}


is_min_var <- function(x, min_var){
  var_epsilon <- 1e-8
  return ((x <= min_var) | (abs(x-min_var) < var_epsilon))
}

#' @importFrom flog appender
start_logger <- function(log_file, debug_mode, verbose){
  config_level <- if (debug_mode) "debug" else "default"
  Sys.setenv(R_CONFIG_ACTIVE = config_level)
  if (is.null(log_file)){
    invisible(flog.appender(appender.console()))
  }else{
    invisible(flog.appender(appender.tee(log_file)))
  }
  invisible(flog.threshold(if(debug_mode) "DEBUG" else "INFO"))
  if (!verbose) (flog.threshold("ERROR"))
}


init_cluster <- function(num_cores = NULL){
  flog.info("Initiate cluster...")
  cl <- makeCluster(get_num_cores(num_cores))
  flog.info("Parallel is on with %s nodes.",get_num_cores(num_cores))
  return(cl)
}

get_num_cores <- function(num_cores){
  if (is.null(num_cores)){
    return (max(1,detectCores() - 1))
  }else{
    return (num_cores)
  }
}

stop_cluster <- function(cl){
  flog.info("Stop cluster")
  stopCluster(cl)
}

# Calculates two-sided p-values under standard normal distribution
calc_pvalues <- function(x, calc_fdr = TRUE){
  pvalues <- 2*pnorm(-abs(x))
  if (calc_fdr){
    return(list("pvalues" = pvalues, "pvalues.adj" = p.adjust(pvalues, method = "fdr")  ))
  }else{
    return(pvalues)
  }
}

get_all_pairs <- function(m){
  pairs <- matrix(0,choose(m,2),2)
  counter <- 1
  for (j in 1:(m-1)){
    for (l in (j+1):m){
      pairs[counter,] <- c(j,l)
      counter <- counter + 1
    }
  }
  return(pairs)
}
