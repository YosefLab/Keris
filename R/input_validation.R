
keris.validate_input <- function(X, W, y, C1, C2, covar_pairs, num_perms, fdr, max_iters, max_sds, penalty, learn_tau, constrain_mus, scale_data, means_ref, vars_ref, rand_seed, parallel, num_cores, debug){

  critical_warnings <- list()
  warnings <- list()

  # critical warnings
  if (!is.null(covar_pairs)){
    if (length(union(paste0(covar_pairs[,1],"@",covar_pairs[,2]),paste0(covar_pairs[,2],"@",covar_pairs[,1]))) != nrow(covar_pairs)*2){
      critical_warnings <- list.append(critical_warnings,
                                       "covar_pairs should not include the same pair of features more than once; note that a pair (x,y) is the same as (y,x). Also, make sure there are no pairs that include the same feature twice.")
    }
    if (!(dim(covar_pairs)[1]>0 & dim(covar_pairs)[2]==2)){
      critical_warnings <- list.append(critical_warnings, "The dimensions of covar_pairs are wrong.")
    }
  }
  if (!(is.null(covar_pairs) | is.matrix(covar_pairs))){
    critical_warnings <- list.append(critical_warnings,
                                     "covar_pairs must be either NULL or a matrix.")
  }
  if (!is.null(y)){
    #if (!is.matrix(y)) critical_warnings <- list.append(critical_warnings, "y must be of class matrix.")
    if (!all(sort(unique(y)) == c(0,1))) critical_warnings <- list.append(critical_warnings, "y must be encoded as 0/1.")
  }
  if (!(is.numeric(num_perms) && (num_perms>=0))){
    critical_warnings <- list.append(critical_warnings,
                                     "num_perms must be non-negative integer.")
  }

  if (is.null(rownames(X)) | is.null(colnames(X))) critical_warnings <- list.append(critical_warnings,"X must include row names and column names.")
  if (is.null(rownames(W)) | is.null(colnames(W))) critical_warnings <- list.append(critical_warnings,"W must include row names and column names.")
  if (!is.null(C1) ){
    if ((ncol(C1) > 0) & (is.null(rownames(C1)) | is.null(colnames(C1)))) critical_warnings <- list.append(critical_warnings,"C1 must include row names and column names.")
  }
  if (!is.null(C2)){
    if ((ncol(C2) > 0) & (is.null(rownames(C2)) | is.null(colnames(C2)))) critical_warnings <- list.append(critical_warnings,"C2 must include row names and column names.")
  }
  if (num_perms > 0){
    if (num_perms < 1/fdr) critical_warnings <- list.append(critical_warnings,"num_perms can be set either as 0 or as a number that is greater or equal to 1/fdr.")
  }

  if (length(critical_warnings)){
    for (i in 1:length(critical_warnings)) flog.warn(critical_warnings[[i]])
    return (FALSE)
  }

  # # non-critical warnings
  # if (is.null(covar_pairs) && (ncol(X)>=100)){
  #   warnings <- list.append(warnings,
  #                           sprintf("Warning: %d covariate pairs will be evaluated; this may take a while.",choose(ncol(X),2)))
  # }

  if (length(warnings)){
    for (i in 1:length(warnings)) flog.warn(warnings[[i]])
  }
  return (TRUE)

}
