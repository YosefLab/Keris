
keris.fit_means <- function(X, W, W_squared, C2, C1_, y, W.inter.0, W.inter.y, valid_data, mus_range, max_iters, penalty, cl, globals){

  flog.debug("Running keris.fit_means...")

  m <- ncol(X)
  k <- ncol(W)
  p1 <- ncol(C1_)/k
  p2 <- ncol(C2)

  mus_hat <- matrix(0, nrow=m, ncol=k)
  rownames(mus_hat) <- colnames(X)
  colnames(mus_hat) <- colnames(W)
  betas_hat <- matrix(0, nrow=m, ncol=p2)
  rownames(betas_hat) <- colnames(X)
  colnames(betas_hat) <- colnames(C2)
  gammas_hat <- matrix(0, nrow=m, ncol=(p1+!is.null(y))*k)
  rownames(gammas_hat) <- colnames(X)
  colnames(gammas_hat) <- colnames(C1_)
  X_means_hat <- NULL

  if (penalty){
    ols_solver <- function(z,x) crossprod(z, t(tcrossprod(inv(crossprod(x,x)+eye(ncol(x))+penalty*sd(z)),x)))
  }else{
    ols_solver <- function(z,x) crossprod(z, t(tcrossprod(inv(crossprod(x,x)),x)))
  }

  flog.debug("Start optimization of means...")
  for (iter in 1:max_iters){
    flog.debug(sprintf("Iteration %d...", iter))
    U_sqrt_diag <- calc_U_sqrt_diag(X, W, W_squared, X_means_hat, iter, globals)
    U_sqrt_diag[which(U_sqrt_diag < globals[["MIN_SQRT_INV_WEIGHT"]])] = globals[["MIN_SQRT_INV_WEIGHT"]]
    if (!is.null(cl)){
      var_list <- c("X","valid_data","U_sqrt_diag","W","C2","C1_","W.inter.0","W.inter.y","k","p1","p2","mus_range","penalty", "globals","ols_solver","repmat","eye","pcls")
      clusterExport(cl, varlist = var_list, envir=environment())
    }
    flog.debug("Learning parameters...")
    res <- pblapply(1:m,function(j) {
      # Estimate all the mean parameters for feature j
      r <- X[valid_data[,j],j]/U_sqrt_diag[valid_data[,j],j]
      x <- cbind(W.inter.0[valid_data[,j],],C2[valid_data[,j],],cbind(C1_[valid_data[,j],],W.inter.y[valid_data[,j],]))/t(repmat(U_sqrt_diag[valid_data[,j],j],2*k+p1*k+p2,1))
      # For numeric stability, normalize the design matrix and adjust the final estimates accordingly; note that using the sqrt of the norms turns out to be more stable in pcls.
      norms <- (colSums(x**2))**0.5
      norms <- sqrt(norms)
      x <- x/repmat(norms,nrow(x),1)
      if (is.null(mus_range)){
        return( tryCatch({
          ols_solver(r, x)/norms
        }, warning = function(w) {
          if (grepl("Matrix appears to be singular", as.character(w), fixed = TRUE)){
            numeric(2*k+p2+k*p1) + NA
          }else{
            stop();
          }
        }) )
      }else{
        lb_means <- numeric(2*k+p2+k*p1) + globals[["MIN_MEAN"]]
        lb_means[c(1:k,(k+p2+k*p1+1):length(lb_means))] <- mus_range[1]
        ub_means <- numeric(2*k+p2+k*p1) + globals[["MAX_MEAN"]]
        ub_means[c(1:k,(k+p2+k*p1+1):length(lb_means))] <- mus_range[2]
        range.0 <- 1:k
        range.1 <- (k+p2+p1*k+1):length(norms)
        x0.0 <- 0.99*(norms[range.0]*lb_means[range.0])+0.01*(norms[range.0]*ub_means[range.0])
        x0.1 <- 0.99*(norms[range.1]*lb_means[range.1])+0.01*(norms[range.1]*ub_means[range.1])
        M <- list(y = r,
                  w = numeric(length(r))+1,
                  X = x,
                  C = matrix(0,0,0),
                  S = list(eye(ncol(x))),
                  sp = penalty*sd(r),
                  off = 0,
                  p = c(x0.0,numeric(ncol(x)-2*k),x0.1),
                  Ain = rbind(eye(ncol(x)),-eye(ncol(x))),
                  bin = c(norms*lb_means, -norms*ub_means) )
        return(pcls(M)/norms)
      }
    }, cl = cl)
    flog.debug("Updating estimates...")
    for (j in 1:m){
      if (!any(is.na(res[[j]]))){
        mus_hat[j,] <- res[[j]][1:k]
        betas_hat[j,] <- res[[j]][seq(k+1,k+p2,length=p2)]
        gammas_hat[j,c((1*(p1>0)):(p1*k),p1*k+1:k)] <- res[[j]][seq(k+p2+1,k*2+p2+p1*k,length=p1*k+k)]
      }
    }
    flog.debug("Updating X_means_hat...")
    X_means_hat <- tcrossprod(W.inter.0,mus_hat)+tcrossprod(C2,betas_hat)+tcrossprod(cbind(C1_,W.inter.y),gammas_hat)
    if (iter > 1){
      flog.debug("Evaluating convergence...")
      diff <- sum(((mus_hat-mus_hat.prev)**2)/(nrow(mus_hat)*ncol(mus_hat)))
      if (p1) diff <- diff + sum(((gammas_hat-gammas_hat.prev)**2)/(nrow(gammas_hat)*ncol(gammas_hat)))
      if (p2) diff <- diff + sum(((betas_hat-betas_hat.prev)**2)/(nrow(betas_hat)*ncol(betas_hat)))
      flog.debug(paste("Average F norm with estimates from previous iteration: ",diff,sep=""))
      if (diff < globals[["EPSILON"]]) break
    }
    mus_hat.prev <- mus_hat
    betas_hat.prev <- betas_hat
    gammas_hat.prev <- gammas_hat
  }
  return(list("mus_hat" = mus_hat, "betas_hat" = betas_hat, "gammas_hat" = gammas_hat, "X_means_hat" = X_means_hat))
}



## TODO add logs
# Finds the best critical values in terms of the most number of hits while controlling for fdr; considered only a subset of the extreme values in S.perm (controlled by window_size).
get_critical_values <- function(S, S.perm, fdr, window_size = 10){
  flog.debug(sprintf("Running get_critical_values with window_size=%f", window_size))
  m <- length(S)
  num_perms <- nrow(S.perm)
  S.perm.all <- numeric(num_perms*m)
  for (p in 1:num_perms){
    S.perm.all[((p-1)*m+1):(p*m)] <- S.perm[p,]
  }
  S.perm.all <- c(S.perm.all,max(S.perm.all)+1e-8, min(S.perm.all)-1e-8) # consider max(s.perm.all)+epsilon and min(s.perm.all)-epsilon since it may be the case that FDR at level fdr may only be achieved if the threshold is set to a value larger (smaller) than the highest (lowest) permutation value.
  values.sorted <- sort(unique(S.perm.all))
  if (window_size > (length(values.sorted)%/%2)){
    flog.debug(sprintf("get_critical_values: changing window_size from %f to %f", window_size, length(values.sorted)%/%2))
    return(get_critical_values(S, S.perm, fdr, window_size = length(values.sorted)%/%2))
  }
  left_values <- values.sorted[1:window_size]
  right_values <- tail(values.sorted,n=window_size)
  S.trimmed <- S[(S <= max(left_values)) | (S >= min(right_values))]
  if (length(S.trimmed) == 0) return (c(-Inf,Inf))
  S.perm.trimmed <- vector("list", length = num_perms)
  for (p in 1:num_perms){
    S.perm.trimmed[[p]] <-S.perm[p,(S.perm[p,]<=max(left_values)) | (S.perm[p,]>=min(right_values))]
  }
  fdr_array <- matrix(0,window_size,window_size)
  num_hits <- matrix(0,window_size,window_size)
  for (i in 1:window_size){
    for (j in 1:window_size){
      num_hits[i,j] <- sum((S.trimmed <= left_values[i]) | (S.trimmed >= right_values[j]))
      fps_frac <- numeric(p)
      for (p in 1:num_perms){
        fps_frac[p] <- sum((S.perm.trimmed[[p]] <= left_values[i]) | (S.perm.trimmed[[p]] >= right_values[j]))
        if (fps_frac[p]+num_hits[i,j]) fps_frac[p] <- fps_frac[p]/(fps_frac[p]+num_hits[i,j])
      }
      fdr_array[i,j] <- mean(fps_frac)
    }
  }
  if (any(fdr_array <= fdr)){
    control_fdr <- which(fdr_array <= fdr)
    best_position <- (1:(window_size**2))[control_fdr][order(-num_hits[control_fdr])[1]]
    edge_positions <- c(seq(10,window_size**2,window_size),1:10) # these positions in the fdr_array matrix represent at least one value from c(max(left_values),min(right_values)), so we want to try increasing the window in case we can get more hits while controlling for fdr.
    if (is.element(best_position,edge_positions) && (2*window_size<=length(values.sorted)%/%2)){
      return(get_critical_values(S, S.perm, fdr, window_size = min(2*window_size,length(values.sorted)%/%2)))
    }
    # return the best critical values
    return (c(left_values[(best_position-1)%%window_size+1],right_values[(best_position-1)%/%window_size+1]))
  }else{
    return (c(-Inf, Inf))
  }
}


calc_t_statistics <- function(X, W, S, U_hat_sqrt_inv, theta_hat, valid_data, globals, cl, covar_pairs.indices = NULL){
  k <- ncol(theta_hat)
  d <- ncol(S)
  m <- nrow(theta_hat)
  SU_hat_starS_j <- matrix(0,d,d)
  SU_hat_tildeS_j <- matrix(0,d,d)
  lower.tri.indices.diag <- lower.tri(SU_hat_starS_j, diag = TRUE)
  upper.tri.indices <- upper.tri(SU_hat_starS_j, diag = FALSE)
  tstat <- matrix(1,m,k)
  rownames(tstat) <- rownames(theta_hat)
  colnames(tstat) <- colnames(W)
  U_hat_tilde_sqrt_inv <- U_hat_sqrt_inv
  U_hat_tilde_sqrt_inv[which(U_hat_tilde_sqrt_inv < globals[["MIN_SQRT_INV_WEIGHT"]])] = globals[["MIN_SQRT_INV_WEIGHT"]]
  U_hat_tilde <- 1/(U_hat_tilde_sqrt_inv**2)
  U_hat <- U_hat_sqrt_inv**(-2)
  U_hat_star <- U_hat
  indices <- U_hat>=(1/sqrt(globals[["MIN_SQRT_INV_WEIGHT"]]))
  U_hat_star[indices] <- 1/U_hat_star[indices]
  if (!is.null(cl)){
    var_list <- c("calc_S_cov_raw_flat", "S", "valid_data", "U_hat_tilde", "valid_data", "U_hat_star", "lower.tri.indices.diag",
                  "upper.tri.indices", "chol2inv", "theta_hat", "SU_hat_tildeS_j","SU_hat_starS_j", "k")
    clusterExport(cl, varlist = var_list, envir=environment())
  }
  res <- pblapply(1:m,function(j){
    if (is.null(covar_pairs.indices)){
      valid_data_j <- valid_data[,j]
    }else{
      valid_data_j <- valid_data[,covar_pairs.indices[j,1]] & valid_data[,covar_pairs.indices[j,2]]
    }
    S_cov_raw_flat_j <- calc_S_cov_raw_flat(S[valid_data_j,])
    SU_hat_tildeS_j_tri <- S_cov_raw_flat_j %*% U_hat_tilde[valid_data_j,j]
    SU_hat_starS_j_tri <- S_cov_raw_flat_j %*% U_hat_star[valid_data_j,j]
    SU_hat_tildeS_j[lower.tri.indices.diag] <- SU_hat_tildeS_j_tri
    SU_hat_tildeS_j[upper.tri.indices] <- t(SU_hat_tildeS_j)[upper.tri.indices]
    #SU_hat_tildeS_j.inv <- chol2inv(chol(SU_hat_tildeS_j+diag(X_tilde.penalties[j],d,d)))
    SU_hat_tildeS_j.inv <- chol2inv(chol(SU_hat_tildeS_j))
    SU_hat_starS_j[lower.tri.indices.diag] <- SU_hat_starS_j_tri
    SU_hat_starS_j[upper.tri.indices] <- t(SU_hat_starS_j)[upper.tri.indices]
    theta_hat[j,]/sqrt(diag(SU_hat_tildeS_j.inv %*% SU_hat_starS_j %*% SU_hat_tildeS_j.inv)[(k+1):(2*k)])
  }, cl = cl)
  for (j in 1:m){
    tstat[j,] <- res[[j]]
  }
  pvalues <- calc_pvalues(tstat, calc_fdr = FALSE)
  pvalues.adj <- p.adjust(Reshape(pvalues,nrow(pvalues)*ncol(pvalues),1), method="fdr")
  pvalues.adj <- Reshape(pvalues.adj,nrow(pvalues),ncol(pvalues))
  rownames(pvalues) <- rownames(theta_hat)
  rownames(pvalues.adj) <- rownames(theta_hat)
  colnames(pvalues) <- colnames(W)
  colnames(pvalues.adj) <- colnames(W)
  return(list("pvalues" = pvalues, "pvalues.adj" = pvalues.adj, "tstat" = tstat))
}


# each row is S_{1j}*S_{1l},...,S_{nj}*S_{nl} for a pair of features j,l
calc_S_cov_raw_flat <- function(S){
  S_cov_raw_flat <- matrix(0,choose(ncol(S),2)+ncol(S),nrow(S))
  counter <- 1
  for (j in 1:ncol(S)){
    for (l in j:ncol(S)){
      S_cov_raw_flat[counter,] <- S[,j]*S[,l]
      counter <- counter + 1
    }
  }
  return(S_cov_raw_flat)
}


keris.fit_vars <- function(X, X_means_hat, W, W_squared, y, W_.inter.0, W_.inter.y, W_cross_prods, squared_residuals, valid_data, max_iters, penalty, learn_tau, globals, cl){

  m <- ncol(X)
  k <- ncol(W)

  sigmas_squared_hat <- matrix(0, nrow=m, ncol=k)
  sigmas_squared_hat.delta <- matrix(0, nrow=m, ncol=k*!(is.null(y)))
  tau_squared_hat <- matrix(0, nrow=m, ncol=1)
  if (learn_tau) tau_squared_hat <- tau_squared_hat + globals[["MIN_VAR"]]

  for (iter in 1:max_iters){
    V_sqrt_diag <- calc_V_sqrt_diag(X, W_.inter.0, W_.inter.y, sigmas_squared_hat, sigmas_squared_hat.delta, tau_squared_hat, squared_residuals, iter, globals)
    V_sqrt_diag[which(V_sqrt_diag < globals[["MIN_SQRT_INV_WEIGHT"]])] = globals[["MIN_SQRT_INV_WEIGHT"]]
    if (!is.null(cl)) clusterExport(cl, varlist = c("V_sqrt_diag","valid_data","W","sigmas_squared_hat","sigmas_squared_hat.delta",
                                                    "y","tau_squared_hat","W_cross_prods","squared_residuals","k","W_squared",
                                                    "W_.inter.0","W_.inter.y","globals","calc_V_sqrt_diag","eye","pcls"), envir=environment())
    res <- pblapply(1:m,function(j){
      V_jj_sqrt_diag <- V_sqrt_diag[valid_data[,j],j]
      x <- cbind(W_.inter.0[valid_data[,j],],W_.inter.y[valid_data[,j],])/t(repmat(V_jj_sqrt_diag,2*k,1))
      if (learn_tau) x <- cbind(x, (numeric(sum(valid_data[,j])) + 1)/V_jj_sqrt_diag)
      r <- squared_residuals[valid_data[,j],j]/V_jj_sqrt_diag
      # For numeric stability, normalize the design matrix and adjust the final estimates accordingly; note that using the sqrt of the norms turns out to be more stable in pcls.
      norms <- (colSums(x**2))**0.5
      norms <- sqrt(norms)
      x <- x/repmat(norms,nrow(x),1)

      x0 <- c(norms*globals[["MIN_VAR"]]+1) # initial point
      if (learn_tau) x0 <- c(x0, norms[2*k+1]*globals[["MIN_VAR"]]+1)
      Ain <- eye(2*k+learn_tau,2*k+learn_tau)
      bin <- c(norms*globals[["MIN_VAR"]])

      M <- list(y = r,
                w = numeric(length(r))+1,
                X = x,
                C = matrix(0,0,0),
                S = list(eye(2*k+learn_tau,2*k+learn_tau)),
                sp = penalty*sd(r),
                off = 0,
                p = x0,
                Ain = Ain,
                bin = bin )
      return(pcls(M)/norms)
    }, cl = cl)
    # update estimates
    for (j in 1:m){
      sigmas_squared_hat[j,] <- res[[j]][1:k]
      diag_j <- res[[j]][(k+1):(2*k)]
      sigmas_squared_hat.delta[j,] <- diag_j-res[[j]][1:k]
      if (learn_tau) tau_squared_hat[j] <- res[[j]][length(res[[j]])]
    }
    if (iter > 1){
      diff <- sum(((sigmas_squared_hat-sigmas_squared_hat.prev)**2)/nrow(sigmas_squared_hat)*ncol(sigmas_squared_hat))
      if (!is.null(sigmas_squared_hat.delta)){
        diff <- diff + sum(((sigmas_squared_hat.delta-sigmas_squared_hat.delta.prev)**2)/nrow(sigmas_squared_hat.delta)*ncol(sigmas_squared_hat.delta))
      }
      flog.debug(paste("Average F norm with estimates from previous iteration: ",diff,sep=""))
      if (diff < globals[["EPSILON"]]) break
    }
    sigmas_squared_hat.prev <- sigmas_squared_hat
    sigmas_squared_hat.delta.prev <- sigmas_squared_hat.delta
    tau_squared_hat.prev <- tau_squared_hat
  }
  return(list("sigmas_squared_hat" = sigmas_squared_hat, "sigmas_squared_hat.delta" = sigmas_squared_hat.delta, "tau_squared_hat" = tau_squared_hat))
}

keris.fit_covars <- function(X, W, W_squared, covar_pairs.indices, y, W_.inter.0, W_.inter.y, W_cross_prods, squared_residuals.covars, valid_data, sigmas_squared_hat, sigmas_squared_hat.delta, tau_squared_hat, max_iters, penalty, globals, cl){

  m <- ncol(X)
  k <- ncol(W)

  sigmas_squared_hat.covars <- matrix(0,nrow(covar_pairs.indices),k)
  sigmas_squared_hat.covars.delta <- matrix(0,nrow(covar_pairs.indices),k*!is.null(y))

  for (iter in 1:max_iters){
    V_sqrt_diag <- calc_V_sqrt_diag(X, W_.inter.0, W_.inter.y, sigmas_squared_hat.covars, sigmas_squared_hat.covars.delta, tau_squared_hat = 0, squared_residuals = squared_residuals.covars, iter = iter, globals = globals)
    V_sqrt_diag[which(V_sqrt_diag < globals[["MIN_SQRT_INV_WEIGHT"]])] = globals[["MIN_SQRT_INV_WEIGHT"]]
    if (!is.null(cl)) clusterExport(cl, varlist = c("calc_V_sqrt_diag","iter","sigmas_squared_hat.covars","sigmas_squared_hat.covars.delta",
                                                    "W_.inter.y","k","sigmas_squared_hat","sigmas_squared_hat.delta","covar_pairs.indices","X","W","y",
                                                    "tau_squared_hat", "W_cross_prods","squared_residuals.covars","valid_data","is_min_var"), envir=environment())
    res <- pblapply(1:nrow(covar_pairs.indices),function(j){
      diag_j <- sigmas_squared_hat[covar_pairs.indices[j,1],]
      diag_l <- sigmas_squared_hat[covar_pairs.indices[j,2],]
      coefs_to_learn.0 <- (1-is_min_var(diag_j, globals[["MIN_VAR"]])) & (1-is_min_var(diag_l, globals[["MIN_VAR"]]))
      if (is.null(y) & !sum(coefs_to_learn.0)) return(list(numeric(k)))
      bnd.0 <- sqrt(diag_j*diag_l)[coefs_to_learn.0]
      A <- NULL
      b <- NULL
      if (!is.null(y)){
        diag_j.1 <- sigmas_squared_hat[covar_pairs.indices[j,1],]+sigmas_squared_hat.delta[covar_pairs.indices[j,1],]
        diag_l.1 <- sigmas_squared_hat[covar_pairs.indices[j,2],]+sigmas_squared_hat.delta[covar_pairs.indices[j,2],]
        coefs_to_learn.1 <- (!is_min_var(diag_j.1, globals[["MIN_VAR"]])) & (!is_min_var(diag_l.1, globals[["MIN_VAR"]]))
        if(!(sum(coefs_to_learn.0) + sum(coefs_to_learn.1))) return(list(numeric(k),numeric(k)))
        coefs_to_learn <- c(coefs_to_learn.0, coefs_to_learn.1)
        diags_prod <- diag_j.1*diag_l.1
        diags_prod[diags_prod<0] <- 0
        bnd.1 <- sqrt(diags_prod)[coefs_to_learn.1]
        bnd <- c(bnd.0, bnd.1)
      }else{
        coefs_to_learn <- coefs_to_learn.0
        coefs_to_learn.1 <- c()
        bnd <- bnd.0
      }
      valid_data_j <- valid_data[,covar_pairs.indices[j,1]] & valid_data[,covar_pairs.indices[j,2]]
      n_filtered_j <- sum(valid_data_j)
      V_jl_sqrt_diag <- V_sqrt_diag[valid_data_j,j]
      x <- cbind(W_.inter.0[valid_data_j,coefs_to_learn.0,drop=F],W_.inter.y[valid_data_j,coefs_to_learn.1,drop=F])/t(repmat(V_jl_sqrt_diag,sum(coefs_to_learn.0)+sum(coefs_to_learn.1),1))
      r <- squared_residuals.covars[valid_data_j,j]/V_jl_sqrt_diag
      # For numeric stability, normalize the design matrix and adjust the final estimates accordingly; note that using the sqrt of the norms turns out to be more stable in pcls.
      norms <- (colSums(x**2))**0.5
      norms <- sqrt(norms)
      x <- x/repmat(norms,n_filtered_j,1);
      M <- list(y = r,
                w = numeric(length(r))+1,
                X = x,
                C = matrix(0,0,0),
                S = list(eye(ncol(x))),
                sp = penalty*sd(r),
                off = 0,
                p = numeric(length(bnd)),
                Ain = rbind(eye(ncol(x)),-eye(ncol(x))),
                bin = c(-bnd*norms, -bnd*norms) )
      coefs_to_learn.coefs <- pcls(M)/norms
      coefs0 <- numeric(k)
      coefs0[coefs_to_learn.0] <- coefs_to_learn.coefs[1:sum(coefs_to_learn.0)]
      if (is.null(y)){
        return(list(coefs0))
      }else{
        coefs1 <- numeric(k)
        coefs1[coefs_to_learn.1] <- coefs_to_learn.coefs[(sum(coefs_to_learn.0)+1):length(coefs_to_learn.coefs)]
        return(list(coefs0,coefs1))
      }
    }, cl = cl)

    for (j in 1:nrow(covar_pairs.indices)){
      sigmas_squared_hat.covars[j,] <- res[[j]][[1]]
      if(length(res[[j]])>1){
        sigmas_squared_hat.covars.delta[j,] <- res[[j]][[2]]-res[[j]][[1]]
      }
    }
    if (iter > 1){
      diff <- sum(((sigmas_squared_hat.covars-sigmas_squared_hat.covars.prev)**2)/(nrow(sigmas_squared_hat.covars)*ncol(sigmas_squared_hat.covars))) +
        sum(((sigmas_squared_hat.covars.delta-sigmas_squared_hat.covars.delta.prev)**2)/(nrow(sigmas_squared_hat.covars)*ncol(sigmas_squared_hat.covars)))
      flog.debug(paste("Average F norm with estimates from previous iteration: ",diff,sep=""))
      if (diff < globals[["EPSILON"]]) break
    }
    sigmas_squared_hat.covars.prev <- sigmas_squared_hat.covars
    sigmas_squared_hat.covars.delta.prev <- sigmas_squared_hat.covars.delta
  }
  return(list("sigmas_squared_hat.covars" = sigmas_squared_hat.covars, "sigmas_squared_hat.covars.delta" = sigmas_squared_hat.covars.delta))
}


calc_U_sqrt_diag <- function(X, W, W_squared, X_means_hat, iter, globals){
  n <- nrow(W)
  m <- ncol(X)
  k <- ncol(W)
  if (iter == 1){
    U_sqrt_diag <- t(repmat(rowSums(W_squared)**0.5,m,1))
  }else{
    U_sqrt_diag <- abs(X-X_means_hat)
  }
  return(U_sqrt_diag)
}


calc_V_sqrt_diag <- function(X, W_.inter.0, W_.inter.y, sigmas_squared_hat, sigmas_squared_hat.delta, tau_squared_hat, squared_residuals, iter, globals){
  if (length(tau_squared_hat) == 1) assert(tau_squared_hat == 0)

  k <- ncol(W_.inter.0)
  n <- nrow(X)
  m <- nrow(sigmas_squared_hat)
  if (iter == 1){
    V_sqrt_diag <- matrix(1,n,m)
  }else{
    V_sqrt_diag <- squared_residuals - tcrossprod(W_.inter.0,sigmas_squared_hat) - tcrossprod(W_.inter.y,sigmas_squared_hat+sigmas_squared_hat.delta)
  }
  if (length(tau_squared_hat)>1){
    V_sqrt_diag <- V_sqrt_diag-t(repmat(tau_squared_hat,1,n))
  }
  V_sqrt_diag <- abs(V_sqrt_diag)
  return(V_sqrt_diag)
}
