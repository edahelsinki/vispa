#' Compute p-value for a test
#'
#' @param test list(t0, t_null) where 
#' `t0` is the observed test statistic given as a numeric value, a vector of n test statistics, or a (n,m) matrix, and 
#' `t_null` are samples of test statistics from the null, given respectively as a numeric vector, a (n,r) matrix, or a (n,m,r) array
#' @param type "one_sided" or "two_sided" (not implemented)
#' @param na.rm if true, sets NAs equal to -Inf in t_null and t0 (assuming NAs are insignicant. -Inf is insignificant in a one-sided test)
#' @param use_minp if TRUE, minP algorithm is used
#'
#' @return p-value, same dimensions as t0
#' @export
#'
#' @examples
test_to_pvalue <- function(test, type="one_sided", na.rm=T, use_minp=TRUE) {
  t0 = test$t0
  t_null = test$t_null
  if (na.rm) {
    t0[is.na(t0)] = -Inf
    t_null[is.na(t_null)] = -Inf
  }
    
  if (length(t0)==1) {
    p = pvalue(t0, t_null)
  } else if (is.matrix(t_null)) {
    if (use_minp) {
      p = minP(t0, t_null)
    } else {
      p = pvalue_vec(t0, t_null)
    }
  } else if (is.matrix(t0)) {
    if (use_minp) {
      p = minP_array(t0, t_null)
    } else {
      p = pval_array(t0, t_null)
    }
  }
  p
}

#' Empirical p-value (one-sided)
#' 
#' @param t0 observed test statistic
#' @param t_null vector of test statistic samples from the null distribution
#'
#' @return p-value
pvalue <- function(t0, t_null) {
  n <- length(t_null)
  r <- sum(t_null >= t0)
  (r + 1) / (n + 1)
}

#' Empirical p-value (one sided) for multiple test statistics
#'
#' @param t0 vector of `n` observed test statistics
#' @param t_null (n,m) matrix of `m` test statistic samples from the null distribution
#'
#' @return vector of `n` p-values
pvalue_vec <- function(t0, t_null) {
  m <- ncol(t_null) 
  if (!is.null(dim(t0)) && dim(t0)[2]==1)
    t0 = t0[,1]
  r <- rowSums(t_null >= t0)
  (r + 1) / (m + 1)
} 

#' Calculate empirical p-value for array of test statistics
#'
#' @param t0 (n,m) array of observed test statistics
#' @param t_null (n,m,r) array of `r` test statistic samples from the null distribution
#'
#' @return nXm array of empirical pvalues
#' @export
pval_array <- function(t0, t_null) {
  pv = matrix(0, nrow=nrow(t0), ncol=ncol(t0), dimnames=dimnames(t0))
  for (i in 1:dim(t_null)[3]) 
    pv = pv + (t0 <= t_null[,,i])
  pv = (pv + 1) / (dim(t_null)[3] + 1)
  pv
}

#' Compute minP adjusted FWER-controlled p-values
#' 
#' Following algorithm in Box 4 of https://doi.org/10.1007/BF02595811
#'
#' @param x vector of `n` observed test statistics
#' @param Y (n,m) matrix of `m` test statistic samples from the null distribution
#' @return Vector of `n` FWER adjusted p-values
#' @export
minP <- function(x, Y) {
  n <- length(x)
  if(!is.matrix(Y)) stop("maxT: !is.matrix(Y)")
  if(n!=dim(Y)[1]) stop("maxT: n!=dim(Y)[1]")
  m <- dim(Y)[2]
  Z <- matrix(c(x, Y), n, m + 1)
  P <- t(apply(Z, 1, to_pvalues))
  
  p_raw = P[,1]
  i <- order(p_raw, decreasing=T)
  Q = P[i,-1, drop=F]
  if (n > 1) 
    Q = apply(Q, 2, cummin)
  p_adj <- (1 + rowSums(Q <= matrix(rep(p_raw[i], m), n, m))) / (1 + m)
  # Revert to original ordering
  p_adj = p_adj[order(i)]
  names(p_adj) = names(x)
  # Enforce p-value monotonicity
  j = rev(i)
  cummax(p_adj[j])[order(j)]
}

#' Compute p-values using order statistics
#' 
#' Following the procedure in section 4.4.1 of https://doi.org/10.1007/BF02595811
#'
#' @param x vector of test statistics
#' @return vector of pvalues
#' @export
#'
#' @examples
to_pvalues <- function(x) {
  n = length(x)
  i_na = is.na(x)
  if (any(i_na)) {
    p <- rep(NA,n)
    p[!i_na] <- to_pvalues(x[!i_na]) 
  } else {
    r = order(x, decreasing = T)
    x_sort = x[r]
    p_sort = (1:n) / n
    if (n>1) {
      # Ties are handled as follows (x=x_sort, p=p_sort).
      # Normally p_j = j/n, but: 
      # If x1=x2, then p1=p2=2/n
      # If x5=x6=x7, then p5=p6=p7=7/n
      for(j in n:2) {
        if(x_sort[j-1] == x_sort[j]) 
          p_sort[j-1] <- p_sort[j]
      }
    }  
    p = p_sort[order(r)]
  }
  p
}

#' minP adjusted FWER-controlled p-values
#'
#' @param t0 (n,m) array of observed test statistics
#' @param t_null (n,m,r) array of `r` test statistic samples from the null distribution
#'
#' @return (n,m) array of minP adjusted p-values
#' @export
minP_array <- function(t0, t_null) {
  t0_vec = as.vector(t0) # (n,m) array to n*m vector
  t_null_mat = apply(t_null, 3, as.vector) # (n,m,r) array to (n*m,r) matrix
  minP_vec = minP(t0_vec, t_null_mat)
  matrix(minP_vec, nrow=nrow(t0), ncol=ncol(t0), dimnames = dimnames(t0)) # n*m vector to (n,m) array
}
