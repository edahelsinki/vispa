#' Compute maxT adjusted FWER-controlled p-values
#' 
#' @param x vector of length n, test statistic in the original data
#' @param y matrix nXm, columns are the test statistic in m samples from null hypothesis
#' @return Vector of length n, maxT adjusted mid-P p-values for the test statistic
#' 
#' @export
maxT <- function(x,y) {
  n <- length(x)
  if(!is.matrix(y)) stop("maxT: !is.matrix(y)")
  if(n!=dim(y)[1]) stop("maxT: n!=dim(y)[1]")
  m <- dim(y)[2]
  
  i <- order(x)
  j <- rev(i)
  xm <- matrix(rep(x,m),n,m)
  ym <- if(n>1) apply(y[i,,drop=FALSE],2,cummax)[order(i),,drop=FALSE] else y
  p <- (1+rowSums(xm<=ym))/(1+m)
  names(p) <- names(x)
  cummax(p[j])[order(j)]
}

#' Compute minP adjusted FWER-controlled p-values
#' 
#' @param x vector of length n, test statistic in the original data
#' @param y matrix nXm, columns are the test statistic in m samples from null hypothesis
#' @return Vector of length n, maxT adjusted mid-P p-values for the test statistic
#' 
#' @export
minP <- function(x,y) {
  n <- length(x)
  if(!is.matrix(y)) stop("maxT: !is.matrix(y)")
  if(n!=dim(y)[1]) stop("maxT: n!=dim(y)[1]")
  m <- dim(y)[2]
  
  z <- matrix(c(x,y),n,m+1)
  p <- t(apply(z,1,topvalues))
  rownames(p) <- names(x)
  maxT(-p[,1],-p[,-1,drop=FALSE])
}

#' Compute unadjusted p-values
#' 
#' @param x vector of length n, test statistic in the original data
#' @param y matrix nXm, columns are the test statistic in m samples from null hypothesis
#' @return Vector of length n, unadjusted mid-P p-values for the test statistic
#' 
#' @export
rawP <- function(x,y) {
  n <- length(x)
  if(!is.matrix(y)) stop("rawP: !is.matrix(y)")
  if(n!=dim(y)[1]) stop("rawP: n!=dim(y)[1]")
  m <- dim(y)[2]
  
  p <- (1+rowSums(x<=y))/(1+m)
  names(p) <- names(x)
  p
}

#' Convert vector of test statistic to the vector of p-values
#' 
#' @param x vector of length n, test statistic
#' @return Vector of length n, p-values.
#' 
#' @export
topvalues <- function(x) {
  n <- length(x)
  idx <- is.na(x)
  if(any(idx)) {
    p <- rep(NA,n)
    p[!idx] <- topvalues(x[!idx])
  } else {
    i <- order(x)
    p <- (n+1-order(i))/n
    if(n>1) {
      ## If there are equally sized test statistic...
      for(j in 2:n) {
        if(x[i[j]]==x[i[j-1]]) p[i[j]] <- p[i[j-1]]
      }
    }
  }
  p
}
