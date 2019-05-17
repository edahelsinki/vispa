#' Compute test_stat on `data` and its surrogates (generated from `sample_from_null`). 
#'
#' @param data              nXm matrix
#' @param nresample         number of surrogates to sample from `sample_from_null`
#' @param test_stat         function that accepts `data` and outputs a vector of k test statistics
#' @param sample_from_null  function that outputs matrix nXm (surrogate of `data`)
#' @param ...               additional parameters for `test_stat`
#'
#' @return
tester <- function(data, nresample = 1000, test_stat, sample_from_null, ...) { 
  test_stat0 <- function(x) test_stat(x, ...) 
  t0 <- test_stat0(data)
  tt <- replicate(nresample, test_stat0(sample_from_null()))
  if (is.list(tt)) stop("Replicate output is a list. Check output of test_stat (e.g. NULL causes replicate to return a list).")
  if (is.matrix(tt)) tt <- t(tt) # Replicate gives repetitions in columns.
  list(t0 = t0, 
       t = tt, 
       params = list(...))
}

#' Simulated user experiments. 
#'
#' @param P           vector of user expertise values
#' @param N           vector of numbers of test statistics
#' @param nresample   number of surrogates
#' @param beta 
#' @param mu          mean of outlier (x1)
#'
#' @return list of results $res and $input. Output has attribute "filename" based on input.
user_exp <- function(P, N, nresample, beta, mu) {
  k_iter <- function(n, p, b = beta) {
    max(ceiling(log(1-b) / (log((n-1)/n) + log(1-p))), 1)
  }
  
  resu <- matrix(1, nrow = length(P), ncol = length(N), 
                 dimnames = list(paste0("p=", round(P,2)), paste0("n=", N)))
  
  for (i in seq_along(P)) {
    p <- P[i]
    for (j in seq_along(N)) {
      n <- N[j]
      x0 <- c(rnorm(n = 1, mean = mu, sd = 1), 
              rnorm(n = n-1, mean = 0, sd = 1))
      k <- k_iter(n, p, beta)
      
      S <- integer(length = k)
      for (ii in 1:k) {
        success <- sample(0:1, 1, prob = c(1-p, p))
        S[ii] <- if (success) 1 else sample(n, 1) 
      }
      S <- unique(S)
      
      ind_one <- which(S == 1)
      
      if (length(ind_one)>0) {
        surrogates <- replicate(nresample, rnorm(length(S)))
        
        if (length(S) > 1) {
          pv_adj <- minP(x0[S], surrogates)[ind_one]
        } else {
          pv_adj <- (1 + (sum(x0[S] <= surrogates))) / (1 + nresample)
        }
        resu[i,j] <- pv_adj
      }
    }
  }
  out <- list(res = resu, 
              input = list(P = P, N = N, nresample = nresample, beta = beta))
  attr(out, "filename") <- paste0("user-model-experiments-", "P", length(P), "N", length(N), "Nr", nresample, "b", beta, ".Rdata")
  out
}

#' Test statistic: scagnostics
#'
#' @param data    nXm matrix
#' @param pvec      mX2 matrix with projection vectors (vectors as columns)
#'
#' @return        1X9 vector of scagnostics for `data` in view defined by `pvec`
scagnosticss <- function(data, pvec) {
  xy <- data %*% pvec
  as.matrix(scagnostics(xy))[,]  
}

#' Test statistic: number of points inside selected region.
#'
#' @param data    nXm matrix
#' @param region  kX2 matrix with x,y coordinates of a polygon selection
#' @param pvec      mX2 matrix with projection vectors (vectors as columns)
#'
#' @return        number of points inside selected region.
num_points_in_selection <- function(data, pvec, region) {
  length(points_in_selection(data %*% pvec, region))
}

#' Test statistic: Apply function f to intervals [a, b] of x.
#'
#' @param x matrix, array or data.frame nXm
#' @param f function(x, a, b) of interval [a, b]
#' @param L vector of interval sizes. f is applied only on intervals with lengths in L.
#' @param returnNamed if TRUE, return named vector. Name="a_b" for value at [a,b].
#'
#' @return array kXp where k = output size of f, p = # of intervals in x.
#' @examples apply_to_intervals(x = array(1:100), f = function(x, a, b) x[b] - x[a], L = c(4,5)) returns vector of x[b]-x[a] for all intervals [a,b] of length 4 or 5. 
apply_to_intervals <- function(x, f, L=NULL, returnNamed = F) {
  n <- dim(x)[1]
  x_intervals <- combn(1:n, 2)
  if (!is.null(L)) {
    if (L >= n) stop("L>=n: interval size L cannot be longer than dim(x)[1]")
    colDiffs <- function(x) colSums(x * c(-1,1))
    x_intervals <- x_intervals[ , colDiffs(x_intervals) %in% L, drop=F]
  }
  y <- apply(x_intervals, 2, FUN=function(interval) f(x, interval[1], interval[2]))
  if (returnNamed) names(y) <- apply(x_intervals, 2, function(interval) paste0(interval[1], "_", interval[2]))
  y
}

#' Squared exponential kernel (for Gaussian process)
#'
#' @param x 
#' @param y 
#' @param l length scale
kernel_sqexp <- function(x,y=x,l=1) { 
  exp(-((x %o% rep(1,length(y)))-(rep(1,length(x)) %o% y))^2/(2*l^2))
}

#' Null distribution: Gaussian process regression
#'
#' @param x_test 
#' @param x_train 
#' @param y_train 
#' @param kern 
#' @param eps 
GP <- function(x_test, x_train = c(), y_train = c(), 
               kern = kernel_sqexp, eps = .Machine$double.eps^0.5) {
  ch <- function(x) t(chol(x + eps * diag(dim(x)[1])))
  diag0 <- function(x) if(dim(x)[1] > 1) diag(x) else c(x)
  lsq <- function(x,y) lsfit(x, y, intercept = FALSE)$coefficients
  K_ss <- kern(x_test, x_test)
  L_ss <- ch(K_ss)
  mu0 <- rep(0,length(x_test))
  ss0 <- sqrt(diag0(K_ss))
  if(length(x_train)>0) {
    K <- kern(x_train,x_train)
    L <- ch(K)
    K_s <- kern(x_train,x_test)
    L_k <- lsq(L,K_s)
    mu <- c(t(L_k) %*% lsq(L,y_train))
    a <- diag0(K_ss)-colSums(L_k^2)
    if(min(a)< -eps) warning(sprintf("GP: imaginary sd (%g)\n",min(a)))
    ss <- sqrt(pmax(0,a))
    L_sample <- ch(K_ss-t(L_k) %*% L_k)
  } else {
    mu <- mu0
    ss <- ss0
    L_sample <- L_ss
  }
  list(
    v=function() list(mu=mu,sd=ss,x=x_test,mu0=mu0,sd0=ss0),
    prior=function(n=1) {
      L_ss %*% matrix(rnorm(dim(L_ss)[1]*n),dim(L_ss)[1],n)
    },
    post=function(n=1) {
      mu %o% rep(1,n)+L_sample %*% matrix(rnorm(dim(L_sample)[1]*n),dim(L_sample)[1],n)
    }
  )
}

#' Create function to sample days from X. 
#'
#' @param X data.frame. First column is timestamp.
#' @param returnMatrixWithoutTimestamp boolean
#'
#' @return function that subsets one random day from X as a data frame (or matrix)
make_sample_days <- function(X, returnMatrixWithoutTimestamp=F) {
  dates <- format(X[,1], "%Y-%m-%d")
  if (returnMatrixWithoutTimestamp) X <- as.matrix(X[ ,-1])
  function(n=1,d=NULL) {
    if (is.null(d)) d <- sample(dates, n) 
    if (any(!d %in% dates)) stop(paste0("Day ", d[!d %in% dates], " is not in X."))
    i <- which(dates %in% d)
    X[i, ]
  }
}

#' Create function to sample FULL days using sample.days().
#'
#' @param sample.days function that outputs matrix or data.frame
#' @param day_length integer, number of data samples in one day
#'
#' @return function
make_sample_full_days <- function(sample.days, day_length) {
  function() {
    d <- sample.days()
    while (nrow(d) < day_length) d <- sample.days()
    d
  }
}

#' Null distribution: independent column permutations (with possibility to add constraints)
#'
#' @param n nrow(data)
#' @param m ncol(data)
#'
#' @return tile
tiling <- function(n,m,Tm=NULL,Tl=NULL,count=NULL,df=function(x) paste(x,collapse="_")) {
  ## Initial tiling
  if(is.null(Tm)) Tm <- matrix(as.character(rep(1:m,each=n)),nrow=n,ncol=m)
  ## Additionally, we maintain a hash table of tiles (R,C).
  if(is.null(Tl)) {
    Tl <- new.env(hash=TRUE)
    for(j in 1:m) Tl[[as.character(j)]] <- list(R=1:n,C=j)
  }
  if(is.null(count)) count <- m
  
  ## Output unique id number
  newid <- function() {
    count <<- count+1
    as.character(count)
  }
  
  ## Tiling that freezes everything within tile.
  add0tile <- function(R=1:n,C=1:m) for(i in R) addtile(i,C)
  
  ## Add tile to a structure
  addtile <- function(R=1:n,C=1:m) {
    if(length(R)>0 && length(C)>0) {
      S <- new.env(hash=TRUE)
      ids <- NULL
      for(i in R) {
        raw <- Tm[i,C]
        K <- df(raw) # hash of permutation ids
        if(exists(K,envir=S)) {
          ##if(any(raw!=S[[K]]$raw)) stop("addtile: hash collision.") # This is probably over-cautious, but can however catch some programming errors (e.g., NAs)
          S[[K]]$rows <- c(S[[K]]$rows,i)
        } else {
          id <- unique(raw) 
          ids <- union(ids,id)
          S[[K]] <- list(rows=i,id=id,raw=raw) 
        } 
      }
      for(K in ls(S)) {
        Cp <- NULL
        for(id in S[[K]]$id) Cp <- union(Cp,Tl[[id]]$C)
        id <- newid()
        Tm[S[[K]]$rows,Cp] <<- id
        Tl[[id]] <- list(R=S[[K]]$rows,C=Cp)
      }
      for(id in ids) { # remove overlapping rows from existing tiles
        Rnew <- setdiff(Tl[[id]]$R,R)
        if(length(Rnew)>0)
          Tl[[id]]$R <- Rnew
        else 
          rm(list=id,envir=Tl) # if there are no rows left remove empty tile
      }
    }
  }
  
  copy <- function() {
    newTl <- new.env(hash=TRUE)
    for(x in ls(Tl,all.names=TRUE)) assign(x,get(x,Tl),newTl)
    tiling(n,m,Tm=Tm,Tl=newTl,count=count)
  }
  
  permute <- function() {
    p <- matrix(1:(n*m),nrow=n,ncol=m)
    for(key in ls(Tl)) {
      R <- Tl[[key]]$R
      C <- Tl[[key]]$C
      if(length(R)>1) p[R,C] <- p[sample(R),C]
    }
    c(p)
  }
  
  permutedata_matrix <- function(x,p) {
    matrix(x[p],nrow=n,ncol=m,dimnames=dimnames(x))
  }
  
  permutedata_df <- function(x,p) {
    pp <- p-n*rep(0:(m-1),each=n)
    if(any(pp<1 | n<pp)) stop("permutedata_df: permutations are not within column.")
    as.data.frame(lapply(1:m,function(j) x[pp[((j-1)*n+1):(j*n)],j]),
                  row.names=rownames(x),col.names=colnames(x))
  }
  
  permutedata_list <- function(x,p) {
    pp <- p-n*rep(0:(m-1),each=n)
    if(any(pp<1 | n<pp)) stop("permutedata_list: permutations are not within column.")
    a <- lapply(1:m,function(j) x[[j]][pp[((j-1)*n+1):(j*n)]])
    names(a) <- names(x)
    a
  }
  
  permutedata <- function(x,p=permute(),nmin=n) {
    if(is.matrix(x)) {
      if(nmin>n) {
        k <- 1+(nmin-1)%/%dim(x)[1]
        y <- matrix(NA,k*dim(x)[1],dim(x)[2],dimnames=dimnames(x))
        for(i in 1:k) y[((i-1)*dim(x)[1]+1):(i*dim(x)[1]),] <- permutedata_matrix(x,p=permute())
        y
      } else {
        permutedata_matrix(x,p)
      }
    } else if(is.data.frame(x)) {
      permutedata_df(x,p)
    } else {
      permutedata_list(x,p)
    }
  }
  
  
  cov_local <- function(x) {
    if(!is.matrix(x)) stop("needs matrix")
    x <- scale(x,center=TRUE,scale=FALSE)
    r <- matrix(NA,m,m,dimnames=list(colnames(x),colnames(x)))
    x0 <- matrix(NA,n,m)
    for(i in 1:m) x0[,i] <- ave(x[,i],Tm[,i])
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        a <- x[,i]
        b <- x[,j]
        idx <- which(Tm[,i]!=Tm[,j])
        if(length(idx)>0) {
          a[idx] <- x0[idx,i]
          b[idx] <- x0[idx,j]
        }
        r[i,j] <- r[j,i] <- mean(a*b)
      }
    }
    diag(r) <- apply(x,2,function(y) mean(y^2))
    r
  }
  
  cor_local <- function(x) {
    cov2cor(cov_local(x))
  }
  
  list(add0tile=add0tile,
       addtile=addtile,
       copy=copy,
       permute=permute,
       permutedata=permutedata,
       cov=cov_local,
       cor=cor_local,
       status=function() list(n=n,m=m,Tm=Tm,Tl=Tl,count=count))
}