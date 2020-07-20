#' Null distribution: Gaussian process regression
#' @author Kai Puolamäki
#' @param x_test 
#' @param x_train 
#' @param y_train 
#' @param kern 
#' @param eps 
#' @examples
#' l = 6
#' gp = train_gp(x_test = 1:100, x_train = 50:55, y_train=rep(5,6), kern = function(x,y) kernel_sqexp(x, y, l))
#' plot(gp$prior(), type="l")
#' for (i in 1:10)
#'   lines(gp$prior(), col="grey")
#' plot(gp$post(), type="l")
#' for (i in 1:10)
#'   lines(gp$post(), col="grey")
train_gp <- function(x_test, x_train = c(), y_train = c(), 
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

#' Squared exponential kernel (for Gaussian process)
#' @author Kai Puolamäki
#' @param x 
#' @param y 
#' @param l length scale
kernel_sqexp <- function(x,y=x,l=1) { 
  exp(-((x %o% rep(1,length(y)))-(rep(1,length(x)) %o% y))^2/(2*l^2))
}

#' Create function to sample days from X
#'
#' @param X data frame. First column is POSIX timestamp.
#' @param day_length if not null, then only samples days of length day_length
#' @param return_matrix if TRUE returns matrix (without timestamp column)
#'
#' @return function that subsets one random day from X as a data frame (or matrix)
#' @export
make_sample_days <- function(X, day_length=NULL, return_matrix=F) {
  dates <- format(X[,1], "%Y-%m-%d")
  if (return_matrix) 
    X <- as.matrix(X[ ,-1, drop=F])
  
  sample_days = function(n=1, d=NULL) {
    if (is.null(d)) 
      d <- sample(dates, n)
    if (any(!d %in% dates)) 
      stop(paste0("Day ", d[!d %in% dates], " is not in X."))
    i <- which(dates %in% d)
    X[i, , drop=F]
  }
  
  if (!is.null(day_length)) {
    function() {
      d <- sample_days()
      while (nrow(d) != day_length) 
        d <- sample_days()
      d
    } 
  } else {
    sample_days
  }
}

#' Null distribution: independent column permutations (based on tile constraints)
#'
#' Source: Kai Puolamäki, 2019. CORAND library. \url{https://github.com/edahelsinki/corand}
#' As described in: Kai Puolamäki, Emilia Oikarinen, Andreas Henelius, 2019. Guided Visual Exploration of Relations in Data Sets. \url{https://arxiv.org/abs/1905.02515}
#'
#' @param n nrow(data)
#' @param m ncol(data)
#'
#' @return tile
#' @export
tiling <- function(n,m,Tm=NULL,Tl=NULL,count=NULL,df=function(x) paste(x,collapse="_")) {
  ## Initial tiling
  ## Tiling matrix as in Alg. 1
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
  
  permute <- function(keys=ls(Tl)) {
    p <- matrix(1:(n*m),nrow=n,ncol=m)
    for(key in keys) {
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
  
  permutedata <- function(x,keys=ls(Tl),p=permute(keys=keys),nmin=n) {
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
  
  findtiles <- function(R=1:n,C=1:m) {
    ## find all tiles which interserct R and C
    ls(Tl)[
      vapply(ls(Tl),
             function(key) {
               ( length(intersect(R,Tl[[key]]$R))>0 &&
                   length(intersect(C,Tl[[key]]$C))>0 )
             },
             TRUE)
    ]
  }
  
  twocolumn <- function(R=1:n,C=1:2,size=1) {
    if(length(C)!=2) stop("twocolumn: length(R)!=2")
    ## Find the relevant tiles for columns 1 and 2
    k1 <- unique(Tm[R,C[1]])
    k2 <- unique(Tm[R,C[2]])
    ## We need only to permute tiles which have at least two rows
    k1 <- k1[vapply(k1,function(key) length(Tl[[key]]$R)>1,TRUE)]
    k2 <- k2[vapply(k2,function(key) length(Tl[[key]]$R)>1,TRUE)]
    ## Separate tiles that occur in only 1 and 2 or in both
    k12 <- intersect(k1,k2)
    k1 <- setdiff(k1,k12)
    k2 <- setdiff(k2,k12)
    ## Find the rows corresponding to the tiles
    l1 <- lapply(k1,function(key) Tl[[key]]$R)
    l2 <- lapply(k2,function(key) Tl[[key]]$R)
    l12 <- lapply(k12,function(key) Tl[[key]]$R)
    ## Create size random permutations for the two columns
    ## x[,a,b] is the ath permutation for column C[b]
    x <- array(rep(1:n,2*size),dim=c(n,size,2))
    for(i in l1) {
      x[i,,1] <- replicate(size,sample(i))
    }
    for(i in l2) {
      x[i,,2] <- replicate(size,sample(i))
    }
    for(i in l12) {
      x[i,,] <- rep(replicate(size,sample(i)),2)
    }
    x
  }
  
  
  list(add0tile=add0tile,
       addtile=addtile,
       copy=copy,
       permute=permute,
       permutedata=permutedata,
       cov=cov_local,
       cor=cor_local,
       findtiles=findtiles,
       twocolumn=twocolumn,
       status=function() list(n=n,m=m,Tm=Tm,Tl=Tl,count=count))
}

#' Wrapper for faster use of tiling with data frames and matrices
#'
#' @param X data frame or matrix
#'
#' @return see [tiling()]
#' @export
#'
#' @examples
#' tt = tiling_df(iris)
#' head(tt$status()$Tm)
tiling_df <- function(X) {
  n = nrow(X)
  m = ncol(X)
  Tm = matrix(as.character(rep(1:m, each = n)), nrow = n, ncol = m)
  colnames(Tm) = colnames(X)
  rownames(Tm) = rownames(X)
  tiling(n = n, m = m, Tm = Tm)
}
