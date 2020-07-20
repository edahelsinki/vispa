#' Contingency table permutation test
#'
#' @param X data frame. Columns X[,C[1]] and X[,C[2]] are factors.
#' @param null_tiling tiling of same size as X
#' @param nresample number of samples to draw from tiling
#' @param R integer vector. Subset of data rows. 
#' @param C character or integer vector of length 2. Two columns in X to compute contingency table.
#' @return array tab where tab[,,i] is the contingency table from the i-th sample
#' @export
tester_conttab <- function(X, null_tiling, nresample, R=NULL, C=NULL) {
  if (is.null(R)) 
    R = 1:nrow(X)
  if (is.null(C))
    C = 1:2
  stopifnot(length(C) == 2, (is.integer(C) | is.character(C)), length(R) > 1, is.integer(R))
  if (is.character(C)) {
    c1 = which(names(X) == C[1])
    c2 = which(names(X) == C[2])
  } else {
    c1 = C[1]
    c2 = C[2]
  }
  x1 = X[, c1]
  x2 = X[, c2]
  if (!is.factor(x1) | !is.factor(x2))
    stop("X[ , C[1]] and X[ , C[2]] should be factors.")
  null_tables = sample_contingency_table(size=nresample, tt=null_tiling, x1=x1, x2=x2, c1=c1, c2=c2, R=R)
  if (is.integer(C))
    C = paste0("C", C)
  names(dimnames(null_tables)) <- c(C[1], C[2], "Sample")
  
  list(t0 = table(x1[R], x2[R], dnn = C), t_null = null_tables)
}

#' Sample contingency tables
#'
#' @author Kai Puolamaki
#' @param size number of samples to obtain
#' @param tt tiling structure
#' @param x1 column 1 as a vector of factors
#' @param x2 column 2 as a vector of factors
#' @param c1 position of column 1 in tt
#' @param c2 position of column 2 in tt
#' @param R integer vector. Subset of rows in x1, x2 in which to compute contingency table
#' @param use_only_R if TRUE, then permutes only within subset R by adding a tile on rows outside R.
#' @return Array tab where tab[,,i] is the contingency table from the ith sample.
#' @export
sample_contingency_table <- function(size = 1, tt, x1, x2, c1, c2, R = 1:length(x1), use_only_R=F) {
  tab <- array(0,
               dim = c(length(levels(x1)), length(levels(x2)), size),
               dimnames = list(levels(x1), levels(x2), 1:size))
  stopifnot(is.numeric(R), length(R) > 1)
  tt2 = tt$copy()
  if (use_only_R)
    tt2$addtile(R = setdiff(1:length(x1), R))
  perm <- tt2$twocolumn(C = c(c1, c2), size = size)
  for (i in 1:size)
    tab[,,i] <- table(x1[perm[R, i, 1]], x2[perm[R, i, 2]])
  tab
}

#' Convert contingency table to vector
#'
#' @param tab contingency table
#'
#' @return numeric vector
#' @export
#'
#' @examples
#' tab = table(mtcars[,c("cyl", "gear")])
#' tab
#' conttab_to_vec(tab)
conttab_to_vec <- function(tab) {
  if (!is.null(dimnames(tab)) & !is.null(names(dimnames(tab)))) {
    dimnames_df = expand.grid(dimnames(tab))
    vec_names = paste0(paste0(names(dimnames_df)[1], "=", dimnames_df[[1]]), ",",
                       paste0(names(dimnames_df)[2], "=", dimnames_df[[2]]))
  } else {
    vec_names = NULL
  }
  vec = as.vector(tab)
  names(vec) = vec_names
  vec
}

#' Wrapper for converting contingency tables to vectors inside tests
#' Converts test with contingency tables as test statistics to test with vectors as test statistics.
#'
#' @param test list(t0, t_null) where t0 is a `(n,m)` table and t_null is an `(n,m,r)` array of tables. See [test_conttab()].
#' @return test: list(t0, t_null) where t0 is a `n*m` vector and t_null is a `(n*m,r)` matrix
#' @export
conttab_to_vec_test <- function(test) {
  test_vec = list()
  test_vec$t0 = conttab_to_vec(test$t0)
  test_vec$t_null = apply(test$t_null, 3, as.vector)
  dimnames(test_vec$t_null) <- list(names(test_vec$t0), NULL)
  test_vec
}

