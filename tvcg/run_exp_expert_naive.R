source("funs/pvalue.R", chdir = TRUE)
overwrite=FALSE
filename = "results/exp_expert_naive.rds"
if (file.exists(filename) & overwrite==FALSE) {
  exp_expert_naive = readRDS(filename)
} else {
  #' Simulated user experiment
  #'
  #' @param K           vector of user expertise values
  #' @param N           vector of numbers of test statistics
  #' @param nresample   number of surrogates
  #' @param beta 
  #' @param mu          mean of outlier (x1)
  #' @param n_rep       number of repetitions over which to average
  #'
  #' @return list of results $res and $input. Output has attribute "filename" based on input.
  do_user_experiment <- function(K, N, nresample, beta, mu, n_rep, verbose=FALSE) {
    m_star <- function(n, k, b = beta) {
      max(ceiling(log(1-b) / (log((n-1)/n) + log(1-k))), 1)
    }

    result <- array(1, dim=c(length(K), length(N), n_rep), 
                    dimnames = list(paste0("k=", K), 
                                    paste0("n=", N), 
                                    NULL))
    for (r in 1:n_rep) {
      if (r%%10==0 || r==1)
        cat(paste0(r, "/", n_rep), "\n")
      for (i in seq_along(K)) {
        k <- K[i]
        for (j in seq_along(N)) {
          n <- N[j]
          x0 <- c(rnorm(n = 1, mean = mu, sd = 1), 
                  rnorm(n = n-1, mean = 0, sd = 1))
          m <- m_star(n, k, beta)
          
          S <- integer(length = m)
          for (ii in 1:m) {
            success <- sample(0:1, 1, prob = c(1-k, k))
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
            result[i, j, r] <- pv_adj
          }
        }
      }
    }
    out <- list(res = result, 
                input = list(K = K, N = N, nresample = nresample, beta = beta))
    out
  }
  
  set.seed(2020)

  exp_expert_naive = do_user_experiment(
    K = seq(from = 0, to = 1, length = 101),
    N = c(1, 2, 5, 10, 25, 50, 100, 250, 500, 1000),
    nresample = 1000,
    mu = 3,
    beta = 0.5,
    n_rep = 1000, 
    verbose = TRUE
  )
  
  saveRDS(exp_expert_naive, filename)
}
