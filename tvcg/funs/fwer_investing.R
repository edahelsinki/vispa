#' FWER investing
#'
#' @param a_t 
#' @param a_b 
#' @param R_t 
#'
#' @return
#' @export
#'
#' @examples
fwer_investing = function(a_t, a_b, R_t) {
  (a_b - a_t * (1-R_t)) / (1 - a_t * (1-R_t))
}


#' Apply the FWER alpha investing algorithm on a sequence of p-values
#'
#' @param pv vector of pvalues
#' @param a_fwer desired upper bound for FWER
#' @param update_alpha_budget function(a_t, a_b, R_t) that updates the alpha budget 
#' for the next step, depending on R_t = p_t <= a_t
#' @param update_alpha_investment function(a_b) that determines a_t (how much budget is used at each step)
#'
#' @return data frame
#' @export
#'
#' @examples
do_fwer_investing <- function(pv, a_fwer=0.1, 
                              update_alpha_budget=fwer_investing, 
                              update_alpha_investment=function(a_b) a_b/2) {
  n_steps = length(pv)
  a_b = numeric(n_steps)
  a_t = numeric(n_steps)
  a_b[1] = a_fwer
  for (i in 1:n_steps) {
    a_t[i] = update_alpha_investment(a_b[i])
    if (i < n_steps)
      a_b[i+1] = update_alpha_budget(a_t[i], a_b[i], pv[i] <= a_t[i])
  }
  data.frame(pv=pv, a_t=a_t, a_b=a_b)
}


#' Apply a list of alpha investing functions on a data frame with experiments
#'
#' @param exp_df data frame from [do_experiment_num_size()]
#' @param list_of_funs 
#'
#' @return
#' @export
#'
#' @examples
apply_alpha_investing_funs <- function(exp_df, 
                                       list_of_funs = 
                                         list(fwer_investing = do_fwer_investing)) {
  for (j in seq_along(list_of_funs)) {
    fun_name = names(list_of_funs)[j]
    fun = list_of_funs[[j]]
    tmp = vector("list", nrow(exp_df))
    for (i in 1:nrow(exp_df)) {
      pv_min = apply(exp_df$experiment[[i]]$pv, 1, min)
      tmp[[i]] = fun(pv_min)
    }  
    exp_df[[fun_name]] = I(tmp)
  }
  exp_df
}

