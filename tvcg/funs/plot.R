library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)

#' Formatted time series plot
#'
#' @param timestamp,y numeric vectors for x and y axes 
#' @param surrogates (n,r) matrix of r surrogates where n=length(y)
#'
#' @return
#' @export
#'
#' @examples
#' plot_ts_gg(
#' airquality[airquality$Month==5, ]$Day, 
#' airquality[airquality$Month==5, ]$Temp, 
#' replicate(10, airquality[airquality$Month==5, ]$Temp[sample(sum(airquality$Month==5))]))
plot_ts_gg <- function(timestamp, y, surrogates = NULL) {
  df = data.frame(timestamp=timestamp, 
                  y = y)
  gg_surrogates = NULL
  if (!is.null(surrogates)) {
    df_surrogates = surrogates %>% 
      as.data.frame %>% 
      mutate(timestamp=df$timestamp) %>% 
      pivot_longer(-timestamp)
    gg_surrogates = geom_line(data=df_surrogates, aes(timestamp, value, group=name), col="grey", alpha=0.5)
  }
  ggplot(df) + 
    gg_surrogates + 
    geom_line(aes(timestamp, y), lwd=1) + 
    scale_x_datetime(date_labels = "%b %d\n%H:%M") + 
    theme_classic() + 
    theme(axis.line.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(y="", x="") 
}

#' Formatted scatterplot
#'
#' @param mat matrix
#' @param surrogates list of surrogates as matrices, same size as mat
#'
#' @return
#' @export
#'
#' @examples
#' scatterplot_gg(project(german_viz, pvec_t1))
#' scatterplot_gg(project(german_viz, pvec_t1), lapply(surrogates_perm[1:10], function(X) X[[1]] %*% pvec_t1))
scatterplot_gg <- function(mat, surrogates=NULL, polygon=NULL) {
  gg_surrogates = NULL
  if (!is.null(surrogates)) {
    surrogates_df =  tibble(surr = lapply(surrogates, function(x) as.data.frame(x)), 
                            surr_id = paste0("Surrogate_", 1:length(surrogates))) %>% 
      unnest_longer(surr) %>% 
      do.call(data.frame, .)
    gg_surrogates = geom_point(data=surrogates_df, aes_string(colnames(surrogates_df)[1], colnames(surrogates_df)[2]), col="grey", size=0.3)
  }
  
  gg_region = NULL
  if (!is.null(polygon)) {
    polygon_df = as.data.frame(polygon)
    polygon_df = rbind(polygon_df, head(polygon_df, 1)) # To connect last and first point
    gg_region = geom_path(data=polygon_df, 
                          aes_string(x = colnames(polygon_df)[1], 
                                     y = colnames(polygon_df)[2]), 
                          col="darkred", 
                          alpha=0.5,
                          lwd=1)
  }

  df = as.data.frame(mat)
  x_var = colnames(df)[1]
  y_var = colnames(df)[2]
  ggplot(df) + theme_classic() +
    labs(x=x_var, y=y_var) + 
    gg_surrogates + 
    geom_point(aes(x=.data[[x_var]], 
                   y=.data[[y_var]]), 
               size=0.8) + 
    theme(axis.title = element_text(size=8))+
    gg_region
}

#' Banana plot
#'
#' @param data data frame 
#' @param day string, date in format YYYY-MM-DD
#'
#' @return plot
#' @importFrom ggplot2 ggplot aes theme_light scale_fill_viridis_c theme margin labs scale_y_log10 geom_raster
#' @export
#'
#' @examples
plot_banana <- function(data, day = min(data$timestamp)) {
  stopifnot("timestamp" %in% names(data))
  j_particles = grepl("^d\\d+", names(data)) # Columns with particle concentration
  dat = data %>% 
    dplyr::select(timestamp, which(j_particles)) %>% 
    dplyr::filter(lubridate::date(timestamp) == lubridate::date(day)) %>% 
    tidyr::gather(particle, concentration, starts_with("d")) %>% 
    dplyr::mutate(particle_size_nm = as.factor(particle_colname_to_nanometers(particle))) %>% 
    dplyr::select(timestamp, particle_size_nm, concentration)

  y_breaks = unique(dat$particle_size_nm)
  y_breaks = y_breaks[floor(seq(1, length(y_breaks), length.out = 10))]
  y_breaks = 10^(0:5)
  y_labels = y_breaks
  z_breaks = 10^(1:4)
  
  ggplot(dat, aes(x = timestamp)) + 
    theme_light() + 
    scale_fill_viridis_c(trans = "log1p", breaks = z_breaks, na.value = "transparent", 
                         limits = range(z_breaks), oob = scales::squish) + 
    scale_y_discrete(labels = y_labels, breaks = y_breaks) + 
    scale_x_datetime(date_labels = "%b %d\n%H:%M") + 
    theme(legend.margin = margin(), 
          legend.title = element_text(size=8), 
          rect = element_blank(), 
          axis.title.y = element_text(size=8)) + 
    labs(x = "", y = "Particle Size (nm)", fill = "Concentration \n(cm-3)") + 
    geom_raster(aes(y = particle_size_nm, fill = concentration))
}

#' Plot line segments over a gg plot
#'
#' @param gg_heatmap ggplot object, where gg_heatmap$data[[1]] contains the 
#' x-axis values, and gg_heatmap$data[[2]] contains the y-axis values
#' @param lines (n, 4) matrix of n lines to be drawn
#' @param alpha,col,... graphical parameters to pass to geom_segment()
#'
#' @return
#' @export
#'
#' @examples
plot_lines <- function(gg, lines, alpha = 0.2, color = "white", ...) {
  x_values = gg$data[[1]]
  y_values = unique(gg$data[[2]])
  gg + 
    geom_segment(data = as.data.frame(lines), 
                 aes(x=x_values[x0], xend=x_values[x1], 
                     y=y_values[y0], yend=y_values[y1]), 
                 alpha = alpha, 
                 color = color, ...)
}

#' Line plot of p-values
#' 
#' Shows raw, minP-adjusted, Bonferroni adjusted p-values
#'
#' @param test see [tester()]
#'
#' @export
#'
#' @examples
plot_pvalues <- function(test, show_x_labels=FALSE) {
  x_axis = 1:length(test$t0)
  x_labels = if (is.null(names(test$t0))) x_axis else names(test$t0)
  gg_x_axis = NULL
  gg_x_axis_rotate_labels = NULL
  if (show_x_labels) {
    gg_x_axis = scale_x_continuous(breaks=x_axis, labels=x_labels)  
    gg_x_axis_rotate_labels = theme(axis.text.x = element_text(angle=90, size=8))   
  }
  
  p_raw = test_to_pvalue(test, use_minp = F)
  p_minp = test_to_pvalue(test, use_minp = T)
  df_pv = data.frame(test_stat = x_axis, 
                     raw = p_raw, 
                     bonf = p.adjust(p_raw), 
                     minp = p_minp) 
  legend_labels = c("bonf"="Bonferroni", "raw"="None", "minp"="minP")
  df_pv %>% 
    tidyr::pivot_longer(-test_stat, "method") %>% 
    ggplot() + theme_classic() + 
    geom_line(aes(x=test_stat, 
                  y=value, 
                  linetype=method, 
                  color=method, 
                  size=method)) + 
    labs(x= "Test statistic", 
         y = "p-value", 
         color="p-value adjustment", 
         linetype="p-value adjustment", 
         size="p-value adjustment") + 
    scale_linetype_manual(labels=legend_labels, 
                          values = c("bonf"=3, "raw"=2, "minp"=1)) + 
    scale_size_manual(labels=legend_labels, 
                      values = c("bonf"=0.3, "raw"=0.3, "minp"=0.6)) + 
    scale_color_manual(labels=legend_labels, 
                       values = c("bonf"="black", "raw"="grey", "minp"="black")) + 
    ylim(c(0, 1)) + 
    gg_x_axis + 
    gg_x_axis_rotate_labels
}

#' Plot test statistics and null distributions
#'
#' @param t_obs vector (n) observed test statistics
#' @param t_null matrix (n,m) of m samples from null distribution for each test 
#' statistic. Only the first 100 samples are plotted for efficiency.
#'
#' @export
#'
#' @examples
plot_test_stat = function(t_obs, t_null, jitter=F, show_x_labels=F) {
  x_axis = 1:length(t_obs)
  x_labels = if (is.null(names(t_obs))) x_axis else names(t_obs)
  t_null = t_null[,1:min(100, ncol(t_null))]
  gg_x_axis = NULL
  gg_x_axis_rotate_labels = NULL
  if (show_x_labels) {
    gg_x_axis = scale_x_continuous(breaks=x_axis, labels=x_labels)  
    gg_x_axis_rotate_labels = theme(axis.text.x = element_text(angle=90, size=8))   
  }
    
  df = t_null %>% 
    as.data.frame %>% 
    dplyr::mutate(test_stat = x_axis) %>% 
    tidyr::pivot_longer(-test_stat)
  if (jitter) {
    position = position_jitter(height=0, width=0.1)
  } else {
    position = position_identity()
  }
  ggplot() + 
    theme_classic() + 
    labs(x="Test Statistic", y="Value") + 
    geom_point(aes(test_stat, value), data=df, col="grey", cex=0.3, alpha=0.6, 
               position=position) +
    geom_point(aes(x=1:length(t_obs), y=t_obs)) + 
    gg_x_axis + 
    gg_x_axis_rotate_labels
}

#' Plot spectrum
#'
#' @param f vector of n frequencies
#' @param y vector length n. Spectrum values for f.
#' @param surrogates nXm matrix 
#'
#' @return ggplot
#' @export
#'
#' @examples
plot_spectrum <- function(f, y, surrogates = NULL) {
  df = data.frame(f = f, 
                  y = y)
  gg_surrogates = NULL
  if (!is.null(surrogates)) {
    df_surrogates = surrogates %>% 
      as.data.frame() %>% 
      mutate(f=df$f) %>% 
      pivot_longer(-f)
    gg_surrogates = geom_line(data=df_surrogates, aes(f, value, group=name), col="grey", alpha=0.5)
  }
  ggplot(df) + 
    gg_surrogates + 
    geom_line(aes(f, y), lwd=1) + 
    scale_y_log10() + 
    scale_x_log10() + 
    theme_classic() + 
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(y="", x="") 
}

#' Histogram
#'
#' @param obs 
#' @param null_sample 
#' @param breaks 
#'
#' @return
#' @export
#'
#' @examples
plot_hist = function(obs, null_samples, breaks) {
  # df = data.frame(id=1:length(obs), obs=obs)
  df_samples = data.frame(id=1:length(obs), null=null_samples)
  df_samples_long = pivot_longer(df_samples, -id)
  ggplot(data.frame(obs=obs)) + 
    # geom_histogram(aes(value, after_stat(density), fill=name, color=name), 
    #                alpha=0.7, position = "identity", breaks=breaks) +
    geom_histogram(aes(value, after_stat(density), group=name), data=df_samples_long, 
                   alpha=0.1, fill="lightgrey", breaks=breaks, position="identity") + 
    geom_histogram(aes(obs, after_stat(density)), color="black", fill=NA, breaks=breaks) +
    theme_minimal() + 
    theme(legend.position = "none") + 
    labs(x="", y="")
}

# Utilities ####
#' Project matrix X onto vector v
#' 
#' Wrapper for `X %*% v` that adds column names 
#' based on largest absolute weights in v.
#'
#' @param X matrix (n,m)
#' @param v vector (m,k)
#'
#' @return matrix (n,k)
#' @examples 
#' project(german_viz, pvec_t1)
project = function(X, v, topm=3) {
  stopifnot(ncol(X) == nrow(v))
  X_proj = X %*% v
  new_colnames = character(ncol(v))
  for (j in 1:ncol(v)) {
    i_highest = rank(-abs(v[,j]), ties.method = "first") <= topm
    new_colnames[j]  = paste0(paste0(round(v[,j],2), "*", colnames(X))[i_highest], collapse=" ")
  }  
  colnames(X_proj) = new_colnames
  X_proj
}
