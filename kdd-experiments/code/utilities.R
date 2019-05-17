check_if_packages_are_present <- function(required_packages) {
  missing_packages  <- setdiff(required_packages, installed.packages())
  
  if (length(missing_packages) > 0) {
    cat("Error! Cannot continue due to missing packages\n")
    cat("Please install the following packages:\n")
    for (i in missing_packages)
      cat(paste0("\t- ", i, "\n"))
    stop("Install packages and try again.")
  }
}

#' Returns rownames of data points inside a polygon selection. Utilizes sp library. 
#'
#' @param data      nX2 matrix with xy coordinates
#' @param selection kX2 matrix of xy coordinates for polygon vertices 
#'
#' @return
points_in_selection <- function(data, selection) {
  data <- data.frame(data)
  selection <- data.frame(selection)
  idx <- sp::point.in.polygon(point.x = data[ , 1], point.y = data[ , 2], 
                              pol.x = selection[ , 1], pol.y = selection[ , 2])
  if (length(idx) < 2) warning(length(idx), " points in selection.")
  as.integer(rownames(data)[as.logical(idx)])
}


#' Empirical p-value (one-sided)
#'
#' @param test: output from tester. test$t0: observed test statistic. test$t: vector of replicate test statistics.
#'
#' @return p-value
pvalue <- function(test, na.rm=T) {
  tt <- test$t
  t0 <- test$t0
  
  if (na.rm) tt <- tt[!is.na(tt)]
  n <- length(tt)
  
  r <- sum(tt >= t0)
  (r + 1) / (n + 1)
}

#' Empirical p-value for multiple test statistics. One sided.
#'
#' @param test tester2 output. 
#' test$t0: vector of m test statistics on the observed data. 
#' test$t: nXm matrix of m test statistics on n replicate data.
#'
#' @return vector of m p-values (one for each test statistic)
pvalue2 <- function(test, na.rm=T) {
  n <- nrow(test$t)  
  r <- rowSums(t(test$t) >= test$t0, na.rm = na.rm)
  (r + 1) / (n + 1)
} 

#' Iteration adjusted p-value (Equation 1)
#'
#' @param pval  p-value
#' @param t     iteration step
#'
pv_adj_iter <- function(pval, t) {
  min(1, 2^t * pval)
}

#' Read and preprocess air quality data
#'
#' @return data.frame
readAQ <- function() {
  AQ <- read.csv("data/AirQualityUCI.csv",
                 header=TRUE,sep=";",dec=",")
  AQ <- AQ[!is.na(AQ[,"CO.GT."]),]
  TI <- mapply(function(x,y) as.POSIXct(strptime(sprintf("%s %s",x,y),
                                                 format="%d/%m/%Y %H.%M.%S")),
               AQ[,"Date"],AQ[,"Time"])
  for(i in 3:15) AQ[AQ[,i]==-200,i] <- NA
  AQ <- AQ[!is.na(TI),]
  TI <- TI[!is.na(TI)]
  data.frame(time=mapply(function(x,y)
    as.POSIXct(strptime(sprintf("%s %s",x,y),format="%d/%m/%Y %H.%M.%S"),origin="1970-01-01"),
    AQ[,"Date"],AQ[,"Time"]),
    AQ[,3:15])
}

#' Formatted time series plot
#'
#' @param x 
#' @param y 
#' @param x_test 
#' @param f 
#' @param mar 
#' @param ... 
plottime <- function(x, y, x_test=c(), f=c(), mar=c(2, 0, 0, 0) + 0.1, bty="n",...) {
  par(mar=mar)
  plot(
    range(c(x_test, x), na.rm = T),
    range(c(f, y), na.rm = T),
    type = "n", bty = bty, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
    ...
  )
  axis.POSIXct(side=1,
               at=as.POSIXct(c(seq(from=min(x),to=max(x),by=4*3600), max(x)), origin="1970-01-01"), 
               format = "%H:%M")
  if (!is.null(x) & !is.null(f)) {
    for(i in 1:dim(f)[2]) lines(x_test,f[,i],col="gray")
  }
  
  lines(x,y, lwd=3)
}

#' Compute polygon of convex hull of points. 
#'
#' @param xy  matrix nX2
#' @param R   vector of m indices (for a subset of xy)
#'
#' @return    matrix kX2 with k points on convex hull of xy[R,]
convhull_subset <- function(xy, R) {
  xyR <- xy[R,]
  xyR[chull(xyR),]
}

#' Formatted scatterplot
#'
#' @param x nX2 matrix or data frame
plott <- function(x, xlab="", ylab="", bg="darkgray", pch=20, xaxt="n", yaxt="n", bty="n", cex = 0.3, 
                  mar=c(1,1,1,1)+0.1, asp=1, tick.num=c(4,4), ...) {
  par(mar=mar)
  plot(x, xlab=xlab, ylab=ylab, bg=bg, pch=pch, xaxt=xaxt, yaxt=yaxt, cex=cex, asp=asp, bty=bty, ...)
  axis(side=1,
       at=seq(from=min(x[,1]), to=max(x[,1]), length.out=tick.num[1]),
       tck=-0.02
  )
  axis(side=2,
       at=seq(from=min(x[,2]), to=max(x[,2]), length.out=tick.num[2]),
       tck=-0.02
  )
}

#' Formatted polygon
#'
#' @param pol 
plot_pol <- function(pol, lwd=2, border=rgb(0.8, 0, 0, 0.5), ...) {
  polygon(pol, lwd=lwd, border=border, ...)
}

#' Save a figure as PDF
#'
#' @param filename 
#' @param fig function that outputs plot
#' @param width 
#' @param height 
savePDF <- function(filename, fig, width, height) {
  pdf(file = filename, width = width, height = height)
  fig
  invisible(dev.off())
}

#' Return table with test statistics and p-values (unadjusted, Bonferroni, minP)
#'
#' @param Test 
#' @param includeBaselines 
#'
#' @return
results <- function(Test, includeBaselines = F) {
  
  if (includeBaselines) { 
    round(t(rbind(t0 = Test$t0,
                  pval.adj = minP(Test$t0, t(Test$t)),
                  pval.unadj = pvalue2(Test),
                  pval.bonf = p.adjust(pvalue2(Test), method="bonferroni")
    )), 2)
  } else {
    round(t(rbind(t0 = Test$t0,
                  pval.adj = minP(Test$t0, t(Test$t)))), 2)    
  }
}
