plot.summary.lmw <- function(x, stats, abs = TRUE, var.order = "data", threshold = NULL, layout = "vertical", ...) {
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  un <- !is.null(x[["bal.un"]])
  base.weighted <- !is.null(x[["bal.base.weighted"]])
  weighted <- !is.null(x[["bal.weighted"]])

  if (!un && !base.weighted && !weighted) {
    stop("plot() can only be used on summary.lmw objects when stat = \"balance\" was specified in the call to summary().", call. = FALSE)
  }

  standard.sum <- if (un) x[["bal.un"]] else x[["bal.weighted"]]

  if (!any(startsWith(colnames(standard.sum), "TSMD"))) {
    stop("Not appropriate for unstandardized summary. Run summary() with the standardize = TRUE option, and then plot.", call. = FALSE)
  }

  if (missing(stats)) {
    stats <- colnames(standard.sum)[startsWith(colnames(standard.sum), "TSMD")]
  }
  stats <- match_arg(stats, colnames(standard.sum), several.ok = TRUE)

  if (un) {
    stats.un <- as.data.frame(x[["bal.un"]][,stats, drop = FALSE])
  }
  if (base.weighted) {
    stats.base.weighted <- as.data.frame(x[["bal.base.weighted"]][,stats, drop = FALSE])
  }
  if (weighted) {
    stats.weighted <- as.data.frame(x[["bal.weighted"]][,stats, drop = FALSE])
  }

  var.names <- rownames(standard.sum)

  var.order <- match_arg(var.order, c("data", "alphabetical", "unadjusted"[un]))

  layout <- match_arg(layout, c("vertical", "horizontal"))

  if (abs) {
    if (un) {
      for (i in seq_along(stats.un))
        stats.un[[i]] <- abs(stats.un[[i]])
    }
    if (base.weighted) {
      for (i in seq_along(stats.base.weighted))
        stats.base.weighted[[i]] <- abs(stats.base.weighted[[i]])
    }
    if (weighted) {
      for (i in seq_along(stats.weighted))
        stats.weighted[[i]] <- abs(stats.weighted[[i]])
    }
  }

  xlab <- vapply(stats, rename_summary_stat, character(1L), abs = abs)

  ord <- switch(var.order,
                "data" = rev(seq_along(var.names)),
                "alphabetical" = order(var.names, decreasing = TRUE),
                "unadjusted" = order(stats.un[[1]]))

  minx <- min(if (un) unlist(stats.un), if (base.weighted) unlist(stats.base.weighted), if (weighted) unlist(stats.weighted), if (abs) 0 else -.01)
  maxx <- max(if (un) unlist(stats.un), if (base.weighted) unlist(stats.base.weighted), if (weighted) unlist(stats.weighted), .01)

  xlim <- c(minx, 1.75*maxx - minx)

  if (layout == "vertical") layout(do.call("rbind", as.list(seq_along(stats))))
  else if (layout == "horizontal") layout(do.call("cbind", as.list(seq_along(stats))))
  par(mar=c(3.75, 6.5, 1.25, 0.5),
      mgp=c(1.5, 0.5, 0))

  if (base.weighted) {
    if (identical(x$base.weights.origin, "MatchIt")) {
      legend.text <- sprintf(c("Before matching"[un],
                               "After matching",
                               "After matching + \n%s regression"[weighted]), x$method)
    }
    else if (identical(x$base.weights.origin, "WeightIt")) {
      legend.text <- sprintf(c("Before weighting"[un],
                               "After weighting",
                               "After weighting + \n%s regression"[weighted]), x$method)
    }
    else {
      legend.text <- sprintf(c("Before base weighting"[un],
                               "After base weighting",
                               "After base weighting + \n%s regression"[weighted]), x$method)
    }
  }
  else {
    legend.text <- sprintf(c("Before regression"[un],
                             "After %s regression"[weighted]), x$method)
  }


  for (i in seq_along(stats)) {
    dotchart(if (un) stats.un[[i]][ord] else stats.weighted[[i]][ord],
             labels = var.names[ord], xlab = xlab[[i]],
             xlim = xlim, lcolor = NA,
             cex=.8,
             bg = NA, color = NA, ...)
    abline(v = 0)

    if (un) {
      points(x = stats.un[[i]][ord], y = seq_along(stats.un[[i]]),
             pch = 4, cex=1.5)
    }
    if (base.weighted) {
      points(x = stats.base.weighted[[i]][ord], y = seq_along(stats.base.weighted[[i]]),
             pch = 1, cex=1.5)
    }
    if (weighted) {
      points(x = stats.weighted[[i]][ord], y = seq_along(stats.weighted[[i]]),
             pch = 19, cex=1.5)
    }

    if (!is.null(threshold)) {
      if (abs) {
        abline(v = threshold, lty = seq_along(threshold))
      }
      else {
        abline(v = threshold, lty = seq_along(threshold))
        abline(v = -threshold, lty = seq_along(threshold))
      }
    }

    # title(ylab = "Covariate", mgp = c(5.25, 0.5, 0), cex = .8)

    legend(x=maxx + .1*(maxx-minx),
           y=length(var.names)+0.65,
           legend=legend.text,
           pch=c(4, 1, 19)[c(un, base.weighted, weighted)],
           bty="n",
           cex = .8)
  }

  invisible(x)
}

rename_summary_stat <- function(x, abs = FALSE) {
  splitted <- strsplit(x, " ", fixed = TRUE)[[1]]
  stat <- splitted[1]
  if (endsWith(stat, "KS")) stat <- paste(stat, "statistic")
  else if (endsWith(stat, "SMD") && abs) stat <- sub("SMD", "ASMD", stat, fixed = TRUE)

  if (length(splitted) > 1) {
    group <- paste(splitted[-1], collapse = " ")
    if (tolower(group) %in% c("treated", "control")) {
      group <- paste(tolower(group), "group")
    }
    sprintf("%s (%s)", stat, group)
  }
  else stat
}
