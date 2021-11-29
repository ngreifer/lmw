plot.summary.lmw <- function(x, abs = TRUE, var.order = "data", threshold = NULL, layout = "vertical", ...) {
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  weighted <- !is.null(x[["sum.weighted"]])
  un <- !is.null(x[["sum.un"]])
  standard.sum <- if (un) x[["sum.un"]] else x[["sum.weighted"]]

  if (!all(c("TSMD Treated", "TSMD Control") %in% colnames(standard.sum))) {
    stop("Not appropriate for unstandardized summary. Run summary() with the standardize = TRUE option, and then plot.", call. = FALSE)
  }

  if (un) {
    sd.all <- list(`1` = x[["sum.un"]][,"TSMD Treated"],
                   `0` = x[["sum.un"]][,"TSMD Control"])
  }
  if (weighted) {
    sd.weighted <- list(`1` = x[["sum.weighted"]][,"TSMD Treated"],
                        `0` = x[["sum.weighted"]][,"TSMD Control"])
  }

  var.names <- rownames(standard.sum)

  var.order <- match_arg(var.order, c("data", "alphabetical"))

  layout <- match_arg(layout, c("vertical", "horizontal"))

  if (abs) {
    if (un) {
      for (i in c("1", "0"))
        sd.all[[i]] <- abs(sd.all[[i]])
    }
    if (weighted) {
      for (i in c("1", "0"))
        sd.weighted[[i]] <- abs(sd.weighted[[i]])
    }
    xlab <- list(`1` = "TASMD (treated group)",
                 `0` = "TASMD (control group)")
  }
  else {
    xlab <- list(`1` = "TSMD (treated group)",
                 `0` = "TSMD (control group)")
  }

  ord <- switch(var.order,
                "data" = rev(seq_along(var.names)),
                "alphabetical" = order(var.names, decreasing = TRUE))

  minx <- min(if (un) c(sd.all[["1"]], sd.all[["0"]]), if (weighted) c(sd.weighted[["1"]], sd.weighted[["0"]]), if (abs) 0 else -.01)
  maxx <- max(if (un) c(sd.all[["1"]], sd.all[["0"]]), if (weighted) c(sd.weighted[["1"]], sd.weighted[["0"]]), .01)

  xlim <- c(minx, 1.75*maxx - minx)

  if (layout == "vertical") layout(rbind(1,2))
  else if (layout == "horizontal") layout(cbind(1,2))
  par(mar=c(3.75, 6.5, 1.25, 0.5),
      mgp=c(1.5, 0.5, 0))

  for (i in c("1", "0")) {
    dotchart(if (un) sd.all[[i]][ord] else sd.weighted[[i]][ord],
             labels = var.names[ord], xlab = xlab[[i]],
             xlim = xlim, lcolor = NA,
             cex=.8,
             bg = NA, color = NA, ...)
    abline(v = 0)

    if (un) {
      points(x = sd.all[[i]][ord], y = seq_along(sd.all[[i]]),
             pch = 4, cex=1.5)
    }
    if (weighted) {
      points(x = sd.weighted[[i]][ord], y = seq_along(sd.weighted[[i]]),
             pch = 20, cex=1.5)
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

    title(ylab = "Covariate", mgp = c(5.25, 0.5, 0), cex = .8)

      legend(x=maxx + .1*(maxx-minx),
             y=length(var.names)+0.65,
             legend=c("Before weighting",
                      "After weighting")[c(un, weighted)],
             pch=c(4, 20)[c(un, weighted)],
             bty="n",
             cex = .8)
  }

  invisible(x)
}
