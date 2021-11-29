plot.lmw <- function(x, type = "weights", ...) {
  type <- match_arg(type, c("weights", "influence", "extrapolation"))

  if (type == "weights") {
    weights_plot(x, ...)
  }
  else if (type == "extrapolation") {
    extrapolation_plot1(x, ...)
  }

  invisible(x)
}

weights_plot <- function(x, rug = TRUE, mean = TRUE, ...) {
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  t <- x$treat
  w <- x$weights

  layout(rbind(1,2))

  par(mar = c(2.75, 3, 1.75, 1),
      mgp = c(1.5, 0.5, 0))

  for (i in c(1,0)) {
    wi <- w[t==i]

    n <- length(wi)

    mw <- sum(wi)/n

    ess <- ESS(wi)
    ess_un <- if (is.null(x$s.weights)) n else ESS(x$s.weights)

    dens <- density(wi)

    xlim <- range(dens$x)

    ylim <- c(0, max(dens$y))

    plot(dens,
         ylab = "Density",
         xlab = "Weights",
         yaxt = "n",
         cex.lab = 0.8,
         cex.axis = 0.8,
         main = "")
    title(main = if (i == 1) "Treated Weights"
          else "Control Weights",
          cex.main = 0.9,
          font.main = 1,
          line = 0.5)
    # title(main = "Dispersion and effective sample size", outer=TRUE, adj=0.7125, line=-1, cex.main=1, font.main=1)
    # axis(side = 2,
    #      at = ylim,
    #      cex.axis = 0.8)
    if (rug) {
      rug(wi,
          col = "black",
          ticksize = 0.1,
          lwd = 0.3,
          side = 1)
    }
    if (mean) {
      rug(mw, col = "red", ticksize = 0.3, lwd = 1)
    }

    legend("topright",
           legend = paste0("ESS = ", round(ess, 1), "\n",
                           "N = ", round(ess_un, 1)),
           bty = "n",
           cex = 0.8)
  }
}

# #Influence plot needs outcome regression model; maybe use a different method
# influence_plot <- function(x, ...) {
#   .pardefault <- par(no.readonly = TRUE)
#   on.exit(par(.pardefault))
#
#   par(mar=c(2.5, 3, 2, 2),
#       mgp=c(1.5, 1, 0))
#
#   # Sample influence curves under MRI and URI
#   a_SIC_MRI = rep(0,N)
#   a_SIC_URI = rep(0,N)
#
#   for(i in 1:N)
#   {
#     a_SIC_MRI[i] = (nt*Z[i] + (1-Z[i])*nc)*abs(e_MRI[i]*w_trial_MRI2[i])/(1-h_MRI[i])
#     a_SIC_URI[i] = (N-1)*abs(e_URI[i]*w_trial_URI2[i])/(1-H_tilde[i,i])
#   }
#
#
#   ## Scaled SIC
#   a_SIC_MRI_std = a_SIC_MRI/max(a_SIC_MRI)
#   a_SIC_URI_std = a_SIC_URI/max(a_SIC_URI)
#
#   matplot(1:N,
#           cbind(-a_SIC_MRI_std,a_SIC_URI_std),
#           type = c("h","h"),
#           lty = c("solid","solid"),
#           xlab = "Index",
#           ylab = "Scaled SIC",
#           cex.lab = 0.8,
#           col = grey(c(0,0.6)),
#           axes = FALSE)
#   title(main="Influence",
#         outer=TRUE,
#         adj=0.555,
#         line=-35.5, cex.main=1, font.main=1)
#   axis(side = 1,
#        at = c(1,seq(200,2500,200),N),
#        cex.axis = 0.8,
#        mgp=c(0, 0.5, 0),
#        cex.axis=0.7)
#   axis(side = 2,
#        ylim = range(a_SIC_URI_std),
#        at = seq(0,1,0.2),
#        cex.axis = 0.7,
#        pos=-50)
#   axis(side = 4, ylim = -range(a_SIC_MRI_std),
#        at = -seq(0,1,0.2),
#        labels = seq(0,1,0.2),
#        cex.axis = 0.7,
#        pos=2750)
#   legend(x=2200,
#          y=1,
#          legend=c("URI", "MRI"),
#          lty=c(1,1),
#          pch=c(NA, NA),
#          col=grey(c(0.6,0)),
#          cex = 0.6,
#          bty="n")
#   par(xpd=NA)
#   abline(v=-500, lwd=0.1)
#   abline(h=1.2, lwd=0.1)
# }

extrapolation_plot1 <- function(x, var, ...) {
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  # set.seed(1234)

  v <- x$covs[[var]]
  t <- x$treat
  w <- x$weights

  layout(rbind(1,2))

  par(mar = c(2.75, 1, 1.75, 1),
      mgp = c(1.5, 0.5, 0))

  xlim <- range(v)

  if (length(unique(v)) <= 2) v <- v + rnorm(length(v), 0,.01)

  for (i in c(1,0)) {

    plot(x = v[t==i],
         y = (1 + 1*(w[t==i] >= 0)) + rnorm(sum(t==i), 0,.05),
         pch = 20,
         col = ifelse(w[t==i] >= 0, adjustcolor("black", alpha.f = .4), adjustcolor("red", alpha.f = .5)),
         ylim = c(0.5, 2.5),
         xlim = xlim,
         xlab = var,
         yaxt = "n",
         cex = 300*abs(w[t==i])/sum(abs(w[t==i])),
         ylab = "",
         main = "",
         frame.plot = TRUE,
         cex.axis = 0.8,
         cex.lab = 0.8)
    title(main = if (i == 1) "Treated Group"
          else "Control Group",
          cex.main = 0.9,
          font.main = 1,
          line = 0.5)
    abline(v = mean_w(v, w, t==i))
    points(x = switch(x$estimand, "ATT" = mean_w(v, x$s.weights, t==1),
                      "ATC" = mean_w(v, x$s.weights, t==0),
                      "ATE" = mean_w(v, x$s.weights)),
           y = 1.5, pch = 4, cex = 4, lwd = 0.5)
  }

}

extrapolation_plot2 <- function(x, var, ...) {
  .pardefault <- par(no.readonly = TRUE)
  on.exit(par(.pardefault))

  # set.seed(1234)

  v <- x$covs[[var]]
  t <- x$treat
  w <- x$weights

  par(mar = c(2.75, 3, 1.75, 1),
      mgp = c(1.5, 0.5, 0))

  xlim <- range(v)

  if (length(unique(v)) <= 2) v <- v + rnorm(length(v), 0,.01)

  cex <- numeric(length(w))
  for (i in c(1,0)) cex[t==i] <- 300*abs(w[t==i])/sum(abs(w[t==i]))
  plot(x = v,
       y = .5 + 2*(t==1) + 1*(w >= 0), #+ rnorm(sum(t==i), 0,.05),
       pch = 20,
       col = ifelse(w >= 0, adjustcolor("black", alpha.f = .4), adjustcolor("red", alpha.f = .5)),
       ylim = c(0, 4),
       xlim = xlim,
       xlab = var,
       yaxt = "n",
       cex = cex,
       ylab = "",
       main = "",
       frame.plot = TRUE,
       cex.axis = 0.8,
       cex.lab = 0.8)
  # title(main = if (i == 1) "Treatment Group"
  #       else "Control Group",
  #       cex.main = 0.9,
  #       font.main = 1,
  #       line = 0.5)
  # abline(v = mean_w(v, w, t==i))
  abline(h = 2)
  segments(x0 = mean_w(v, w, t==0), y0 = -1, y1 = 2)
  segments(x0 = mean_w(v, w, t==1), y0 = 2, y1 = 5)
  points(x = rep(switch(x$estimand, "ATT" = mean_w(v, x$s.weights, t==1),
                    "ATC" = mean_w(v, x$s.weights, t==0),
                    "ATE" = mean_w(v, x$s.weights)), 2),
         y = c(1, 3), pch = 4, cex = 4, lwd = 0.5)
  axis(2, at = c(1,3), labels = c("Control", "Treated"),
       tick = FALSE)

}
