instant_pkgs <- function(pkgs) { 
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss)
  }
  
  if (length(pkgs_miss) == 0) {
    message("\n ...Packages were already installed!\n")
  }
  
  # install packages not already loaded:
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss)
  }
  
  # load packages not already loaded:
  attached <- search()
  attached_pkgs <- attached[grepl("package", attached)]
  need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
  
  if (length(need_to_attach) > 0) {
    for (i in 1:length(need_to_attach)) require(need_to_attach[i], character.only = TRUE)
  }
  
  if (length(need_to_attach) == 0) {
    message("\n ...Packages were already loaded!\n")
  }
}

import_data <- function(obj, types){
  for (i in 1:length(obj)){
    FUN = switch(types[i],
                 numeric = as.numeric, 
                 integer = as.integer,
                 factor = as.factor
    )
    
    obj[,i] = FUN(obj[,i])
  }
  obj
}

pres <- function(raw_resids, trts, merMod) {
  pres <- rep(0, length(raw_resids))
  for (i in 1:length(raw_resids)) {
    res <- raw_resids[i]
    trt <- trts[i]
    j <- switch(trt, 
                control = 3,
                human = 2,
                raptor = 1)
    
    var <- summary(modId)$varcor[[j]][1] + sigma(merMod)^2
    pres[i] <- res / sqrt(var)
  }
  pres
}


pairs.panels <- function (x, science = c("biol", "phys"), smooth = TRUE, scale = TRUE, 
                          density = TRUE, ellipses = FALSE, digits = 2, method = "spearman", pch = 20, lm = FALSE, 
                          cor = TRUE, jiggle = FALSE, factor = 2, hist.col = "lightblue3", show.points = TRUE, 
                          breaks = "Sturges", cex.cor = 10, ...) 
{
  science <- match.arg(science)
  
  "panel.hist.density" <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    h <- hist(x, breaks = breaks, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = hist.col)
    if (density) {
      tryd <- try(d <- density(x, na.rm = TRUE, bw = "nrd", 
                               adjust = 1.2), silent = TRUE)
      if (class(tryd) != "try-error") {
        d$y <- d$y/max(d$y)
        lines(d)
      }
    }
  }
  
  "panel.cor" <- function(x, y, digits = 2, prefix = "", ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = "pairwise", method = method)
    # Add background fill
    colors <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", 
                                 "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", 
                                 "#D6604D", "#B2182B", "#67001F"))(200)
    col.fill <- colors[findInterval(r, seq(-1, 1, length.out = 200))]
    # Scale box size
    r.range <- seq(0, 1, length.out = 200)
    if (science == "biol") {
      fancy.size <- 4/3 - 4/3*exp(-log(2)/0.5 * r.range)
    } else {
      fancy.size <- scale_vec(exp(4*r.range) - 1, c(0, 1))
    }
    adj <- fancy.size[findInterval(abs(r), r.range)] * 0.5
    polygon(x = c(0.5-adj, 0.5+adj, 0.5+adj, 0.5-adj),
            y = c(0.5+adj, 0.5+adj, 0.5-adj, 0.5-adj),
            border = col.fill, col = col.fill)
    txt <- format(c(round(r, digits), 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    cex <- cex.cor * 0.8/strwidth(txt)
    if (scale) {
      cex1 <- cex * adj * 2
      if (cex1 < 0.25) 
        cex1 <- 0.25
      text(0.5, 0.5, txt, cex = cex1)
    } else {
      text(0.5, 0.5, txt, cex = cex)
    }
  }
  "panel.smoother" <- function(x, y, pch = par("pch"), col.smooth = "red", 
                               span = 2/3, iter = 3, ...) {
    xm <- mean(x, na.rm = TRUE)
    ym <- mean(y, na.rm = TRUE)
    xs <- sd(x, na.rm = TRUE)
    ys <- sd(y, na.rm = TRUE)
    r = cor(x, y, use = "pairwise", method = method)
    if (jiggle) {
      x <- jitter(x, factor = factor)
      y <- jitter(y, factor = factor)
    }
    if (show.points) 
      points(x, y, pch = pch, ...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
    panel.ellipse1(xm, ym, xs, ys, r, col.smooth = col.smooth, 
                   ...)
  }
  "panel.smoother.no.noellipse" <- function(x, y, pch = par("pch"), 
                                            col.smooth = "red", span = 2/3, iter = 3, ...) {
    xm <- mean(x, na.rm = TRUE)
    ym <- mean(y, na.rm = TRUE)
    xs <- sd(x, na.rm = TRUE)
    ys <- sd(y, na.rm = TRUE)
    r = cor(x, y, use = "pairwise", method = method)
    if (jiggle) {
      x <- jitter(x, factor = factor)
      y <- jitter(y, factor = factor)
    }
    if (show.points) 
      points(x, y, pch = pch, ...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
  }
  "panel.lm" <- function(x, y, pch = par("pch"), col.lm = "red", 
                         ...) {
    ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin, xmin), max(ymax, xmax))
    xlim <- ylim
    if (jiggle) {
      x <- jitter(x, factor = factor)
      y <- jitter(y, factor = factor)
    }
    points(x, y, pch = pch, ylim = ylim, xlim = xlim, ...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      abline(lm(y[ok] ~ x[ok]), col = col.lm, ...)
  }
  "panel.lm.ellipse" <- function(x, y, pch = par("pch"), col.lm = "red", 
                                 ...) {
    ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin, xmin), max(ymax, xmax))
    xlim <- ylim
    if (jiggle) {
      x <- jitter(x, factor = factor)
      y <- jitter(y, factor = factor)
    }
    points(x, y, pch = pch, ylim = ylim, xlim = xlim, ...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      abline(lm(y[ok] ~ x[ok]), col = col.lm, ...)
    xm <- mean(x, na.rm = TRUE)
    ym <- mean(y, na.rm = TRUE)
    xs <- sd(x, na.rm = TRUE)
    ys <- sd(y, na.rm = TRUE)
    r = cor(x, y, use = "pairwise", method = method)
    panel.ellipse1(xm, ym, xs, ys, r, col.smooth = col.lm, 
                   ...)
  }
  "panel.ellipse1" <- function(x = 0, y = 0, xs = 1, ys = 1, 
                               r = 0, col.smooth, add = TRUE, segments = 51, ...) {
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    if (!is.na(r)) {
      if (abs(r) > 0) 
        theta <- sign(r)/sqrt(2)
      else theta = 1/sqrt(2)
      shape <- diag(c(sqrt(1 + r), sqrt(1 - r))) %*% matrix(c(theta, 
                                                              theta, -theta, theta), ncol = 2, byrow = TRUE)
      ellipse <- unit.circle %*% shape
      ellipse[, 1] <- ellipse[, 1] * xs + x
      ellipse[, 2] <- ellipse[, 2] * ys + y
      points(x, y, pch = 19, col = col.smooth, cex = 1.5)
      lines(ellipse, ...)
    }
  }
  "panel.ellipse" <- function(x, y, pch = par("pch"), col.smooth = "red", 
                              ...) {
    segments = 51
    xm <- mean(x, na.rm = TRUE)
    ym <- mean(y, na.rm = TRUE)
    xs <- sd(x, na.rm = TRUE)
    ys <- sd(y, na.rm = TRUE)
    r = cor(x, y, use = "pairwise", method = method)
    if (jiggle) {
      x <- jitter(x, factor = factor)
      y <- jitter(y, factor = factor)
    }
    if (show.points) 
      points(x, y, pch = pch, ...)
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    if (!is.na(r)) {
      if (abs(r) > 0) 
        theta <- sign(r)/sqrt(2)
      else theta = 1/sqrt(2)
      shape <- diag(c(sqrt(1 + r), sqrt(1 - r))) %*% matrix(c(theta, 
                                                              theta, -theta, theta), ncol = 2, byrow = TRUE)
      ellipse <- unit.circle %*% shape
      ellipse[, 1] <- ellipse[, 1] * xs + xm
      ellipse[, 2] <- ellipse[, 2] * ys + ym
      points(xm, ym, pch = 19, col = col.smooth, cex = 1.5)
      if (ellipses) 
        lines(ellipse, ...)
    }
  }
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  if (missing(cex.cor)) 
    cex.cor <- 1
  for (i in 1:ncol(x)) {
    if (is.character(x[[i]])) {
      x[[i]] <- as.numeric(as.factor(x[[i]]))
      colnames(x)[i] <- paste(colnames(x)[i], "*", sep = "")
    }
  }
  if (!lm) {
    if (smooth) {
      if (ellipses) {
        pairs(x, diag.panel = panel.hist.density, lower.panel = panel.cor, 
              upper.panel = panel.smoother, pch = pch, ...)
      }
      else {
        pairs(x, diag.panel = panel.hist.density, lower.panel = panel.cor, 
              upper.panel = panel.smoother.no.noellipse, 
              pch = pch, ...)
      }
    }
    else {
      pairs(x, diag.panel = panel.hist.density, lower.panel = panel.cor, 
            upper.panel = panel.ellipse, pch = pch, ...)
    }
  }
  else {
    if (!cor) {
      if (ellipses) {
        pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm.ellipse, 
              lower.panel = panel.lm.ellipse, pch = pch, 
              ...)
      }
      else {
        pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm, 
              lower.panel = panel.lm, pch = pch, ...)
      }
    }
    else {
      if (ellipses) {
        pairs(x, diag.panel = panel.hist.density, lower.panel = panel.cor, 
              upper.panel = panel.lm.ellipse, pch = pch, 
              ...)
      }
      else {
        pairs(x, diag.panel = panel.hist.density, lower.panel = panel.cor, 
              upper.panel = panel.lm, pch = pch, ...)
      }
    }
  }
}