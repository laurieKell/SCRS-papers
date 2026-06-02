# Multivariate log-normal (MVLN) sampling from SS3 covariance and derived quantities.
# Adapted from mvln/R/ssMVLNFunction.R (Mosqueira / WUR); no FLCore dependency.

cor2cov <- function(cor, var) {
  sd <- sqrt(var)
  outer(sd, sd) * cor
}

sshat <- function(hat) {
  names(hat) <- tolower(names(hat))
  y <- rbind(
    hat[grep("^Recr_", hat$label), , drop = FALSE],
    hat[grep("^Bratio_", hat$label), , drop = FALSE],
    hat[grep("^F_", hat$label), , drop = FALSE]
  )
  y <- y[substr(y$label, 1, 2) %in% c("F_", "Br", "Re"), , drop = FALSE]
  y <- y[!y$label %in% c(
    "Recr_Initial", "Recr_unfished", "Recr_Virgin", "Ret_Catch_MSY"
  ), , drop = FALSE]
  y$cv2 <- (y$stddev / y$value)^2
  y$var <- log(1 + (y$stddev / y$value)^2)
  y <- y[, c("label", "value", "stddev", "cv2", "var"), drop = FALSE]
  names(y) <- c("label", "hat", "std", "cv2", "var")
  y
}

sscor <- function(covar) {
  names(covar) <- tolower(names(covar))
  flag <- unique(sort(c(
    grep("^Recr_", covar$label.i),
    grep("^Recr_", covar$label.j),
    grep("^Bratio_", covar$label.i),
    grep("^Bratio_", covar$label.j),
    grep("^F_", covar$label.i),
    grep("^F_", covar$label.j)
  )))
  if (!length(flag)) {
    return(NULL)
  }
  cor <- covar[flag, c("label.i", "label.j", "corr"), drop = FALSE]
  keep <- substr(cor$label.i, 1, 2) %in% c("Re", "F_", "Br") &
    substr(cor$label.j, 1, 2) %in% c("Re", "F_", "Br")
  cor <- cor[keep, , drop = FALSE]
  cor <- rbind(
    cor,
    data.frame(
      label.i = cor$label.j,
      label.j = cor$label.i,
      corr = cor$corr,
      stringsAsFactors = FALSE
    )
  )
  cor <- cor[!(is.na(cor$label.i) | is.na(cor$label.j)), , drop = FALSE]
  labs <- sort(unique(c(cor$label.i, cor$label.j)))
  mat <- matrix(0, nrow = length(labs), ncol = length(labs), dimnames = list(labs, labs))
  for (i in seq_len(nrow(cor))) {
    mat[cor$label.i[i], cor$label.j[i]] <- cor$corr[i]
  }
  diag(mat) <- 1
  mat[!is.finite(mat)] <- 0
  mat
}

ssmvln <- function(covar, hat, mc = 500L, seed = NULL) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package 'mvtnorm' is required for MVLN sampling.", call. = FALSE)
  }
  if (is.null(covar) || !NROW(covar)) {
    return(NULL)
  }
  if (is.null(hat) || !NROW(hat)) {
    return(NULL)
  }
  cor <- sscor(covar)
  if (is.null(cor) || !length(cor)) {
    return(NULL)
  }
  if (mc <= 1L) {
    return(cor)
  }

  names(hat) <- tolower(names(hat))
  hat <- hat[hat$label %in% rownames(cor), , drop = FALSE]
  if (!NROW(hat)) {
    return(NULL)
  }
  hat <- sshat(hat)
  cor <- cor[hat$label, hat$label, drop = FALSE]
  cvr <- cor2cov(cor, hat$std^2)
  cvr[!is.finite(cvr)] <- 0
  mvlmu <- log(hat$hat)
  names(mvlmu) <- hat$label

  if (!is.null(seed)) {
    set.seed(seed)
  }
  rtn <- mvtnorm::rmvnorm(mc, mean = mvlmu, sigma = cvr, method = "svd")
  rtn <- exp(rtn)
  colnames(rtn) <- hat$label

  dat <- data.frame(iter = seq_len(mc), rtn, check.names = FALSE)
  long <- reshape2::melt(dat, id.vars = "iter", variable.name = "label", value.name = "value")
  split_lab <- strsplit(as.character(long$label), "_", fixed = TRUE)
  long$metric <- vapply(split_lab, function(x) x[[1]], character(1))
  long$year <- vapply(split_lab, function(x) as.integer(x[[2]]), integer(1))
  long[, c("iter", "metric", "year", "label", "value")]
}

fmsy_from_rep <- function(rep) {
  dq <- rep[["derived_quants"]]
  if (is.null(dq) || !NROW(dq)) {
    return(NA_real_)
  }
  hit <- dq[["Label"]] %in% c("annF_MSY", "Fmsy", "F_MSY")
  if (!any(hit)) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(dq[hit, "Value"][[1]]))
}

mvln_long_to_kobe <- function(mvln_long, fmsy) {
  if (is.null(mvln_long) || !NROW(mvln_long) || !is.finite(fmsy) || fmsy <= 0) {
    return(NULL)
  }
  stock <- mvln_long[mvln_long$metric == "Bratio", , drop = FALSE]
  harvest <- mvln_long[mvln_long$metric == "F", , drop = FALSE]
  if (!NROW(stock) || !NROW(harvest)) {
    return(NULL)
  }
  names(stock)[names(stock) == "value"] <- "stock"
  names(harvest)[names(harvest) == "value"] <- "harvest"
  harvest$harvest <- harvest$harvest / fmsy
  merge(
    stock[, c("iter", "year", "stock"), drop = FALSE],
    harvest[, c("iter", "year", "harvest"), drop = FALSE],
    by = c("iter", "year"),
    all = FALSE
  )
}
