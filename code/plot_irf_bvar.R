#!/usr/bin/env Rscript
# code/plot_bvar_irf_from_draws.R
# Compute and plot orthogonalized BVAR IRFs from posterior draws saved in output/bvar_results.RData

# Minimal dependencies (base R)

datfile <- "output/bvar_results.RData"
if (!file.exists(datfile)) stop("File not found: ", datfile, "  -- run code/wdi_clean.R to produce it.")

load(datfile)   # expects B_draws (ns x K x M), Sigma_draws (ns x M x M), optionally p

if (!exists("B_draws") || !exists("Sigma_draws")) stop("B_draws or Sigma_draws not found in RData. Re-run the sampler.")
ns <- dim(B_draws)[1]
K  <- dim(B_draws)[2]
M  <- dim(B_draws)[3]

if (!exists("p")) {
  p_guess <- floor((K - 1) / M)
  p <- p_guess
  message("`p` not found in .RData; guessing p = ", p)
}

H <- 10  # horizons 0..H

build_A_list <- function(Bmat, M, p){
  A_list <- vector("list", p)
  for (lag in 1:p) {
    A <- matrix(0, M, M)
    for (r in 1:M) for (c in 1:M) {
      rowpos <- 1 + (lag - 1) * M + c
      A[r, c] <- Bmat[rowpos, r]
    }
    A_list[[lag]] <- A
  }
  A_list
}

compute_Theta <- function(A_list, H){
  p <- length(A_list)
  M <- nrow(A_list[[1]])
  Theta <- vector("list", H+1)
  Theta[[1]] <- diag(M) # Theta_0
  for (h in 1:H) {
    Th <- matrix(0, M, M)
    for (l in 1:p) if ((h - l) >= 0) Th <- Th + A_list[[l]] %*% Theta[[h - l + 1]]
    Theta[[h+1]] <- Th
  }
  Theta
}

IRF_draws <- array(NA_real_, dim = c(ns, H+1, M, M))
for (s in 1:ns){
  Bmat <- B_draws[s, , ]
  Sig  <- Sigma_draws[s, , ]
  if (any(is.na(Bmat)) || any(is.na(Sig))) next
  A_list <- build_A_list(Bmat, M, p)
  Theta  <- compute_Theta(A_list, H)
  L <- tryCatch(t(chol(Sig)), error = function(e) t(chol(Sig + diag(1e-8, nrow(Sig)))))
  for (h in 0:H) IRF_draws[s, h+1, , ] <- Theta[[h+1]] %*% L
}

IRF_med <- apply(IRF_draws, c(2,3,4), median, na.rm = TRUE)
IRF_lo  <- apply(IRF_draws, c(2,3,4), quantile, probs = 0.025, na.rm = TRUE)
IRF_hi  <- apply(IRF_draws, c(2,3,4), quantile, probs = 0.975, na.rm = TRUE)

if (exists("dep_names")) resp_names <- dep_names else if (exists("Y") && !is.null(colnames(Y))) resp_names <- colnames(Y) else resp_names <- paste0("y", 1:M)
shock_names <- resp_names

out_png <- "output/plot_irf_bvar.png"
# make a wider, page-width-friendly canvas (wide and not too tall)
png(out_png, width = 11, height = 7, units = "in", res = 300)

# try to create short labels similar to the previous plot: drop a leading 'g' (growth) if present
if (all(grepl('^g', resp_names))) {
  var_labels <- setNames(sub('^g', '', resp_names), resp_names)
} else {
  var_labels <- setNames(resp_names, resp_names)
}

hseq <- 0:H

# Use the same base-R multi-panel, margins and styling as the earlier IRF plot
par(
  mfrow = c(M, M),
  mar   = c(1.5, 3.2, 1.8, 0.5),   # increase left margin so left labels sit further from axis
  oma   = c(1.75, 4.5, 2, 2),           # more space bottom & left
  mgp   = c(1.4, 0.35, 0),         # axis label spacing
  xpd   = FALSE                    # clip everything inside each panel
)


for (i in 1:M) {
  for (j in 1:M) {
    med <- IRF_med[, i, j]
    lo  <- IRF_lo[, i, j]
    hi  <- IRF_hi[, i, j]

    x <- hseq
    ylim <- range(c(lo, hi, med), na.rm = TRUE)
    plot(x, med, type = "n", ylim = ylim,
         main = paste("coeff on", var_labels[[ shock_names[j] ]]), cex.main = 1.3, font.main = 2,
         xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")

    # light grid for readability
    grid_col <- "#f0f0f0"
    abline(v = x, col = grid_col, lty = 3)

    # shaded 95% credible band (grey) and solid blue median line
    band_col <- adjustcolor("grey60", alpha.f = 0.35)
    polygon(c(x, rev(x)), c(hi, rev(lo)), col = band_col, border = NA)
    lines(x, med, lwd = 2, col = "#1f77b4")

    # zero line
    abline(h = 0, col = "darkgrey", lty = 2)

    # axes: show x labels only on bottom row, y labels only on left column
  # show x-axis (horizon) only on the bottom row
  if (i == M) axis(1, at = x)
    if (j == 1) axis(2) else axis(2, labels = FALSE)

    box()
    # left-side row labels: only on the left-most column
    if (j == 1) {
      mtext(side = 2, text = paste("pred", var_labels[[ resp_names[i] ]]), font = 2, line = 2, cex = 0.9, las = 1)
    }
  }
}

dev.off()

cat("Wrote:", out_png, "\n")

png_irf_bvar <- image_read(out_png)
png_irf_bvar <- image_trim(png_irf_bvar)
image_write(png_irf_bvar, path = out_png, format = "png")


