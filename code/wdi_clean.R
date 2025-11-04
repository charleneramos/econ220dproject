# =============================================================================#
#   Setup R for World Bank World Development Indicators (WDI)                  # 
#   https://databank.worldbank.org/source/world-development-indicators         #
# =============================================================================#

cat("\014")
rm(list = ls())
setwd('/Users/charleneramos/Documents/Git/econ220dproject')

# Install and load packages
install.packages("broom")
install.packages("BVAR")
install.packages("bvartools")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("magick")
install.packages("MCMCpack")
install.packages("purrr")
install.packages("quarto")
install.packages("rmarkdown")
install.packages("scales")
install.packages("strucchange")
install.packages("tidyr")
install.packages("tidyverse")
install.packages("tseries")
install.packages("WDI") # https://github.com/vincentarelbundock/WDI
install.packages("wbstats") # https://econandrew.github.io/wdi-api-R-gettingstarted/using-r-to-access-the-world-bank-api.html
install.packages("vars")

library(BVAR)
library(bvartools)
library(broom)
library(dplyr)
library(ggplot2)
library(magick)
library(MCMCpack)
library(purrr)
library(quarto)
library(rmarkdown)
library(strucchange)
library(scales)
library(tidyr)
library(tidyverse)
library(tseries)
library(WDI)
library(wbstats) 
library(vars)
library(MASS) # for mvrnorm

# Load functions

# Function to build V0 (K x K) as diagonal Minnesota-type matrix given lambda
build_V0 <- function(lambda_val) {
  Vdiag <- numeric(K)
  Vdiag[1] <- 1e8
  pos <- 2
  for (lag in 1:p) for (j in 1:M) {
    Vdiag[pos] <- (lambda_val^2) / (lag^2)
    pos <- pos + 1
  }
  diag(Vdiag)
}

# robust LaTeX escaping
escape_tex <- function(s) {
  if (is.null(s)) return(s)
  s <- as.character(s)
  s <- gsub('\\\\', '\\\\\\\\', s)   # backslash
  s <- gsub('\\$', '\\\\$', s)
  s <- gsub('%', '\\\\%', s)
  s <- gsub('&', '\\\\&', s)
  s <- gsub('#', '\\\\#', s)
  s <- gsub('_', '\\\\_', s)
  s <- gsub('\\{', '\\\\{', s)
  s <- gsub('\\}', '\\\\}', s)
  s <- gsub('\\^', '\\\\^', s)
  s <- gsub('~', '\\\\textasciitilde{}', s)
  s
}

glog <- function(x, lag_period = 1, scale = 100) {
  lx  <- log(x)
  llx <- dplyr::lag(lx, lag_period)
  scale * (lx - llx)
}

gham <- function(x, lag_period = 1, scale = 100) { # See Baumeister and Hamilton (2023) at pp. 5 and 39
  lag_x <- dplyr::lag(x, lag_period)
  num <- x - lag_x
  den <- 0.5 * (x + lag_x)
  out <- num / den
  # avoid Inf when den == 0
  out[den == 0] <- NA
  scale * out
}

# helper: map regressor names to table-1 style labels (const -> c, ggdp -> gdp, gfdi -> fdi, lags -> math subscripts)
make_reg_label <- function(s) {
  s <- as.character(s)
  if (is.na(s) || s == '') return('')
  s0 <- tolower(s)
  if (s0 %in% c('const','intercept')) return('$c$')
  m <- regexpr('^(.*)_lag(\\d+)$', s, perl = TRUE)
  if (m[1] != -1) {
    varname <- sub('^(.*)_lag(\\d+)$', '\\1', s, perl = TRUE)
    lagnum  <- sub('^(.*)_lag(\\d+)$', '\\2', s, perl = TRUE)
    varshort <- sub('^g', '', varname)
    return(paste0('$', varshort, '_{t-', lagnum, '}$'))
  }
  # no lag, remove leading g if present and show contemporaneous variable
  varshort <- sub('^g', '', s)
  return(paste0('$', varshort, '_t$'))
}

# =============================================================================#
#   II(a) Search for relevant macroeconomic variables                                # 
# =============================================================================#

# Get a list of countries and identify argument entry for Philippines
countries <- WDI_data$country # PH

# GDP (current US$)
gdp <- WDIsearch(string = 'gdp', field = 'name', cache = NULL)
ph_gdp <- WDI(country = 'PH', indicator = 'NY.GDP.MKTP.CD')
sum(!is.na(ph_gdp$NY.GDP.MKTP.CD)) # 65

# Unemployment, total (% of total labor force) (modeled ILO estimate)
unrate <- WDIsearch(string = 'unemployment', field = 'name', cache = NULL)
ph_unrate <- WDI(country = 'PH', indicator = 'SL.UEM.TOTL.ZS')
sum(!is.na(ph_unrate$SL.UEM.TOTL.ZS)) # 34

# Inflation, consumer prices (annual %)
cpi <- WDIsearch(string = 'inflation', field = 'name', cache = NULL)
ph_cpi <- WDI(country = 'PH', indicator = 'FP.CPI.TOTL.ZG')
sum(!is.na(ph_cpi$FP.CPI.TOTL.ZG)) # 14

# Population, total
pop <- WDIsearch(string = 'population', field = 'name', cache = NULL)
ph_pop <- WDI(country = 'PH', indicator = 'SP.POP.TOTL')
sum(!is.na(ph_pop$SP.POP.TOTL)) # 65

# Labor force, total
lf <- WDIsearch(string = 'labor force', field = 'name', cache = NULL)
ph_lf <- WDI(country = 'PH', indicator = 'SL.TLF.TOTL.IN')
sum(!is.na(ph_lf$SL.TLF.TOTL.IN)) # 35

# Foreign direct investment, net inflows (BoP, current US$)
fdi <- WDIsearch(string = 'foreign direct investment', field = 'name', cache = NULL)
ph_fdi <- WDI(country = 'PH', indicator = 'BX.KLT.DINV.CD.WD')
sum(!is.na(ph_fdi$BX.KLT.DINV.CD.WD)) # 55 

# Net official development assistance received (current US$)
oda <- WDIsearch(string = 'Net official development assistance received', field = 'name', cache = NULL)
ph_oda <- WDI(country = 'PH', indicator = 'DT.ODA.ODAT.CD')
sum(!is.na(ph_oda$DT.ODA.ODAT.CD)) # 64 

# Personal remittances, received (current US$)
remit <- WDIsearch(string = 'personal remittances', field = 'name', cache = NULL)
ph_remit <- WDI(country = 'PH', indicator = 'BX.TRF.PWKR.CD.DT')
sum(!is.na(ph_remit$BX.TRF.PWKR.CD.DT)) # 48

rm(list = c('cpi', 'fdi', 'gdp', 'lf', 'oda', 'pop', 'remit', 'unrate'))
rm(list = c('ph_unrate','ph_cpi','ph_lf'))

# =============================================================================#
#   II(b) Create merged dataset                                                      # 
# =============================================================================#

# Remove NA values
ph_fdi <- na.omit(ph_fdi)
ph_gdp <- na.omit(ph_gdp)
ph_oda <- na.omit(ph_oda)
ph_pop <- na.omit(ph_pop)
ph_remit <- na.omit(ph_remit)

# Reorder by ascending year
ph_fdi <- ph_fdi[order(ph_fdi$year), ]
ph_gdp <- ph_gdp[order(ph_gdp$year), ]
ph_oda <- ph_oda[order(ph_oda$year), ]
ph_pop <- ph_pop[order(ph_pop$year), ]
ph_remit <- ph_remit[order(ph_remit$year), ]

# Drop irrelevant columns
ph_fdi <- ph_fdi[, !(names(ph_fdi) %in% c('country','iso3c'))]
ph_gdp <- ph_gdp[, !(names(ph_gdp) %in% c('country','iso3c'))]
ph_oda <- ph_oda[, !(names(ph_oda) %in% c('country','iso3c'))]
ph_pop <- ph_pop[, !(names(ph_pop) %in% c('country','iso3c'))]
ph_remit <- ph_remit[, !(names(ph_remit) %in% c('country','iso3c'))]

# Rename columns
names(ph_fdi) <- c('country','year','fdi')
names(ph_gdp) <- c('country','year','gdp')
names(ph_oda) <- c('country','year','oda')
names(ph_pop) <- c('country','year','pop')
names(ph_remit) <- c('country','year','remit')

# inner-join all series on country + year and keep desired columns
df <- Reduce(function(x, y) merge(x, y, by = c('country', 'year')), 
             list(ph_fdi, ph_gdp, ph_oda, ph_pop, ph_remit))

df <- df[, c('year','pop','gdp','fdi','oda','remit')]

# long format with tidyverse
df_long <- df %>%
    pivot_longer(cols = c(pop, gdp, fdi, oda, remit),
                             names_to = "variable",
                             values_to = "value")

rm('countries','ph_fdi','ph_gdp','ph_oda','ph_pop','ph_remit')

# Save datasets
write.csv(df, file = 'data/wdi_philippines.csv', row.names = FALSE)
write.csv(df_long, file = 'data/wdi_philippines_long.csv', row.names = FALSE)
 
df <- read.csv('data/wdi_philippines.csv', stringsAsFactors = FALSE) 
df_long <- read.csv('data/wdi_philippines_long.csv', stringsAsFactors = FALSE) 

# =============================================================================#
#   II(c) Plot assembled data                                                        # 
# =============================================================================#

# Principal component analysis
pca <- prcomp(df[, c('pop','gdp','fdi','oda','remit')],
                 center = TRUE, scale. = TRUE)
summary(pca)

# Population, total
plot_pop <- df_long[df_long$variable == 'pop', ]
plot_pop$value_mil <- plot_pop$value / 1e6 
plot_pop$value_bil <- plot_pop$value / 1e9 

png('output/plot_pop.png', width = 800, height = 600, res = 120)

plot(plot_pop$year, plot_pop$value_mil, type = 'l',
     col = 'black', lwd = 2, bty = 'l',
     xlab = 'Year', ylab = 'Population Levels (Millions)',
     col.axis = 'black', col.lab = 'black',
     col.main = 'black', fg = 'black',
     las = 1)

box(lwd = 1.2, col = "black")

dev.off()

png_pop <- image_read('output/plot_pop.png')
png_pop <- image_trim(png_pop)
image_write(png_pop, path = 'output/plot_pop.png')

# GDP (current US$)
plot_gdp <- df_long[df_long$variable == 'gdp', ]
plot_gdp$value_mil <- plot_gdp$value / 1e6 
plot_gdp$value_bil <- plot_gdp$value / 1e9 

png('output/plot_gdp.png', width = 800, height = 600, res = 120)

plot(plot_gdp$year, plot_gdp$value_bil, type = 'l',
     col = 'black', lwd = 2, bty = 'l',
     xlab = 'Year', ylab = 'Current US$ (Billions)',
     col.axis = 'black', col.lab = 'black',
     col.main = 'black', fg = 'black',
     las = 1)

box(lwd = 1.2, col = "black")

dev.off()

png_gdp <- image_read('output/plot_gdp.png')
png_gdp <- image_trim(png_gdp)
image_write(png_gdp, path = 'output/plot_gdp.png')

# Foreign direct investment, net inflows (BoP, current US$)
plot_fdi <- df_long[df_long$variable == 'fdi', ]
plot_fdi$value_mil <- plot_fdi$value / 1e6 
plot_fdi$value_bil <- plot_fdi$value / 1e9 

png('output/plot_fdi.png', width = 800, height = 600, res = 120)

plot(plot_fdi$year, plot_fdi$value_bil, type = 'l',
     col = 'black', lwd = 2, bty = 'l',
     xlab = 'Year', ylab = 'Current US$ (Billions)',
     col.axis = 'black', col.lab = 'black',
     col.main = 'black', fg = 'black',
     las = 1)

box(lwd = 1.2, col = "black")

dev.off()

png_fdi <- image_read('output/plot_fdi.png')
png_fdi <- image_trim(png_fdi)
image_write(png_fdi, path = 'output/plot_fdi.png')

# =============================================================================#
#   III VAR and IRF analysis                                                   # 
# =============================================================================#

# Transform data for VAR
dfg <- df %>%
  mutate(
    gpop   = glog(pop),
    ggdp   = glog(gdp),
    gfdi = gham(fdi),
    goda = gham(oda),
    gremit = glog(remit)
  ) %>% 
  dplyr::select(year, ggdp, gpop, gfdi, goda, gremit) %>% 
  na.omit()

X <- dfg %>% dplyr::select(ggdp, gfdi, goda, gremit) %>% as.matrix()

# Lag selection
lag_select <- VARselect(X, lag.max = 8, type = "const")
print(lag_select)
p <- as.integer(lag_select$selection["SC(n)"])

# Estimate VAR and IRF
var <- VAR(X, p = p, type = "const")
irf95 <- irf(var,  n.ahead = 10, boot = TRUE , ci = 0.95, ortho = FALSE)

# Format VAR table output
var_tidy <- broom::tidy(var)

var_tidy <- var_tidy %>%
  mutate(
    star = case_when(
      p.value < 0.01 ~ "\\textsuperscript{***}",
      p.value < 0.05 ~ "\\textsuperscript{**}",
      p.value < 0.10 ~ "\\textsuperscript{*}",
      TRUE ~ ""
    ),
    est_star = paste0(round(estimate, 4), star)
  )

var_table <- var_tidy %>%
  dplyr::select(group, term, est_star) %>%
  pivot_wider(names_from = group, values_from = est_star)

# Format IRFs table output
response <- "ggdp"
impulse <- c("gfdi","goda","gremit")

irf_df <- purrr::map_dfr(impulse, function(imp) {
  tibble(
    horizon = seq_len(nrow(irf95$irf[[response]])),
    impulse = imp,
    irf     = irf95$irf[[response]][, imp],
    lower   = irf95$Lower[[response]][, imp],
    upper   = irf95$Upper[[response]][, imp]
  )
})

irf_table <- irf_df %>%
  mutate(
    # Determine which significance level each CI excludes zero at
    star = case_when(
      (lower * upper > 0) & (abs(irf) > 0.01) ~ "$\\textsuperscript{***}$",  # strongest
      (lower * upper > 0) & (abs(irf) > 0.005) ~ "$\\textsuperscript{**}$",
      (lower * upper > 0) ~ "$\\textsuperscript{*}$",
      TRUE ~ ""
    ),
    cell = paste0(sprintf("%.4f", irf), star)
  ) %>%
  dplyr::select(horizon, impulse, cell) %>%
  pivot_wider(names_from = impulse, values_from = cell)

# Plot IRFs for GDP response
vars <- names(irf95$irf) # Always take names from the IRF object to avoid mismatches
# Rename them for labeling only
var_labels <- c(
  ggdp   = "gdp",
  gfdi   = "fdi",
  goda   = "oda",
  gremit = "remit"
)

k <- length(vars)
h <- nrow(irf95$irf[[1]]) # Use the actual number of horizons in the matrices

png("output/plot_irf.png", width = 11, height = 9, units = "in", res = 300)

par(
  mfrow = c(k, k),
  mar   = c(1.5, 3.2, 1.8, 0.5),   # increase left margin so left labels sit further from axis
  oma   = c(8, 9, 4, 2),           # more space bottom & left
  mgp   = c(1.4, 0.35, 0),         # axis label spacing
  xpd   = FALSE                    # clip everything inside each panel
)

for (i in vars) {
  for (j in vars) {
    mat <- irf95$irf[[i]]
    low <- irf95$Lower[[i]]
    up  <- irf95$Upper[[i]]
    
    if (is.null(dim(mat))) {
      y <- mat; lower <- low; upper <- up
    } else {
      y <- mat[, j]; lower <- low[, j]; upper <- up[, j]
    }
    
    x <- seq_along(y)
    ylim <- range(y, lower, upper, na.rm = TRUE)
    plot(x, y, type = "n", ylim = ylim,
         main = paste("coeff on", var_labels[j]), cex.main = 1.2, font.main = 2,
         xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")

    # subtle vertical grid lines
    grid_col <- "#f0f0f0"
    abline(v = x, col = grid_col, lty = 3)

    # shaded 95% CI band and blue median line
    band_col <- grDevices::adjustcolor("grey60", alpha.f = 0.35)
    polygon(c(x, rev(x)), c(upper, rev(lower)), col = band_col, border = NA)
    lines(x, y, lwd = 2, col = "#1f77b4")

    # horizontal zero line
    abline(h = 0, col = "darkgrey", lty = 2)

    # axes: x only on bottom row, y only on left column
  # show x-axis (horizon) ticks and labels only on the bottom row
  if (i == vars[length(vars)]) axis(1, at = x)
    if (j == vars[1]) axis(2) else axis(2, labels = FALSE)

    box()
    if (j == vars[1])
      mtext(side = 2, text = paste("pred", var_labels[i]), font = 2, line = 2, cex = 0.9, las = 1)
  }
}

dev.off()

png_irf <- image_read("output/plot_irf.png")
image_write(image_trim(png_irf), path = "output/plot_irf.png", format = "png")

# =============================================================================#
#   IV Bayesian analysis                                                       # 
# =============================================================================#

if (!exists("p")) p <- 1
Tfull <- nrow(X)
rows <- (p + 1):Tfull
Y <- X[rows, , drop = FALSE]
Z <- matrix(1, nrow = length(rows), ncol = 1)
for (lag in 1:p) Z <- cbind(Z, X[rows - lag, , drop = FALSE])

Tn <- nrow(Z); K <- ncol(Z); M <- ncol(Y)

# OLS quantities
XtX <- t(Z) %*% Z
XtX_inv <- solve(XtX)
B_ols <- XtX_inv %*% t(Z) %*% Y
E_ols <- Y - Z %*% B_ols
S <- t(E_ols) %*% E_ols    # sum of squared residuals matrix (MxM)

# Posterior with proper NIW prior on Sigma and Normal prior on B with Minnesota-style V0(lambda).
# We'll run Gibbs for B | Sigma and Sigma | B, and Metropolis-Hastings for lambda (overall tightness).

# Prior hyperparameters
lambda <- 0.2
lambda_prior_shape <- 2
lambda_prior_rate  <- 2

# NIW prior for Sigma: Sigma ~ Inv-Wishart(nu0, S0)
nu0 <- M + 2
S0 <- diag(rep(0.1, M))

# Prior mean for B (K x M)
B0 <- matrix(0, nrow = K, ncol = M)

set.seed(123)
niter <- 3000; burn <- 1000; thin <- 1
nsave <- floor((niter - burn) / thin)

# storage
B_draws <- array(NA, dim = c(nsave, K, M))
Sigma_draws <- array(NA, dim = c(nsave, M, M))
lambda_draws <- numeric(nsave)

# initialize
V0 <- build_V0(lambda)
V0_inv <- solve(V0)
Sigma_cur <- MCMCpack::riwish(nu0 + Tn, S0 + S) # initialize from posterior-ish draw
B_cur <- B_ols

save_idx <- 0
for (it in 1:niter) {
  # --- Gibbs step: sample B | Sigma_cur, lambda ---
  Vn_inv <- V0_inv + XtX
  Vn <- solve(Vn_inv)
  Bn <- Vn %*% (V0_inv %*% B0 + t(Z) %*% Y)
  # draw matrix-normal: B = Bn + A * Zm * t(C)
  A <- tryCatch(chol(Vn), error = function(e) { chol(Vn + diag(1e-8, nrow(Vn))) })
  C <- tryCatch(chol(Sigma_cur), error = function(e) { chol(Sigma_cur + diag(1e-8, nrow(Sigma_cur))) })
  Zm <- matrix(rnorm(K * M), nrow = K, ncol = M)
  B_cur <- Bn + A %*% Zm %*% t(C)

  # --- Gibbs step: sample Sigma | B_cur, Y ---
  Ecur <- Y - Z %*% B_cur
  S_post <- S0 + t(Ecur) %*% Ecur
  nu_post <- nu0 + Tn
  Sigma_cur <- MCMCpack::riwish(nu_post, S_post)

  # --- MH step: update lambda conditional on B_cur and Sigma_cur ---
  lam_prop_log <- log(lambda) + rnorm(1, 0, 0.05)
  lam_prop <- exp(lam_prop_log)
  if (lam_prop > 1e-6 && lam_prop < 5) {
    V0_prop <- build_V0(lam_prop)
    # compute log p(B_cur | Sigma_cur, lambda)
    # logdet term: det(Sigma_cur  V0_prop) = det(Sigma_cur)^K * det(V0_prop)^M
    logdet_prop <- K * as.numeric(determinant(Sigma_cur, logarithm = TRUE)$modulus) + M * as.numeric(determinant(V0_prop, logarithm = TRUE)$modulus)
    logdet_cur  <- K * as.numeric(determinant(Sigma_cur, logarithm = TRUE)$modulus) + M * as.numeric(determinant(V0, logarithm = TRUE)$modulus)
    # quadratic form: tr(Sigma^{-1} t(Bdiff) V0^{-1} Bdiff)
    Bdiff <- B_cur - B0
    V0_prop_inv <- solve(V0_prop)
    V0_inv <- V0_inv # from current
    quad_prop_mat <- t(Bdiff) %*% V0_prop_inv %*% Bdiff
    quad_cur_mat  <- t(Bdiff) %*% V0_inv %*% Bdiff
    Sigma_inv <- solve(Sigma_cur)
    quad_prop <- sum(Sigma_inv * quad_prop_mat)
    quad_cur  <- sum(Sigma_inv * quad_cur_mat)

    loglik_prop <- -0.5 * (logdet_prop + quad_prop)
    loglik_cur  <- -0.5 * (logdet_cur  + quad_cur)

    logprior_prop <- dgamma(lam_prop, shape = lambda_prior_shape, rate = lambda_prior_rate, log = TRUE)
    logprior_cur  <- dgamma(lambda,   shape = lambda_prior_shape, rate = lambda_prior_rate, log = TRUE)

    log_accept <- (loglik_prop + logprior_prop) - (loglik_cur + logprior_cur)
    if (log(runif(1)) < log_accept) {
      lambda <- lam_prop
      V0 <- V0_prop
      V0_inv <- solve(V0)
    }
  }

  # save draws
  if (it > burn && ((it - burn) %% thin == 0)) {
    save_idx <- save_idx + 1
    B_draws[save_idx, , ] <- B_cur
    Sigma_draws[save_idx, , ] <- Sigma_cur
    lambda_draws[save_idx] <- lambda
  }

  if (it %% 500 == 0) message('NIW+MH iter: ', it, ' lambda=', round(lambda,3))
}

# posterior summaries
B_post_mean <- apply(B_draws, c(2,3), mean)
Sigma_post_mean <- apply(Sigma_draws, c(2,3), mean)
lambda_post_mean <- mean(lambda_draws)

save(B_draws, Sigma_draws, lambda_draws, B_post_mean, Sigma_post_mean, lambda_post_mean, file = 'output/bvar_results.RData')

# -------------------------------
# Bayesian VAR Estimates Summary
# -------------------------------
dir.create('output', showWarnings = FALSE)

# Names for regressors in Z: intercept + lags
colnames_Z <- c('const')
vars_Z <- colnames(X)
if (is.null(vars_Z)) vars_Z <- paste0('var', seq_len(ncol(X)))
for (lag in 1:p) for (v in vars_Z) colnames_Z <- c(colnames_Z, paste0(v, '_lag', lag))

# Dependent variable names (equations)
dep_names <- c('gdp', 'fdi', 'oda', 'remit')
if (is.null(dep_names)) dep_names <- paste0('eq', seq_len(M))

res_rows <- list()
idx <- 1
for (m in 1:M) {
  for (k_i in 1:K) {
    draws <- B_draws[, k_i, m]
    mean_v <- mean(draws)
    median_v <- median(draws)
    sd_v <- sd(draws)
    q025 <- quantile(draws, 0.025)
    q975 <- quantile(draws, 0.975)
    pr_pos <- mean(draws > 0)
    t_like <- ifelse(sd_v > 0, mean_v / sd_v, NA)

    res_rows[[idx]] <- data.frame(
      equation = dep_names[m],
      regressor = colnames_Z[k_i],
      mean = mean_v,
      median = median_v,
      sd = sd_v,
      q025 = as.numeric(q025),
      q975 = as.numeric(q975),
      Pr_pos = pr_pos,
      t_like = t_like,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}

bvar_sum <- do.call(rbind, res_rows)
write.csv(bvar_sum, file = 'output/bvar_sum.csv', row.names = FALSE)

## Write a plain LaTeX fragment (Rmd child will include raw LaTeX)
tex_lines <- c()
# begin condensed table fragment: smaller font, reduced spacing, and auto-scale to text width
# Use [H] to force placement at the include point (requires placeins or float package/preamble)
tex_lines <- c(tex_lines, '\\begin{table}[H]')
tex_lines <- c(tex_lines, '\\centering')
tex_lines <- c(tex_lines, '\\caption{Bayesian VAR posterior estimates}')
tex_lines <- c(tex_lines, '\\label{tab:bvar_coef_post}')
tex_lines <- c(tex_lines, '{\\setlength{\\tabcolsep}{5pt}\\renewcommand{\\arraystretch}{0.95}\\small\\begin{adjustbox}{max width=\\textwidth}')
 # columns: Equation | Regressor | Mean | Median | SD | 95% CI | Pr(>0)
 # Use llrrrcr so textual columns (first two) are left aligned and numeric columns right aligned,
 # but center the header text using \multicolumn{1}{c}{...}
tex_lines <- c(tex_lines, '\\begin{tabular}{llrrrcr}')
tex_lines <- c(tex_lines, '\\hline')
tex_lines <- c(tex_lines, ' & Regressor & Mean & Median & SD & 95\\% CI & Pr($>$0) \\\\')
tex_lines <- c(tex_lines, '\\hline')

for (m in unique(bvar_sum$equation)) {
  # grouped by equation, but do not print an explicit equation header row
  subset_m <- bvar_sum[bvar_sum$equation == m, ]
  for (i in seq_len(nrow(subset_m))) {
    r <- subset_m[i, ]
    # show equation name only on the first row for this equation, else blank
    eqn  <- if (i == 1) escape_tex(as.character(r$equation)) else ''
    # regressor label matching Table 1 style (math mode)
    reg_label <- make_reg_label(as.character(r$regressor))
  tex_lines <- c(tex_lines, sprintf('%s & %s & %.4f & %.4f & %.4f & (%.4f, %.4f) & %.3f \\\\',
                    eqn, reg_label, r$mean, r$median, r$sd, r$q025, r$q975, r$Pr_pos))
  }
  tex_lines <- c(tex_lines, '\\hline')
}

tex_lines <- c(tex_lines, '\\end{tabular}')
tex_lines <- c(tex_lines, '\\end{adjustbox}}')
tex_lines <- c(tex_lines, '\\end{table}')

writeLines(tex_lines, con = 'output/bvar_table.Rmd')

# -------------------------------
# Stationarity and Structural Breaks
# -------------------------------

vars_to_test <- colnames(X)
if (is.null(vars_to_test)) vars_to_test <- paste0('var', seq_len(ncol(X)))

test_sum <- list()
for (v in vars_to_test) {
  series <- dfg[[v]]
  # drop NA
  ok <- !is.na(series)
  series <- series[ok]
  years <- dfg$year[ok]

  # Augmented Dickey-Fuller (null = unit root)
  adf_out <- tryCatch(tseries::adf.test(series, alternative = "stationary"), error = function(e) NULL)
  adf_stat <- if (!is.null(adf_out)) adf_out$statistic else NA
  adf_pval <- if (!is.null(adf_out)) adf_out$p.value else NA

  # KPSS (null = stationarity)
  kpss_out <- tryCatch(tseries::kpss.test(series, null = "Level"), error = function(e) NULL)
  kpss_stat <- if (!is.null(kpss_out)) kpss_out$statistic else NA
  kpss_pval <- if (!is.null(kpss_out)) kpss_out$p.value else NA

  # Structural breaks (Bai-Perron style) for mean shifts using strucchange::breakpoints
  bp <- tryCatch(strucchange::breakpoints(series ~ 1), error = function(e) NULL)
  bps <- if (!is.null(bp) && !is.null(bp$breakpoints)) bp$breakpoints else integer(0)
  # remove NA entries (sometimes breakpoints contains NA) and keep only valid indices
  bps <- bps[!is.na(bps)]
  bp_years <- if (length(bps) > 0) years[bps] else integer(0)

  # save small plot of breaks
  pngfile <- file.path('output', sprintf('break_%s.png', v))
  png(pngfile, width = 900, height = 600, res = 150)
  plot(years, series, type = 'l', main = paste('Series:', v), xlab = 'Year', ylab = v)
  if (length(bp_years) > 0) abline(v = bp_years, col = 'red', lwd = 2)
  dev.off()

  test_sum[[v]] <- data.frame(
    variable = v,
    adf_stat = as.numeric(adf_stat),
    adf_pval = as.numeric(adf_pval),
    kpss_stat = as.numeric(kpss_stat),
    kpss_pval = as.numeric(kpss_pval),
    n_breaks = ifelse(length(bp_years) > 0, length(bp_years), 0),
    break_years = ifelse(length(bp_years) > 0, paste(bp_years, collapse = '; '), NA),
    stringsAsFactors = FALSE
  )
}

test_sum <- do.call(rbind, test_sum)
write.csv(test_sum, file = 'output/test_sum.csv', row.names = FALSE)

# -------------------------------
# Stationarity & Breaks LaTeX table (Rmd fragment)
# -------------------------------
tex_lines_ts <- c()
tex_lines_ts <- c(tex_lines_ts, '\\begin{table}[H]')
tex_lines_ts <- c(tex_lines_ts, '\\centering')
tex_lines_ts <- c(tex_lines_ts, '\\caption{Testing of stationarity and structural breaks}')
tex_lines_ts <- c(tex_lines_ts, '\\label{tab:stationarity_breaks}')
tex_lines_ts <- c(tex_lines_ts, '{\\setlength{\\tabcolsep}{5pt}\\renewcommand{\\arraystretch}{0.95}\\small\\begin{adjustbox}{max width=\\textwidth}')
tex_lines_ts <- c(tex_lines_ts, '\\begin{tabular}{lrrrrcr}')
tex_lines_ts <- c(tex_lines_ts, '\\hline')
tex_lines_ts <- c(tex_lines_ts, 'Variable & ADF stat & ADF p-val & KPSS stat & KPSS p-val & \# breaks & Break years \\\\')
tex_lines_ts <- c(tex_lines_ts, '\\hline')

for (i in seq_len(nrow(test_sum))) {
  r <- test_sum[i, ]
  var_label <- escape_tex(as.character(r$variable))
  adf_stat <- ifelse(is.na(r$adf_stat), 'NA', sprintf('%.4f', as.numeric(r$adf_stat)))
  adf_pval <- ifelse(is.na(r$adf_pval), 'NA', sprintf('%.3f', as.numeric(r$adf_pval)))
  kpss_stat <- ifelse(is.na(r$kpss_stat), 'NA', sprintf('%.4f', as.numeric(r$kpss_stat)))
  kpss_pval <- ifelse(is.na(r$kpss_pval), 'NA', sprintf('%.3f', as.numeric(r$kpss_pval)))
  n_breaks <- ifelse(is.na(r$n_breaks), '0', as.character(r$n_breaks))
  breaks_raw <- as.character(r$break_years)
  breaks_fmt <- ifelse(is.na(breaks_raw) || breaks_raw == 'NA', 'NA', escape_tex(breaks_raw))

  tex_lines_ts <- c(tex_lines_ts, sprintf('%s & %s & %s & %s & %s & %s & %s \\\\',
                                           var_label, adf_stat, adf_pval, kpss_stat, kpss_pval, n_breaks, breaks_fmt))
}

tex_lines_ts <- c(tex_lines_ts, '\\hline')
tex_lines_ts <- c(tex_lines_ts, '\\end{tabular}')
tex_lines_ts <- c(tex_lines_ts, '\\end{adjustbox}}')
tex_lines_ts <- c(tex_lines_ts, '\\end{table}')

writeLines(tex_lines_ts, con = 'output/test_sum_table.Rmd')



