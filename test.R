library(ECTools)
library(purrr)
library(dplyr)
library(tibble)

x = read.csv("/Users/n87557/Documents/Metodologia de capital/skewt/vectortest.csv") %>%
  as_tibble()

# fit skewt ---------------------------------------------------------------------

stx = x %>%
  map(fit_skewt)

stx = x %>%
  apply(2, fit_skewt)

st = fit_skewt(x$serie1, n_days = 1)

c = stx %>% lapply(get_param)
stx %>% lapply(ks_test)

plot(x = st)
summary(st)
qst(0.9995, dp = get_param(st))
get_param(st)
ks_test(st)
qqPlot(st)

# fit saddle-point --------------------------------------------------------------

sp = fit_saddlepoint(x = x$serie1, N_SIMUL = 1e4, n_days = 1, nm = 4)
ks_test(sp, alternative = "two.sided", exact = NULL)
plot(sp)
qqPlot(sp)
sp

qqPlot(x$serie1, y = sp)

# transition matrix -------------------------------------------------------------

mat = read_mat(file = "transition_matrix.rda")

tran(mat, horizon =  10)
abs_pd(mat, horizon =  10)
con_pd(mat, horizon =  10)
cum_pd(mat, horizon =  10)


# Anderson-Darling test ---------------------------------------------------

ad.test(x$serie1, "pst", dp = get_param(st))


# cor_optim ---------------------------------------------------------------

library(readxl)
library(tibble)
library(magrittr)

# carga de matriz de correlaciones historicas --------------------------------

hist = read_excel("/Users/n87557/Documents/Metodologia de capital/Correlations/data/IRBB_3f.xlsx",
                  sheet = "matrix") %>%
  as.data.frame() %>%
  column_to_rownames("X__1") %>%
  apply(c(1,2), function(x) max(x, 0)) %>%
  as.matrix()

# carga de la matriz de correlaciones de los credit drivers ------------------

CD = read_excel("/Users/n87557/Documents/Metodologia de capital/Correlations/data/cor_FL_2018.xlsx") %>%
  as.data.frame() %>%
  column_to_rownames("X__1") %>%
  abs() %>%
  as.matrix()


# mapeo ----------------------------------------------------------------------

map = read_excel("/Users/n87557/Documents/Metodologia de capital/Correlations/data/IRBB_3f.xlsx",
                 sheet = "map")

library(ECTools)

x = cor_optim(map = map[,-3],hist = hist, CD = CD, lim = 1, maxiter = 1000, minFG = 0.4)

summary(x)

x = cor_optim(map = map,hist = hist, CD = CD, lim = 1, maxiter = 500, run = 10, suggestions = get_suggestions(x))

fitted_cor(x)

cbind(map[,1],
      as.data.frame(x))

R_2 = function(FG = FG, FL = FL, RU = map, hist = hist, CD = CD, col = col) {

  mat = mat(FG = FG, FL = FL, RU = map, col = col)
  mat1 = mat %*% CD %*% t(mat)
  diag(mat1) = 1

  return(sum((hist - mat) ^ 2))
}

xx = get_sensitivities(x)

R_2(FG = xx[,2], FL = as.matrix(xx[,-(1:2)]), RU = map, hist = hist, CD = CD, col = colnames(CD))


xx %>% write.csv2("/Users/n87557/Desktop/sensi.txt")


# ALM ---------------------------------------------------------------------

library(ECTools)
library(readxl)
library(tibble)

hist = read_excel("/Users/n87557/Desktop/cor_hist_ALM_max_1.xlsx",
           sheet = "matriz") %>%
  as.data.frame() %>%
  column_to_rownames("X__1") %>%
  apply(c(1,2), function(x) max(x, 0)) %>%
  as.matrix()

map = read_excel("/Users/n87557/Desktop/cor_hist_ALM_max_1.xlsx",
                 sheet = "mapeo")

CD = read_excel("/Users/n87557/Desktop/matriz_prueba18_02112017.xlsx") %>%
  as.data.frame() %>%
  column_to_rownames("X__1") %>%
  abs() %>%
  as.matrix()

x = cor_optim(map = map,
              hist = hist,
              CD = CD,
              lim = 0.6,
              maxiter = 5000,
              minFG = 0.27,
              run = 1000)

summary(x)
maxFL  = sqrt(lim - (minFG ^ 2))

c(rep(minFG, nrow(map)), rep(0, nrow(map) * (ncol(map) - 1)))
c(rep(1, nrow(map)), rep(maxFL, nrow(map) * (ncol(map) - 1)))

0.7331439 ^ 2 + minFG ^ 2

qqPlot(rst(100, omega = 1, alpha = 0.6, nu = 2),
       cex = 0.7,
       pch = 16,
       ylab = "Sample Quantiles",
       xlab = "Theorical Quantiles", cex.lab = 1, cex.axis = 0.1)
