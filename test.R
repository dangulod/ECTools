library(ECTools)
library(purrr)

x = read.csv("/Users/n87557/Documents/Metodologia de capital/skewt/vectortest.csv")

# fit skewt ---------------------------------------------------------------------

stx = x %>% map(fit_skewt)
stx = x %>% apply(2, fit_skewt)

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

mat = read_mat(file = "transition_matrix.xlsx")

tran(mat, horizon =  10)
abs_pd(mat, horizon =  10)
con_pd(mat, horizon =  10)
cum_pd(mat, horizon =  10)


# Anderson-Darling test ---------------------------------------------------

ad.test(x$serie1, "pst", dp = get_param(st))

# C++ classes R6 ----------------------------------------------------------------

pepe = new(Persona, 39, 0)
juan = new(Persona, 27, 0)

pepe$print()

pepe$ano(1)
juan$ano(1)

