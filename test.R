library(fda)
library(Rcpp)

x <- seq(0, 1, length.out = 20)
knots <- seq(0, 1, length.out = 400)
sourceCpp("main.cpp")

start_time <- proc.time()
B <- calcB(x, knots, 8);
end_time <- proc.time()
end_time - start_time

start_time2 <- proc.time()
basis <- create.bspline.basis(breaks = knots, norder = 8)
B1 <- eval.basis(x, basis)
end_time2 <- proc.time()
end_time2 - start_time2

start_time3 <- proc.time()
Omega_B <- calc_Omega(knots, 4);
end_time3 <- proc.time()
end_time3 - start_time3

start_time4 <- proc.time()
basis <- create.bspline.basis(breaks = knots, norder = 4)
Omega_B1 <- fda::bsplinepen(basis)
end_time4 <- proc.time()
end_time4 - start_time4

