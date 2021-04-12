set.seed(12345)
n <- 10; p <- 2; d <- 5;
A <- matrix(rnorm(n*p), n, p); R1 <- matrix(rnorm(p * p), p, p) ; R2 <- matrix(rnorm(p * p), p, p)
X1 <- A %*% R1 %*% t(A) ; X2 <- A %*% R2 %*% t(A)
X  <- list(X1,X2)

result1 <- RESCAL_ALS(X, p, lam.A=0, lam.R=0, compute_fit=TRUE, conv=1e-16, maxIter=100)
Ahat <- result1$A ; R1hat <- result1$R[[1]]; R2hat <- result1$R[[2]];


test_that("RESCAL works", {
  expect_equal(norm(X1 - Ahat %*% R1hat %*%  t(Ahat)), 0, tolerance=1e-4)
})

test_that("RESCAL works", {
  expect_equal(norm(X2 - Ahat %*% R2hat %*%  t(Ahat)), 0, tolerance=1e-4)
})

result2 <- PRESCAL_ALS(X, p, P=diag(n), lam.A=0, lam.R=0, compute_fit=TRUE, conv=1e-19, maxIter=1000)
Ahat <- result2$A ; R1hat <- result2$R[[1]]; R2hat <- result2$R[[2]];


test_that("PRESCAL works", {
  expect_equal(norm(X1 - Ahat %*% R1hat %*%  t(Ahat)), 0, tolerance=1e-4)
})

test_that("PRESCAL works", {
  expect_equal(norm(X2 - Ahat %*% R2hat %*%  t(Ahat)), 0, tolerance=1e-4)
})
