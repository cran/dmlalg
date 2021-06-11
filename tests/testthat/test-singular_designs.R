n <- 100
K <- 2
cond_method <- rep("forest", 3)
S <- 10
gamma <- exp(seq(-4, 10, length.out = 20))


####################
# generate data
####################

discrete_1 <- function(n, beta0 = 1) {
  W <- as.numeric(rnorm(n, 0, 1) >= 0)
  A <- cbind(W + 0.0005 * (rnorm(n, 0, 1) >= 1.5),
             W + 0.05 * rnorm(n, 0, 1),
             W + rnorm(n, 0, 1))
  H <- 1 * W + 0.25 * rnorm(n, 0, 1)
  X <- cbind(A[, 1] + W, W + 0.05 * rnorm(n, 0, 1), W + rnorm(n, 0, 1))
  Y <- X %*% rep(beta0, 3) + H - W + 0.25 * rnorm(n, 0, 1)
  return(list(aa = as.matrix(A), xx = as.matrix(X),
              ww = data.frame(w = W), yy = as.matrix(Y)))
}

discrete_2 <- function(n, beta0 = 1) {
  W <- cbind(as.numeric(rnorm(n) <= 0.2533471),
             as.numeric(rnorm(n) <= -0.2533471))
  A <- 1 * rnorm(n, 0, 1)
  H <- 1 * W[, 1] + 0.25 * rnorm(n, 0, 1)
  X <- W
  Y <- X %*% rep(beta0, ncol(X)) + H - tanh(W[, 2]) + 0.25 * rnorm(n, 0, 1)
  return(list(aa = as.matrix(A), xx = as.matrix(X),
              ww = data.frame(w = W), yy = as.matrix(Y)))
}

data_3 <- function(n, beta0 = 1) {
  W <- rnorm(n, 0, 1)
  A <- cbind(W, rnorm(n, 0, 1))
  index <- 50
  A[1:index, 1] <- rnorm(index, 0, 1)
  H <- rnorm(n, 0, 1)
  X <- W
  Y <- beta0 * X + H
  return(list(aa = as.matrix(A), xx = as.matrix(X),
              ww = data.frame(w = W), yy = as.matrix(Y)))
}

data_4 <- function(n, beta0 = 1) {
  W <- rnorm(n, 0, 1)
  A <- W
  index <- 1
  A[index] <- A[index] + 0.000001 * rnorm(length(index), 0, 1)
  H <- rnorm(n, 0, 1)
  X <- W
  Y <- beta0 * X + H
  return(list(aa = as.matrix(A), xx = as.matrix(X),
              ww = data.frame(w = W), yy = as.matrix(Y)))
}


####################
# perform tests
####################

test_that("error and warnings work with (partly) singular design", {
  set.seed(5)
  data <- discrete_1(n)
  res <- tryCatch_W_E(regsdml(a = data$aa, w = data$ww, x = data$xx, y = data$yy,
                              gamma = gamma, S = S,
                              DML = "DML1",
                              cond_method = cond_method), "error-found")
  expect_null(res$error)

  set.seed(5)
  data <- discrete_2(n)
  res <- tryCatch_W_E(regsdml(a = data$aa, w = data$ww, x = data$xx, y = data$yy,
                              gamma = gamma, S = S,
                              DML = "DML1",
                              cond_method = cond_method), "error-found")
  expect_equal(res$value, "error-found")

  set.seed(50)
  data <- data_3(n)
  res <- tryCatch_W_E(regsdml(a = data$aa, w = data$ww, x = data$xx, y = data$yy,
                              gamma = gamma, S = S,
                              DML = "DML1",
                              cond_method = rep("ols", 3)), "error-found")
  expect_equal(res$warning,
               "\nWarning messages:\nEssentially perfect fit: do S more repetitions.")
  expect_null(res$error)

  set.seed(5)
  data <- data_4(n)
  res <- tryCatch_W_E(regsdml(a = data$aa, w = data$ww, x = data$xx, y = data$yy,
                              gamma = gamma, S = S,
                              DML = "DML2",
                              cond_method = rep("ols", 3)), "error-occurred")
  expect_equal(res$warning, "\nWarning messages:\nEssentially perfect fit: DML summary may be unreliable.")
})
