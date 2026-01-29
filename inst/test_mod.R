# library(cmdstanr)


# simulate data
n <- 100
x <- runif(n, 0, 8)
b0 <- 2
b1 <- 2
b2 <- 0.5
b3 <- 80

a0 <- b0 - b2 * b1^2
a1 <- 2 * b1 * b2
a2 <- -b2

z <- a0 + a1 * x + a2 * x^2

plot(x, z)

p <- 1 / (1 + exp(-z))

plot(x, p)

m <- b3 * p

plot(x, m)

# m <- exp(b1 + b2 * x)
y <- rnbinom(n, mu = m, size = 1)

plot(x, y)

dat <- list(N = n,
            y = y,
            x = x)


mod <- cmdstan_model("inst/mod.stan")

fit <- mod$sample(
    data = dat,
    seed = 123,
    chains = 4,
    parallel_chains = 4,
    refresh = 500 # print update every 500 iters
)

fit$summary()

draws_df <- fit$draws(c("a_0", "a_1", "a_2"), format = "df")

hist(draws_df$a_3)
abline(v = a2)
