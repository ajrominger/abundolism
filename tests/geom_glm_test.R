
x <- sample(1:100, size = 100, replace = TRUE)
fx <- 2 - 0.01 * (x - 50)^2
p <- 1 / (1 + exp(-fx))
plot(x, p)
y <- rgeom(length(x), 1 - p) + 1
dat <- data.frame(x, y)
plot(dat, log = "y")

boo <- geom_glm(y ~ x + I(x^2), data = dat)
coef(boo)
points(1:100, pred_geom(boo, data.frame(x = 1:100)), col = "blue")

plot(x, pred_geom(boo, type = "link"))
points(x, p, col = "red")

lm(y ~ x + I(x^2), data = dat)

lines(1:100, pred_geom(boo, data.frame(x = 1:100)))

plot(p, pred_geom(boo, data.frame(x = x)), log = "xy")


