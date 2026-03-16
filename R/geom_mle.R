#' @title Maximum likelihood estimation for geometric glm
#'
#' @description
#' Find maximum likelihood parameters for glm assuming geometric
#' distribution and logit link
#'
#' @param formula an object of class "formula" that describes the model
#' @param data an optional data.frame where the data live
#'
#' @returns an object of class `glmGeom` inheriting from `glm` and `lm`
#'
#' @export

geomglm <- function(formula, data) {
    if(missing(data)) {
        m <- model.frame(formula)
    } else {
        m <- model.frame(formula, data)
    }

    # model matrix and reponse
    X <- model.matrix(m)
    y <- model.response(m)

    # use linear model with gaussian errors for starting
    # param values
    p0 <- as.formula(m) |> lm() |> coef()

    # mle
    fit <- optim(p0, fn = geo_ll, mx = X, my = y, method = "BFGS",
                 hessian = TRUE,
                 control = list(fnscale = -1))

    if(fit$convergence != 0) {
        stop("bad convergence when maximizing log likelihood function")
    }

    # fitted values
    beta_hat <- fit$par
    lin_pred <- as.vector(X %*% beta_hat)
    p_hat <- 1 / (1 + exp(-lin_pred))
    fitted_vals <- as.vector((1 - p_hat) / p_hat)
    err <- y - fitted_vals
    ll <- structure(fit$value, df = ncol(X), class = "logLik")

    # output needed for `glm` and `lm`
    res <- list(
        coefficients = beta_hat,
        fitted.values = fitted_vals,
        residuals = err,
        linear.predictors = lin_pred,
        df.residual = nrow(X) - ncol(X),
        df.null = nrow(X) - 1,
        rank = ncol(X),
        call = match.call(),
        terms = attr(m, "terms"),
        model = m,
        x = X,
        y = y,
        qr = qr(X),
        family = list(family = "binomial", link = "logit", dispersion = 1),
        deviance = -2 * fit$value,
        null.deviance = -2 * null_geo_ll(y),
        converged = fit$convergence == 0,
        aic = 2 * (ncol(X) - fit$value),
        logLik = ll,
        R = solve(fit$hessian)
    )

    class(res) <- c("glmGeom", "glm", "lm")

    return(res)



}

boo <- geomglm(y ~ x + I(x^2))
boo0 <- geomglm(y ~ x)
anova(boo, boo0)

summary(boo)

foo0 <- glm(y ~ 1, family = "poisson")

anova(foo, foo0)

# non-exported function for calculating log likelihood
#' @param pars vector of parameter values
#' @param mx model matrix
#' @param my model response

geo_ll <- function(pars, mx, my) {
    # linear prediction
    fx <- mx %*% pars

    # logit-transformed
    p <- 1 / (1 + exp(-fx))

    # avoid really tiny p
    p <- pmax(pmin(p, 1 - 1e-10), 1e-10)

    # log probabilities for geo dist with support >= 1
    logd <- log(1 - p) + (my - 1) * log(p)

    return(sum(logd))
}


# non-exported function for calculating null log likelihood
#' @param y response data

null_geo_ll <- function(y) {
    p <- 1 - 1 / mean(y)

    sum(log(1 - p) + (y - 1) * log(p))
}


logLik.glmGeom <- function(object) {
    object$logLik
}

geometric <- function() {
    structure(
        list(
            family = "geometric",
            link = "logit",
            linkfun = qlogis,
            linkinv = plogis,
            variance = function(mu) mu + mu^2,
            dispersion = 1
        ),
        class = "family"
    )
}
