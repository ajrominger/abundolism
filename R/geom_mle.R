#' @title Maximum likelihood estimation for geometric glm
#'
#' @description
#'     Find maximum likelihood parameters for glm assuming geometric
#'     distribution and logit link
#'
#' @param formula an object of class "formula" that describes the model
#' @param data an optional data.frame where the data live
#' @param ... arguments passed to `control` of `optim`
#'
#' @details
#'     Some methods that work for objects of class `glm` and `lm` will
#'     also work for the output of this function, specifically:
#'     `print`, `coef`, `logLik`, and `AIC`. A custom likelihood ratio test
#'     function is provided by `lrt_geom`
#'
#' @returns an object of class `glmGeom` inspired by, but not strictly
#'     inheriting from, `glm` and `lm`
#'
#' @export

geom_glm <- function(formula, data, ...) {
    if(missing(data)) {
        m <- model.frame(formula)
    } else {
        m <- model.frame(formula, data)
    }

    # model matrix and reponse
    X <- model.matrix(formula, m)
    y <- model.response(m)

    # use linear model with gaussian errors for starting
    # param values
    p0 <- lm(formula, data = data) |> coef()

    # mle
    fit <- optim(p0, fn = geo_ll, mx = X, my = y,
                 control = list(fnscale = -1, reltol = .Machine$double.eps^0.75,
                                ...))

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

    # subset of output from `glm` and `lm`
    res <- list(
        coefficients = beta_hat,
        # fitted.values = fitted_vals,
        # residuals = err,
        # linear.predictors = lin_pred,
        # rank = ncol(X),
        call = match.call(),
        x = X,
        y = y,
        deviance = -2 * fit$value,
        aic = 2 * (ncol(X) - fit$value),
        logLik = ll,
        terms = attr(m, "terms")
    )

    class(res) <- "glmGeom"

    return(res)



}


#' @export
logLik.glmGeom <- function(object) {
    object$logLik
}


#' @export
print.glmGeom <- function(object) {
    cat("Call:  ")
    print(object$call)

    cat("\nCoefficients:", "\n")
    print(object$coefficients)

    cat("\nLog likelihood:", "\n")
    print(logLik(object))

}


#' @title Likelihood ratio test for geometric glm
#'
#' @description
#'     Preform likelihood ratio test on two models fitted with `geom_glm`
#'     distribution and logit link
#'
#' @param mod0 reduced model
#' @param mod1 full model
#'
#'
#' @returns an ANOVA table like object inheriting from `anova` and `data.frame`.
#'     The LR test statistic can be extracted with `$Chisq[2]` and the P-value
#'     extracted with \code{$`Pr(>Chisq)`}.
#'
#' @export

lrt_geom <- function(mod0, mod1) {
    ll0 <- logLik(mod0)
    ll1 <- logLik(mod1)
    df0 <- attr(ll0, "df")
    df1 <- attr(ll1, "df")


    lr <- -2 * (ll0 - ll1)
    df <- df1 - df0

    res <- data.frame(df = c(df0, df1),
                      LogLik = c(ll0, ll1),
                      Df = c(NA, df1 - df0),
                      Chisq = c(NA, lr),
                      Pr = c(NA, pchisq(lr, df, lower.tail = FALSE)))
    names(res)[1] <- "#Df"
    names(res)[5] <- "Pr(>Chisq)"


    rownames(res) <- c("Model 0", "Model 1")

    attr(res, "heading") <- c("Likelihood ratio test\n",
                              sprintf("Model 0: %s\nModel 1: %s\n",
                                      as.character(mod0$call)[2],
                                      as.character(mod1$call)[2]))

    class(res) <- c("anova", "data.frame")

    res
}


#' @title Simple prediction for geometric glm's
#'
#' @description
#'     Predict response or link function
#'
#' @param mod fitted `geomGLM` object
#' @param newdata optional `data.frame` of new predictor variable values
#' @param type whether to predict response or link function values
#'
#' @returns column vector of predictions
#'
#' @export

pred_geom <- function(mod, newdata, type = c("response", "link")) {
    if(missing(newdata)) {
        mm <- model.matrix(mod)
    } else {
        tt <- delete.response(mod$terms)
        m <- model.frame(tt, newdata)
        mm <- model.matrix(tt, m)
    }

    lin <- mm %*% coef(mod)
    pr <- 1 / (1 + exp(-lin))

    # return either prob (=link) or response
    type <- match.arg(type, c("response", "link"))

    if(type == "response") {
        return(1 / (1 - pr))
    } else {
        return(pr)
    }
}


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
#
# coef(arth_mod)
# logLik(arth_mod)
