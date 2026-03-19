#' @title Maximum likelihood estimation for beta geometric glm
#'
#' @description
#'     Find maximum likelihood parameters for glm assuming geometric
#'     distribution and logit link
#'
#' @param formula an object of class "formula" that describes the model
#' @param data a data.frame where the data live (unlike `glm` and `lm`,
#'             the data object is *required*)
#' @param ... arguments passed to `control` of `optim`
#'
#' @details
#'     Some methods that work for objects of class `glm` and `lm` will
#'     also work for the output of this function, specifically:
#'     `print`, `coef`, `logLik`, and `AIC`. A custom likelihood ratio test
#'     function is provided by `lrt_geom`
#'
#' @returns an object of class `glmBGeom` inspired by, but not strictly
#'     inheriting from, `glm` and `lm`
#'
#' @export

glm_bgeom <- function(formula, data, ...) {
    m <- model.frame(formula, data)

    # model matrix and reponse
    X <- model.matrix(formula, m)
    y <- model.response(m)

    # use yule-simons model for starting param values
    m0 <- suppressWarnings(
        VGAM::vglm(formula, data = data,
                   family = VGAM::yulesimon)
    )

    p0 <- c(coef(m0), 0)

    # mle
    fit <- optim(p0, fn = bgeo_ll, mx = X, my = y,
                 hessian = TRUE,
                 control = list(fnscale = -1, reltol = .Machine$double.eps^0.75,
                                ...))

    if(fit$convergence != 0) {
        stop("bad convergence when maximizing log likelihood function")
    }

    beta_hat <- fit$par
    names(beta_hat)[length(beta_hat)] <- "(Intercept):beta"

    ll <- structure(fit$value, df = ncol(X) + 1, class = "logLik")

    # subset of output from `glm` and `lm`
    res <- list(
        coefficients = beta_hat,
        call = match.call(),
        x = X,
        y = y,
        hess = fit$hessian,
        logLik = ll,
        terms = attr(m, "terms")
    )

    class(res) <- "glmBGeom"

    return(res)



}


#' @export
model.matrix.glmBGeom <- function(object) {
    m <- object$x
    m <- cbind(m, 1)
    colnames(m)[ncol(m)] <- "(Intercept):beta"

    m
}

#' @export
logLik.glmBGeom <- function(object) {
    object$logLik
}


#' @export
print.glmBGeom <- function(object) {
    cat("Call:  ")
    print(object$call)

    cat("\nCoefficients:", "\n")
    print(object$coefficients)

    cat("\nLog likelihood:", "\n")
    print(logLik(object))

}


#' @title Likelihood ratio test for beta-geometric glm
#'
#' @description
#'     Preform likelihood ratio test on two models fitted with `glm_bgeom`
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

lrt_bgeom <- function(mod0, mod1) {
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


#' @title Confidence intervals and prediction for beta-geometric glm's
#'
#' @description
#'     Calculate confidence intervals and prediction
#'
#' @param mod fitted `geomGLM` object
#' @param newdata optional `data.frame` of new predictor variable values
#' @param level confidence level
#' @param type whether to predict response or probability
#'
#' @returns column vector of predictions
#'
#' @export

ci_bgeom <- function(mod, newdata, level = 0.05,
                     type = c("prob", "response")) {
    type <- match.arg(type)

    X <- model.matrix(mod)

    # model matrix for just alpha terms
    Xa <- X[, -ncol(X)]

    # model matrix for just beta terms
    Xb <- X[, ncol(X), drop = FALSE]

    cc <- coef(mod)
    cca <- cc[-ncol(X)] # just alpha coef
    ccb <- cc[ncol(X)] # just beta coef

    # variance co-variance matrix
    vcov_mat <- solve(-mod$hess)

    # linear predictions for alpha and beta
    eta_a <- Xa %*% cca
    eta_b <- Xb %*% ccb
    ab <- eta_a - eta_b

    # prediction for p
    p_hat <- plogis(ab)

    # gradient for delta method se calculation
    grad <- as.vector(p_hat * (1 - p_hat)) * cbind(Xa, -Xb)

    # se of fitted values for p
    # equivalent to `diag(grad %*% vcov_mat %*% t(grad))`
    se_p <- rowSums((grad %*% vcov_mat) * grad) |>
        sqrt()

    # signif cuttoff
    crit <- qnorm(level / 2, lower.tail = FALSE)

    if(type == "prob") {
        p_lwr <- plogis(ab - crit * se_p / (p_hat * (1 - p_hat)))
        p_upr <- plogis(ab + crit * se_p / (p_hat * (1 - p_hat)))

        res <- data.frame(p_hat = 1 - p_hat,
                          p_lwr = 1 - p_lwr,
                          p_upr = 1 - p_upr)
    } else {
        se_y <- se_p / p_hat^2
        y_hat <- 1 / p_hat
        y_lwr <- y_hat - crit * se_y
        y_upr <- y_hat + crit * se_y

        res <- data.frame(y_hat = y_hat,
                          y_lwr = y_lwr,
                          y_upr = y_upr)
    }


    return(res)
}


# non-exported function for calculating log likelihood
#' @param pars vector of parameter values
#' @param mx model matrix
#' @param my model response

bgeo_ll <- function(pars, mx, my) {
    lin_pars <- pars[1:ncol(mx)]
    bet_pars <- pars[ncol(mx) + 1]

    # linear prediction
    fx <- mx %*% lin_pars

    # log-transformed
    a <- exp(fx)
    b <- exp(bet_pars)

    # log probabilities for beta geo dist with support >= 1
    logd <- VGAM::dbetageom(my - 1, a, b, log = TRUE)

    return(sum(logd))
}

