#' Custom predictions from a binomial glmmTMB model object
#'
#' Makes prediction on response scale across an evenly spaced sequence of
#' fixed effect in `object`, then generates response prediction and SE
#' without random effects (i.e. population-level sensu `?glmmTMB`).
#'
#' @param object A fitted `glmmTMB` model object of binomial family.
#' @param length_out The length of the evenly spaced sequence
#'   generated for each non-fixed predictor. Defaults to 50.
#'
#' @return A data frame combining the prediction grid (`newdata`) with
#'         the predicted probability (`fit`) and its standard error
#'         (`se_fit`)
#'
#' @export

pred_binom <- function(object, length_out = 50) {
    mf <- model.frame(object)

    pred_vars <- formula(object, fixed.only = TRUE) |>
        terms() |>
        delete.response() |>
        all.vars()

    var_seqs <- lapply(pred_vars, function(v) {
        seq(min(mf[[v]]), max(mf[[v]]), length.out = length_out)
    })
    names(var_seqs) <- pred_vars

    newdata <- do.call(expand.grid, var_seqs)

    pred <- predict(object, newdata = newdata, type = "response",
                    se.fit = TRUE, re.form = NA)

    cbind(newdata, fit = pred$fit, se_fit = pred$se.fit)
}
