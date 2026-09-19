#' @title Plot logistic data as paired histograms
#' @param data simulatin output `data.frame`; must have a `speciation` column
#' @param x name (unquoted) of predictor variable
#' @param xlab label for x axis
#' @param xlog logical, should x-axis be log-transformed (default FALSE)
#' @param ymin optional minimum for y-axis
#' @param ymax optional maximum for y-axis
#' @param formula optional formula to pass to `geom_smooth`
#' @export

logistic_hist_plot <- function(data, x, xlab, xlog, ymin, ymax, formula) {
    # data$fitted <- glm(formula, data = data, family = "binomial") |>
    #     predict(type = "response")

    p <- ggplot(data, aes({{ x }})) +
        scale_x_continuous(transform = ifelse(xlog, "log10", "identity")) +
        geom_smooth(mapping = aes({{ x }}, speciation),
                    method = "glm",
                    formula = formula,
                    method.args = list(family = "binomial"),
                    color = "#006B99", fill = "#006B99", alpha = 0.5) +
        xlab(xlab) +
        ylab("Prob. of Full Speciation") +
        cowplot::theme_cowplot()

    if(missing(ymin)) ymin <- layer_scales(p)$y$get_limits()[1]
    if(missing(ymax)) ymax <- layer_scales(p)$y$get_limits()[2]

    p <- p +
        hist_helper(filter(data, speciation == 0), {{ x }},
                    ymin = ymin,
                    ymax = ymax) +
        hist_helper(filter(data, speciation == 1), {{ x }},
                    ymin = ymin,
                    ymax = ymax,
                    top = TRUE)
    p$layers <- p$layers[c(2, 3, 1)]

    p
}

hist_helper <- function(data, x, ymin = 0, ymax = 1, top = FALSE, bins = 30) {
    geom_histogram(
        data = data,
        aes(
            x = {{ x }},
            y = ifelse(top, -1, 1) * 0.3 * abs(ymax - ymin) *
                after_stat(count) / max(after_stat(count))
        ),
        bins = bins,
        fill = "gray80",
        position = position_nudge(y = ifelse(top, ymax, ymin))
    )
}


