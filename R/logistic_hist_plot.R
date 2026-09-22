#' @title Plot logistic data as paired histograms
#' @param data simulatin output `data.frame`; must have a `speciation` column
#' @param x name (unquoted) of predictor variable
#' @param xlab label for x axis
#' @param xlog logical, should x-axis be log-transformed (default FALSE)
#' @param ymin optional minimum for y-axis
#' @param ymax optional maximum for y-axis
#' @export

logistic_hist_plot <- function(data, x, xlab, xlog, ymin, ymax) {
    # data$fitted <- glm(formula, data = data, family = "binomial") |>
    #     predict(type = "response")

    p <- ggplot2::ggplot(data, aes({{ x }})) +
        ggplot2::scale_x_continuous(transform = ifelse(xlog,
                                                       "log10",
                                                       "identity")) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab("Prob. of Full Speciation") +
        cowplot::theme_cowplot()

    if(missing(ymin)) ymin <- ggplot2::layer_scales(p)$y$get_limits()[1]
    if(missing(ymax)) ymax <- ggplot2::layer_scales(p)$y$get_limits()[2]

    p +
        hist_helper(dplyr::filter(data, speciation == 0), {{ x }},
                    ymin = ymin,
                    ymax = ymax) +
        hist_helper(dplyr::filter(data, speciation == 1), {{ x }},
                    ymin = ymin,
                    ymax = ymax,
                    top = TRUE)
}

hist_helper <- function(data, x, ymin = 0, ymax = 1, top = FALSE, bins = 30) {
    ggplot2::geom_histogram(
        data = data,
        aes(
            x = {{ x }},
            y = ifelse(top, -1, 1) * 0.3 * abs(ymax - ymin) *
                ggplot2::after_stat(count) / max(ggplot2::after_stat(count))
        ),
        bins = bins,
        fill = "gray80",
        position = ggplot2::position_nudge(y = ifelse(top, ymax, ymin))
    )
}


