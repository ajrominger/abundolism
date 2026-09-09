#' @title Plot logistic data as paired histograms
#' @param data simulatin output `data.frame`; must have a `speciation` column
#' @param x name (unquoted) of predictor variable
#' @param xlab label for x axis
#' @param xlog logical, should x-axis be log-transformed (default FALSE)
#' @param formula optional formula to pass to `geom_smooth`

logistic_hist_plot <- function(data, x, xlab, xlog, formula) {
    # browser()
    data$fitted <- glm(formula, data = data, family = "binomial") |>
        predict(type = "response")

    ggplot(data, aes({{ x }})) +
        scale_x_continuous(transform = ifelse(xlog, "log10", "identity")) +
        hist_helper(filter(data, speciation == 0), {{ x }},
                    ymin = min(data$fitted),
                    ymax = max(data$fitted)) +
        hist_helper(filter(data, speciation == 1), {{ x }},
                    ymin = min(data$fitted),
                    ymax = max(data$fitted),
                    top = TRUE) +
        # geom_line(aes(y = fitted)) +
        xlab(xlab) +
        ylab("Speciation")
}

hist_helper <- function(data, x, ymin = 0, ymax = 1, top = FALSE, bins = 30) {
    geom_histogram(
        data = data,
        aes(
            x = {{ x }},
            y = ifelse(top, -1, 1) * 0.3 * abs(ymax - ymin) * stat(count) /
                max(stat(count))
        ),
        bins = bins,
        position = position_nudge(y = ifelse(top, ymax, ymin))
    )
}


