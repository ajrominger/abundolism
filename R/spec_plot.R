#' @title Custom plotting routine for speciation as response variable
#'
#' @description
#'     Plot speciation against a predictor variable with glm prediction
#'
#' @param dat simulatin output `data.frame`
#' @param dat output from simulation
#' @param x name (unquoted) of predictor variable
#' @param xlab label for x axis
#'
#'
#' @details
#'     Not a general-purpose function, use only for intended scenarios
#'
#' @returns a `ggplot` object
#'
#' @export

spec_plot <- function(dat, x, xlab,
                      xlog = FALSE, formula = NULL) {
    ggplot(dat,
           aes({{ x }},
               speciation + ifelse(speciation == 1, -1, 1) *
                   runif(nrow(dat), 0, 0.05))) +
        geom_pointdensity(method = "kde2d") +
        scale_shape_binned() +
        scale_x_continuous(transform = ifelse(xlog, "log10", "identity")) +
        scale_color_viridis_c(
            transform = "log",
            breaks = \(x) {
                log10(x) |>
                    (\(x) {
                        seq(floor(x[1]), ceiling(x[2]), by = 0.2)
                    })() |>
                    round(1) |>
                    (\(x) 10^x)() |>
                    round()
            }
        ) +
        xlab(xlab) +
        ylab("Probability of speciation") +
        geom_smooth(mapping = aes({{ x }}, speciation),
                    method = "glm",
                    formula = formula,
                    method.args = list(family = "binomial"),
                    color = "black")
}
