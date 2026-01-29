library(dplyr)
library(ggplot2)
library(cmdstanr)

# read gruner data
arth <- read.csv("inst/raw_data/arth.csv")

# read checklist data
arth_check <- read.csv("inst/raw_data/hawaii_arthropod_checklist.csv")

gen_rich <- filter(arth_check, Genus != "") |>
    group_by(Genus) |>
    summarize(nspp = length(unique(Genus.and.species))) |>
    rename(genus = Genus)


gen_abund <- group_by(arth, genus, species_code) |>
    summarize(abund = n()) |>
    ungroup() |>
    group_by(genus) |>
    summarize(hmean_abund = n() / sum(1 / abund),
              mean_abund = mean(abund))

gen <- inner_join(gen_rich, gen_abund)

ggplot(gen, aes(hmean_abund, nspp)) +
    geom_point() +
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10")


dat <- list(N = nrow(gen),
            y = gen$nspp,
            x = log(gen$hmean_abund) + 1)


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


plot(1 + log(gen$hmean_abund), gen$nspp, log = "y")
b0 <- 40
b1 <- 3
b2 <- -1
curve(1 / (1 + exp(-(b0 + b1 * x + b2 * x^2))), add = TRUE)




