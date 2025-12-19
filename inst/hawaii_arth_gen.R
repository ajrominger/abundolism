library(dplyr)
library(ggplot2)
library(rstan)

# read gruner data
arth <- read.csv("inst/raw_data/arth.csv")

# read checklist dadta
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

