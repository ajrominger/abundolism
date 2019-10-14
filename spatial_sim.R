library(pika)
library(plyr)
library(socorro)
library(MASS)

bbs <- read.csv('data/bbs/bbs2011.csv', as.is = TRUE)
bbsSpInfo <- read.csv('data/bbs/speciesTableBody.csv', as.is = TRUE)
bbsSpInfo <- bbsSpInfo[, c('sppKey', 'mass')]
names(bbsSpInfo) <- c('spp', 'mass')

x <- merge(bbs, bbsSpInfo, by = 'spp')
x$metab <- x$mass^0.75

summ <- ddply(x, 'spp', function(dat) {
    c(abund = sum(dat$abund), metab = mean(dat$metab))
})
b <- 1 / mean(summ$metab[summ$abund == 1])
S <- 10000
aprop <- 0.01

nn <- rfish(S, 0.01)
rr <- rexp(S, b * nn)
# plot(log(nn), log(rr))
# mod <- lm(log(rr) ~ log(nn))
# abline(mod, col = 'red')
# mod$coeff

nhat <- rnbinom(S, size = 0.01, mu = nn * aprop)
plot(nhat[nhat > 0], rr[nhat > 0], log = 'xy')
modSamp <- lm(log(rr[nhat > 0]) ~ log(nhat[nhat > 0]))
modSamp$coefficients
