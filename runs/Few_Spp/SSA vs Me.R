# PARAMETERS AND INITIAL CONDITIONS
require('gruyere')

run_ODE <- FALSE
immigration <- TRUE
N0 <- c(1e7, 1e2)
plotfile <- 'SSA_vs_MeA1'

realtf <- 5
data(TL84)
TL84_pm <- PredationMatrix(TL84)
dimnames(TL84_pm) <- list(NULL, NULL)

vertebrate_indx <- TL84$nodes$category == 'vert.ecto'
vertebrates <- TL84$nodes$node[vertebrate_indx]
invert_indx <- TL84$nodes$category == 'invertebrate'

spp_indx <- c(4, 32)
spp <- TL84$nodes$node[spp_indx]


community <- TL84
class(community) <- 'list'
community$nodes <- community$nodes[spp_indx, ]
valid_links <- apply(community$trophic.links, 1, function(link) {
                 2 == sum(spp %in% link)
             })
community$trophic.links <- community$trophic.links[valid_links, ]

community$title <- paste('A few species', Sys.time())
class(community) <- c('Community', 'list')

M <- community$nodes$M
N <- community$nodes$N
pm <- PredationMatrix(community)

x0 <- M * N0

preparams <- list(
  a.constants = YodzisInnes92AConstants(),
  f.constants = AllFConstantsEqual(),
  e.producer = 0.45,
  e.consumer = 0.85,
  fe = 1,
  W = 10^(-3.5), # from Lawrence Paper (means kind of modal values taken from graphs)
  d = 0,
  q = 1,
  K = 10^(-0.2),       # from Lawrence paper
  a = 1,

  # I_rates params
  immigration.c = 1e-14, #1e-12,   # less immigration
  immigration.e = -0.75 + 0.2)


# RUN SCRIPTS
cat('Running SSA.R\n')
source('SSA.R')
cat('SSA.R completed\nRunning mine\n')
source('RunSimpleGillespieWithMyModel.R', chdir = TRUE)
cat('My code completed\n')


# PLOTS AND PRINTS

if (!exists('.newdevs')) assign('.newdevs', FALSE, envir = .GlobalEnv) 
if (!.newdevs && (length(dev.list()) > 1)) dev.off()
if (.newdevs) dev.new()

mysim <- sim$TX
ssasim <- a.sim$data

sim <- cbind(sim$TX[, 1], log10(sim$TX[, -1]))
a.sim$data <- cbind(a.sim$data[, 1], log10(a.sim$data[, -1]))

MinOfMax <- function(a, b, na.rm = TRUE) {
    min(max(a, na.rm = na.rm), max(b, na.rm = na.rm))
}
minmax.t <- MinOfMax(a.sim$data[, 1][a.sim$data[, 1] != Inf],
    sim[, 1][sim[, 1] != Inf])
minmax.B <- MinOfMax(a.sim$data[, -1][a.sim$data[, -1] != Inf],
    sim[, -1][sim[, -1] != Inf])
minld <- min(sim[, -1][sim[, -1] != -Inf],
    a.sim$data[, -1][a.sim$data[, -1] != -Inf], na.rm = TRUE)
zeroval <- minld - abs(minld) * 0.1
sim[which(sim == -Inf)] <- zeroval
a.sim$data[which(a.sim$data == -Inf)] <- zeroval
rant <- c(0, minmax.t)

YLIMS <- c(zeroval, minmax.B)
max.indx <- min(nrow(a.sim$data), nrow(sim))
indices <- 1 : max.indx

pdf(plotfile)
par(mfrow = c(2, 2))                    
matplot(sim[, 1],  sim[, -1],
    type = 'l', xlab = 't', ylab = 'Log10(B)',
    ylim = YLIMS, xlim = rant, col = 1:2, lty = 2)
matplot(a.sim$data[, 1], a.sim$data[, -1], type = 'l', lty = 3,
    col = 3:4, add = TRUE)
points(seq(rant[1], rant[2], length.out = 20), rep(zeroval, 20))


plot(a.sim$data[indices, 1], sim[indices, 1], type = 'l',
     main = 'T.Mine against T.SSA',
     xlab = 'T.SSA', ylab = 'T.Mine')
abline(a = 0, b = 1, lty = 2, col = 'red')

plot(sim[indices, 1] / a.sim$data[indices, 1],
     main = 'Ratio of times (Me/SSA)', type = 'l', xlab = 'event #')
abline(h = 1, col = 'blue', lty = 2)

matplot(sim[indices, -c(1)] / a.sim$data[indices, -c(1)],
        main = 'Ratio of B at each event (Me/SSA)', type = 'l',
        ylab = 'ratio (Me/SSA)', xlab = 'event #')
abline(h = 1, col = 'blue', lty = 2)
dev.off()

print(cbind(sim[(max.indx-10) : max.indx, 1],
            a.sim$data[(max.indx-10) : max.indx, 1],
            sim[(max.indx-10) : max.indx, 1] /
            a.sim$data[(max.indx-10) : max.indx, 1]))
print(c(nrow(mysim), nrow(ssasim)))

# USING GRUYERE TO RUN ODE SIM TO minmax.t


