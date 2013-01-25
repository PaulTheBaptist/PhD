if(!('MyYodzisInnesFlux' %in% ls())) {
    source('../../event_models/YodzisInnes92Style/MyModel.R', chdir = TRUE)
}

if(!('Gillespie' %in% ls())) {
    source('../../gillespie/SimpleGillespie.R', chdir = TRUE)
}

data(TL84)

vertebrate_indx <- TL84$nodes$category == 'vert.ecto'
vertebrates <- TL84$nodes$node[vertebrate_indx]
community <- TL84
class(community) <- 'list'
community$nodes <- community$nodes[!vertebrate_indx, ]
tmp <- apply(community$trophic.links, 1, function(link) {
                 ifelse(any(vertebrates %in% link), FALSE, TRUE)
             })
community$trophic.links <- community$trophic.links[tmp, ]
rm(tmp)
community$title <- paste(community$title, 'without vertebrates.', Sys.time())
class(community) <- c('Community', 'list')
invert_indx <- TL84$nodes$category == 'invertebrate'
M <- community$nodes$M
N <- community$nodes$N
pm <- PredationMatrix(community)

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
  immigration.c = 0.5,
  immigration.e = -0.75 + 0.2)

pars <- do.call(ModelParamsSpec, preparams[-(11:12)])
pars <- IntermediateModelParams(community, pars)
params <- FinalModelParams(community, pars, 1/4)
params <- append(params, list(immigration = preparams[[11]] *
                                    M ^ preparams[[12]])) # Puts it on same scale
params$t.prime <- 1 / params$one.t.prime

a_args <- expression(list(args = list(B = state[-1], params = params)))
a_func <- function(args) {
    with(args, YodzisInnesPropensityFunc(B, params))
}

nu <- YodzisInnesNu(M)

###################################################
# Modifications to SimpleGillespie.R functions for this purpose

# There are rounding errors
# Evidence: by species 35 going down to -1.29247e-26 despite the lowest mass of any
# species being 3.03e-14
# This could possibly be due to adding the wrong stuff... but my tests on my code
# would suggest that u is always chosen correctly and that the propensity function
# works
minM <- min(M) * 0.9
# minM made slightly smaller, in case we accidentally wipe out the smallest species
updateState <- cmpfun(function(x, t, a, a0, nu, r1, r2) {

    # Event choosing
    u <- getU(r1, a, a0)
    # Update state vector
    x <- x + nu[, u]
    x[x < minM] <- 0

    # Update system time
    t <- t + getTau(r2, a0)

    return (c(t, x))
})

updateState <- cmpfun(function(x, t, a, a0, nu, r1, r2) {

    # Event choosing
    u <- getU(r1, a, a0)
    # Update state vector
    x <- x + nu[, u]

    # Update system time
    t <- t + getTau(r2, a0)

    return (c(t, x))
})

# Modified Gillespie() function for my odd storage conditions
OddGillespie <- cmpfun(function(x0, t0, a_func, a_args, nu,
                      chunk_size,
                      store_every,
                      # Store_time = list(how long you want to store for,
                      #                do that storage every so many...,
                      #                time unit ('seconds', 'minutes', 'hours'))
                      store_time = list(length = 30, every = 40, unit = 'seconds'),
                      filename = 'out',
                      plot = TRUE) {
    if (store_time$length > store_time$every) {
        stop('The length stored for must be shorter than how often it stores')
    }
    smh <- switch(store_time$unit, 'seconds' = 1, 'minutes' = 60,
                  'hours' = 60^2, 1)
    store_time$every <- store_time$every * smh
    store_time$length <- store_time$length * smh

    TXEMPTYSIZE <- chunk_size * 10 # Rows

    time_start <- Sys.time()
    print(paste('Starting run at:', time_start))

    init_state <- state <- c(t0, x0)
    # system state vector always in c(t, x) format
    TX_empty <- matrix(NA, byrow = TRUE, nrow = TXEMPTYSIZE,
        ncol = length(t0) + length(x0),
        dimnames = list(NULL, names(c(t0, x0))))
    # having loads of space to fill in avoids "growing" TX, see R inferno for details
    
    events_so_far <- 0
    chunk <- 0
    n <- 0 # last event number stored

    outnum <- 0
    store <- FALSE
    while (TRUE) {
        time_in <- as.numeric(Sys.time()) -  as.numeric(time_start)
        if (time_in %% store_time$every < store_time$length) {
            if (!store) {
                outnum <- outnum + 1
                TX <- TX_empty
                TX_cur_row <- 0
                store_t <- Sys.time()
                print(paste('Storing from:', store_t))
                outfile <- paste(filename, "_", outnum, "_",
                             format(store_t, "%y%m%d"),
                             ".csv", sep = '')
                store <- TRUE
            }
        } else {
            if (store) {
                write.csv(TX[1:TX_cur_row, ], outfile)
                if (any(is.na(TX[1:TX_cur_row, ]))) {
                       stop('Something became NA in TX')
                }
                print(paste('Store complete at:', Sys.time()))
                print(paste(TX_cur_row * store_every, 'events occurred'))
                print(TX[TX_cur_row, ], digits = 2)
                if (plot) {
                    pdf(file = 'CurrentStoreState.pdf')
                    par(mfrow = c(1, 2))
                    logdat <- log10(TX[1:TX_cur_row, -1])
                    minld <- min(logdat[logdat != -Inf])
                    zeroval <- minld - abs(minld) * 0.1
                    logdat[which(logdat == -Inf)] <- zeroval
                    rant <- range(TX[1:TX_cur_row, 1])
                    TL84B <- log10(M * N)
                    matplot(TX[1:TX_cur_row, 1],  logdat,
                            type = 'l', xlab = 't', ylab = 'Log10(B)',
                            ylim = range(logdat, TL84B))
                    points(seq(rant[1], rant[2], length.out = 20),
                           rep(zeroval, 20))
                    matplot(rant, matrix(rep(TL84B, 2), nrow = 2, byrow = TRUE),
                            add = FALSE, ylim = range(logdat, TL84B),
                            type = 'l', ylab = '', xlab = 't')
                    dev.off()
                }
            }
            store <- FALSE
        }

        chunk <- chunk + 1

        r1 <- runif(chunk_size)
        r2 <- runif(chunk_size)
        if (events_so_far == 0) {
            tau1 <- getTau(r2[1],
                         a0 = sum(do.call(what = a_func, args = eval(a_args))))
            print(paste('The first event occurred at',
                        t0 + tau1))
            #state[1] <- t0 - tau1
            #print('t of event 1 has been altered to equal t0')
        }

        for (event in seq.int(chunk_size)) {
            a <- do.call(what = a_func, args = eval(a_args))
            a0 <- sum(a)

            state <- updateState(x = state[-1], # where x is in state vector
                                 t = state[1], # where t is in state vector
                                 a = a,
                                 a0 = a0,
                                 nu = nu,
                                 r1 = r1[event],
                                 r2 = r2[event])
            events_so_far <- events_so_far + 1

            if (store) {
                if ((events_so_far %% store_every) == 0) {
                    TX_cur_row <- TX_cur_row + 1
                    TX[TX_cur_row, ] <- state

                    if ((TX_cur_row %% TXEMPTYSIZE) == 0) {
                        TX <- rbind(TX, TX_empty)
                    }
                }
            }
        }
    }
})


new <- TRUE
if (new) {
OddGillespie(x0 = rep(0, params$n.species),
          t0 = 0,
          a_func = a_func,
          a_args = a_args,
          nu = nu,
          chunk_size = 5e4,
          store_every = 1e2,
          store_time = list(length = 15, every = 60, unit = 'seconds'),
          filename = 'out')
} else {print('nuthin')}







