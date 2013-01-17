if(!(exists('MyYodzisInnesFlux'))) source('MyModel.R')
if(!(exists('Gillespie'))) source('../../gillespie/SimpleGillespie.R')

preparams <- preparams
pars <- do.call(ModelParamsSpec, preparams[-(11:12)])
pars <- IntermediateModelParams(community, pars)
params <- BuildModelParams(community, pars, 1/4)
params <- append(params, list(immigration = preparams[[11]] *
                                    M ^ preparams[[12]] * 
                                    params$one.t.prime)) # Puts it on same scale
params$t.prime <- (params$one.t.prime)^(-1)

a_args <- expression(list(args = list(B = state[-1], params = params)))
a_func <- function(args) {
    with(args, YodzisInnesPropensityFunc(B, params))
}

nu <- YodzisInnesNu(M)

###################################################
# Modifications to SimpleGillespie.R functions

# There are rounding errors
# Evidence: species 35 going down to -1.29247e-26 despite the lowest mass of any
# species being 3.03e-14
# This could possibly be due to adding the wrong stuff... but tests on my code
# would suggest that u is always chosen correctly and that the propensity function
# works
#minM <- min(M) * 0.9
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
####################################################

termination_conditions <- list(
    list(
        expression(difftime(Sys.time(), time_start, units = 'secs') > realtf),
        'Max real time limit reached'))
# a max. real time limit

sim <- Gillespie(x0, t0 = 0, a_func, a_args, nu, termination_conditions,
                     chunk_size = 1e5, store_every = 1)






