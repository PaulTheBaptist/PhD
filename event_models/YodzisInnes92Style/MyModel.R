# MyModel in R

# Script to run my modified Lawrence code which just computes what we want
require('gruyere')
require('compiler')

if (!is.loaded("MyYodzisInnesState")) dyn.load("myGruyere_model.so")
if (!is.loaded("MyYodzisInnesState")) stop("myGruyere_model.so didn't load...")

MyYodzisInnesFlux <- function (B, params) 
cmpfun({
    res <- .C("MyYodzisInnesState", params$n.species, params$K, 
        params$a, params$q, params$d, params$W, params$producers.c, 
        params$n.producers, params$consumers.c, params$n.consumers, 
        params$rho, params$x, params$y, params$e, params$fe, 
        B,
        #dydt = numeric(params$n.species), My function does not compute dydt
        # whereas Lawrence's YodzisInnesState does
        growth = numeric(params$n.species), 
        respiration = numeric(params$n.species),
        assimilation = numeric(params$n.species * params$n.species),
        consumption = numeric(params$n.species * params$n.species),
        NAOK = TRUE, DUP = FALSE)
    dim(res$assimilation) <- dim(res$consumption) <- c(params$n.species, 
        params$n.species)
    return(res[c("growth", "respiration", "assimilation", "consumption")])
})


YodzisInnesPropensityFunc <- cmpfun(function(B, params)
{
    rates <- MyYodzisInnesFlux(B, params) 
    births <- (rates$growth > 0) * rates$growth + # Producer Births
              colSums(rates$assimilation) # Consumer Births
    deaths <- rowSums(rates$consumption) - # Producer Deaths
            rates$respiration #Consumer Deaths

    # Transformation of propensities
    out <- c(births, deaths, params$immigration) * params$t.prime
    return(out)
})


YodzisInnesNu <- function(M, immigration = TRUE) {
    numspecies <- length(M)
    M_birth <- M_immigration <- diag(numspecies) * M # Creates Matrix with M on diagonal
    M_death <- M_birth * -1
    # So that nu is rows as species, columns as events
    out <- ifelse(immigration, yes = cbind(M_birth, M_death, M_immigration),
        no = cbind(M_birth, M_death))
    
    return (out)
}








