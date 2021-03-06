require('GillespieSSA')

SSAParams <- function(preparams) {
    stopifnot(preparams$fe == 1)
    with(preparams,
        {
            c(ar = unname(a.constants$ar['producer']),
                aT = unname(a.constants$aT['invertebrate']),
                aJ = unname(a.constants$aJ['invertebrate']),
                fr = unname(f.constants$fr['producer']),
                fT = unname(f.constants$fT['invertebrate']),
                fJ = unname(f.constants$fJ['invertebrate']),
                e.producer = e.producer,
                e.consumer = e.consumer,
                fe = fe,
                W = W,
                d = d,
                q = q,
                K = K,
                a = a,
                i1 = immigration.c * M[1] ^ immigration.e,
                i2 = immigration.c * M[2] ^ immigration.e,
                M1 = M[1],
                M2 = M[2])
        }
    )
}

parms <- SSAParams(preparams)
names(x0) <- c('B1', 'B2')
a <- c('max(fr * ar * M1 ^ (-0.25) * B1 * (1 - B1 / K) / M1, 0)', # P birth
        'B2 * fJ * aJ *M2 ^ (-0.25) * (B1/W)^(1+q)/(1+d*B2+(B1/W)^(1+q)) / M1', # P death
        'i1 / M1', # P immigration
        'B2 * (B1/W)^(1+q)/(1+d*B2+(B1/W)^(1+q)) / M2', # C birth
        'B2 * fT * aT * M2 ^ (-0.25) / M2', # C death
        'i2 / M2') # C immigration

if (immigration) {
    nu <- matrix(0, nrow = 2, ncol = 6,
        dimnames = list(c('sp1', 'sp2'), c('PB', 'PD', 'PI', 'CB', 'CD', 'CI')))
    nu[1, 1:3] <- rep(M[1], 3) * c(1, -1, 1)
    nu[2, 4:6] <- rep(M[2], 3) * c(1, -1, 1)
} else {
    a <- a[c(1, 2, 4, 5)]
    nu <- matrix(0, nrow = 2, ncol = 4,
        dimnames = list(c('sp1', 'sp2'), c('PB', 'PD', 'CB', 'CD')))
    nu[1, 1:2] <- rep(M[1], 2) * c(1, -1)
    nu[2, 3:4] <- rep(M[2], 2) * c(1, -1)
}

a.sim <- ssa(x0, a, nu, parms, tf = Inf, maxWallTime = realtf)
