library(compiler)

# A simple version of the Gillespie Algorithm
# No fuss, no long functions
# Unit testable and sexy

# a: Propensity vector- same order as state change columns
# nu: State change matrix- rows, state variables (same order as state vector);
#+cols, events
# x: State vector
# t: Time

# r1, r2: vectors of random numbers for event choosing and time choosing

# .TEST: if TRUE, provides variables for testing and runs tests
if (!exists('.TEST')) .TEST <- FALSE
if (.TEST) {
    t0 <- c(t = 0)
    x0 <- c(A = 0, B = 0, C = 0)
    a <- c(0, 0, 1, 0, 3, 0)
    a0 <- sum(a)
    nu <- matrix(byrow = TRUE, nrow = length(x0),
                 data = c(
                 -1,  1,  0,  0,  0,  0,
                  0,  0, 10,-10,  0,  0,
                  0,  0,  0,  0, 10, -10))
    r1 <- c(0, 0.5, 1, 2, 4) / a0
    r2 <- runif(length(r1))
}


# getTau:
#     Inputs:
#         r: ~Unif(0, 1)
#         a0: sum of the propensity vector
#
#     Outputs:
#         tau = -log(r)/a0
getTau <- function(r, a0) return (-log(r) / a0)

if (.TEST) {
    t1 <- apply(X = t(r2), MARGIN = 2,
               FUN = function(r) getTau(r, a0))
    cat('\nTau test:\n    t =\n')
    print(t1)
}


# getU:
#     Inputs:
#         r: ~Unif(0, 1)
#         a: propensity vector
#
#     Outputs:
#         u: such that,
#            sum(a[i-1]) < r * a0 <= sum(a[i])
getU <- cmpfun(function(r, a, a0) {
    i <- 0
    sum <- 0
    decider <- r * a0

    #if (.TEST) {
        #checkPrint <- function() {
        #    print(c(i = i, sum = sum, decider = decider))
        #}

        #checkPrint()
    #}

    # Increment i until a[i] != 0
    while (sum == 0) {
        i <- i + 1
        sum <- sum + a[i]
        
        #if (.TEST) checkPrint()
    }

    # Now choose event
    while (sum < decider) {
        i <- i + 1
        sum <- sum + a[i]

        #if (.TEST) checkPrint()
    }

    # i is equivalent to u
    return (i)
})

if (.TEST) {
    u <- apply(X = t(r1), MARGIN = 2,
               FUN = function(r) {
                   #print(r)
                   u = getU(r, a, a0)
                   #print(c(u = u))
                   return (u)
               })
    cat('\ngetU test:\n    u =\n')
    print(u)
}


# updateState:
#     Inputs:
#         x: state vector
#         t: the system time
#         a: propensity vector
#         nu: state change matrix
#         r1: ~Unif(0, 1) for determining u
#         r2: ~Unif(0, 1) for determining tau
#
#     Outputs:
#         system state: vector of updated time and state c(t, x)
updateState <- cmpfun(function(x, t, a, a0, nu, r1, r2) {

    # Event choosing
    u <- getU(r1, a, a0)
    # Update state vector
    x <- x + nu[, u]

    # Update system time
    t <- t + getTau(r2, a0)

    return (c(t, x))
})

if (.TEST) {
    state_1 <- matrix(NA, ncol = length(t0) + length(x0), nrow = length(r1))
    for (test in seq(length(r1))) {
        #print(test)
        state_1[test, ] <- updateState(x0, t0, a, a0, nu, r1[test], r2[test])
    }
    cat('\nupdateState test:\n    state_1=\n')
    print(state_1)
}

# A toy function, useful in concept for when I formulate more complex versions
# with more advanced stopping conditions
# createTerminator:
#     Relies on termination_conditions existing as a list where each element is
#    +in order of importance to the user list(expression(condition), message),
#    +where condition is made a boolean expression and message is a relevant
#    +descriptive string
#
#     Outputs:
#         Terminator: an expression to be evaluated that
#                   +if eval(termination_conditions)
#                   +super assigns reason for termination to termination_reason
#                   +returns FALSE to terminate the simulation
#                   +else conditions hold and it
#                   +returns TRUE
createTerminator <- function() {
    out <- expression({
        any_broken <- ifelse(length(termination_reason) > 0,
                             yes = TRUE, no = FALSE)

        for (term_test in seq(length(termination_conditions))) {
            broken <- eval(termination_conditions[[term_test]][[1]])

            if (broken) {
                termination_reason <- append(termination_reason,
                    termination_conditions[[term_test]][[2]])
                any_broken <- TRUE
            }
        }
        return (!any_broken)
        })
    return (out)
}


# Gillespie:
#     Inputs:
#         x0
#         t0
#         a_func: takes a_args as arguments to compute the propensity,
#                +must take a_args as arguments
#         a_args: must be whatever a_func takes
#         nu
#         termination_conditions: expression(condition), should be in order of
#                                +importance for the user to know about
#                                +should also super assign (<<-) a descriptive
#                                +string to termination_reason
#                                +so that an appropriate reason is given
#                             +NB checking for 0 propensities is done by default
#         chunk_size: how many events should pass before termination_conditions
#                    +are checked again.
#                    +If negative states was a termination condition,
#                    +e.g. all(x >= 0),
#                    +then the last post-negative states must be removed
#                    +from the last chunk
#         store_every: so many events 
#
#     Outputs:
#         list of the following-
#         TX: Assuming n is the number of the last event stored, it is a
#            +matrix of ( t0 , x0 )
#                       ( |  , |  )
#                       ( tn , xn )
#         info: gives list of-
#              +time_start = time simulation started,
#              +time_end = time simulation ended,
#              +TX1 = inital system state vector
#              +TXn = the last stored system state
#              +n = number of last stored event
#              +termination_reason = description of why it terminated
#TODO add store_conditions, using new knowledge of eval(parse(text = ...))
Gillespie <- cmpfun(function(x0, t0, a_func, a_args, nu,
                      termination_conditions, chunk_size,
                      store_every) {
    TXEMPTYSIZE <- chunk_size * 20 # Rows

    time_start <- Sys.time()
    cat('Starting run at:', time_start, '\n')

    termination_reason <- NULL
    # somewhere for eval(Terminator) to assign a string to
    Terminator <- createTerminator()
    # see createTerminator to understand

    init_state <- state <- c(t0, x0)
    # system state vector always in c(t, x) format
    TX <- TX_empty <- matrix(NA, byrow = TRUE, nrow = TXEMPTYSIZE,
        ncol = length(t0) + length(x0),
        dimnames <- list(NULL, names(c(t0, x0))))
    # having loads of space to fill in avoids growing TX
    TX[1, ] <- init_state
    TX_cur_row <- 1
    
    events_so_far <- 0
    chunk <- 0
    n <- 0 # last event number stored
    while (eval(Terminator)) {
        chunk <- chunk + 1

        r1 <- runif(chunk_size)
        r2 <- runif(chunk_size)

        for (event in seq.int(chunk_size)) {
            a <- do.call(what = a_func, args = eval(a_args))
            a0 <- sum(a)
            if (a0 == 0) {
                termination_reason <- append(termination_reason,
                    'All propensities == 0')
                break
            }
            state <- updateState(x = state[-1], # where x is in state vector
                                 t = state[1], # where t is in state vector
                                 a = a,
                                 a0 = a0,
                                 nu = nu,
                                 r1 = r1[event],
                                 r2 = r2[event])
            events_so_far <- events_so_far + 1

#TODO Use the power of eval(parse(text = ...)) here, as with the termination
#    +conditions above
            if ((events_so_far %% store_every) == 0) {
                TX_cur_row <- TX_cur_row + 1
                TX[TX_cur_row, ] <- state

                if ((TX_cur_row %% TXEMPTYSIZE) == 0) {
                    TX <- rbind(TX, TX_empty)
                }
            }

        }
    }

    # Calculate extra information to be stored
    time_end <- Sys.time()
    rownames(TX) <- NULL
    # Next line concerns evaluating store_conditions
    if (n == 0) n <- events_so_far - (events_so_far %% store_every)

    info <- list(TX1 = TX[1, ],
                 TXn = TX[TX_cur_row, ],
                 n = n,
                 num_chunks = chunk,
                 time_start = time_start,
                 time_end = time_end,
                 time_elapsed = as.numeric(
                     difftime(time_end, time_start, units = "secs"),
                     units = "secs"),
                 termination_reason = termination_reason)

    out <- list(TX = TX[seq(TX_cur_row), ], info = info)
    print(time_end)
    return (out)
})


# Gillespie tests
if (.TEST) {
    x0 <- c(N = 500, P = 100)
    t0 <- c(t = 0)

    a_func <- function(args) {
        # Calculation of the new a would usually occur in the with statement
        #+below

        # # Example of how to adapt to my previous code
        # a <- with(args, {
        #               rates <- ratefunction(B, params)
        #               # so a_args would get
        #               # expression(list(B = state[-1], params = params))
        #               a <- {< modify rates so it becomes a >}
        #               a
        #           }

        # Demonstration of Lotka-Volterra Pred. Prey
        a <- with(args, with(params, {
            unname(c(r * x[1] * (1 - x[1] / K),
              a * x[1] * x[2],
              e * a * x[1] * x[2],
              d * x[2]))
            }))
        return (a)
    }

    a_args <-expression(list(args = list(
        t = state[1], x = state[-1],
        params = list(r = 1.3, K = 1000, e = 0.5, a = 0.01, d = 1.1))))
    # Example of useful a_args, as a may depend on system state

    nu <- matrix(byrow = TRUE, nrow = length(x0),
                 data = c(1, -1, 0, 0,
                          0,  0, 1, -1))

    realtf <- 30

    termination_conditions <- list(
        list(
        # Examples of what you can do, in order of importance to user
        expression(events_so_far >= 1e10), 'Max event limit reached'),
        # an event limit OR
        list(
        expression(difftime(Sys.time(), time_start, units = 'secs') > realtf),
        'Max real time limit reached'),
        # a max. real time limit OR
        list(
        expression(state[1] > 1e6), 'Simulation time limit reached'),
        # a max. simulation time limit OR
        list(
        expression(sum(state[-1]) > 1e6), 'Max system size reached'),
        # a max. system size OR
        list(
        expression(any(state[-1] > 8e4)), 'Max pop size for a species reached'))
        # a max. size one population can get to
    #TODO figure out a good implementation of checking for negative states-

    store_every <- 1
    chunk_size <- 1e5

    sim <- Gillespie(x0, t0, a_func, a_args, nu, termination_conditions,
                     chunk_size, store_every)

    cat("\nGillespie test:\n    sim$info =\n")
    print(sim$info)

    source('FiguringOutSSA.r')

    par(mfrow = c(2, 2))
    poop <- min(nrow(sim.SSA$data), nrow(sim$TX))
    pindx <- seq(poop)

    matplot(x = sim$TX[pindx, 1], y = sim$TX[pindx, -1], type = 'l',
        xlab = 'Time',
        ylab = 'State',
        main = 'My Stochastic LV Pred-Prey')
    legend('topright', legend = c('Prey', 'Predator'), lty = 1, col = 1:2)

    matplot(x = sim.SSA$data[pindx, 1], y = sim.SSA$data[pindx, -1], type = 'l',
        xlab = 'Time',
        ylab = 'State',
        main = 'SSA Stochastic LV Pred-Prey')
    legend('topright', legend = c('Prey', 'Predator'), lty = 1, col = 1:2)

    plot(x = sim$TX[pindx, 'N'], y = sim$TX[pindx, 'P'],
        type = 'l',
        xlab = 'N',
        ylab = 'P',
        main = 'Phase plane')

    plot(x = sim.SSA$data[pindx, 'N'], y = sim.SSA$data[pindx, 'P'],
        type = 'l',
        xlab = 'N',
        ylab = 'P',
        main = 'Phase plane')

    print(cbind(head(sim$TX), head(sim.SSA$data)))
    print('.................................')
    print(cbind(tail(sim$TX), tail(sim.SSA$data)))
    print(c(sim$info$time_elapsed, sim.SSA$stats$elapsedWallTime))
    print(c(sim$info$n, sim.SSA$stats$nSteps))
    eps <- c(Me = sim$info$n / sim$info$time_elapsed,
        SSA = unname(sim.SSA$stats$nSteps / sim.SSA$stats$elapsedWallTime))
    print(c(eps, ratio = eps[1]/eps[2]))
}
























